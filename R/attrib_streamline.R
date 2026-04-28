#' Generate and attribute stream lines and points
#'
#' @noRd
attrib_streamline <- function(
  input,
  extra_attr = c(
    "link_slope",
    "cont_slope",
    "USChnLn_To",
    "Elevation",
    "StOrd_Hack",
    "StOrd_Str",
    "StOrd_Hort",
    "StOrd_Shr"
  ),
  points = NULL,
  site_id_col = NULL,
  snap_distance = NULL,
  break_on_noSnap = TRUE,
  return_products = FALSE,
  temp_dir = NULL,
  compress = FALSE,
  verbose = FALSE
) {
  check_ihydro(input)
  extra_attr <- match.arg(extra_attr, several.ok = TRUE)

  # Validate point-related args

  if (!is.null(points)) {
    validate_site_id_col(site_id_col)
    if (is.null(snap_distance)) {
      snap_distance <- Inf
    }
    if (!is.numeric(snap_distance)) {
      cli::cli_abort("{.arg snap_distance} must be numeric.")
    }
    stopifnot(is.logical(break_on_noSnap))
  } else {
    site_id_col <- "link_id"
  }
  stopifnot(
    is.logical(return_products),
    is.logical(verbose),
    is.logical(compress)
  )

  temp_dir <- ensure_temp_dir(temp_dir)
  db_fp <- input$outfile

  # Configure tools

  whitebox::wbt_options(
    exe_path = whitebox::wbt_exe_path(),
    verbose = verbose > 2,
    wd = temp_dir,
    compress_rasters = compress
  )
  terra::terraOptions(verbose = verbose > 3, tempdir = temp_dir)
  gdal_arg <- if (compress) "COMPRESS=NONE" else NULL

  # Extract rasters from gpkg to temp
  for (lyr in c("dem_d8", "dem_streams_d8_sub", "dem_final")) {
    terra::writeRaster(
      read_ihydro(input, lyr),
      file.path(temp_dir, paste0(lyr, ".tif")),
      overwrite = TRUE,
      gdal = gdal_arg,
      datatype = "FLT4S"
    )
  }

  dem_final <- terra::rast(file.path(temp_dir, "dem_final.tif"))
  names(dem_final) <- "Elevation"
  terra::writeRaster(
    dem_final,
    file.path(temp_dir, "Elevation.tif"),
    overwrite = TRUE,
    gdal = gdal_arg,
    datatype = "FLT4S"
  )

  # ── Compute stream link attributes ──────────────────────────────────────
  if (verbose) {
    message("Calculating stream link attributes")
  }
  compute_wbt_stream_attrs(temp_dir, extra_attr)

  # ── Vectorise stream links ──────────────────────────────────────────────
  whitebox::wbt_raster_streams_to_vector(
    streams = "link_id.tif",
    d8_pntr = "dem_d8.tif",
    output = "strm_link_id.shp"
  )
  strm <- sf::read_sf(file.path(temp_dir, "strm_link_id.shp")) |>
    dplyr::select(STRM_VAL)
  sf::st_crs(strm) <- terra::crs(dem_final)
  colnames(strm)[1] <- "link_id"

  # ── Extract point-level attributes ──────────────────────────────────────
  if (verbose) {
    message("Extracting stream link attributes")
  }
  final_points <- extract_stream_point_attrs(temp_dir, extra_attr)

  # ── Insert sampling points ──────────────────────────────────────────────
  snapped_points <- NULL
  original_points <- NULL

  if (!is.null(points)) {
    if (verbose) {
      message("Snapping points")
    }

    snap_result <- snap_points_to_network(
      points,
      final_points,
      site_id_col,
      snap_distance,
      break_on_noSnap,
      dem_final,
      temp_dir
    )
    final_points <- snap_result$final_points
    snapped_points <- snap_result$snapped_points
    original_points <- snap_result$original_points
    strm <- snap_result$strm
  }

  sf::write_sf(strm, file.path(temp_dir, "stream_lines.shp"))

  # ── Identify upstream/downstream links ──────────────────────────────────
  if (verbose) {
    message("Identifying upstream links")
  }
  final_us <- identify_adjacent_links(
    final_points,
    temp_dir,
    direction = "upstream"
  )

  if (verbose) {
    message("Identifying downstream links")
  }
  final_ds <- identify_adjacent_links(
    final_points,
    temp_dir,
    direction = "downstream"
  )

  # ── Assemble final link table ───────────────────────────────────────────
  links <- final_points |>
    dplyr::group_by(link_id) |>
    dplyr::filter(USChnLn_Fr == max(USChnLn_Fr)) |>
    dplyr::ungroup()

  final_links <- links |>
    dplyr::full_join(final_us, by = c("link_id", "trib_id")) |>
    dplyr::full_join(final_ds, by = c("link_id", "trib_id"))

  if (any(duplicated(final_links$link_id)) || any(duplicated(final_links$ID))) {
    warning("Duplicated link_id/ID detected; check upstream/downstream IDs.")
  }

  # ── Write to GeoPackage ─────────────────────────────────────────────────
  sf::write_sf(
    final_links |> dplyr::select(link_id),
    file.path(temp_dir, "stream_links.shp")
  )
  sf::write_sf(
    final_points |> dplyr::select(ID),
    file.path(temp_dir, "stream_points.shp")
  )

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  dplyr::copy_to(
    con,
    final_links |> tibble::as_tibble() |> dplyr::select(-geometry),
    "stream_links_attr",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = c("link_id", "trib_id"),
    analyze = TRUE,
    in_transaction = TRUE
  )

  dplyr::copy_to(
    con,
    final_points |> tibble::as_tibble() |> dplyr::select(-geometry),
    "stream_points_attr",
    overwrite = TRUE,
    temporary = FALSE,
    unique_indexes = "ID",
    indexes = c("link_id", "trib_id"),
    analyze = TRUE,
    in_transaction = TRUE
  )

  dplyr::copy_to(
    con,
    tibble::tibble(site_id_col = site_id_col),
    "site_id_col",
    overwrite = TRUE,
    temporary = FALSE,
    analyze = TRUE,
    in_transaction = TRUE
  )

  if (verbose) {
    message("Writing vector layers to GeoPackage")
  }
  for (lyr in c(
    "snapped_points",
    "original_points",
    "stream_links",
    "stream_lines",
    "stream_points"
  )) {
    shp <- file.path(temp_dir, paste0(lyr, ".shp"))
    if (file.exists(shp)) {
      sf::write_sf(
        sf::read_sf(shp),
        db_fp,
        layer = lyr,
        delete_layer = TRUE,
        append = TRUE
      )
    }
  }

  # ── Build return object ─────────────────────────────────────────────────
  output <- input[
    !names(input) %in%
      c("stream_lines", "links", "points", "snapped_points", "original_points")
  ]

  if (return_products) {
    output <- c(
      list(
        stream_lines = strm,
        links = final_links,
        points = final_points,
        snapped_points = snapped_points,
        original_points = original_points
      ),
      output
    )
  }

  clean_temp_files(temp_dir)
  structure(output, class = "ihydro")
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Validate site_id_col
#' @noRd
validate_site_id_col <- function(site_id_col) {
  if (!is.character(site_id_col) || length(site_id_col) != 1) {
    cli::cli_abort("{.arg site_id_col} must be a single character string.")
  }
  if (site_id_col == "link_id") {
    cli::cli_abort("{.arg site_id_col} cannot be {.val link_id}.")
  }
}

#' Run whitebox stream attribution tools
#' @noRd
compute_wbt_stream_attrs <- function(temp_dir, extra_attr) {
  # Essential attributes
  whitebox::wbt_stream_link_identifier(
    d8_pntr = "dem_d8.tif",
    streams = "dem_streams_d8_sub.tif",
    output = "link_id.tif"
  )
  whitebox::wbt_stream_link_class(
    d8_pntr = "dem_d8.tif",
    streams = "dem_streams_d8_sub.tif",
    output = "link_class.tif"
  )
  whitebox::wbt_stream_link_length(
    d8_pntr = "dem_d8.tif",
    linkid = "link_id.tif",
    output = "link_lngth.tif"
  )
  whitebox::wbt_tributary_identifier(
    d8_pntr = "dem_d8.tif",
    streams = "dem_streams_d8_sub.tif",
    output = "trib_id.tif"
  )
  whitebox::wbt_farthest_channel_head(
    d8_pntr = "dem_d8.tif",
    streams = "dem_streams_d8_sub.tif",
    output = "USChnLn_Fr.tif"
  )

  # Optional attributes
  attr_map <- list(
    link_slope = list(
      fn = whitebox::wbt_stream_link_slope,
      args = list(
        d8_pntr = "dem_d8.tif",
        linkid = "link_id.tif",
        dem = "Elevation.tif",
        output = "link_slope.tif"
      )
    ),
    cont_slope = list(
      fn = whitebox::wbt_stream_slope_continuous,
      args = list(
        d8_pntr = "dem_d8.tif",
        streams = "dem_streams_d8_sub.tif",
        dem = "Elevation.tif",
        output = "cont_slope.tif"
      )
    ),
    StOrd_Hack = list(
      fn = whitebox::wbt_hack_stream_order,
      args = list(
        streams = "dem_streams_d8_sub.tif",
        d8_pntr = "dem_d8.tif",
        output = "StOrd_Hack.tif"
      )
    ),
    StOrd_Str = list(
      fn = whitebox::wbt_strahler_stream_order,
      args = list(
        streams = "dem_streams_d8_sub.tif",
        d8_pntr = "dem_d8.tif",
        output = "StOrd_Str.tif"
      )
    ),
    StOrd_Hort = list(
      fn = whitebox::wbt_horton_stream_order,
      args = list(
        streams = "dem_streams_d8_sub.tif",
        d8_pntr = "dem_d8.tif",
        output = "StOrd_Hort.tif"
      )
    ),
    StOrd_Shr = list(
      fn = whitebox::wbt_shreve_stream_magnitude,
      args = list(
        streams = "dem_streams_d8_sub.tif",
        d8_pntr = "dem_d8.tif",
        output = "StOrd_Shr.tif"
      )
    ),
    USChnLn_To = list(
      fn = whitebox::wbt_length_of_upstream_channels,
      args = list(
        d8_pntr = "dem_d8.tif",
        streams = "dem_streams_d8_sub.tif",
        output = "USChnLn_To.tif"
      )
    )
  )

  for (attr_name in extra_attr) {
    if (attr_name %in% names(attr_map)) {
      entry <- attr_map[[attr_name]]
      do.call(entry$fn, entry$args)
    }
  }
}

#' Extract point-level stream attributes from raster layers
#' @noRd
extract_stream_point_attrs <- function(temp_dir, extra_attr) {
  id_rast <- terra::rast(file.path(temp_dir, "link_id.tif"))
  id_pts <- terra::as.points(id_rast)
  id_sf <- sf::st_as_sf(id_pts)

  # Core attribute rasters
  core_files <- c(
    "link_class.tif",
    "link_lngth.tif",
    "trib_id.tif",
    "USChnLn_Fr.tif"
  )
  extra_files <- paste0(extra_attr, ".tif")
  extra_files <- extra_files[file.exists(file.path(temp_dir, extra_files))]
  all_files <- c(core_files, extra_files)

  attr_rasts <- lapply(file.path(temp_dir, all_files), terra::rast)
  names(attr_rasts) <- gsub("\\.tif$", "", all_files)

  extracted <- suppressWarnings(suppressMessages(
    terra::extract(do.call(c, unname(attr_rasts)), id_pts)
  ))

  extracted <- extracted |>
    dplyr::mutate(
      link_type = dplyr::case_when(
        link_class == 1 ~ "Exterior Link",
        link_class == 2 ~ "Interior Link",
        link_class == 3 ~ "Source Node (head water)",
        link_class == 4 ~ "Link Node",
        link_class == 5 ~ "Sink Node"
      )
    ) |>
    tibble::tibble() |>
    dplyr::select(ID, link_class, link_type, tidyselect::everything())

  names(extracted) <- abbreviate(names(extracted), 10)
  dplyr::bind_cols(id_sf, extracted)
}

#' Snap sampling points to stream network and split reaches
#' @noRd
snap_points_to_network <- function(
  points,
  final_points,
  site_id_col,
  snap_distance,
  break_on_noSnap,
  dem_final,
  temp_dir
) {
  points <- process_input(
    points,
    align_to = terra::vect(utils::head(final_points[, 1]))
  )

  original_sf <- sf::st_as_sf(points) |>
    dplyr::select(tidyselect::any_of(site_id_col), tidyselect::everything())
  sf::write_sf(original_sf, file.path(temp_dir, "original_points.shp"))

  # Collapse duplicates to centroids
  points_sf <- sf::st_as_sf(points) |>
    dplyr::filter(!sf::st_is_empty(sf::st_as_sf(points))) |>
    dplyr::group_by(!!rlang::sym(site_id_col)) |>
    dplyr::summarize(
      dplyr::across(tidyselect::contains("geom"), sf::st_union)
    ) |>
    sf::st_centroid() |>
    dplyr::ungroup()

  if (!site_id_col %in% names(points_sf)) {
    cli::cli_abort("{.arg site_id_col} must be a column in {.arg points}.")
  }

  # Snap to nearest interior/exterior link point
  snap_targets <- final_points |>
    dplyr::filter(
      !is.na(link_type),
      !link_type %in%
        c(
          "Link Node",
          "Sink Node",
          "Source Node (head water)"
        )
    ) |>
    dplyr::group_by(link_id) |>
    dplyr::filter(USChnLn_Fr != max(USChnLn_Fr)) |>
    dplyr::select(ID, link_id)

  snapped <- suppressWarnings(
    suppressMessages(
      sf::st_join(
        points_sf,
        snap_targets,
        join = nngeo::st_nn,
        k = 1,
        maxdist = snap_distance,
        progress = TRUE
      )
    )
  )

  snapped <- snapped |>
    tibble::as_tibble() |>
    dplyr::select(-geometry) |>
    dplyr::left_join(
      final_points |> dplyr::select(ID, link_id),
      by = c("ID", "link_id")
    ) |>
    sf::st_as_sf()

  # Handle un-snappable points
  not_snapped <- is.na(snapped$link_id)
  if (any(not_snapped)) {
    failed_ids <- snapped[[site_id_col]][not_snapped]
    msg <- paste0("Could not snap: ", paste(failed_ids, collapse = ", "))
    if (break_on_noSnap) cli::cli_abort(msg) else warning(msg)
  }

  snapped <- snapped |>
    dplyr::filter(!is.na(link_id)) |>
    dplyr::select(tidyselect::any_of(site_id_col), tidyselect::everything()) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(dplyr::across(
      tidyselect::everything(),
      ~ utils::head(., 1)
    )) |>
    dplyr::ungroup()

  # Create new sub-link IDs at snap locations
  new_pts <- final_points
  new_pts$link_class[new_pts$ID %in% snapped$ID] <- 6
  new_pts$link_type[new_pts$ID %in% snapped$ID] <- "Sample Point"

  new_pts <- new_pts |>
    dplyr::left_join(
      snapped |>
        tibble::as_tibble() |>
        dplyr::select(ID, tidyselect::any_of(site_id_col)),
      by = "ID"
    ) |>
    dplyr::filter(!is.na(!!rlang::sym(site_id_col))) |>
    dplyr::group_by(link_id) |>
    dplyr::arrange(link_id, dplyr::desc(USChnLn_Fr)) |>
    dplyr::mutate(
      link_id_new = dplyr::row_number(),
      link_id_new = formatC(
        link_id_new,
        width = nchar(max(link_id_new)),
        format = "d",
        flag = 0
      ),
      link_id_new = as.numeric(paste0(link_id, ".", link_id_new))
    )

  final_points <- final_points |>
    dplyr::filter(!ID %in% new_pts$ID) |>
    dplyr::bind_rows(new_pts) |>
    dplyr::group_by(link_id) |>
    dplyr::arrange(link_id, dplyr::desc(USChnLn_Fr)) |>
    dplyr::mutate(
      link_id_new = dplyr::case_when(
        dplyr::row_number() == 1 & is.na(link_id_new) ~ link_id,
        TRUE ~ link_id_new
      )
    ) |>
    tidyr::fill(link_id_new, .direction = "down") |>
    dplyr::select(-link_id) |>
    dplyr::rename(link_id = link_id_new) |>
    dplyr::select(link_id, tidyselect::everything()) |>
    dplyr::ungroup() |>
    dplyr::arrange(ID)

  snapped_out <- final_points |>
    dplyr::select(link_id, tidyselect::any_of(site_id_col)) |>
    dplyr::filter(
      !is.na(link_id),
      !dplyr::if_any(tidyselect::any_of(site_id_col), is.na)
    )
  sf::write_sf(snapped_out, file.path(temp_dir, "snapped_points.shp"))

  # Rebuild stream raster with new link IDs
  strm_vect <- final_points |>
    dplyr::select(link_id) |>
    dplyr::arrange(link_id) |>
    terra::vect()
  terra::writeVector(strm_vect, file.path(temp_dir, "all_points.shp"))

  whitebox::wbt_vector_points_to_raster(
    input = file.path(temp_dir, "all_points.shp"),
    output = file.path(temp_dir, "new_stream_layer.tif"),
    field = "link_id",
    assign = "last",
    nodata = TRUE,
    base = file.path(temp_dir, "dem_final.tif")
  )
  whitebox::wbt_raster_streams_to_vector(
    streams = file.path(temp_dir, "new_stream_layer.tif"),
    d8_pntr = file.path(temp_dir, "dem_d8.tif"),
    output = file.path(temp_dir, "new_stream_layer.shp")
  )

  # Fix rounding from wbt_raster_streams_to_vector
  un_id <- unique(final_points$link_id)
  strm <- sf::read_sf(file.path(temp_dir, "new_stream_layer.shp")) |>
    dplyr::select(STRM_VAL) |>
    dplyr::rename(link_id = STRM_VAL) |>
    dplyr::rowwise() |>
    dplyr::mutate(link_id = un_id[which.min(abs(link_id - un_id))]) |>
    dplyr::ungroup()
  sf::st_crs(strm) <- terra::crs(dem_final)

  list(
    final_points = final_points,
    snapped_points = snapped_out,
    original_points = points_sf,
    strm = strm
  )
}

#' Identify upstream or downstream adjacent links via D8 flow direction
#' @noRd
identify_adjacent_links <- function(
  final_points,
  temp_dir,
  direction = c("upstream", "downstream")
) {
  direction <- match.arg(direction)

  st_r <- terra::rast(file.path(temp_dir, "dem_streams_d8_sub.tif"))
  d8_pntr <- terra::rast(file.path(temp_dir, "dem_d8.tif"))

  # Select nodes: upstream = min USChnLn_Fr, downstream = max
  if (direction == "upstream") {
    nodes <- final_points |>
      dplyr::group_by(link_id) |>
      dplyr::filter(USChnLn_Fr == min(USChnLn_Fr)) |>
      dplyr::ungroup()
  } else {
    nodes <- final_points |>
      dplyr::group_by(link_id) |>
      dplyr::filter(USChnLn_Fr == max(USChnLn_Fr)) |>
      dplyr::ungroup()
  }

  node_vect <- nodes |>
    dplyr::select(ID) |>
    dplyr::arrange(ID) |>
    terra::vect()
  cv <- terra::cells(d8_pntr, node_vect)

  all_cell <- terra::cells(d8_pntr, terra::vect(final_points))
  all_cell <- tibble::as_tibble(
    dplyr::mutate(final_points, cell = all_cell[, 2])
  )

  # Get queen neighbours
  adj_names <- c("TL", "TM", "TR", "ML", "MR", "BL", "BM", "BR")
  cv1 <- tibble::tibble(target = cv[, 2]) |>
    dplyr::bind_cols(
      terra::adjacent(d8_pntr, cells = cv[, 2], directions = "queen") |>
        as.data.frame() |>
        stats::setNames(adj_names)
    ) |>
    tidyr::gather("direction", "cell_num", -target) |>
    dplyr::arrange(target) |>
    dplyr::filter(!is.na(cell_num), !is.nan(cell_num))

  # Filter to cells on stream
  cv1 <- cv1 |>
    dplyr::mutate(
      on_stream = data.frame(
        terra::extract(st_r, terra::xyFromCell(st_r, cell_num))
      )$dem_streams_d8_sub
    ) |>
    dplyr::filter(on_stream == 1)

  # Join link/trib IDs
  cv1 <- cv1 |>
    dplyr::left_join(
      all_cell |> dplyr::select(link_id, cell),
      by = c("cell_num" = "cell")
    ) |>
    dplyr::left_join(
      all_cell |> dplyr::select(trib_id, cell),
      by = c("cell_num" = "cell")
    ) |>
    dplyr::left_join(
      all_cell |> dplyr::select(target_link_id = link_id, cell),
      by = c("target" = "cell")
    ) |>
    dplyr::left_join(
      all_cell |> dplyr::select(target_trib_id = trib_id, cell),
      by = c("target" = "cell")
    )

  # Determine flow direction
  if (direction == "upstream") {
    # Neighbours whose flow points INTO the target cell
    flow_dir_col <- "cell_num"
    flow_map <- c(
      TL = 4,
      TM = 8,
      TR = 16,
      ML = 2,
      MR = 32,
      BL = 1,
      BM = 128,
      BR = 64
    )
  } else {
    # Target cell flow direction points to neighbour
    flow_dir_col <- "target"
    flow_map <- c(
      TL = 64,
      TM = 128,
      TR = 1,
      ML = 32,
      MR = 2,
      BL = 16,
      BM = 8,
      BR = 4
    )
  }

  cv2 <- cv1 |>
    dplyr::mutate(
      flow_dir = terra::extract(
        d8_pntr,
        terra::xyFromCell(d8_pntr, .data[[flow_dir_col]])
      )$dem_d8
    ) |>
    dplyr::mutate(
      flow_in = purrr::map2_lgl(direction, flow_dir, function(dir, fd) {
        !is.na(fd) && !is.na(flow_map[dir]) && fd == flow_map[dir]
      })
    ) |>
    dplyr::group_by(target) |>
    dplyr::mutate(n_inflows = sum(flow_in)) |>
    dplyr::ungroup()

  # Build result table
  prefix <- if (direction == "upstream") "us" else "ds"

  result <- cv2 |>
    dplyr::filter(flow_in) |>
    dplyr::select(
      n_inflows,
      link_id,
      trib_id,
      target_link_id,
      target_trib_id
    ) |>
    dplyr::group_by(target_link_id, target_trib_id) |>
    dplyr::mutate(nm = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::rename(
      !!paste0(prefix, "trib_id") := trib_id,
      !!paste0(prefix, "link_id") := link_id
    ) |>
    tidyr::pivot_wider(
      id_cols = c(target_link_id, target_trib_id),
      names_from = nm,
      names_sep = "",
      values_from = c(
        !!paste0(prefix, "trib_id"),
        !!paste0(prefix, "link_id")
      )
    ) |>
    dplyr::rename(link_id = target_link_id, trib_id = target_trib_id)

  names(result) <- abbreviate(names(result), 10)
  result
}
