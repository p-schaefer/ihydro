#' Extract polygon subbasins from flow direction products
#'
#' @noRd
generate_subbasins <- function(
  input,
  points,
  return_products = FALSE,
  temp_dir = NULL,
  verbose = FALSE
) {
  check_ihydro(input)
  stopifnot(is.logical(return_products), is.logical(verbose))

  temp_dir <- ensure_temp_dir(temp_dir)
  db_fp <- input$outfile

  whitebox::wbt_options(
    exe_path = whitebox::wbt_exe_path(),
    verbose = verbose > 2,
    wd = temp_dir
  )
  terra::terraOptions(verbose = verbose > 3, tempdir = temp_dir)

  # Extract required rasters

  for (lyr in c("dem_d8", "dem_streams_d8_sub")) {
    terra::writeRaster(
      read_ihydro(input, lyr),
      file.path(temp_dir, paste0(lyr, ".tif")),
      overwrite = TRUE,
      datatype = "FLT4S"
    )
  }

  site_id_col <- read_site_id_col(db_fp)

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # ── Generate subbasin polygons ──────────────────────────────────────────
  whitebox::wbt_subbasins(
    d8_pntr = "dem_d8.tif",
    streams = "dem_streams_d8_sub.tif",
    output = "Subbasins.tif"
  )

  if (verbose) {
    message("Converting subbasins to polygons")
  }
  subb_rast <- terra::rast(file.path(temp_dir, "Subbasins.tif"))
  # subb_rast_cells <- subb_rast
  # subb_rast_cells[] <- 1
  # names(subb_rast_cells) <- "sum_nCells"
  #
  # subb_cells <- terra::zonal(subb_rast_cells, subb_rast, fun = "sum")

  subb <- subb_rast |>
    terra::as.polygons(dissolve = TRUE) |>
    sf::st_as_sf() |>
    dplyr::mutate(sbbsn_area = sf::st_area(geometry)) #|>
  #dplyr::left_join(subb_cells, by = "Subbasins")
  names(subb)[1] <- "link_id"

  sf::write_sf(subb, file.path(temp_dir, "Subbasins_poly.shp"))

  # ── Read existing link attributes ───────────────────────────────────────
  read_links <- function(con, input, site_id_col) {
    attrs <- read_ihydro(input, "stream_links_attr") |>
      dplyr::mutate(
        dplyr::across(
          c(link_id, tidyselect::any_of(site_id_col)),
          as.character
        ),
        dplyr::across(
          tidyselect::any_of(site_id_col),
          ~ dplyr::na_if(., "")
        )
      )
    geom <- sf::read_sf(input$outfile, layer = "stream_links") |>
      dplyr::mutate(dplyr::across(
        c(link_id, tidyselect::any_of(site_id_col)),
        as.character
      ))
    dplyr::left_join(geom, attrs, by = "link_id")
  }

  # ── Split subbasins at sampling points ──────────────────────────────────
  if (!is.null(points)) {
    if (verbose) {
      message("Splitting subbasins at sampling points")
    }
    stream_links <- read_links(con, input, site_id_col) |>
      dplyr::mutate(link_id = as.numeric(link_id))

    # Identify links with sample points
    links_with_pts <- stream_links |>
      dplyr::mutate(target_crs = list(sf::st_crs(subb))) |>
      dplyr::filter(
        floor(link_id) %in%
          floor(link_id[!is.na(!!rlang::sym(site_id_col))])
      )

    new_data <- links_with_pts |>
      dplyr::mutate(link_id_base = floor(link_id)) |>
      dplyr::select(
        link_id_base,
        link_id,
        tidyselect::any_of(site_id_col),
        target_crs
      ) |>
      tibble::as_tibble() |>
      dplyr::rename(point = geom) |>
      dplyr::group_by(link_id_base) |>
      tidyr::nest() |>
      dplyr::ungroup() |>
      dplyr::left_join(
        subb |>
          dplyr::select(-sbbsn_area) |> #, -sum_nCells
          tibble::as_tibble() |>
          dplyr::rename(subb_poly = geometry),
        by = c("link_id_base" = "link_id")
      )

    p <- progressr::progressor(steps = nrow(new_data))
    progressr::with_progress(enable = verbose, {
      new_data$new_subb <- furrr::future_pmap(
        list(
          data = new_data$data,
          link_id = new_data$link_id_base,
          subb_poly = new_data$subb_poly,
          temp_dir = rep(list(temp_dir), nrow(new_data)),
          p = rep(list(p), nrow(new_data))
        ),
        .options = furrr::furrr_options(
          globals = FALSE,
          seed = NULL,
          scheduling = 1L
        ),
        split_subbasin_worker
      )
    })

    subb2 <- new_data |>
      dplyr::select(new_subb) |>
      tidyr::unnest(cols = new_subb) |>
      sf::st_as_sf() |>
      dplyr::mutate(sbbsn_area = sf::st_area(geometry))

    subb <- subb |>
      dplyr::filter(!link_id %in% new_data$link_id_base) |>
      dplyr::bind_rows(subb2) |>
      dplyr::arrange(link_id)

    sf::write_sf(subb, file.path(temp_dir, "Subbasins_poly.shp"))

    # Re-read and join area
    all_links <- read_links(con, input, site_id_col)
    final_links <- all_links |>
      dplyr::left_join(
        subb |>
          tibble::as_tibble() |>
          dplyr::select(link_id, sbbsn_area) |>
          dplyr::mutate(link_id = as.character(link_id)),
        by = "link_id"
      )
  } else {
    stream_links <- read_links(con, input, site_id_col)
    final_links <- stream_links |>
      dplyr::left_join(
        subb |>
          tibble::as_tibble() |>
          dplyr::mutate(dplyr::across(
            c(link_id, tidyselect::any_of(site_id_col)),
            as.character
          )) |>
          dplyr::select(link_id, sbbsn_area),
        by = "link_id"
      )
  }

  # ── Write outputs ───────────────────────────────────────────────────────
  sf::write_sf(
    final_links |> dplyr::select(link_id),
    file.path(temp_dir, "stream_links.shp")
  )

  dplyr::copy_to(
    con,
    final_links |>
      tibble::as_tibble() |>
      dplyr::select(-tidyselect::any_of(c("geom", "geometry"))),
    "stream_links_attr",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = c("link_id", "trib_id"),
    analyze = TRUE,
    in_transaction = TRUE
  )

  for (lyr in c("Subbasins_poly", "stream_links")) {
    shp <- file.path(temp_dir, paste0(lyr, ".shp"))
    if (file.exists(shp)) {
      sf::write_sf(
        sf::read_sf(shp),
        db_fp,
        layer = lyr,
        append = TRUE,
        delete_layer = TRUE
      )
    }
  }

  # ── Return ──────────────────────────────────────────────────────────────
  output <- input[!names(input) %in% c("subbasins", "links")]
  if (return_products) {
    output <- c(list(subbasins = subb, links = final_links), output)
  }

  clean_temp_files(temp_dir)
  structure(output, class = "ihydro")
}


# ── Worker function for parallel subbasin splitting ───────────────────────────

#' @noRd
split_subbasin_worker <- carrier::crate(
  function(data, link_id, subb_poly, temp_dir, p) {
    `%>%` <- magrittr::`%>%`

    # Single-point case: no splitting needed

    if (nrow(data) == 1) {
      catch_poly <- data %>%
        dplyr::select(link_id) %>%
        dplyr::mutate(geom = sf::st_geometry(subb_poly)) %>%
        sf::st_as_sf(crs = data$target_crs[[1]]) %>%
        dplyr::mutate(sbbsn_area = sf::st_area(geom)) %>%
        dplyr::select(link_id, sbbsn_area, geom)
      p()
      return(catch_poly)
    }

    # Multi-point case: use unnest basins
    pnt_file <- file.path(temp_dir, paste0("Tempsite_", link_id, ".shp"))
    sf::write_sf(
      data %>% dplyr::select(link_id, point) %>% sf::st_as_sf(),
      pnt_file
    )

    cr <- terra::vect(subb_poly)
    d8_cropped <- terra::rast(file.path(temp_dir, "dem_d8.tif"))
    terra::crs(cr) <- terra::crs(d8_cropped)

    d8_cropped <- terra::crop(
      d8_cropped,
      y = cr,
      mask = TRUE,
      snap = "near",
      touches = FALSE,
      filename = file.path(
        temp_dir,
        paste0("d8_", link_id, ".tif")
      ),
      overwrite = TRUE
    )

    whitebox::wbt_unnest_basins(
      d8_pntr = paste0("d8_", link_id, ".tif"),
      pour_pts = paste0("Tempsite_", link_id, ".shp"),
      output = paste0("Catch_", link_id, ".tif")
    )

    catch_fls <- list.files(temp_dir, pattern = paste0("Catch_", link_id, "_"))
    catch_rast <- terra::rast(file.path(temp_dir, catch_fls))
    catch_rast[is.na(catch_rast)] <- 0
    catch_rast[catch_rast > 0] <- 1
    catch_rast <- terra::app(catch_rast, sum)
    catch_rast[catch_rast == 0] <- NA

    catch_poly <- catch_rast %>%
      terra::as.polygons() %>%
      sf::st_as_sf() %>%
      sf::st_join(data %>% sf::st_as_sf()) %>%
      dplyr::mutate(sbbsn_area = sf::st_area(geometry)) %>%
      dplyr::select(link_id, sbbsn_area, geometry)

    p()
    catch_poly
  }
)
