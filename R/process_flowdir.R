#' Process flow direction/accumulation and extract streams from DEM
#'
#' Produces D8 flow direction, flow accumulation, and stream extraction rasters
#' from a DEM, storing everything in a GeoPackage.
#'
#' @details
#' If `burn_streams` is specified:
#' - with `burn_slope_depth` and `burn_slope_dist`, then a gradual slope of
#' `burn_slope_depth` (map depth units) over `burn_slope_dist` meters is applied
#' to directing surface flows to the nearest `burn_streams`.
#' - with `burn_depth`, then the dem covered by the supplied stream lines is depressed
#' into the landscape by `burn_depth` meters.
#' - [whitebox::wbt_feature_preserving_smoothing] is applied
#' - if `stream_keep_thresh` is provided, then the overlap between `burn_streams`
#'  and `dem_streams_d8` is calculated and only stream segments >= `stream_keep_thresh`
#'  will be kept in the final product (note that this option should only be specified
#'  if the `burn_streams` is considered highly accurate, and is burned into the DEM
#'  deeply enough that extracted streams closely align)
#'
#' `depression_corr` can apply
#' filling [whitebox::wbt_fill_depressions()] or breaching
#' [whitebox::wbt_breach_depressions()] to enforce flow from internal pit cells.
#'
#' @param dem Character path, `SpatRaster`, or `PackedSpatRaster` for a DEM.
#' @param threshold Integer. Flow accumulation threshold for stream initiation.
#' @param burn_streams Optional stream layer to burn into the DEM.
#' @param burn_depth Numeric burn depth in map elevation units (default `NULL`).
#' @param stream_keep_thresh Numeric minimum proportion of overlap between `burn_streams`
#' and stream segments extracted from DEM to keep in final stream layer (see details)
#' @param burn_slope_dist Numeric burn slope width in meters (default `NULL`).
#' @param burn_slope_depth Numeric burn slope depth in map elevation units (default `NULL`).
#' @param min_length Numeric. Minimum tributary length in meters; shorter 1st-order
#'   tributaries are removed.
#' @param depression_corr `NULL`, `"fill"` [whitebox::wbt_fill_depressions()], or `"breach"` [whitebox::wbt_breach_depressions()].
#' @param output_filename Character path ending in `.gpkg`.
#' @param return_products Logical. Return raster products in addition to the
#'   file path?
#' @param temp_dir Character directory for temporary files.
#' @param compress Logical. Compress output rasters?
#' @param verbose Logical.
#'
#' @return An object of class `ihydro`.
#' The contained rasters are:
#' - `dem_raw`: the input DEM.
#' - `dem_final`: the final DEM after burning and depression correction.
#' - `dem_d8`: D8 flow direction pointer raster.
#' - `dem_accum_d8`: D8 flow accumulation in number of upstream cells.
#' - `dem_accum_d8_sca`: D8 flow accumulation in upstream contributing area (SCA).
#' - `dem_streams_d8`: D8-derived stream raster based on the specified flow accumulation threshold.
#' - `dem_streams_d8_sub`: Subset of dem_streams_d8 which meet minimum stream_keep_thresh.
#' If `return_products` is `TRUE`, the returned `ihydro` object will also include the above rasters as `SpatRaster` objects for immediate use in R. Otherwise, only file paths are included and temporary files are cleaned up.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
#' tdir <- tempdir()
#'
#' ex_dem <- ex_data("elev_ned_30m.tif")
#' ex_streams <- ex_data("streams.shp")
#' output_gpkg <- file.path(tdir, "output.gpkg")
#'
#' ihydro_obj <- ihydro::process_flowdir(
#'    dem = ex_dem,
#'    threshold = 1000,
#'    burn_streams = ex_streams,
#'    burn_depth = 5,
#'    burn_slope_dist = 250,
#'    burn_slope_depth = 5,
#'    min_length = 3,
#'    depression_corr = "breach",
#'    output_filename = output_gpkg,
#'    return_products = TRUE,
#'    verbose = TRUE
#' )
#'
#' ihydro::ihydro_layers(ihydro_obj, n = Inf)
#'
#' flow_accumulation <- ihydro::read_ihydro(ihydro_obj, "dem_accum_d8")
#' terra::plot(log10(flow_accumulation), main = "D8 Flow Accumulation")
#' }
#'
process_flowdir <- function(
    dem,
    threshold,
    burn_streams = NULL,
    stream_keep_thresh = NULL,
    chunk_size = 250,
    stream_keep_burnbuff = 25,
    stream_keep_dembuff = 1,
    burn_depth = NULL,
    burn_slope_dist = NULL,
    burn_slope_depth = NULL,
    min_length = NULL,
    depression_corr = c(NULL, "fill", "breach"),
    output_filename,
    return_products = FALSE,
    temp_dir = NULL,
    compress = FALSE,
    verbose = FALSE
) {
  # ── Validation ────────────────────────────────────────────────────────────
  if (!is.integer(threshold)) {
    if (as.integer(threshold) != threshold) {
      cli::cli_alert_info("{.arg threshold} was converted to an integer.")
    }
    threshold <- as.integer(threshold)
  }
  if (
    !is.null(burn_streams) && (is.null(burn_depth) && is.null(burn_slope_depth))
  ) {
    cli::cli_abort(
      "{.arg burn_depth} and/or {.arg burn_slope_depth} must be provided when {.arg burn_streams} is present"
    )
  }
  if (!is.null(burn_streams)) {
    if (!is.null(burn_depth) && !is.numeric(burn_depth)) {
      cli::cli_abort("{.arg burn_depth} must be numeric.")
    }
    if (!is.null(burn_slope_dist) && !is.numeric(burn_slope_dist)) {
      cli::cli_abort("{.arg burn_slope_dist} must be numeric.")
    }
    if (!is.null(burn_slope_depth) && !is.numeric(burn_slope_depth)) {
      cli::cli_abort("{.arg burn_slope_depth} must be numeric.")
    }
    if (!is.null(stream_keep_thresh) && !is.numeric(stream_keep_thresh)) {
      cli::cli_abort("{.arg stream_keep_thresh} must be numeric.")
    }
    if (!is.null(stream_keep_thresh) && (stream_keep_thresh <=0 |stream_keep_thresh > 1)) {
      cli::cli_abort("{.arg stream_keep_thresh} must be > 0 and <= 1")
    }
    if (!is.null(stream_keep_thresh) && (!is.numeric(stream_keep_burnbuff) | !is.numeric(stream_keep_dembuff))) {
      cli::cli_abort("{.arg stream_keep_burnbuff} and {.arg stream_keep_dembuff} must be numeric")
    }
  }

  if (!is.null(min_length) && !is.numeric(min_length)) {
    cli::cli_abort("{.arg min_length} must be numeric.")
  }
  stopifnot(
    is.logical(return_products),
    is.logical(compress),
    is.logical(verbose)
  )

  output_filename <- normalizePath(output_filename, mustWork = FALSE)
  if (!grepl("\\.gpkg$", output_filename)) {
    cli::cli_abort("{.arg output_filename} must end in {.val .gpkg}.")
  }
  if (file.exists(output_filename)) {
    cli::cli_abort("{.arg output_filename} already exists.")
  }

  temp_dir <- ensure_temp_dir(temp_dir)
  depression_corr <- match.arg(depression_corr)

  out_dir <- dirname(output_filename)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  if (normalizePath(out_dir) == temp_dir) {
    cli::cli_abort("{.arg output_filename} must not be inside {.arg temp_dir}.")
  }

  # ── Configure tools ──────────────────────────────────────────────────────
  whitebox::wbt_options(
    exe_path = whitebox::wbt_exe_path(),
    verbose = verbose > 2,
    wd = temp_dir,
    compress_rasters = compress
  )
  terra::terraOptions(verbose = verbose > 3, tempdir = temp_dir)

  gdal_arg <- if (compress) "COMPRESS=NONE" else NULL

  # ── Load and write DEM ──────────────────────────────────────────────────
  dem <- process_input(dem, input_name = "dem", working_dir = temp_dir)
  if (!inherits(dem, "SpatRaster")) {
    cli::cli_abort("{.arg dem} must resolve to a {.cls SpatRaster}.")
  }
  target_crs <- terra::crs(dem)
  names(dem) <- "Elevation"
  terra::writeRaster(
    dem,
    file.path(temp_dir, "dem_final.tif"),
    overwrite = TRUE,
    gdal = gdal_arg,
    datatype = "FLT4S"
  )
  terra::writeRaster(
    dem,
    file.path(temp_dir, "dem_raw.tif"),
    overwrite = TRUE,
    gdal = gdal_arg,
    datatype = "FLT4S"
  )

  # ── Optional stream burning ──────────────────────────────────────────────
  if (
    !is.null(burn_streams) && (!is.null(burn_depth) | is.null(burn_slope_depth))
  ) {
    if (verbose) {
      message("Burning Streams into DEM")
    }

    burn_dem(
      dem,
      burn_streams,
      burn_depth,
      burn_slope_dist,
      burn_slope_depth,
      target_crs,
      temp_dir,
      gdal_arg
    )
  }

  # ── Depression correction ────────────────────────────────────────────────
  if (!is.null(depression_corr)) {
    if (verbose) {
      message("Hydrologically conditioning DEM")
    }

    switch(
      depression_corr,
      fill = whitebox::wbt_fill_depressions(
        dem = "dem_final.tif",
        output = "dem_final.tif"
      ),
      breach = whitebox::wbt_breach_depressions(
        dem = "dem_final.tif",
        output = "dem_final.tif",
        fill_pits = TRUE
      )
    )
  }

  # ── D8 processing ───────────────────────────────────────────────────────
  if (verbose) {
    message("Generating D8 pointer")
  }
  whitebox::wbt_d8_pointer(dem = "dem_final.tif", output = "dem_d8.tif")

  if (verbose) {
    message("Generating D8 flow accumulation")
  }
  whitebox::wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
    output = "dem_accum_d8.tif",
    out_type = "cells",
    pntr = TRUE
  )
  whitebox::wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
    output = "dem_accum_d8_sca.tif",
    out_type = "sca",
    pntr = TRUE
  )

  # ── Stream extraction ───────────────────────────────────────────────────
  if (verbose) {
    message("Extracting streams")
  }
  whitebox::wbt_extract_streams(
    flow_accum = "dem_accum_d8.tif",
    output = "dem_streams_d8.tif",
    threshold = threshold
  )
  if (!is.null(min_length)) {
    if (verbose) {
      message("Trimming short streams")
    }
    fr <- file.rename(
      file.path(temp_dir, "dem_streams_d8.tif"),
      file.path(temp_dir, "dem_streams_d8_totrim.tif")
    )
    whitebox::wbt_remove_short_streams(
      d8_pntr = "dem_d8.tif",
      streams = "dem_streams_d8_totrim.tif",
      output = "dem_streams_d8.tif",
      min_length = min_length
    )
  }

  # ── Only keep DEM streams that align with  burn_streams ─────────────────

  if (is.null(stream_keep_thresh)) {
    dem_streams_d8_sub <- terra::rast(file.path(temp_dir, "dem_streams_d8.tif"))
    names(dem_streams_d8_sub) <- "dem_streams_d8_sub"

    terra::writeRaster(
      dem_streams_d8_sub,
      file.path(temp_dir, "dem_streams_d8_sub.tif"),
      overwrite = TRUE
    )

  } else {
    if (verbose) {
      message("Trimming DEM streams to match burn_streams")
    }

    res <- terra::res(dem)[[1]]

    ls <- sf::st_sfc(sf::st_linestring(rbind(c(0,0),c(0,1))))
    sf::st_crs(ls) <- sf::st_crs(target_crs)

    strm_lines_base <- process_input(
      burn_streams,
      align_to = ls
    ) |>
      sf::st_as_sf() |>
      sf::st_buffer(units::as_units(stream_keep_burnbuff,"m"),endCapStyle = "FLAT") |>
      dplyr::summarise()

    whitebox::wbt_stream_link_identifier(
      d8_pntr = "dem_d8.tif",
      streams = "dem_streams_d8.tif",
      output = "dem_streams_d8_link.tif"
    )

    whitebox::wbt_raster_streams_to_vector(
      "dem_streams_d8_link.tif",
      "dem_d8.tif",
      "dem_streams_d8_link.shp"
    )

    dem_lines_base0 <- sf::read_sf(
      file.path(temp_dir,"dem_streams_d8_link.shp")
    ) |>
      sf::st_set_crs(sf::st_crs(target_crs)) |>
      tibble::as_tibble() |>
      dplyr::mutate(
        geometry_split = sf::st_line_sample(geometry,density = units::as_units(chunk_size,"m")),
        geometry = sf::st_snap(geometry,geometry_split,sqrt(res))
      ) |>
      dplyr::mutate(
        geometry = furrr::future_map2(geometry,
                                      geometry_split,
                                      ~ {
                                        if (sf::st_is_empty(.y)) {
                                          parent_line <- sf::st_coordinates(.x)[, 1:2]
                                          parent_line <- sf::st_linestring(parent_line)
                                          parent_line <- tibble::tibble(geometry = sf::st_as_sfc(list(parent_line)))
                                          parent_line <- sf::st_as_sf(parent_line)
                                          sf::st_crs(parent_line) <- sf::st_crs(target_crs)
                                          return(parent_line)
                                        }
                                        .x <- sf::st_snap(.x,.y,sqrt(res))
                                        parent_line <- sf::st_coordinates(.x)[, 1:2]
                                        parent_line <- data.frame(matrix(parent_line,ncol=2,byrow=F))
                                        snap_point <- sf::st_coordinates(.y)[, 1:2]
                                        snap_point <- data.frame(matrix(snap_point,ncol=2,byrow=F))
                                        parent_index <- list()
                                        for (i in 1:nrow(snap_point)) {
                                          parent_index[[length(parent_index) + 1]] <- apply(parent_line, 1, function(x) all(x == snap_point[i,]))
                                        }
                                        parent_index <- sapply(parent_index, which)
                                        parent_index <- rep(parent_index,each = 2)
                                        if (head(parent_index,1) != 1) {
                                          parent_index <- c(1,parent_index)
                                        } else {
                                          parent_index <- c(1,parent_index[parent_index!=1])
                                        }
                                        if (tail(parent_index,1) != nrow(parent_line)) {
                                          parent_index <- c(parent_index,nrow(parent_line))
                                        }else {
                                          parent_index <- c(parent_index[parent_index!=nrow(parent_line)],nrow(parent_line))
                                        }
                                        parent_index <- split(parent_index,rep(1:floor(length(parent_index)/2),each = 2))

                                        parent_line <- lapply(
                                          parent_index,
                                          function(x){
                                            as.matrix(parent_line[x[[1]]:x[[2]], 1:2])
                                          }
                                        )
                                        parent_line <- lapply(parent_line, sf::st_linestring)
                                        parent_line <- tibble::tibble(geometry = sf::st_as_sfc(parent_line))
                                        parent_line <- sf::st_as_sf(parent_line)
                                        sf::st_crs(parent_line) <- sf::st_crs(target_crs)
                                        return(parent_line)
                                      })
      ) |>
      dplyr::select(-geometry_split) |>
      tidyr::unnest(geometry) |>
      dplyr::group_by(STRM_VAL) |>
      dplyr::mutate(sub_grp = 1:dplyr::n()) |>
      dplyr::group_by(STRM_VAL,sub_grp) |>
      sf::st_as_sf()

    dem_lines_base <- dem_lines_base0  |>
      dplyr::summarise() |>
      sf::st_buffer(units::as_units(stream_keep_dembuff,"m"),endCapStyle = "FLAT") |>
      dplyr::mutate(area = sf::st_area(geometry))

    intersect_pct <- sf::st_intersection(
      dem_lines_base,
      strm_lines_base
    )

    intersect_pct <- intersect_pct |>
      dplyr::mutate(intersect_area = as.numeric(sf::st_area(geometry)/area)) |>
      dplyr::select(STRM_VAL,sub_grp,intersect_area) |>
      sf::st_drop_geometry()

    intersect_thres <- intersect_pct |>
      dplyr::group_by(STRM_VAL) |>
      dplyr::mutate(keep = intersect_area > stream_keep_thresh) |>
      dplyr::mutate(max_keep = dplyr::case_when(
        all(!keep) ~ 0,
        T ~ suppressWarnings(min(which(keep)))
      )) |>
      dplyr::mutate(keep = sub_grp >= max_keep) |>
      dplyr::ungroup()

    dem_lines <- dem_lines_base0 |>
      dplyr::left_join(intersect_thres, by = c("STRM_VAL","sub_grp")) |>
      dplyr::filter(keep) |>
      dplyr::select(STRM_VAL) |>
      dplyr::mutate(STRM_VAL = 1)

    tstrm2 <- tempfile(fileext = ".shp")
    sf::write_sf(dem_lines, tstrm2)

    whitebox::wbt_rasterize_streams(
      streams = tstrm2,
      base = "dem_d8.tif",,
      output = "dem_streams_d8_s.tif"
    )

    dem_lines <- terra::rast(file.path(temp_dir, "dem_streams_d8_s.tif"))

    names(dem_lines) <- "dem_streams_d8_sub"
    terra::writeRaster(dem_lines,
                       file.path(temp_dir, "dem_streams_d8_sub.tif"),
                       overwrite = TRUE)
  }


  # ── Write outputs to GeoPackage ─────────────────────────────────────────
  raster_files <- c(
    "dem_raw.tif",
    "dem_final.tif",
    "dem_d8.tif",
    "dem_accum_d8.tif",
    "dem_accum_d8_sca.tif",
    "dem_streams_d8.tif",
    "dem_streams_d8_sub.tif"
  )

  # Write DEM extent polygon
  dem_extent <- terra::rast(file.path(temp_dir, raster_files[1]))
  dem_extent |>
    terra::ext() |>
    terra::as.polygons(crs = terra::crs(dem_extent)) |>
    sf::st_as_sf() |>
    sf::write_sf(
      output_filename,
      layer = "DEM_Extent",
      append = TRUE,
      delete_layer = FALSE,
      delete_dsn = FALSE
    )

  if (verbose) {
    message("Writing rasters to GeoPackage")
  }
  for (f in raster_files) {
    r <- terra::rast(file.path(temp_dir, f))
    write_raster_gpkg(r, output_filename)
  }

  # ── Return ──────────────────────────────────────────────────────────────
  output <- list(outfile = output_filename)

  if (return_products) {
    rast_list <- lapply(
      stats::setNames(raster_files, raster_files),
      function(f) terra::wrap(terra::rast(file.path(temp_dir, f)))
    )
    output <- c(rast_list, output)
  } else {
    clean_temp_files(temp_dir)
  }

  structure(output, class = "ihydro")
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Burn streams into DEM with graduated buffering
#' @noRd
burn_dem <- function(
    dem,
    burn_streams,
    burn_depth,
    burn_slope_dist,
    burn_slope_depth,
    target_crs,
    temp_dir,
    gdal_arg
) {
  bbox <- terra::as.polygons(terra::ext(dem), crs = target_crs)

  resol <- terra::cellSize(dem, unit = "m")
  resol <- unlist(terra::global(resol, "mean"), use.names = F)
  resol <- sqrt(resol)

  if (is.null(burn_slope_dist)) {
    burn_slope_dist <- Inf
  }
  if (is.null(burn_slope_depth)) {
    burn_slope_depth <- 0
  } else {
    burn_slope_depth <- abs(burn_slope_depth)
  }
  if (is.null(burn_depth)) {
    burn_depth <- 0
  } else {
    burn_depth <- abs(burn_depth)
  }

  burn_streams <- process_input(
    burn_streams,
    align_to = terra::as.lines(
      terra::vect(
        "POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",
        crs = target_crs
      )
    ),
    clip_region = bbox,
    input_name = "burn_streams",
    working_dir = temp_dir
  )

  rast_streams <- process_input(
    burn_streams[, 1],
    align_to = dem,
    input_name = "burn_streams",
    working_dir = temp_dir
  )

  rast_streams <- terra::mask(dem, burn_streams, updatevalue = NA)
  rast_streams[!is.na(rast_streams)] <- 1

  dem_final <- dem
  if (burn_slope_depth != 0) {
    rast_streams[is.na(rast_streams)] <- 0

    terra::writeRaster(
      rast_streams,
      file.path(temp_dir, "stream_temp.tif"),
      overwrite = TRUE,
      datatype = "FLT4S"
    )

    whitebox::wbt_euclidean_distance(
      input = file.path(temp_dir, "stream_temp.tif"),
      output = file.path(temp_dir, "stream_ds_dist.tif"),
      verbose_mode = FALSE
    )

    rast_streams[rast_streams != 1] <- NA

    rast_streams <- terra::buffer(rast_streams, width = burn_slope_dist)
    rast_streams <- as.numeric(rast_streams)
    rast_streams[rast_streams == 0] <- NA

    ds_dist <- terra::rast(file.path(temp_dir, "stream_ds_dist.tif"))
    ds_dist <- ds_dist * rast_streams

    max_dist <- unlist(
      terra::global(ds_dist, "max", na.rm = T),
      use.names = FALSE
    )
    if (max_dist == 0) {
      cli::cli_abort(
        "{.arg burn_slope_dist} is too small for the supplied dem."
      )
    }
    ds_dist <- (ds_dist - max_dist) / max_dist
    ds_dist[is.na(ds_dist)] <- 0
    ds_dist <- ds_dist * burn_slope_depth

    dem_final <- dem_final + ds_dist
  }

  if (burn_depth != 0) {
    sr1 <- terra::mask(
      dem_final,
      burn_streams %>%
        sf::st_as_sf() %>%
        sf::st_buffer(resol * 3) %>%
        terra::vect()
    )
    dem_final[!is.na(sr1)] <- sr1 - ceiling(burn_depth / 3)
    sr1 <- terra::mask(
      dem_final,
      burn_streams %>%
        sf::st_as_sf() %>%
        sf::st_buffer(resol * 2) %>%
        terra::vect()
    )
    dem_final[!is.na(sr1)] <- sr1 - ceiling(burn_depth / 3)
    sr1 <- terra::mask(
      dem_final,
      burn_streams %>%
        sf::st_as_sf() %>%
        sf::st_buffer(resol * 1) %>%
        terra::vect()
    )
    dem_final[!is.na(sr1)] <- sr1 - ceiling(burn_depth / 3)

    # rast_streams[rast_streams != 1] <- NA
    #
    # rast_streams <- terra::buffer(
    #   rast_streams,
    #   width = terra::res(rast_streams)[[1]]
    # )
    # rast_streams <- as.numeric(rast_streams)
    #
    # rast_streams[rast_streams == 1] <- (-burn_depth)
    #
    # dem_final <- dem_final + rast_streams
  }

  terra::writeRaster(
    dem_final,
    file.path(temp_dir, "dem_final.tif"),
    overwrite = TRUE,
    gdal = gdal_arg,
    datatype = "FLT4S"
  )

  whitebox::wbt_feature_preserving_smoothing(
    dem = "dem_final.tif",
    output = "dem_final.tif"
  )

  return(NULL)
}

#' Remove temporary files (non-destructively)
#' @noRd
clean_temp_files <- function(temp_dir) {
  files <- list.files(temp_dir, full.names = TRUE, recursive = TRUE)
  suppressWarnings(file.remove(files))
}
