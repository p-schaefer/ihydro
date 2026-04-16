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
#'
#' `depression_corr` can apply
#' filling [whitebox::wbt_fill_depressions()] or breaching
#' [whitebox::wbt_breach_depressions()] to enforce flow from internal pit cells.
#'
#' @param dem Character path, `SpatRaster`, or `PackedSpatRaster` for a DEM.
#' @param threshold Integer. Flow accumulation threshold for stream initiation.
#' @param burn_streams Optional stream layer to burn into the DEM.
#' @param burn_depth Numeric burn depth in map elevation units (default `NULL`).
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

  # ── Write outputs to GeoPackage ─────────────────────────────────────────
  raster_files <- c(
    "dem_raw.tif",
    "dem_final.tif",
    "dem_d8.tif",
    "dem_accum_d8.tif",
    "dem_accum_d8_sca.tif",
    "dem_streams_d8.tif"
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
