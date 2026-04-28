#' Prepare inverse distance weights
#'
#' Computes stream-targeted (iFLS, HAiFLS) and point-targeted (iFLO, HAiFLO)
#' inverse distance weights and stores them in a GeoPackage.
#' Stream-targeted weights are calculated for all stream links, while point-targeted
#' weights are calculated for a subset of target points (e.g., sample points) to save time and storage.
#' The function checks for existing weights in the output GeoPackage and only calculates missing ones.
#' It uses parallel processing to speed up the computation of point-targeted weights, which can be time-consuming.
#'
#' @details
#' The `target_o_type` argument controls the geometry type of the point-targeted weights, which can be
#' "point" (pour points), "segment_whole" (whole
#' stream reaches [ignoring artificial reach beaks introduced with the `points` argument in [generate_vectors()]]),
#' or "segment_point" (segments calculated separately for artificial reach beaks introduced
#' with the `points` argument in [generate_vectors()]).
#'
#' Point-targeted (iFLO, HAiFLO) inverse distance weights are calculated and stored in unnested format,
#' meaning that weights for different target points may be stored in the same raster layer, if they do
#' not overlap. This allows for more efficient storage and retrieval of weights for specific target points
#' without having to load many layers.
#'
#' When run in parallel, the function processes each LOI independently across multiple cores.
#' The `write_strategy` argument controls how processed rasters are persisted to the GeoPackage:
#' * `"sequential"` (default): all layers are processed in parallel first, then written to the
#'   GeoPackage in a single sequential pass. This approach is faster but uses more disk space temporarily.
#'   (note: this may trigger nested-future connection warnings in some environments, but these can be safely ignored).
#' * `"batched"`: layers are processed in chunks (one chunk per batch of `n_cores` layers).
#'   After each chunk finishes, its rasters are written to the GeoPackage before the next
#'   chunk begins. This reduces peak temporary disk usage at the cost of slightly more
#'   scheduling overhead.
#'
#' @param input An `ihydro` object.
#' @param output_filename Optional path for weight GeoPackage.
#' @param sample_points Character vector of site IDs (must exist in the `site_id_col` column provided in [generate_vectors()]), or `NULL`.
#' @param link_id Character vector of link IDs, or `NULL`.
#' @param target_o_type Character. Geometry type for point-targeted weights (see details).
#' @param weighting_scheme Character vector of weighting schemes.
#' @param inv_function Inverse distance function.
#' @param write_strategy Character. How processed weight rasters are written to the
#'   GeoPackage. `"sequential"` (default) waits for all parallel workers to
#'   finish, then writes every raster in one pass. `"batched"` processes unnest
#'   groups in chunks, writing to the GeoPackage between chunks to reduce peak
#'   temporary disk usage.
#' @param temp_dir Temporary file directory.
#' @param verbose Logical.
#'
#' @return An `ihydro` object pointing to the weight file.
#'
#' @export
#' @examples
#' \dontrun{
#' # Example usage
#'
#' future::plan(future::multisession(workers = 3))
#'
#' tdir <- tempdir()
#'
#' ex_dem <- ex_data("elev_ned_30m.tif")
#' ex_streams <- ex_data("streams.shp")
#' ex_points <- ex_data("sites_nc.shp")
#' output_gpkg <- file.path(tdir, "output.gpkg")
#'
#' # Process Hydrology
#' ihydro_obj <- ihydro::process_hydrology(
#'    dem = ex_dem,
#'    threshold = 1000,
#'    burn_streams = ex_streams,
#'    burn_depth = 5,
#'    burn_slope_dist = 250,
#'    burn_slope_depth = 5,
#'    min_length = 3,
#'    depression_corr = "breach",
#'    points = ex_points,
#'    site_id_col = "site_id",
#'    snap_distance = 150,
#'    break_on_noSnap = FALSE,
#'    pwise_dist = TRUE,
#'    pwise_all_links = TRUE,
#'    output_filename = output_gpkg,
#'    return_products = TRUE,
#'    verbose = TRUE
#' )
#'
#' # Prep Weights
#' weight_file <- file.path(tdir, "weights.gpkg")
#' ihydro_weights <- ihydro::prep_weights(
#'    input = ihydro_obj,
#'    output_filename = weight_file,
#'    sample_points = c("1", "25", "80"),
#'    link_id = c("100", "200", "800"),
#'    target_o_type = "segment_point",
#'    weighting_scheme = c("iFLS", "HAiFLS", "iFLO", "HAiFLO"),
#'    inv_function = function(x) (x * 0.001 + 1)^-1,
#'    verbose = TRUE
#' )
#'
#' lookup_tbl <- ihydro::read_ihydro(ihydro_weights,"unnest_catchment")
#' weight_id <- dplyr::filter(lookup_tbl, link_id == "800")$unn_group
#'
#' weight_rast <- ihydro::read_ihydro(ihydro_weights, paste0("iFLO_unn_group", weight_id))
#'
#' terra::plot(weight_rast)
#'
#' }
#'

prep_weights <- function(
  input,
  output_filename,
  weighting_scheme = c(
    "iFLS",
    "HAiFLS",
    "iFLO",
    "HAiFLO"
  ),
  inv_function = function(x) (x * 0.001 + 1)^-1,
  temp_dir = NULL,
  verbose = FALSE
) {
  check_ihydro(input)

  if (!grepl("\\.gpkg$", output_filename)) {
    cli::cli_abort("{.arg output_filename} must end in {.val .gpkg}.")
  }
  if (file.exists(output_filename)) {
    cli::cli_abort("{.arg output_filename} already exists.")
  }

  weighting_scheme <- match.arg(weighting_scheme, several.ok = TRUE)

  temp_dir <- ensure_temp_dir(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # ── Write rasters to temp ───────────────────────────────────────────────
  if (verbose) {
    message("Preparing inverse distance weights")
  }

  lyr_sel <- c("dem_accum_d8", "dem_final")
  target_S <- NULL
  target_O <- NULL

  for (lyr in lyr_sel) {
    targ_rast <- read_ihydro(input, lyr)
    if (lyr == "dem_final") {
      dem_crs <- terra::crs(targ_rast)

      if (grepl("iFLO", weighting_scheme)) {
        target_O <- dplyr::left_join(
          read_ihydro(input, "stream_links"),
          read_ihydro(input, "stream_links_attr"),
          by = "link_id"
        )
        target_O <- dplyr::filter(
          target_O,
          link_type == "Sink Node"
        )
        target_O <- target_O[, "link_id", drop = FALSE]
        target_O <- terra::rasterize(
          terra::vect(target_O),
          targ_rast,
          field = "",
          overwrite = TRUE,
          filename = file.path(temp_dir, paste0("target_O.tif")),
          wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=NONE"))
        )
      }
    }
    if (lyr == "dem_accum_d8") {
      flow_accum <- targ_rast + 1
    }

    terra::writeRaster(
      targ_rast,
      file.path(temp_dir, paste0(lyr, ".tif")),
      overwrite = TRUE,
      datatype = "FLT4S"
    )
  }

  # ── Generating Output ──────────────────────────────
  flow_accum |>
    terra::ext() |>
    terra::as.polygons(crs = terra::crs(flow_accum)) |>
    sf::st_as_sf() |>
    sf::write_sf(
      output_filename,
      layer = "DEM_Extent",
      append = TRUE,
      delete_layer = FALSE,
      delete_dsn = FALSE
    )

  # ── Stream-targeted weights (iFLS, HAiFLS) ──────────────────────────────
  if (grepl("iFLS", weighting_scheme)) {
    target_S <- terra::writeRaster(
      read_ihydro(input, "dem_streams_d8_sub"),
      file.path(temp_dir, "dem_streams_d8_sub.tif"),
      overwrite = TRUE,
      datatype = "FLT4S",
      gdal = c("COMPRESS=NONE")
    )

    whitebox::wbt_downslope_distance_to_stream(
      dem = file.path(temp_dir, "dem_final.tif"),
      streams = file.path(temp_dir, "dem_streams_d8_sub.tif"),
      output = file.path(temp_dir, "wbt_dist_to_stream.tif"),
      verbose_mode = FALSE
    )

    iFLS <- terra::rast(file.path(temp_dir, "wbt_dist_to_stream.tif"))
    iFLS_inv <- apply_inverse(iFLS, "iFLS")
    terra::crs(iFLS_inv) <- dem_crs
    names(iFLS_inv) <- "iFLS"

    if ("iFLS" %in% weighting_scheme) {
      write_raster_gpkg(iFLS_inv, output_filename)
    }

    if ("HAiFLS" %in% weighting_scheme_s) {
      HAiFLS_inv <- iFLS_inv * flow_accum
      HAiFLS_inv <- terra::mask(HAiFLS_inv, target_S, maskvalues = 1)
      names(HAiFLS_inv) <- "HAiFLS"
      write_raster_gpkg(HAiFLS_inv, output_filename)
    }
  }

  # ── Outlet-targeted weights (iFLO, HAiFLO) ──────────────────────────────
  if (grepl("iFLO", weighting_scheme)) {
    whitebox::wbt_downslope_distance_to_stream(
      dem = file.path(temp_dir, "dem_final.tif"),
      streams = file.path(temp_dir, "target_O.tif"),
      output = file.path(temp_dir, "wbt_dist_to_outlet.tif"),
      verbose_mode = FALSE
    )

    iFLO <- terra::rast(file.path(temp_dir, "wbt_dist_to_outlet.tif"))
    iFLO_inv <- apply_inverse(iFLO, "iFLO")
    terra::crs(iFLO_inv) <- dem_crs
    names(iFLO_inv) <- "iFLO"

    if ("iFLO" %in% weighting_scheme) {
      write_raster_gpkg(iFLO_inv, output_filename)
    }

    if ("HAiFLO" %in% weighting_scheme_s) {
      HAiFLO_inv <- iFLO_inv * flow_accum
      names(HAiFLO_inv) <- "HAiFLO"
      write_raster_gpkg(HAiFLO_inv, output_filename)
    }
  }

  sf::write_sf(
    data.frame(
      inv_function = deparse(substitute(inv_function)),
    ),
    output_filename,
    layer = "weight_meta",
    append = TRUE,
    delete_layer = FALSE,
    delete_dsn = FALSE
  )
}
