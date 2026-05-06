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
#' Point-targeted (iFLO, HAiFLO) inverse distance weights are calculated and stored in unnested format,
#' meaning that weights for different target points may be stored in the same raster layer, if they do
#' not overlap. This allows for more efficient storage and retrieval of weights for specific target points
#' without having to load many layers.
#'
#' @param input An `ihydro` object.
#' @param output_filename Optional path for weight GeoPackage.
#' @param sample_points Character vector of site IDs (must exist in the `site_id_col` column provided in [generate_vectors()]), or `NULL`.
#' @param link_id Character vector of link IDs, or `NULL`.
#' @param weighting_scheme Character vector of weighting schemes.
#' @param inv_function Inverse distance function.
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
  weighting_scheme = c("iFLS", "HAiFLS", "iFLO", "HAiFLO"),
  inv_function = function(x) (x * 0.001 + 1)^-1,
  return_products = FALSE,
  mem_fraction = 0.5,
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
    message("Gathering necessary layers")
  }

  lyr_sel <- c("dem_accum_d8", "dem_final")
  target_S <- NULL
  target_O <- NULL

  n_cores <- n_workers()

  for (lyr in lyr_sel) {
    targ_rast <- read_ihydro(input, lyr)
    if (lyr == "dem_final") {
      dem_crs <- terra::crs(targ_rast)

      if (any(grepl("iFLO", weighting_scheme))) {
        target_O <- dplyr::left_join(
          read_ihydro(input, "stream_links"),
          read_ihydro(input, "stream_links_attr"),
          by = "link_id"
        )

        target_O <- dplyr::filter(
          target_O,
          is.na(dstrib_id1)
        )

        dem_outline <- terra::clamp(targ_rast, 1, 1)
        dem_outline <- terra::as.polygons(dem_outline)
        dem_outline <- terra::as.lines(dem_outline)
        dem_outline <- terra::buffer(
          dem_outline,
          width = terra::res(targ_rast)[[1]] * 2
        )
        dem_outline <- sf::st_as_sf(dem_outline)

        target_O <- sf::st_filter(
          target_O,
          dem_outline
        )

        target_O <- target_O[, "link_id", drop = FALSE]
        target_O$link_id <- 1L

        sf::write_sf(
          target_O,
          file.path(temp_dir, paste0("target_O.shp"))
        )

        whitebox::wbt_vector_points_to_raster(
          input = file.path(temp_dir, paste0("target_O.shp")),
          output = file.path(temp_dir, paste0("target_O.tif")),
          base = file.path(temp_dir, paste0(lyr_sel[[1]], ".tif")),
          field = "link_id",
          verbose_mode = FALSE,
          compress_rasters = FALSE
        )
      }
    }
    if (lyr == "dem_accum_d8") {
      flow_accum <- terra::app(targ_rast, fun = function(x) x + 1)
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
  if (any(grepl("iFLS", weighting_scheme))) {
    if (verbose) {
      message("Preparing stream-targeted inverse distance weights")
    }

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
    iFLS_inv <- terra::app(iFLS, fun = inv_function)
    terra::crs(iFLS_inv) <- dem_crs
    names(iFLS_inv) <- "iFLS"

    if ("iFLS" %in% weighting_scheme) {
      write_raster_gpkg(iFLS_inv, output_filename)
    }

    if ("HAiFLS" %in% weighting_scheme) {
      HAiFLS_inv <- iFLS_inv * flow_accum
      HAiFLS_inv <- terra::mask(HAiFLS_inv, target_S, maskvalues = 1)
      names(HAiFLS_inv) <- "HAiFLS"
      write_raster_gpkg(HAiFLS_inv, output_filename)
    }
  }

  # ── Outlet-targeted weights (iFLO, HAiFLO) ──────────────────────────────
  if (any(grepl("iFLO", weighting_scheme))) {
    if (verbose) {
      message("Preparing outlet-targeted inverse distance weights")
    }

    dist_to_outlet_path <- file.path(temp_dir, "wbt_dist_to_outlet.tif")

    whitebox::wbt_downslope_distance_to_stream(
      dem = file.path(temp_dir, "dem_final.tif"),
      streams = file.path(temp_dir, "target_O.tif"),
      output = dist_to_outlet_path,
      verbose_mode = FALSE
    )

    unnest_catchment <- read_ihydro(input, "unnest_catchment")

    ungp <- unique(unnest_catchment$unn_group)

    # Ensure catchments exist
    catch <- get_catchment(
      input = input,
      link_id = unnest_catchment$link_id,
      temp_dir = temp_dir,
      verbose = verbose,
      return = FALSE
    )

    arg_list <- tibble::tibble(
      input_file = input$outfile,
      unn_group = ungp,
      temp_dir = temp_dir,
      dist_to_outlet_path = dist_to_outlet_path,
      inv_function = rep(list(inv_function), length(ungp)),
      weighting_scheme = rep(list(weighting_scheme), length(ungp)),
      mem_fraction = mem_fraction,
      n_cores = n_cores
    )

    arg_list <- split(arg_list, 1:nrow(arg_list))
    arg_list <- lapply(arg_list, as.list)

    progressr::with_progress(enable = verbose, {
      total_tasks <- length(arg_list)
      p <- progressr::progressor(steps = total_tasks)

      arg_list <- lapply(arg_list, function(args) {
        args$progressor <- p
        args
      })

      extract_execute <- lapply(arg_list, function(args) {
        future::futureCall(
          .O_target_worker,
          args = args,
          seed = NULL,
          globals = c("args"),
          packages = c("sf", "terra", "exactextractr", "dplyr", "ihydro")
        )
      })
      extract_value <- lapply(extract_execute, future::value)
    })

    for (i in extract_value) {
      rast_in <- terra::rast(i)
      terra::crs(rast_in) <- dem_crs
      for (ii in names(rast_in)) {
        write_raster_gpkg(rast_in[[ii]], output_filename)
      }
      rm(rast_in)
      file.remove(i)
    }
  }

  sf::write_sf(
    data.frame(
      inv_function = deparse(substitute(inv_function))
    ),
    output_filename,
    layer = "weight_meta",
    append = TRUE,
    delete_layer = FALSE,
    delete_dsn = FALSE
  )

  # ── Return ──────────────────────────────────────────────────────────────
  output <- list(outfile = output_filename)

  if (return_products) {
    rast_list <- lapply(
      setNames(weighting_scheme, weighting_scheme),
      function(x) read_ihydro(as_ihydro(output_filename), x)
    )
    output <- c(output, rast_list)
  }

  structure(output, class = "ihydro")
}

#' Worker for calculating iFLO and HAiFLO weights
#' @noRd
.O_target_worker <- carrier::crate(
  function(
    input_file,
    unn_group,
    temp_dir,
    dist_to_outlet_path,
    inv_function,
    weighting_scheme,
    mem_fraction = 0.5,
    n_cores = 1L,
    progressor = NULL
  ) {
    inv_function <- inv_function[[1]]
    weighting_scheme <- weighting_scheme[[1]]
    # Setup -------------------------------------------------------------------

    input <- ihydro::as_ihydro(input_file)
    temp_dir_sub <- file.path(temp_dir, basename(tempfile()))
    dir.create(temp_dir_sub)

    old_terra_opts <- ihydro:::set_terra_options(
      n_cores = n_cores,
      temp_dir = temp_dir_sub,
      verbose = FALSE
    )
    # NOTE: on.exit runs FIFO, so tmpFiles must be registered BEFORE
    # restore_terra_options to ensure it cleans temp_dir_sub (not the global dir).
    on.exit(
      suppressWarnings(
        terra::tmpFiles(
          current = TRUE,
          orphan = FALSE,
          old = FALSE,
          remove = TRUE
        )
      ),
      add = TRUE
    )
    on.exit(ihydro:::restore_terra_options(old_terra_opts), add = TRUE)
    on.exit(
      suppressWarnings(
        unlink(
          temp_dir_sub,
          recursive = TRUE,
          force = TRUE
        )
      ),
      add = T
    )

    on.exit(gc(verbose = FALSE), add = TRUE)

    max_cells_in_memory <- ihydro:::.max_cells_in_memory_helper(
      mem_fraction = mem_fraction,
      n_cores = n_cores
    )

    # Get subbasins -----------------------------------------------------------
    unnest_catchment <- ihydro::read_ihydro(input, "unnest_catchment")
    unnest_catchment <- unnest_catchment[
      unnest_catchment$unn_group == unn_group,
    ]

    subbasins <- sf::read_sf(
      input,
      query = ihydro:::build_sql_in(
        "Subbasins_poly",
        "link_id",
        unique(unnest_catchment$link_id)
      )
    )

    catchments <- ihydro::get_catchment(
      input = input,
      link_id = unique(unnest_catchment$link_id)
    )

    # Get Distance to outlet raster -------------------------------------------

    dist_to_outlet <- terra::rast(dist_to_outlet_path)

    # Modify Distance to outlet raster to catchments --------------------------
    catchments_rast <- exactextractr::rasterize_polygons(
      catchments,
      dist_to_outlet
    )

    min_vals <- exactextractr::exact_extract(
      dist_to_outlet,
      catchments,
      "min",
      full_colnames = T,
      force_df = T,
      max_cells_in_memory = max_cells_in_memory,
      progress = FALSE
    )

    min_vals <- data.frame(
      link_id = 1:nrow(catchments),
      min_vals = as.numeric(min_vals[[1]])
    )

    catchments_rast <- terra::classify(catchments_rast, as.matrix(min_vals))
    iDW_rasts <- dist_to_outlet - catchments_rast

    out_rasts <- list()
    # iFLO generation ---------------------------------------------------------

    iFLO_inv <- terra::app(iDW_rasts, fun = inv_function)
    if ("iFLO" %in% weighting_scheme) {
      names(iFLO_inv) <- paste0("iFLO_unn_group", unn_group)
      out_rasts[[names(iFLO_inv)]] <- iFLO_inv
    }

    # HAiFLO generation ---------------------------------------------------------

    if ("HAiFLO" %in% weighting_scheme) {
      flow_accum <- terra::rast(input_file, "dem_accum_d8")

      HAiFLO_inv <- iFLO_inv * flow_accum
      names(HAiFLO_inv) <- paste0("HAiFLO_unn_group", unn_group)
      out_rasts[[names(HAiFLO_inv)]] <- HAiFLO_inv
    }

    out_rasts <- terra::rast(out_rasts)
    outfile <- file.path(temp_dir, paste0("OT_", unn_group, ".tif"))
    terra::writeRaster(
      out_rasts,
      filename = outfile,
      overwrite = TRUE,
      datatype = "FLT4S",
      gdal = c("COMPRESS=NONE")
    )

    if (!is.null(progressor)) {
      progressor()
    }

    rm(
      dist_to_outlet,
      catchments_rast,
      iFLO_inv,
      flow_accum,
      HAiFLO_inv,
      iDW_rasts,
      out_rasts
    )

    return(outfile)
  }
)
