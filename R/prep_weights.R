#' #' Prepare inverse distance weights
#' #'
#' #' Computes stream-targeted (iFLS, HAiFLS) and point-targeted (iFLO, HAiFLO)
#' #' inverse distance weights and stores them in a GeoPackage.
#' #' Stream-targeted weights are calculated for all stream links, while point-targeted
#' #' weights are calculated for a subset of target points (e.g., sample points) to save time and storage.
#' #' The function checks for existing weights in the output GeoPackage and only calculates missing ones.
#' #' It uses parallel processing to speed up the computation of point-targeted weights, which can be time-consuming.
#' #'
#' #' @details
#' #' The `target_o_type` argument controls the geometry type of the point-targeted weights, which can be
#' #' "point" (pour points), "segment_whole" (whole
#' #' stream reaches [ignoring artificial reach beaks introduced with the `points` argument in [generate_vectors()]]),
#' #' or "segment_point" (segments calculated separately for artificial reach beaks introduced
#' #' with the `points` argument in [generate_vectors()]).
#' #'
#' #' Point-targeted (iFLO, HAiFLO) inverse distance weights are calculated and stored in unnested format,
#' #' meaning that weights for different target points may be stored in the same raster layer, if they do
#' #' not overlap. This allows for more efficient storage and retrieval of weights for specific target points
#' #' without having to load many layers.
#' #'
#' #' When run in parallel, the function processes each LOI independently across multiple cores.
#' #' The `write_strategy` argument controls how processed rasters are persisted to the GeoPackage:
#' #' * `"sequential"` (default): all layers are processed in parallel first, then written to the
#' #'   GeoPackage in a single sequential pass. This approach is faster but uses more disk space temporarily.
#' #'   (note: this may trigger nested-future connection warnings in some environments, but these can be safely ignored).
#' #' * `"batched"`: layers are processed in chunks (one chunk per batch of `n_cores` layers).
#' #'   After each chunk finishes, its rasters are written to the GeoPackage before the next
#' #'   chunk begins. This reduces peak temporary disk usage at the cost of slightly more
#' #'   scheduling overhead.
#' #'
#' #' @param input An `ihydro` object.
#' #' @param output_filename Optional path for weight GeoPackage.
#' #' @param sample_points Character vector of site IDs (must exist in the `site_id_col` column provided in [generate_vectors()]), or `NULL`.
#' #' @param link_id Character vector of link IDs, or `NULL`.
#' #' @param target_o_type Character. Geometry type for point-targeted weights (see details).
#' #' @param weighting_scheme Character vector of weighting schemes.
#' #' @param inv_function Inverse distance function.
#' #' @param write_strategy Character. How processed weight rasters are written to the
#' #'   GeoPackage. `"sequential"` (default) waits for all parallel workers to
#' #'   finish, then writes every raster in one pass. `"batched"` processes unnest
#' #'   groups in chunks, writing to the GeoPackage between chunks to reduce peak
#' #'   temporary disk usage.
#' #' @param temp_dir Temporary file directory.
#' #' @param verbose Logical.
#' #'
#' #' @return An `ihydro` object pointing to the weight file.
#' #'
#' #' @export
#' #' @examples
#' #' \dontrun{
#' #' # Example usage
#' #'
#' #' future::plan(future::multisession(workers = 3))
#' #'
#' #' tdir <- tempdir()
#' #'
#' #' ex_dem <- ex_data("elev_ned_30m.tif")
#' #' ex_streams <- ex_data("streams.shp")
#' #' ex_points <- ex_data("sites_nc.shp")
#' #' output_gpkg <- file.path(tdir, "output.gpkg")
#' #'
#' #' # Process Hydrology
#' #' ihydro_obj <- ihydro::process_hydrology(
#' #'    dem = ex_dem,
#' #'    threshold = 1000,
#' #'    burn_streams = ex_streams,
#' #'    burn_depth = 5,
#' #'    burn_slope_dist = 250,
#' #'    burn_slope_depth = 5,
#' #'    min_length = 3,
#' #'    depression_corr = "breach",
#' #'    points = ex_points,
#' #'    site_id_col = "site_id",
#' #'    snap_distance = 150,
#' #'    break_on_noSnap = FALSE,
#' #'    pwise_dist = TRUE,
#' #'    pwise_all_links = TRUE,
#' #'    output_filename = output_gpkg,
#' #'    return_products = TRUE,
#' #'    verbose = TRUE
#' #' )
#' #'
#' #' # Prep Weights
#' #' weight_file <- file.path(tdir, "weights.gpkg")
#' #' ihydro_weights <- ihydro::prep_weights(
#' #'    input = ihydro_obj,
#' #'    output_filename = weight_file,
#' #'    sample_points = c("1", "25", "80"),
#' #'    link_id = c("100", "200", "800"),
#' #'    target_o_type = "segment_point",
#' #'    weighting_scheme = c("iFLS", "HAiFLS", "iFLO", "HAiFLO"),
#' #'    inv_function = function(x) (x * 0.001 + 1)^-1,
#' #'    verbose = TRUE
#' #' )
#' #'
#' #' lookup_tbl <- ihydro::read_ihydro(ihydro_weights,"unnest_catchment")
#' #' weight_id <- dplyr::filter(lookup_tbl, link_id == "800")$unn_group
#' #'
#' #' weight_rast <- ihydro::read_ihydro(ihydro_weights, paste0("iFLO_unn_group", weight_id))
#' #'
#' #' terra::plot(weight_rast)
#' #'
#' #' }
#' #'
#'
#' prep_weights_rem <- function(
#'   input,
#'   output_filename = NULL,
#'   sample_points = NULL,
#'   link_id = NULL,
#'   target_o_type = c(
#'     "point",
#'     "segment_point",
#'     "segment_whole"
#'   ),
#'   weighting_scheme = c(
#'     "iFLS",
#'     "HAiFLS",
#'     "iFLO",
#'     "HAiFLO"
#'   ),
#'   inv_function = function(x) (x * 0.001 + 1)^-1,
#'   write_strategy = c("sequential", "batched"),
#'   temp_dir = NULL,
#'   verbose = FALSE
#' ) {
#'   check_ihydro(input)
#'   if (inherits(output_filename, "ihydro")) {
#'     output_filename <- output_filename$outfile
#'   }
#'
#'   write_strategy <- match.arg(write_strategy)
#'
#'   n_cores <- n_workers()
#'   # max_cores_opt <- getOption("parallelly.maxWorkers.localhost")
#'   # on.exit(options(parallelly.maxWorkers.localhost = max_cores_opt), add = TRUE)
#'   # options(parallelly.maxWorkers.localhost = n_cores)
#'
#'   target_o_type <- match.arg(target_o_type)
#'   weighting_scheme <- match.arg(weighting_scheme, several.ok = TRUE)
#'   weighting_scheme <- weighting_scheme[weighting_scheme != "lumped"]
#'
#'   temp_dir <- ensure_temp_dir(temp_dir)
#'   whitebox::wbt_options(
#'     exe_path = whitebox::wbt_exe_path(),
#'     verbose = verbose > 2,
#'     wd = temp_dir
#'   )
#'   terra::terraOptions(verbose = verbose > 3, tempdir = temp_dir)
#'
#'   db_fp <- input$outfile
#'
#'   # ── Resolve output path ─────────────────────────────────────────────────
#'   if (is.null(output_filename) || output_filename == db_fp) {
#'     output_filename <- db_fp
#'     if (!file.exists(output_filename)) {
#'       dir.create(
#'         dirname(output_filename),
#'         showWarnings = FALSE,
#'         recursive = TRUE
#'       )
#'     }
#'   } else {
#'     if (!grepl("\\.gpkg$", output_filename)) {
#'       cli::cli_abort("{.arg output_filename} must end in {.val .gpkg}.")
#'     }
#'     if (!file.exists(output_filename)) {
#'       dir.create(
#'         dirname(output_filename),
#'         showWarnings = FALSE,
#'         recursive = TRUE
#'       )
#'       dem <- read_ihydro(input, "dem_final")
#'       dem |>
#'         terra::ext() |>
#'         terra::as.polygons(crs = terra::crs(dem)) |>
#'         sf::st_as_sf() |>
#'         sf::write_sf(
#'           output_filename,
#'           layer = "DEM_Extent",
#'           append = TRUE,
#'           delete_layer = FALSE,
#'           delete_dsn = FALSE
#'         )
#'     }
#'   }
#'
#'   output_filename <- as.ihydro(output_filename)
#'   lyrs <- ihydro_layers(output_filename)
#'
#'   # ── Check what already exists ───────────────────────────────────────────
#'   weighting_scheme_s <- weighting_scheme[grepl("FLS", weighting_scheme)]
#'   weighting_scheme_s <- weighting_scheme_s[
#'     !weighting_scheme_s %in%
#'       lyrs$layer_name
#'   ]
#'   weighting_scheme_o <- weighting_scheme[!grepl("lumped|FLS", weighting_scheme)]
#'
#'   target_ids <- target_id_fun(
#'     db_fp = db_fp,
#'     sample_points = sample_points,
#'     link_id = link_id
#'   )
#'
#'   # Check available point-targeted weights
#'   avail_weights <- resolve_available_weights(
#'     lyrs,
#'     output_filename,
#'     target_ids,
#'     weighting_scheme_o,
#'     target_o_type = target_o_type
#'   )
#'
#'   if (length(weighting_scheme_s) == 0 && all(avail_weights$layer_name_exists)) {
#'     return(
#'       if (!inherits(output_filename, "ihydro")) {
#'         as.ihydro(output_filename$outfile)
#'       } else {
#'         output_filename
#'       }
#'     )
#'   }
#'
#'   if (verbose) {
#'     message("Preparing inverse distance weights")
#'   }
#'
#'   if (all(weighting_scheme_s %in% lyrs$layer_name)) {
#'     if (verbose) {
#'       message("Stream-targeted weights already exist; won't recalculate.")
#'     }
#'   }
#'
#'   # ── Write rasters to temp ───────────────────────────────────────────────
#'   temp_dir_sub <- file.path(temp_dir, basename(tempfile()))
#'   dir.create(temp_dir_sub)
#'
#'   for (lyr in c("dem_streams_d8_sub", "dem_final", "dem_accum_d8", "dem_d8")) {
#'     terra::writeRaster(
#'       read_ihydro(input, lyr),
#'       file.path(temp_dir_sub, paste0(lyr, ".tif")),
#'       overwrite = TRUE,
#'       datatype = "FLT4S"
#'     )
#'   }
#'
#'   # ── Stream-targeted weights (iFLS, HAiFLS) ──────────────────────────────
#'   if (length(weighting_scheme_s) > 0) {
#'     if (verbose) {
#'       message("Generating stream-targeted weights")
#'     }
#'     temp_hw <- file.path(temp_dir_sub, basename(tempfile()))
#'     dir.create(temp_hw)
#'
#'     suppressMessages(
#'       hw_streams <- hydroweight(
#'         hydroweight_dir = temp_hw,
#'         target_O = NULL,
#'         target_S = file.path(temp_dir_sub, "dem_streams_d8_sub.tif"),
#'         target_uid = "ALL",
#'         OS_combine = FALSE,
#'         dem = file.path(temp_dir_sub, "dem_final.tif"),
#'         flow_accum = file.path(temp_dir_sub, "dem_accum_d8.tif"),
#'         weighting_scheme = weighting_scheme_s,
#'         inv_function = inv_function,
#'         clean_tempfiles = FALSE,
#'         return_products = TRUE,
#'         wrap_return_products = FALSE,
#'         save_output = FALSE
#'       )
#'     )
#'
#'     for (r in hw_streams) {
#'       write_raster_gpkg(r, output_filename$outfile)
#'
#'       # lyr_name <- names(r)
#'       # r2 <- r^2
#'       # lyr_name2 <- lyr_name
#'       # lyr_name2 <- gsub("iFLS", "iFLSSQ", lyr_name2)
#'       # names(r2) <- gsub("iFLS", "iFLSSQ", names(r2))
#'       # terra::varnames(r2) <- lyr_name2
#'       # write_raster_gpkg(r2, output_filename$outfile)
#'     }
#'     rm(hw_streams)
#'     unlink(temp_hw, recursive = TRUE, force = TRUE)
#'   }
#'
#'   # ── Point-targeted weights (iFLO, HAiFLO) ───────────────────────────────
#'   if (length(weighting_scheme_o) > 0) {
#'     compute_point_weights(
#'       db_fp = db_fp,
#'       output_filename = output_filename,
#'       lyrs = lyrs,
#'       target_ids = target_ids,
#'       weighting_scheme_o = weighting_scheme_o,
#'       inv_function = inv_function,
#'       target_o_type = target_o_type,
#'       temp_dir_sub = temp_dir_sub,
#'       n_cores = n_cores,
#'       write_strategy = write_strategy,
#'       verbose = verbose
#'     )
#'   }
#'
#'   unlink(temp_dir, recursive = TRUE, force = TRUE)
#'
#'   output <- input
#'   if (output$outfile != output_filename$outfile) {
#'     output <- list(outfile = output_filename$outfile)
#'   }
#'   if (!inherits(output, "ihydro")) {
#'     output <- as.ihydro(output$outfile)
#'   }
#'   output
#' }
#'
#'
#' # ══════════════════════════════════════════════════════════════════════════════
#' # Internal helpers
#' # ══════════════════════════════════════════════════════════════════════════════
#'
#' #' Resolve which point weights already exist
#' #' @noRd
#' resolve_available_weights <- function(
#'   lyrs,
#'   output_filename,
#'   target_ids,
#'   weighting_scheme_o,
#'   target_o_type = NULL
#' ) {
#'   if ("unnest_catchment" %in% lyrs$layer_name) {
#'     existing_meta <- read_ihydro(output_filename, "unnest_catchment") |>
#'       dplyr::mutate(link_id = as.character(link_id))
#'
#'     # Filter by target_o_type if the column exists and a type was provided
#'     if (
#'       !is.null(target_o_type) && "target_o_type" %in% colnames(existing_meta)
#'     ) {
#'       existing_meta <- dplyr::filter(
#'         existing_meta,
#'         target_o_type == !!target_o_type
#'       )
#'     }
#'
#'     # left_join: keep all target_ids, match existing metadata where available.
#'     # full_join was introducing rows from unnest_catchment for link_ids NOT in
#'     # target_ids, and — more critically — target_ids without a match got
#'     # unn_group = NA, producing "iFLO_unn_groupNA" layer names.
#'     avail <- target_ids |>
#'       dplyr::select(link_id) |>
#'       dplyr::left_join(existing_meta, by = "link_id")
#'
#'     purrr::map_dfr(
#'       stats::setNames(weighting_scheme_o, weighting_scheme_o),
#'       function(ws) {
#'         avail |>
#'           dplyr::mutate(
#'             weight = ws,
#'             layer_name = dplyr::if_else(
#'               is.na(unn_group),
#'               NA_character_,
#'               paste0(ws, "_unn_group", unn_group)
#'             ),
#'             layer_name_exists = !is.na(layer_name) &
#'               layer_name %in% lyrs$layer_name
#'           )
#'       }
#'     )
#'   } else {
#'     purrr::map_dfr(
#'       stats::setNames(weighting_scheme_o, weighting_scheme_o),
#'       function(ws) {
#'         tibble::tibble(link_id = target_ids$link_id) |>
#'           dplyr::mutate(
#'             weight = ws,
#'             layer_name = NA_character_,
#'             layer_name_exists = FALSE
#'           )
#'       }
#'     )
#'   }
#' }
#'
#' #' Generate a short random ID
#' #' @noRd
#' rand_id <- function(n = 1, d = 12) {
#'   paste0(replicate(d, sample(c(LETTERS, letters, 0:9), n, TRUE)), collapse = "")
#' }
#'
#' #' Compute point-targeted inverse distance weights
#' #' @noRd
#' compute_point_weights <- function(
#'   db_fp,
#'   output_filename,
#'   lyrs,
#'   target_ids,
#'   weighting_scheme_o,
#'   inv_function,
#'   target_o_type,
#'   temp_dir_sub,
#'   n_cores,
#'   write_strategy = "sequential",
#'   verbose
#' ) {
#'   temp_hw <- file.path(temp_dir_sub, basename(tempfile()))
#'   dir.create(temp_hw)
#'
#'   # Check for existing metadata
#'   unnest_catchment <- NULL
#'   if ("unnest_catchment" %in% lyrs$layer_name) {
#'     unnest_catchment <- read_ihydro(output_filename, "unnest_catchment")
#'     # Filter by target_o_type if the column exists
#'     if ("target_o_type" %in% colnames(unnest_catchment)) {
#'       unnest_catchment <- dplyr::filter(
#'         unnest_catchment,
#'         target_o_type == !!target_o_type
#'       )
#'     }
#'     unnest_catchment <- unnest_catchment |>
#'       dplyr::filter(
#'         paste0(weighting_scheme_o, "_unn_group", unn_group) %in%
#'           lyrs$layer_name
#'       )
#'   }
#'
#'   target_o <- target_o_fun(
#'     db_fp = db_fp,
#'     target_IDs = target_ids,
#'     target_type = target_o_type
#'   ) |>
#'     dplyr::filter(!link_id %in% unnest_catchment$link_id)
#'
#'   all_points <- target_o_fun(
#'     db_fp = db_fp,
#'     target_IDs = target_ids,
#'     target_type = "point"
#'   ) |>
#'     dplyr::filter(!link_id %in% unnest_catchment$link_id)
#'
#'   if (nrow(target_o) == 0) {
#'     return(invisible(NULL))
#'   }
#'
#'   sf::write_sf(
#'     all_points |> dplyr::select(link_id),
#'     file.path(temp_hw, "pour_points.shp"),
#'     overwrite = TRUE
#'   )
#'
#'   if (verbose) {
#'     message("Generating point-targeted weights (can be slow)")
#'   }
#'
#'   # ── Unnest basins to find non-overlapping groups ────────────────────────
#'   future_unnest <- future::future({
#'     whitebox::wbt_unnest_basins(
#'       d8_pntr = file.path(temp_dir_sub, "dem_d8.tif"),
#'       pour_pts = file.path(temp_hw, "pour_points.shp"),
#'       output = file.path(temp_hw, "unnest.tif")
#'     )
#'   })
#'
#'   rast_out <- collect_unnest_rasters(future_unnest, temp_hw)
#'
#'   link_id_split <- purrr::map(rast_out, ~ all_points$link_id[.[[1]]])
#'
#'   # ── Compute weights in parallel ─────────────────────────────────────────
#'   # if (n_cores > 1) {
#'   #   oplan <- future::plan(future::multisession, workers = n_cores)
#'   #   on.exit(future::plan(oplan), add = TRUE)
#'   # }
#'
#'   target_o_sub <- purrr::map(link_id_split, function(ids) {
#'     dplyr::filter(target_o, link_id %in% ids) |>
#'       dplyr::select(link_id) |>
#'       dplyr::mutate(unn_group = rand_id(1))
#'   })
#'
#'   # Helper to write finished .tif files from temp_hw into the gpkg
#'   write_weight_tifs <- function(pattern, p = NULL) {
#'     fl <- list.files(temp_hw, full.names = TRUE)
#'     fl <- fl[grepl(pattern, fl)]
#'     for (f in fl) {
#'       r <- try(terra::rast(f), silent = TRUE)
#'       if (inherits(r, "try-error")) {
#'         next
#'       }
#'       lyr_name <- gsub("\\.tif$", "", basename(f))
#'       res <- try(
#'         write_raster_gpkg(r, output_filename$outfile, lyr_name),
#'         silent = TRUE
#'       )
#'       if (inherits(res, "try-error")) {
#'         msg <- conditionMessage(attr(res, "condition"))
#'         if (!msg %in% c("stoi", "stol")) stop(msg)
#'       }
#'
#'       # r2 <- r^2
#'       # lyr_name2 <- lyr_name
#'       # lyr_name2 <- gsub("iFLO", "iFLOSQ", lyr_name2)
#'       # names(r2) <- gsub("iFLO", "iFLOSQ", names(r2))
#'       # terra::varnames(r2) <- lyr_name2
#'       # res <- try(
#'       #   write_raster_gpkg(r2, output_filename$outfile, lyr_name2),
#'       #   silent = TRUE
#'       # )
#'       # if (inherits(res, "try-error")) {
#'       #   msg <- conditionMessage(attr(res, "condition"))
#'       #   if (!msg %in% c("stoi", "stol")) stop(msg)
#'       # }
#'
#'       if (!is.null(p)) {
#'         p()
#'       }
#'       try(file.remove(f), silent = TRUE)
#'     }
#'   }
#'
#'   ws_pattern <- paste0(weighting_scheme_o, collapse = "|")
#'
#'   if (write_strategy == "batched") {
#'     # ── Batched: process in chunks, drain to gpkg between chunks ────────
#'     chunks <- split(
#'       seq_along(target_o_sub),
#'       ceiling(seq_along(target_o_sub) / max(n_cores, 1L))
#'     )
#'
#'     progressr::with_progress(enable = verbose, {
#'       p <- progressr::progressor(
#'         steps = length(target_o_sub) * length(weighting_scheme_o)
#'       )
#'
#'       for (chunk in chunks) {
#'         chunk_lst <- target_o_sub[chunk]
#'
#'         furrr::future_pmap(
#'           list(
#'             x = list(chunk_lst),
#'             output_filename = list(output_filename),
#'             temp_dir_sub = list(temp_dir_sub),
#'             temp_dir_sub2 = list(temp_hw),
#'             weighting_scheme_o = list(weighting_scheme_o),
#'             inv_function = list(inv_function),
#'             verbose = list(verbose)
#'           ),
#'           .options = furrr::furrr_options(
#'             globals = FALSE,
#'             seed = NULL
#'           ),
#'           point_weight_worker
#'         )
#'
#'         # Drain this batch's rasters to gpkg
#'         write_weight_tifs(ws_pattern, p)
#'       }
#'     })
#'   } else {
#'     # ── Sequential: process all, then write ─────────────────────────────
#'     splt_lst <- suppressWarnings(terra::split(
#'       target_o_sub,
#'       seq_len(max(n_cores, 1L))
#'     ))
#'     splt_lst <- splt_lst[!sapply(splt_lst, is.null)]
#'
#'     progressr::with_progress(enable = verbose, {
#'       p <- progressr::progressor(
#'         steps = length(target_o_sub) * length(weighting_scheme_o)
#'       )
#'
#'       furrr::future_pmap(
#'         list(
#'           x = splt_lst,
#'           output_filename = list(output_filename),
#'           temp_dir_sub = list(temp_dir_sub),
#'           temp_dir_sub2 = list(temp_hw),
#'           weighting_scheme_o = list(weighting_scheme_o),
#'           inv_function = list(inv_function),
#'           verbose = list(verbose)
#'         ),
#'         .options = furrr::furrr_options(
#'           globals = FALSE,
#'           seed = NULL
#'         ),
#'         point_weight_worker
#'       )
#'     })
#'
#'     if (verbose) {
#'       message("Writing weight outputs")
#'     }
#'     write_weight_tifs(ws_pattern)
#'   }
#'
#'   # Save metadata
#'   if (verbose) {
#'     message("Saving weight metadata")
#'   }
#'   meta <- dplyr::bind_rows(target_o_sub) |>
#'     tibble::as_tibble() |>
#'     dplyr::select(link_id, unn_group) |>
#'     dplyr::mutate(target_o_type = target_o_type)
#'   sf::write_sf(
#'     meta,
#'     output_filename$outfile,
#'     layer = "unnest_catchment",
#'     append = TRUE,
#'     delete_layer = FALSE,
#'     delete_dsn = FALSE
#'   )
#' }
#'
#' #' Collect unnest raster results while the future runs
#' #' @noRd
#' collect_unnest_rasters <- function(future_obj, temp_dir) {
#'   future_status <- future::futureOf(future_obj)
#'   rast_out <- list()
#'
#'   while (!future::resolved(future_status)) {
#'     Sys.sleep(0.2)
#'     fl <- list.files(temp_dir, "unnest_", full.names = TRUE)
#'     if (length(fl) == 0) {
#'       next
#'     }
#'     rasts <- purrr::map(fl, ~ try(terra::rast(.), silent = TRUE))
#'     rasts <- rasts[!sapply(rasts, inherits, "try-error")]
#'     if (length(rasts) > 0) {
#'       rast_out <- c(rast_out, purrr::map(rasts, terra::unique))
#'       suppressWarnings(file.remove(unlist(purrr::map(rasts, terra::sources))))
#'     }
#'   }
#'
#'   # Check for errors
#'   if (length(future_obj$result$conditions) > 0) {
#'     err <- future_obj$result$conditions[[1]]$condition
#'     if (inherits(err, "error")) stop(err)
#'   }
#'
#'   # Final collection
#'   fl <- list.files(temp_dir, "unnest_.*\\.tif$", full.names = TRUE)
#'   rasts <- purrr::map(fl, ~ try(terra::rast(.), silent = TRUE))
#'   rasts <- rasts[!sapply(rasts, inherits, "try-error")]
#'   if (length(rasts) > 0) {
#'     rast_out <- c(rast_out, purrr::map(rasts, terra::unique))
#'     suppressWarnings(file.remove(unlist(purrr::map(rasts, terra::sources))))
#'   }
#'
#'   rast_out
#' }
#'
#' #' Worker for computing O-target weights
#' #' @noRd
#' point_weight_worker <- carrier::crate(
#'   function(
#'     x,
#'     output_filename,
#'     temp_dir_sub,
#'     temp_dir_sub2,
#'     weighting_scheme_o,
#'     inv_function,
#'     verbose
#'   ) {
#'     #suppressPackageStartupMessages(library(sf)) # Not sure why, but this is necessary
#'     #options(dplyr.summarise.inform = FALSE, scipen = 999)
#'     `%>%` <- magrittr::`%>%`
#'
#'     target_S <- terra::rast(file.path(temp_dir_sub, "dem_streams_d8_sub.tif"))
#'     dem <- terra::rast(file.path(temp_dir_sub, "dem_final.tif"))
#'     flow_accum <- terra::rast(file.path(temp_dir_sub, "dem_accum_d8.tif"))
#'
#'     purrr::map(x, function(y) {
#'       if (verbose) {
#'         message("Calculating weights for: ", y$unn_group[[1]])
#'       }
#'
#'       hw_dir <- file.path(temp_dir_sub2, basename(tempfile()))
#'       dir.create(hw_dir)
#'
#'       suppressMessages(
#'         hw_o <- ihydro::hydroweight(
#'           hydroweight_dir = hw_dir,
#'           target_O = y,
#'           target_S = target_S,
#'           target_uid = paste0("unnest_group_", y$unn_group[[1]]),
#'           OS_combine = FALSE,
#'           dem = dem,
#'           flow_accum = flow_accum,
#'           weighting_scheme = weighting_scheme_o,
#'           inv_function = inv_function,
#'           clean_tempfiles = FALSE,
#'           return_products = TRUE,
#'           wrap_return_products = FALSE,
#'           save_output = FALSE
#'         )
#'       )
#'
#'       for (r in hw_o) {
#'         terra::writeRaster(
#'           r,
#'           file.path(
#'             temp_dir_sub2,
#'             paste0(
#'               names(r),
#'               "_unn_group",
#'               y$unn_group[[1]],
#'               ".tif"
#'             )
#'           ),
#'           datatype = "FLT4S",
#'           todisk = TRUE,
#'           overwrite = TRUE,
#'           gdal = "COMPRESS=NONE"
#'         )
#'       }
#'
#'       rm(hw_o)
#'       suppressWarnings(unlink(hw_dir, recursive = TRUE, force = TRUE))
#'       gr <- gc(verbose = FALSE)
#'       NULL
#'     })
#'   }
#' )
#'
#' #' Drain weight raster files to gpkg
#' #' @noRd
#' drain_weight_rasters <- function(
#'   future_proc,
#'   temp_dir,
#'   output_file,
#'   weighting_scheme_o,
#'   p
#' ) {
#'   future_status <- future::futureOf(future_proc)
#'   pattern <- paste0(weighting_scheme_o, collapse = "|")
#'
#'   write_batch <- function(min_age = 5) {
#'     fl <- list.files(temp_dir, full.names = TRUE)
#'     fl <- fl[grepl(pattern, fl)]
#'     if (min_age > 0) {
#'       fl <- fl[file.mtime(fl) < Sys.time() - min_age]
#'     }
#'     for (f in fl) {
#'       Sys.sleep(0.2)
#'       r <- try(terra::rast(f), silent = TRUE)
#'       if (inherits(r, "try-error")) {
#'         next
#'       }
#'       lyr_name <- gsub("\\.tif$", "", basename(f))
#'       res <- try(write_raster_gpkg(r, output_file, lyr_name), silent = TRUE)
#'       if (inherits(res, "try-error")) {
#'         msg <- conditionMessage(attr(res, "condition"))
#'         if (!msg %in% c("stoi", "stol")) stop(msg)
#'       }
#'
#'       # r2 <- r^2
#'       # lyr_name2 <- lyr_name
#'       # lyr_name2 <- gsub("iFLO", "iFLOSQ", lyr_name2)
#'       # names(r2) <- gsub("iFLO", "iFLOSQ", names(r2))
#'       # terra::varnames(r2) <- lyr_name2
#'       # res <- try(
#'       #   write_raster_gpkg(r2, output_filename$outfile, lyr_name2),
#'       #   silent = TRUE
#'       # )
#'       # if (inherits(res, "try-error")) {
#'       #   msg <- conditionMessage(attr(res, "condition"))
#'       #   if (!msg %in% c("stoi", "stol")) stop(msg)
#'       # }
#'
#'       p()
#'       try(file.remove(f), silent = TRUE)
#'     }
#'   }
#'
#'   while (!future::resolved(future_status)) {
#'     Sys.sleep(0.2)
#'     write_batch(5)
#'   }
#'
#'   Sys.sleep(5)
#'
#'   stop_on_future_errors(future_proc)
#'   #if (length(future_proc$result$conditions) > 0) {
#'   #  errs <- purrr::keep(
#'   #    purrr::map(future_proc$result$conditions, "condition"),
#'   #    ~ inherits(., "error")
#'   #  )
#'   #  if (length(errs) > 0) stop(paste(errs, collapse = "\n"))
#'   #}
#'
#'   write_batch(0)
#'   unlink(temp_dir, recursive = TRUE, force = TRUE)
#' }
