#' Attribute sampling points with distance-weighted landscape summaries
#'
#' Extracts weighted summaries of layers of interest (LOI) from upstream
#' catchments and returns a table of attributes per sampling point.
#' Supports numeric and categorical rasters, multiple inverse-distance
#' weighting schemes, and automatic retry on memory errors.
#'
#' @details
#' The function performs the following steps:
#' 1. Normalizes all file inputs to `ihydro` objects for consistent handling.
#' 2. Validates the presence of LOI data and resolves target sampling points or link IDs.
#' 3. Validates the presence of or prepares inverse distance weights if specified in the weighting scheme.
#' 4. Retrieves upstream catchments for the target points.
#' 5. Extracts raster attributes for each catchment polygon, applying the specified weighting schemes.
#' 6. Joins the extracted attributes with the target IDs and writes the final output to a CSV file.
#' 7. If any extractions fail during parallel processing (e.g., due to memory issues), it retries those
#' sequentially to ensure maximum completion.
#'
#' @param input An `ihydro` object (from [process_hydrology()]).
#' @param out_filename Character. CSV output path.
#' @param loi_file Optional `ihydro` object or `.gpkg` path with LOI layers.
#'   Defaults to the layers inside `input`.
#' @param loi_cols Character vector of LOI column names, or `NULL` for all.
#' @param iDW_file Optional `ihydro` object or `.gpkg` path for pre-computed
#'   inverse-distance weights. If `NULL`, weights are computed on the fly.
#' @param store_iDW Logical. If `TRUE` and `iDW_file` is `NULL`, weights are
#'   written into the `input` GeoPackage.
#' @param sample_points Character vector of site IDs, or `NULL` for all.
#' @param link_id Character vector of link IDs, or `NULL`.
#' @param target_o_type One of `"point"`, `"segment_point"`, `"segment_whole"`.
#' @param weighting_scheme Character vector. Any of `"lumped"`, `"iFLS"`,
#'   `"HAiFLS"`, `"iFLO"`, `"HAiFLO"`.
#' @param loi_numeric_stats Character vector. Any of `"mean"`, `"sd"`,
#'   `"median"`, `"min"`, `"max"`, `"sum"`.
#' @param inv_function Inverse-distance function (see [prep_weights()]).
#' @param write_strategy Character. How processed weight rasters are written to the
#'   GeoPackage. `"sequential"` (default) waits for all parallel workers to
#'   finish, then writes every raster in one pass. `"batched"` processes unnest
#'   groups in chunks, writing to the GeoPackage between chunks to reduce peak
#'   temporary disk usage.
#' @param temp_dir Temporary directory for intermediate files.
#' @param verbose Logical.
#'
#' @return A data.frame of weighted attributes, also written to `out_filename`.
#' @export
#'
#' @seealso [process_hydrology()], [process_loi()], [prep_weights()],
#'   [attrib_points()]
#'
#' @rdname fasttrib_points
#'
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
#' # Process LOI
#' ex_loc <- tempdir()
#'
#' ex_dem <- ex_data("elev_ned_30m.tif")
#' ex_data("landuse_r.tif") %>%
#'  setNames("LC") %>%
#'  terra::writeRaster(file.path(ex_loc, "LC.tif"), overwrite = T)
#'
#' terra::writeVector(
#'  ex_data("geology.shp"),
#'  file.path(ex_loc, "geology.shp"),
#'  overwrite = T
#' )
#' terra::writeVector(
#'  ex_data("pointsources.shp"),
#'  file.path(ex_loc, "pointsources.shp"),
#'  overwrite = T
#' )
#'
#' landuse_r_path <- file.path(ex_loc, "LC.tif")
#' geology_path <- file.path(ex_loc, "geology.shp")
#' pointsources_path <- file.path(ex_loc, "pointsources.shp")
#'
#' sf::read_sf(pointsources_path) %>%
#'  mutate(pointsource = "pontsrc") %>%
#'  st_buffer(60) %>%
#'  sf::write_sf(file.path(ex_loc, "pointsources.shp"), overwrite = T)
#'
#' pointsources_path <- file.path(ex_loc, "pointsources.shp")
#'
#' whitebox::wbt_slope(
#'  dem = file.path(ex_loc, "toy_dem.tif"),
#'  output = file.path(ex_loc, "slope.tif")
#' )
#'
#' output_filename_loi <- file.path(ex_loc, "Processed_loi.gpkg")
#'
#' loi_combined <- ihydro::process_loi(
#'  dem = toy_dem,
#'  num_inputs = list(
#'    # Can be given as a mixture of input types (file paths, or any sf or terra format)
#'    slope = file.path(ex_loc, "slope.tif")
#'  ),
#'  cat_inputs = list(
#'    # Can be given as a mixture of input types (file paths, or any sf or terra format)
#'    landcover = landuse_r_path,
#'    geology = geology_path,
#'    pointsources = pointsources_path
#'  ),
#'  variable_names = list(
#'    # any unlisted inputs will be used in their entirety
#'    geology = "GEO_NAME", # names listed here will subset those attributes or layers from the inputs
#'    pointsources = "pontsrc"
#'  ),
#'  output_filename = output_filename_loi,
#'  return_products = T,
#'  temp_dir = NULL,
#'  verbose = T
#' )
#'
#' # Attribute sample points
#' output_csv <- file.path(ex_loc, "sample_points_attributes.csv")
#' ihydro::fasttrib_points(
#'  input = output_gpkg,
#'  out_filename = output_csv,
#'  loi_file = output_filename_loi,
#'  loi_cols = NULL,
#'  iDW_file = NULL,
#'  store_iDW = FALSE,
#'  sample_points = c("1", "25", "80"),
#'  link_id = c("100", "200", "800"),
#'  target_o_type = "segment_point",
#'  weighting_scheme = c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO"),
#'  loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
#'  inv_function = function(x) (x * 0.001 + 1)^-1,
#'  temp_dir = NULL,
#'  verbose = TRUE,
#'  backend = c("tibble", "SQLite", "data.table"),
#'  SQLdb_path = ":memory:"
#' )
#'
#'
#' }
#'

fasttrib_points <- function(
  input,
  out_filename,
  loi_file = NULL,
  loi_cols = NULL,
  iDW_file = NULL,
  store_iDW = FALSE,
  sample_points = NULL,
  link_id = NULL,
  target_o_type = c("point", "segment_point", "segment_whole"),
  weighting_scheme = c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO"),
  loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
  inv_function = function(x) (x * 0.001 + 1)^-1,
  write_strategy = c("sequential", "batched"),
  temp_dir = NULL,
  verbose = FALSE
) {
  # ─ Validate inputs ───────────────────────────
  .check_ihydro(input)
  stopifnot(is.logical(store_iDW), length(store_iDW) == 1L)

  target_o_type <- match.arg(target_o_type)
  weighting_scheme <- match.arg(weighting_scheme, several.ok = TRUE)
  loi_numeric_stats <- match.arg(loi_numeric_stats, several.ok = TRUE)
  loi_numeric_stats <- stats::setNames(loi_numeric_stats, loi_numeric_stats)
  write_strategy <- match.arg(write_strategy)

  input_path <- input$outfile
  temp_dir <- .ensure_temp_dir(temp_dir)
  n_cores <- .n_workers()

  # ─ Configure external tools ───────────────────────
  whitebox::wbt_options(
    exe_path = whitebox::wbt_exe_path(),
    verbose = verbose > 2,
    wd = temp_dir
  )
  old_terra <- .set_terra_opts(temp_dir = temp_dir, verbose = verbose > 3)
  on.exit(.restore_terra_opts(old_terra), add = TRUE)

  old_opts <- options(
    scipen = 999,
    dplyr.summarise.inform = FALSE,
    future.rng.onMisuse = "ignore"
  )
  on.exit(options(old_opts), add = TRUE)

  # ─ Resolve LOI file ──────────────────────────
  loi_file <- .resolve_loi(input, loi_file)
  loi_path <- loi_file$outfile

  loi_meta <- sf::read_sf(loi_path, "loi_meta")
  if (is.null(loi_cols)) {
    loi_cols <- loi_meta$loi_var_nms
  }
  bad_cols <- loi_cols[!loi_cols %in% loi_meta$loi_var_nms]
  if (length(bad_cols) > 0L) {
    cli::cli_abort("LOI columns not found: {.val {bad_cols}}")
  }
  loi_meta <- dplyr::filter(loi_meta, loi_var_nms %in% loi_cols)

  # ─ Resolve target IDs ─────────────────────────
  target_ids <- target_id_fun(
    db_fp = input_path,
    sample_points = sample_points,
    link_id = link_id,
    segment_whole = target_o_type == "segment_whole",
    target_o_type = target_o_type
  )

  # ─ Resolve iDW file and prepare weights ─────────────────
  iDW_path <- .resolve_idw_path(input, iDW_file, store_iDW)
  target_o_meta <- dplyr::mutate(target_ids, unn_group = "1")

  if (!all(weighting_scheme == "lumped")) {
    if (verbose) {
      message("Preparing inverse distance weights")
    }

    # ─ Check if existing weights match the requested target_o_type ────────
    needs_recalc <- FALSE
    if (file.exists(iDW_path)) {
      iDW_lyrs <- tryCatch(
        ihydro_layers(as.ihydro(iDW_path)),
        error = function(cnd) tibble::tibble(layer_name = character(0))
      )
      if ("target_o_meta" %in% iDW_lyrs$layer_name) {
        stored_meta <- tryCatch(
          sf::read_sf(iDW_path, "target_o_meta"),
          error = function(cnd) NULL
        )
        if (
          !is.null(stored_meta) &&
            "target_o_type" %in% colnames(stored_meta) &&
            nrow(stored_meta) > 0L
        ) {
          stored_types <- unique(stored_meta$target_o_type)
          if (!target_o_type %in% stored_types) {
            cli::cli_warn(c(
              "!" = "Stored weights in {.path {iDW_path}} were computed with
                     {.val {stored_types}}, not {.val {target_o_type}}.",
              "i" = "Recalculating weights for {.val {target_o_type}}."
            ))
            needs_recalc <- TRUE
          }
        }
      }
    }

    iDW_out <- prep_weights(
      input = input,
      output_filename = iDW_path,
      sample_points = sample_points,
      link_id = link_id,
      target_o_type = target_o_type,
      weighting_scheme = weighting_scheme[weighting_scheme != "lumped"],
      inv_function = inv_function,
      temp_dir = temp_dir,
      write_strategy = write_strategy,
      verbose = verbose
    )
    if (!inherits(iDW_out, "ihydro")) {
      stop(iDW_out)
    }
    iDW_path <- iDW_out$outfile
    target_o_meta <- tryCatch(
      {
        meta <- sf::read_sf(iDW_path, "target_o_meta")
        if ("target_o_type" %in% colnames(meta)) {
          dplyr::filter(meta, target_o_type == !!target_o_type)
        } else {
          meta
        }
      },
      error = function(cnd) dplyr::mutate(target_ids, unn_group = "1")
    )
  }

  # ─ Build subbasin lookup ────────────────────────
  subb_ids <- ihydro::read_ihydro(input, "us_flowpaths") |>
    dplyr::filter(pour_point_id %in% local(target_ids$link_id)) |>
    dplyr::select(pour_point_id, origin_link_id) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::left_join(
      dplyr::mutate(target_o_meta, link_id = as.character(link_id)),
      by = c("pour_point_id" = "link_id")
    ) |>
    dplyr::arrange(origin_link_id)

  # ─ Pre-compute catchment polygons (main process only) ──────────
  # get_catchment() opens SQLite connections and writes cached results,

  # so it must NOT run inside parallel workers.
  if (verbose) {
    message("Pre-computing catchment polygons")
  }
  all_link_ids <- unique(subb_ids$pour_point_id)
  catchment_polys <- tryCatch(
    get_catchment(
      input,
      link_id = all_link_ids,
      temp_dir = temp_dir,
      verbose = verbose
    ),
    error = function(cnd) {
      cli::cli_warn("get_catchment failed: {conditionMessage(cnd)}")
      sf::st_sf(
        link_id = character(0),
        geom = sf::st_sfc(),
        crs = sf::st_crs(4326)
      )
    }
  )

  # ─ Extract raster attributes ──────────────────────
  temp_dir_sub <- file.path(temp_dir, basename(tempfile()))
  dir.create(temp_dir_sub, recursive = TRUE)
  if (verbose) {
    message("Calculating weighted attributes")
  }

  o1 <- .ft_extract_all(
    catchment_polys = catchment_polys,
    iDW_path = iDW_path,
    loi_path = loi_path,
    weighting_scheme = weighting_scheme,
    loi_numeric_stats = loi_numeric_stats,
    loi_cols = loi_cols,
    loi_meta = loi_meta,
    subb_ids = subb_ids,
    temp_dir_sub = temp_dir_sub,
    verbose = verbose,
    n_cores = n_cores
  )

  # ─ Join results and write ────────────────────────
  target_ids_out <- target_id_fun(
    db_fp = input_path,
    sample_points = sample_points,
    link_id = link_id,
    segment_whole = FALSE,
    target_o_type = target_o_type
  )

  final_out <- dplyr::left_join(
    target_ids_out,
    dplyr::mutate(o1, pour_point_id = as.character(pour_point_id)),
    by = c("link_id" = "pour_point_id"),
    multiple = "all"
  )

  data.table::fwrite(
    x = final_out,
    file = out_filename,
    buffMB = 128L,
    nThread = 1L,
    showProgress = FALSE
  )

  final_out
}


# ────────────────────────────────────────
# Internal helpers ─ NOT exported
# ────────────────────────────────────────

# ─ Input validation & setup ─────────────────────────

#' @noRd
.check_ihydro <- function(x) {
  if (!inherits(x, "ihydro")) {
    cli::cli_abort("{.arg input} must be of class {.cls ihydro}.")
  }
}

#' @noRd
.ensure_temp_dir <- function(temp_dir) {
  if (is.null(temp_dir)) {
    temp_dir <- tempfile()
  }
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  normalizePath(temp_dir)
}

#' @noRd
.n_workers <- function() {
  n <- future::nbrOfWorkers()
  if (is.infinite(n)) {
    n <- future::availableCores(logical = FALSE)
  }
  max(n, 1L)
}

#' Set terra options; return prior values for restoration
#' @noRd
.set_terra_opts <- function(temp_dir, verbose = FALSE) {
  old <- list(
    memfrac = terra::terraOptions(print = F)$memfrac,
    memmin = terra::terraOptions(print = F)$memmin,
    tempdir = terra::terraOptions(print = F)$tempdir,
    verbose = terra::terraOptions(print = F)$verbose
  )
  terra::terraOptions(
    memfrac = 0.8,
    memmin = 1,
    tempdir = temp_dir,
    verbose = verbose
  )
  old
}

#' @noRd
.restore_terra_opts <- function(old) {
  terra::terraOptions(
    memfrac = old$memfrac,
    memmin = old$memmin,
    tempdir = old$tempdir,
    verbose = old$verbose
  )
}

#' Resolve LOI file to an ihydro object
#' @noRd
.resolve_loi <- function(input, loi_file) {
  if (is.null(loi_file)) {
    loi_file <- as.ihydro(input$outfile)
  } else if (inherits(loi_file, "ihydro")) {
    # already good
  } else {
    loi_file <- as.ihydro(loi_file)
  }
  if (!"loi_meta" %in% ihydro_layers(loi_file)$layer_name) {
    cli::cli_abort("No LOI data found in {.arg loi_file}.")
  }
  loi_file
}

#' Resolve iDW path â€” returns a character path to a .gpkg
#' @noRd
.resolve_idw_path <- function(input, iDW_file, store_iDW) {
  if (!is.null(iDW_file)) {
    if (inherits(iDW_file, "ihydro")) {
      return(iDW_file$outfile)
    }
    return(iDW_file)
  }
  if (store_iDW) {
    return(input$outfile)
  }
  tmp <- tempfile()
  dir.create(tmp, recursive = TRUE)
  file.path(tmp, "temp_iDW.gpkg")
}


# ────────────────────────────────────────
# Raster extraction pipeline
# ────────────────────────────────────────

#' Orchestrate parallel extraction across unnest groups
#' @noRd
.ft_extract_all <- function(
  catchment_polys,
  iDW_path,
  loi_path,
  weighting_scheme,
  loi_numeric_stats,
  loi_cols,
  loi_meta,
  subb_ids,
  temp_dir_sub,
  verbose,
  n_cores
) {
  # One task per unnest-group
  tasks <- subb_ids |>
    dplyr::group_by(unn_group) |>
    tidyr::nest() |>
    dplyr::ungroup()

  # ─ First pass (parallel if workers > 1) ─────────────────
  out <- .ft_dispatch_tasks(
    tasks,
    catchment_polys,
    iDW_path,
    loi_path,
    weighting_scheme,
    loi_numeric_stats,
    loi_cols,
    loi_meta,
    temp_dir_sub,
    verbose,
    n_cores,
    parallel = n_cores > 1L
  )

  # ─ Sequential retry for failures ────────────────────
  incomplete <- out$pour_point_id[out$status == "Incomplete"]
  if (length(incomplete) > 0L && n_cores > 1L) {
    if (verbose) {
      message("Retrying ", length(incomplete), " link_id(s) sequentially.")
    }
    retry_tasks <- subb_ids |>
      dplyr::filter(pour_point_id %in% incomplete) |>
      dplyr::group_by(unn_group) |>
      tidyr::nest() |>
      dplyr::ungroup()

    out_retry <- .ft_dispatch_tasks(
      retry_tasks,
      catchment_polys,
      iDW_path,
      loi_path,
      weighting_scheme,
      loi_numeric_stats,
      loi_cols,
      loi_meta,
      temp_dir_sub,
      verbose,
      n_cores = 1L,
      parallel = FALSE
    )

    out <- dplyr::bind_rows(
      dplyr::filter(out, !pour_point_id %in% incomplete),
      out_retry
    )
  }

  out
}


#' Submit tasks via carrier::crate + furrr (parallel) or purrr (sequential)
#'
#' carrier::crate is retained because it serialises closures cleanly for
#' future workers, avoiding accidental capture of large parent-scope objects
#' (consistent with Advanced R Â§10.2.5 on garbage collection in function
#' factories).
#'
#' @noRd
.ft_dispatch_tasks <- function(
  tasks,
  catchment_polys,
  iDW_path,
  loi_path,
  weighting_scheme,
  loi_numeric_stats,
  loi_cols,
  loi_meta,
  temp_dir_sub,
  verbose,
  n_cores,
  parallel = TRUE
) {
  # Build a self-contained worker using carrier::crate.
  # crate captures only explicitly-listed values, preventing the
  # entire parent environment from being serialised to workers.
  worker <- carrier::crate(
    function(
      task_data,
      unn_group,
      catchment_polys,
      iDW_path,
      loi_path,
      loi_meta,
      weighting_scheme,
      loi_numeric_stats,
      loi_cols,
      temp_dir_sub,
      verbose
    ) {
      ihydro:::.ft_process_group(
        task_data = task_data,
        unn_group = unn_group,
        catchment_polys = catchment_polys,
        iDW_path = iDW_path,
        loi_path = loi_path,
        loi_meta = loi_meta,
        weighting_scheme = weighting_scheme,
        loi_numeric_stats = loi_numeric_stats,
        loi_cols = loi_cols,
        temp_dir_sub = temp_dir_sub,
        verbose = verbose
      )
    }
  )

  args <- list(
    task_data = tasks$data,
    unn_group = tasks$unn_group,
    catchment_polys = list(catchment_polys),
    iDW_path = list(iDW_path),
    loi_path = list(loi_path),
    loi_meta = list(loi_meta),
    weighting_scheme = list(weighting_scheme),
    loi_numeric_stats = list(loi_numeric_stats),
    loi_cols = list(loi_cols),
    temp_dir_sub = list(temp_dir_sub),
    verbose = list(verbose)
  )

  if (parallel) {
    results <- furrr::future_pmap(
      args,
      worker,
      .options = furrr::furrr_options(
        globals = FALSE,
        seed = NULL,
        scheduling = 1L # feed one task at a time to free workers
      )
    )
  } else {
    results <- purrr::pmap(args, worker)
  }

  dplyr::bind_rows(results)
}


#' Process one unnest-group: load rasters, iterate over catchments
#' @noRd
.ft_process_group <- function(
  task_data,
  unn_group,
  catchment_polys,
  iDW_path,
  loi_path,
  loi_meta,
  weighting_scheme,
  loi_numeric_stats,
  loi_cols,
  temp_dir_sub,
  verbose
) {
  # Configure terra inside worker
  terra::terraOptions(
    memfrac = 0.8,
    memmin = 1,
    tempdir = temp_dir_sub,
    verbose = FALSE
  )

  # Load LOI rasters
  loi_rasts <- terra::rast(loi_path, loi_cols)

  # Stream-targeted iDW rasters (iFLS / HAiFLS)
  ws_s <- intersect(weighting_scheme, c("iFLS", "HAiFLS"))
  iDWs_rasts <- if (length(ws_s) > 0L) terra::rast(iDW_path, ws_s)

  # Point-targeted iDW rasters (iFLO / HAiFLO) â€” group-specific
  ws_o <- intersect(weighting_scheme, c("iFLO", "HAiFLO"))
  iDWo_rasts <- NULL
  if (length(ws_o) > 0L) {
    o_lyrs <- paste0(
      rep(ws_o, each = length(unique(unn_group))),
      "_unn_group",
      rep(unique(unn_group), times = length(ws_o))
    )
    iDWo_rasts <- terra::rast(iDW_path, o_lyrs)
    names(iDWo_rasts) <- gsub(
      paste0("_unn_group", unn_group),
      "",
      names(iDWo_rasts)
    )
  }

  input_rasts <- c(loi_rasts, iDWs_rasts, iDWo_rasts)

  # Subset pre-computed catchments for this group's pour points
  group_ids <- unique(task_data$pour_point_id)
  input_poly <- catchment_polys[catchment_polys$link_id %in% group_ids, ]

  if (nrow(input_poly) == 0L) {
    return(tibble::tibble(
      pour_point_id = group_ids,
      status = "Incomplete"
    ))
  }

  # Process each catchment
  results <- vector("list", nrow(input_poly))
  for (i in seq_along(results)) {
    results[[i]] <- .ft_one_catchment(
      sub_poly = input_poly[i, ],
      input_rasts = input_rasts,
      loi_cols = loi_cols,
      weighting_scheme = weighting_scheme,
      loi_meta = loi_meta,
      loi_numeric_stats = loi_numeric_stats,
      temp_dir_sub = temp_dir_sub
    )
  }

  dplyr::bind_rows(results)
}


#' Extract and summarise raster data for one catchment polygon
#' @noRd
.ft_one_catchment <- function(
  sub_poly,
  input_rasts,
  loi_cols,
  weighting_scheme,
  loi_meta,
  loi_numeric_stats,
  temp_dir_sub
) {
  sub_id <- sub_poly$link_id
  temp_fl <- tempfile(
    pattern = "ihydro",
    tmpdir = temp_dir_sub,
    fileext = ".tif"
  )
  on.exit(try(file.remove(temp_fl), silent = TRUE), add = TRUE)

  # Rasterize catchment to get cell indices
  sub_rast <- tryCatch(
    terra::rasterize(
      terra::vect(sub_poly),
      input_rasts[[1]],
      fun = sum,
      field = 1,
      filename = temp_fl
    ),
    error = function(cnd) NULL
  )
  if (is.null(sub_rast)) {
    return(tibble::tibble(pour_point_id = sub_id, status = "Incomplete"))
  }
  cells <- terra::cells(sub_rast)
  if (length(cells) == 0L) {
    return(tibble::tibble(pour_point_id = sub_id, status = "Incomplete"))
  }

  # Full extraction with memory-error fallback
  tryCatch(
    .ft_extract_and_summarise(
      input_rasts,
      cells,
      sub_id,
      loi_cols,
      weighting_scheme,
      loi_meta,
      loi_numeric_stats
    ),
    error = function(cnd) {
      # Retry with chunked extraction
      tryCatch(
        .ft_chunked_extract(
          input_rasts,
          cells,
          sub_id,
          loi_cols,
          weighting_scheme,
          loi_meta,
          loi_numeric_stats
        ),
        error = function(cnd2) {
          tibble::tibble(pour_point_id = sub_id, status = "Incomplete")
        }
      )
    }
  )
}


#' Full extraction: all rasters at once
#' @noRd
.ft_extract_and_summarise <- function(
  input_rasts,
  cells,
  point_id,
  loi_cols,
  weighting_scheme,
  loi_meta,
  loi_numeric_stats
) {
  ot <- terra::extract(input_rasts, cells)
  .ft_attr(
    df = ot,
    point_id = point_id,
    weighting_scheme = weighting_scheme,
    loi_meta = loi_meta,
    loi_cols = loi_cols,
    loi_numeric_stats = loi_numeric_stats
  )
}


#' Chunked extraction: IDW first, then LOI in safe-sized chunks
#' @noRd
.ft_chunked_extract <- function(
  input_rasts,
  cells,
  point_id,
  loi_cols,
  weighting_scheme,
  loi_meta,
  loi_numeric_stats
) {
  ws_nonlumped <- weighting_scheme[weighting_scheme != "lumped"]

  if (length(ws_nonlumped) > 0L && any(ws_nonlumped %in% names(input_rasts))) {
    idw_names <- intersect(ws_nonlumped, names(input_rasts))
    ot_idw <- terra::extract(terra::subset(input_rasts, idw_names), cells)
    if (ncol(ot_idw) == 1L) colnames(ot_idw) <- idw_names
  } else {
    ot_idw <- data.frame(lumped = 1)
  }

  # Estimate chunk size from free RAM
  bytes_per_col <- length(cells) * 8L
  free_bytes <- as.numeric(terra::free_RAM()) * 1024
  max_cols <- max(floor((free_bytes * 0.4) / bytes_per_col), 1L)
  chunk_size <- min(max_cols, length(loi_cols))
  loi_chunks <- split(loi_cols, ceiling(seq_along(loi_cols) / chunk_size))

  chunk_results <- lapply(loi_chunks, function(chunk) {
    avail <- intersect(chunk, names(input_rasts))
    if (length(avail) == 0L) {
      return(NULL)
    }
    ot_loi <- terra::extract(terra::subset(input_rasts, avail), cells)
    if (ncol(ot_loi) == 1L) {
      colnames(ot_loi) <- avail
    }

    .ft_attr(
      df = dplyr::bind_cols(ot_loi, ot_idw),
      point_id = point_id,
      weighting_scheme = weighting_scheme,
      loi_meta = loi_meta,
      loi_cols = chunk,
      loi_numeric_stats = loi_numeric_stats
    )
  })

  chunk_results <- Filter(Negate(is.null), chunk_results)
  if (length(chunk_results) == 0L) {
    return(tibble::tibble(pour_point_id = point_id, status = "Incomplete"))
  }
  purrr::reduce(
    chunk_results,
    dplyr::left_join,
    by = c("pour_point_id", "status")
  )
}


# ────────────────────────────────────────
# Attribute computation (data.table backend only via dtplyr)
# ────────────────────────────────────────

#' Compute weighted and lumped summary statistics
#'
#' Pure computation â€” no I/O, no raster work. Uses dtplyr (data.table
#' backend via dplyr syntax) for fast grouped operations.
#'
#' @noRd
.ft_attr <- function(
  df,
  point_id,
  weighting_scheme,
  loi_meta,
  loi_cols,
  loi_numeric_stats
) {
  stopifnot(is.data.frame(df))

  # NaN -> NA
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], function(x) {
    x[is.nan(x)] <- NA_real_
    x
  })

  loi_meta <- dplyr::filter(
    loi_meta,
    loi_var_nms %in% loi_cols,
    loi_var_nms %in% colnames(df)
  )
  weighting_scheme <- weighting_scheme[
    weighting_scheme %in% c("lumped", colnames(df))
  ]

  numb_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "num_rast"]
  cat_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "cat_rast"]

  # dtplyr lazy data.table

  dt <- dtplyr::lazy_dt(df)

  # Pre-compute weighted columns: loi * weight
  for (ws in intersect(
    c("iFLS", "iFLO", "HAiFLS", "HAiFLO"),
    weighting_scheme
  )) {
    dt <- dt |>
      dplyr::mutate(dplyr::across(
        tidyselect::any_of(loi_cols),
        ~ . * (!!rlang::sym(ws)),
        .names = paste0("{.col}_", ws)
      ))
  }
  dt <- dplyr::compute(dt)

  # ─ Collect results ───────────────────────────
  results <- list()

  # Lumped statistics
  if ("lumped" %in% weighting_scheme) {
    results <- c(
      results,
      .ft_lumped_stats(dt, numb_rast, cat_rast, loi_cols, loi_numeric_stats)
    )
  }

  # Weighted mean / proportions
  ws_active <- weighting_scheme[weighting_scheme != "lumped"]
  if (
    length(ws_active) > 0L &&
      (any(loi_numeric_stats == "mean") || length(cat_rast) > 0L)
  ) {
    results <- c(results, list(.ft_weighted_mean(dt, numb_rast, cat_rast)))
  }

  # Weighted SD
  if (
    length(ws_active) > 0L &&
      any(loi_numeric_stats %in% c("sd", "stdev")) &&
      length(numb_rast) > 0L
  ) {
    results <- c(results, list(.ft_weighted_sd(dt, numb_rast, ws_active)))
  }

  # Assemble
  results <- Filter(function(x) !is.null(x) && nrow(x) > 0L, results)

  dplyr::bind_cols(
    tibble::tibble(pour_point_id = point_id, status = "Complete"),
    results
  ) |>
    dplyr::select(
      tidyselect::any_of("pour_point_id"),
      tidyselect::any_of("status"),
      tidyselect::contains(loi_meta$loi_var_nms)
    )
}


#' Lumped (unweighted) summary statistics
#' @noRd
.ft_lumped_stats <- function(
  dt,
  numb_rast,
  cat_rast,
  loi_cols,
  loi_numeric_stats
) {
  results <- list()

  if ("mean" %in% loi_numeric_stats) {
    r <- dt |>
      dplyr::select(tidyselect::any_of(loi_cols)) |>
      dplyr::summarise(
        dplyr::across(
          tidyselect::any_of(numb_rast),
          ~ sum(., na.rm = TRUE) / sum(!is.na(.), na.rm = TRUE)
        ),
        dplyr::across(
          tidyselect::any_of(cat_rast),
          ~ sum(., na.rm = TRUE) / dplyr::n()
        )
      ) |>
      dplyr::collect()

    if (length(numb_rast) > 0L) {
      r <- dplyr::rename_with(
        r,
        .cols = tidyselect::any_of(numb_rast),
        ~ paste0(., "_lumped_mean")
      )
    }
    if (length(cat_rast) > 0L) {
      r <- r |>
        dplyr::rename_with(
          .cols = tidyselect::any_of(cat_rast),
          ~ paste0(., "_lumped_prop")
        ) |>
        dplyr::mutate(dplyr::across(
          tidyselect::ends_with("_prop"),
          ~ ifelse(is.na(.), 0, .)
        ))
    }
    results <- c(results, list(r))
  }

  if ("sd" %in% loi_numeric_stats && length(numb_rast) > 0L) {
    results <- c(
      results,
      list(
        dt |>
          dplyr::select(tidyselect::any_of(loi_cols)) |>
          dplyr::summarise(dplyr::across(
            tidyselect::any_of(numb_rast),
            ~ stats::sd(., na.rm = TRUE)
          )) |>
          dplyr::collect() |>
          dplyr::rename_with(~ paste0(., "_lumped_sd"))
      )
    )
  }

  stat_fns <- list(
    min = function(x) min(x, na.rm = TRUE),
    max = function(x) max(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE)
  )

  for (stat_name in intersect(names(stat_fns), loi_numeric_stats)) {
    if (length(numb_rast) > 0L) {
      collected <- dt |>
        dplyr::select(tidyselect::any_of(numb_rast)) |>
        dplyr::collect() |>
        as.data.frame()
      summ <- vapply(collected, stat_fns[[stat_name]], numeric(1))
      names(summ) <- paste0(names(summ), "_lumped_", stat_name)
      results <- c(results, list(tibble::as_tibble(as.list(summ))))
    }
  }

  if ("count" %in% loi_numeric_stats) {
    results <- c(
      results,
      list(
        dt |>
          dplyr::select(tidyselect::any_of(loi_cols)) |>
          dplyr::summarise(dplyr::across(
            tidyselect::everything(),
            ~ sum(!is.na(.), na.rm = TRUE)
          )) |>
          dplyr::collect() |>
          dplyr::rename_with(~ paste0(., "_lumped_count"))
      )
    )
  }

  results
}


#' Weighted means for all active IDW schemes
#' @noRd
.ft_weighted_mean <- function(dt, numb_rast, cat_rast) {
  r <- dt |>
    dplyr::summarize(
      dplyr::across(
        tidyselect::ends_with("_iFLS"),
        ~ sum(., na.rm = TRUE) / sum(!!rlang::sym("iFLS"), na.rm = TRUE)
      ),
      dplyr::across(
        tidyselect::ends_with("_HAiFLS"),
        ~ sum(., na.rm = TRUE) / sum(!!rlang::sym("HAiFLS"), na.rm = TRUE)
      ),
      dplyr::across(
        tidyselect::ends_with("_iFLO"),
        ~ sum(., na.rm = TRUE) / sum(!!rlang::sym("iFLO"), na.rm = TRUE)
      ),
      dplyr::across(
        tidyselect::ends_with("_HAiFLO"),
        ~ sum(., na.rm = TRUE) / sum(!!rlang::sym("HAiFLO"), na.rm = TRUE)
      )
    ) |>
    dplyr::collect()

  if (length(numb_rast) > 0L) {
    r <- dplyr::rename_with(
      r,
      .cols = tidyselect::starts_with(paste0(numb_rast, "_")),
      ~ paste0(., "_mean")
    )
  }
  if (length(cat_rast) > 0L) {
    r <- r |>
      dplyr::rename_with(
        .cols = tidyselect::starts_with(paste0(cat_rast, "_")),
        ~ paste0(., "_prop")
      ) |>
      dplyr::mutate(dplyr::across(
        tidyselect::ends_with("_prop"),
        ~ ifelse(is.na(.), 0, .)
      ))
  }
  r
}


#' Weighted standard deviations
#' @noRd
.ft_weighted_sd <- function(dt, numb_rast, ws_active) {
  ws_sd <- dt |>
    dplyr::select(
      tidyselect::starts_with(numb_rast),
      tidyselect::any_of(ws_active)
    )

  for (ws in ws_active) {
    sym_ws <- rlang::sym(ws)
    ws_sd <- ws_sd |>
      dplyr::mutate(
        dplyr::across(
          tidyselect::any_of(numb_rast),
          ~ (!!sym_ws) *
            ((. -
              (sum(. * (!!sym_ws), na.rm = TRUE) /
                sum(!!sym_ws, na.rm = TRUE)))^2),
          .names = paste0("{.col}_", ws, "_term1")
        ),
        dplyr::across(
          tidyselect::any_of(numb_rast),
          ~ ((sum(!!sym_ws != 0, na.rm = TRUE) - 1) /
            sum(!!sym_ws != 0, na.rm = TRUE)) *
            sum(!!sym_ws, na.rm = TRUE),
          .names = paste0("{.col}_", ws, "_term2")
        )
      )
  }

  ws_sd |>
    dplyr::summarize(
      dplyr::across(tidyselect::ends_with("_term1"), ~ sum(., na.rm = TRUE)),
      dplyr::across(tidyselect::ends_with("_term2"), ~ .[1])
    ) |>
    dplyr::collect() |>
    tidyr::pivot_longer(cols = tidyselect::everything()) |>
    dplyr::mutate(
      attr = stringr::str_split_fixed(
        name,
        "_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",
        2
      )[, 1],
      term = stringr::str_split_fixed(
        name,
        "_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",
        2
      )[, 2]
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(hw = gsub(paste0(attr, "_", "|", "_", term, ""), "", name)) |>
    dplyr::ungroup() |>
    dplyr::mutate(name = paste0(attr, "_", hw, "_sd")) |>
    dplyr::group_by(name) |>
    dplyr::summarize(
      sd = sqrt(value[term == "term1"] / value[term == "term2"])
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = name, values_from = sd)
}
