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
#' @param weighting_scheme Character vector. Any of `"lumped"`, `"iFLS"`,
#'   `"HAiFLS"`, `"iFLO"`, `"HAiFLO"`.
#' @param loi_numeric_stats Character vector. Any of `"mean"`, `"sd"`,
#'   `"median"`, `"min"`, `"max"`, `"sum"`.
#' @param mem_fraction Numeric. Value between 0.1 and 0.9 (larger values give a warning).
#'   The fraction of RAM that may be used by the program.
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
#' library(ihydro)
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
#' toy_dem <- terra::writeRaster(ex_dem, file.path(ex_loc, "toy_dem.tif"), overwrite = TRUE)
#' ex_data("landuse_r.tif") %>%
#'  setNames("LC") %>%
#'  terra::writeRaster(file.path(ex_loc, "LC.tif"), overwrite = TRUE)
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
#'  dplyr::mutate(pointsource = "pontsrc") %>%
#'  sf::st_buffer(60) %>%
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
#'  input = ihydro::as_ihydro(output_gpkg),
#'  out_filename = output_csv,
#'  loi_file = ihydro::as_ihydro(output_filename_loi),
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
#'  verbose = TRUE
#' )
#'
#'
#' }
#'

fasttrib_points <- function(
    input,
    out_filename = NULL,
    loi_file = NULL,
    loi_cols = NULL,
    iDW_file = NULL,
    sample_points = NULL,
    link_id = NULL,
    weighting_scheme = c("lumped","iFLS","HAiFLS","iFLO","HAiFLO"),
    loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
    mem_fraction = 0.5,
    n_batches = NULL,
    temp_dir = NULL,
    verbose = FALSE,
    ...
) {
  # ─ Validate inputs ───────────────────────────
  check_ihydro(input)

  stopifnot(mem_fraction < 0.9 | mem_fraction > 0.1)
  if (is.null(n_batches)) {
    n_batches <- 1L
  } else {
    stopifnot(is.numeric(n_batches))
  }

  loi_numeric_stats <- match.arg(
    loi_numeric_stats,
    several.ok = TRUE,
    choices = c("mean", "sd", "median", "min", "max", "sum")
  )
  loi_numeric_stats <- stats::setNames(loi_numeric_stats, loi_numeric_stats)

  input_path <- input$outfile
  temp_dir <- ensure_temp_dir(temp_dir)
  n_cores <- n_workers()

  # ─ Process iDW ───────────────────────

  if (any(weighting_scheme != "lumped")){
    if (is.null(iDW_file)) {
      temp_dir_sub <- file.path(temp_dir,basename(tempfile()))
      dir.create(temp_dir_sub,showWarnings = F,recursive = T)
      on.exit(unlink(temp_dir_sub,recursive = T, force = T), add = TRUE)
      args <- list(...)
      inv_function <- args$inv_function
      if (is.null(inv_function)) {
        inv_function <- function(x) (x * 0.001 + 1)^-1
      }

      fun_dep <- trimws(deparse(substitute(inv_function))[3])

      if (verbose) {
        cli::cli_alert_info(c(
          "{.arg iDW_file} is Null. Calculating temporary iDW layers using: {.var {fun_dep}}.\n",
          "To use a different inverse distance function specify {.arg inv_function} in {.fun fasttrib_points}\n",
          "or generate a permanent iDW file with {.fun prep_weights}"
        ))
      }


      iDW_file <- file.path(temp_dir,"temp_idw.gpkg")
      iDW_file <- prep_weights(
        input = input,
        output_filename = iDW_file,
        weighting_scheme = weighting_scheme[weighting_scheme != "lumped"],
        temp_dir = temp_dir_sub,
        verbose = FALSE,
        inv_function
      )
    }

    check_ihydro(iDW_file)
  }

  # ─ Configure external tools ───────────────────────
  wbt_opt_orig <- whitebox::wbt_options()
  names(wbt_opt_orig) <- gsub("^whitebox\\.", "", names(wbt_opt_orig))
  on.exit(do.call(whitebox::wbt_options, wbt_opt_orig), add = TRUE)

  whitebox::wbt_options(
    exe_path = whitebox::wbt_exe_path(),
    verbose = verbose > 2,
    wd = temp_dir
  )

  old_terra_opts <- set_terra_options(
    n_cores = 1L,
    temp_dir = temp_dir,
    verbose = verbose > 3
  )
  on.exit(restore_terra_options(old_terra_opts), add = TRUE)

  max_mem <- memuse::Sys.meminfo()$totalram * mem_fraction
  max_square <- memuse::howmany(max_mem)
  max_cells_in_memory <- max_square[[1]] * max_square[[2]]

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
    link_id = link_id
  )

  # ─ Resolve iDW file and prepare weights ─────────────────
  idw_layers <- ihydro_layers(iDW_file)$layer_name
  missing_layers <- setdiff(weighting_scheme, idw_layers)
  missing_layers <- missing_layers[missing_layers != "lumped"]
  if (length(missing_layers) > 0L) {
    cli::cli_abort(c(
      "The following weighting layers are missing from iDW_file: ",
      "x" = "{.val {missing_layers}}"
    ))
  }

  # ─ Build subbasin lookup ────────────────────────
  #unnest_catchment <- sf::read_sf(input, "unnest_catchment")
  site_id_col <- read_ihydro(input, "site_id_col")[[1]]

  subb_ids <- read_ihydro(input, "stream_links_attr") |>
    dplyr::filter(link_lngth > 0) |> # small catchments without streams
    dplyr::filter(link_id %in% target_ids$link_id) |> # small catchments without streams
    dplyr::select(link_id, tidyselect::any_of(site_id_col), USChnLn_To) #|>
  #dplyr::left_join(
  #  dplyr::mutate(unnest_catchment, link_id = as.character(link_id)),
  #  by = c("link_id")
  #)

  #subb_ids <- subb_ids[!is.na(subb_ids$unn_group), ]

  subb_lookup <- ihydro::read_ihydro(
    input,
    "us_flowpaths"
  ) |>
    dplyr::select(
      final_link_id = pour_point_id,
      link_id = origin_link_id #subbasin ID
    ) |>
    dplyr::filter(
      final_link_id %in% subb_ids$link_id
    )

  # Ensure catchments exist
  catch <- get_catchment(
    input = input,
    link_id = subb_ids$link_id,
    temp_dir = temp_dir,
    verbose = verbose,
    return = FALSE
  )

  # Identify numeric and categorical rasters from LOI metadata
  numb_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "num_rast"]
  cat_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "cat_rast"]

  # Setup and execute parallel extraction of attributes for each subbasin
  fun_sel <- unique(c("sum",loi_numeric_stats))

  sub_summ <- .extract_planner(
    input = input,
    subb_ids = subb_ids,
    loi_rast_input = loi_file,
    numb_rast = numb_rast,
    cat_rast = cat_rast,
    iDW_rast_input = iDW_file,
    iDW_cols = weighting_scheme,
    max_cells_in_memory = max_cells_in_memory,
    n_cores = n_cores,
    chunks_per_worker = n_batches,
    fun = fun_sel,
    quantiles = NULL,
    temp_dir = temp_dir,
    verbose = verbose
  )
  sub_summ <- unlist(sub_summ, recursive = FALSE)

  if (verbose) {
    message("Extracting attributes for each subbasin...")
  }

  gc(verbose = FALSE)

  progressr::with_progress(enable = verbose, {
    total_tasks <- length(sub_summ)
    p <- progressr::progressor(steps = total_tasks)

    sub_summ <- lapply(sub_summ, function(args) {
      args$progressor <- p
      args
    })

    extract_execute <- lapply(sub_summ, function(args) {
      future::futureCall(
        .extract_subbasins,
        args = args,
        seed = NULL,
        globals = c("args"),
        packages = c("sf", "terra", "exactextractr", "dplyr")
      )
    })
    extract_value <- lapply(extract_execute, future::value)
  })

  if (verbose) {
    message("Combining subbasin attributes across catchments.")
  }

  # Combine extracted attributes into a single table
  extract_value_comb <- list()
  extract_value_nms <- names(extract_value)
  ws_s_nms <- extract_value_nms[grepl("^ws_s\\.", extract_value_nms)]
  #ws_o_nms <- extract_value_nms[grepl("^ws_o\\.", extract_value_nms)]
  ws_lump_nms <- extract_value_nms[grepl("^ws_lump\\.", extract_value_nms)]

  extract_value_comb$ws_s <- dplyr::bind_rows(extract_value[ws_s_nms])
  #extract_value_comb$ws_o <- dplyr::bind_rows(extract_value[ws_o_nms])
  extract_value_comb$ws_lump <- dplyr::bind_rows(extract_value[ws_lump_nms])

  if (length(ws_s_nms) > 0) {
    extract_value_comb$ws_s <- dplyr::left_join(
      dplyr::rename(
        subb_lookup,
        subbasin_link_id = link_id,
        link_id_otarget = final_link_id
      ),
      dplyr::select(
        extract_value_comb$ws_s,
        -link_id_otarget
      ),
      by = "subbasin_link_id"
    )
  }

  # if (length(ws_o_nms) > 0) {
  #   extract_value_comb$ws_o <- extract_value_comb$ws_o |>
  #     dplyr::left_join(
  #       dplyr::select(subb_ids, link_id), #unn_group
  #       by = c("link_id_otarget" = "link_id")
  #     )
  # }

  # Summarize catchment attributes and join with target IDs
  result <- .summarize_catchment(
    extract_value = extract_value_comb,
    weighting_scheme = weighting_scheme,
    loi_numeric_stats = loi_numeric_stats,
    numeric_vars = numb_rast,
    cat_vars = cat_rast
  )

  result <- subb_ids |>
    dplyr::select(link_id, tidyselect::any_of(site_id_col)) |>
    dplyr::left_join(
      result,
      by = "link_id"
    )

  if (!is.null(out_filename)){
    write.csv(
      result,
      out_filename
    )
  }


  return(result)
}


# ────────────────────────────────────────
# Internal helpers ─ NOT exported
# ────────────────────────────────────────

# ─ Input validation & setup ─────────────────────────

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

#' Extract raster attributes
#' @noRd
.extract_fun <- function(
    all_rasts,
    subbasins,
    x_cols = NULL,
    weight_cols = NULL,
    #sqweight_cols = NULL,
    default_value = NA_real_,
    default_weight = NA_real_,
    fun,
    quantiles = NULL,
    max_cells_in_memory = 3e+07,
    include_count = FALSE,
    n_cores = 1L,
    temp_dir_sub = NULL
) {
  if (is.null(x_cols) && is.null(weight_cols)) {
    return(NULL)
  }
  all_rasts <- terra::crop(
    all_rasts,
    terra::vect(subbasins)
  )

  weights_resolve <- function(all_rasts, weight_cols = NULL) {
    if (is.null(weight_cols)) {
      return(NULL)
    }
    all_rasts[[weight_cols]]
  }

  # Per-worker terra memory management
  if (is.null(temp_dir_sub)) {
    temp_dir_sub <- tempfile()
  }
  if (!dir.exists(temp_dir_sub)) {
    dir.create(temp_dir_sub)
  }
  old_terra_opts <- ihydro:::set_terra_options(
    n_cores = n_cores,
    temp_dir = temp_dir_sub,
    verbose = FALSE
  )
  on.exit(ihydro:::restore_terra_options(old_terra_opts), add = TRUE)
  on.exit(
    suppressWarnings(unlink(temp_dir_sub, recursive = TRUE, force = TRUE)),
    add = TRUE
  )
  on.exit(
    terra::tmpFiles(
      current=TRUE,
      orphan=FALSE,
      old=FALSE,
      remove=TRUE
    ),
    add=T
  )
  on.exit(gc(verbose = FALSE), add = TRUE)

  if (include_count) {
    all_rasts[["count_internal"]] <- all_rasts[[1]]
    all_rasts[["count_internal"]][] <- 1
    x_cols <- c(x_cols, "count_internal")
    fun <- unique(c(fun, "count"))
  }

  # Ensure x_cols and weight_cols exist in all_rasts
  if (!is.null(x_cols)) {
    missing_x <- setdiff(x_cols, names(all_rasts))
    if (length(missing_x) > 0) {
      stop(
        "The following x_cols are not present in all_rasts: ",
        paste(missing_x, collapse = ", ")
      )
    }
  }
  if (!is.null(weight_cols)) {
    missing_w <- setdiff(weight_cols, names(all_rasts))
    if (length(missing_w) > 0) {
      stop(
        "The following weight_cols are not present in all_rasts: ",
        paste(missing_w, collapse = ", ")
      )
    }
  }

  out <- list()
  if (is.null(weight_cols)) {
    out[[1]] <- exactextractr::exact_extract(
      x = all_rasts[[x_cols]],
      y = subbasins,
      fun = fun,
      quantiles = quantiles,
      default_value = default_value,
      max_cells_in_memory = max_cells_in_memory,
      progress = FALSE,
      force_df = TRUE,
      full_colnames = TRUE
    )
  } else {
    # Duplicate loi layers for each ws
    loi_idx <- rep(x_cols, each = length(weight_cols))
    loi_idxnm <- paste0(
      loi_idx,
      "_",
      rep(weight_cols, length.out = length(loi_idx))
    )

    all_rast2 <- list()
    for (i in loi_idx) {
      idx <- length(all_rast2) + 1
      all_rast2[[idx]] <- all_rasts[[i]]
      names(all_rast2[[idx]]) <- loi_idxnm[[idx]]
      terra::varnames(all_rast2[[idx]]) <- loi_idxnm[[idx]]
    }

    # Duplicate iDW layers for each lio
    iDW_idx <- rep(
      weight_cols,
      length.out = length(x_cols) * length(weight_cols)
    )
    iDW_idxnm <- paste0(iDW_idx, "_", rep(x_cols, each = length(weight_cols)))

    cnt <- 0
    for (i in iDW_idx) {
      cnt <- cnt + 1
      idx <- length(all_rast2) + 1
      all_rast2[[idx]] <- all_rasts[[i]]
      names(all_rast2[[idx]]) <- iDW_idxnm[[cnt]]
      terra::varnames(all_rast2[[idx]]) <- iDW_idxnm[[cnt]]
    }

    # combine new rasters
    all_rast2 <- terra::rast(all_rast2)

    out[[1]] <- exactextractr::exact_extract(
      x = all_rast2[[loi_idxnm]],
      y = subbasins,
      fun = fun,
      quantiles = quantiles,
      weights = all_rast2[[iDW_idxnm]],
      default_value = default_value,
      default_weight = default_weight,
      max_cells_in_memory = max_cells_in_memory,
      progress = FALSE,
      force_df = TRUE,
      full_colnames = TRUE
    )

    # rename columns
    final_names <- colnames(out[[1]])
    for (i in seq_along(loi_idxnm)) {
      final_names <- gsub(loi_idxnm[[i]], loi_idx[[i]], final_names)
    }
    for (i in seq_along(iDW_idxnm)) {
      final_names <- gsub(iDW_idxnm[[i]], iDW_idx[[i]], final_names)
    }

    colnames(out[[1]]) <- final_names

    # Get sum of weights for each loi (to account for NAs)
    all_rasts[[x_cols]] <- terra::clamp(all_rasts[[x_cols]], 1, 1)
    # all_rasts[[x_cols]] <- terra::classify(all_rasts[[x_cols]], c(-Inf, Inf, 1))

    # duplicated iDW cols
    idw_comb <- c(weight_cols) # sqweight_cols
    iDW_idx <- rep(idw_comb, each = length(x_cols))
    iDW_idxnm <- paste0(iDW_idx, "_", rep(x_cols, length.out = length(iDW_idx)))

    all_rast2 <- list()

    gc(verbose = FALSE)
    for (i in iDW_idx) {
      idx <- length(all_rast2) + 1
      all_rast2[[idx]] <- all_rasts[[i]]
      names(all_rast2[[idx]]) <- iDW_idxnm[[idx]]
      terra::varnames(all_rast2[[idx]]) <- iDW_idxnm[[idx]]
    }

    # Duplicate loi layers for each iDW
    loi_idx <- rep(x_cols, length.out = length(iDW_idx))
    loi_idxnm <- paste0(loi_idx, "_", rep(idw_comb, each = length(x_cols)))

    cnt <- 0
    for (i in loi_idx) {
      cnt <- cnt + 1
      idx <- length(all_rast2) + 1
      all_rast2[[idx]] <- all_rasts[[i]]
      names(all_rast2[[idx]]) <- loi_idxnm[[cnt]]
      terra::varnames(all_rast2[[idx]]) <- loi_idxnm[[cnt]]
    }

    # combine new rasters
    all_rast2 <- terra::rast(all_rast2)

    out[[2]] <- exactextractr::exact_extract(
      x = all_rast2[[iDW_idxnm]],
      y = subbasins,
      fun = "weighted_sum",
      weights = all_rast2[[loi_idxnm]],
      default_value = NA_real_,
      default_weight = 0,
      max_cells_in_memory = max_cells_in_memory,
      progress = FALSE,
      force_df = TRUE,
      full_colnames = TRUE
    )

    # rename columns
    final_names <- colnames(out[[2]])
    for (i in seq_along(loi_idxnm)) {
      final_names <- gsub(loi_idxnm[[i]], loi_idx[[i]], final_names)
    }
    for (i in seq_along(iDW_idxnm)) {
      final_names <- gsub(iDW_idxnm[[i]], iDW_idx[[i]], final_names)
    }
    #for (i in x_cols) {
    #  final_names <- gsub(paste0(".", i, "SQ$"), paste0(".", i), final_names)
    #}

    colnames(out[[2]]) <- final_names
  }

  dplyr::bind_cols(out)
}


#' Extract raster attributes for a set of subbassins (for parallel)
#' @noRd
.extract_subbasins <- #carrier::crate(
  function(
    input_file,
    link_id,
    link_id_otarget = NA_character_,
    loi_rast_input,
    loi_summary = TRUE,
    numb_rast = NULL,
    cat_rast = NULL,
    iDW_rast_input = NULL,
    iDW_cols = NULL,
    #iDWSQ_rast_input = NULL,
    #iDWSQ_cols = NULL,
    catch_source = c("Subbasins_poly", "Catchment_poly"),
    median = FALSE,
    quantiles = NULL,
    max_cells_in_memory = 3e+07,
    n_cores = 1L,
    include_count = FALSE,
    progressor = NULL
  ) {
    catch_source <- match.arg(catch_source)

    subbasins <- sf::read_sf(
      input_file,
      query = ihydro:::build_sql_in(catch_source, "link_id", unique(link_id))
    )

    iDW_rasts <- NULL
    #iDWSQ_rasts <- NULL
    loi_rasts <- NULL
    iDW_out <- NULL
    #iDWSQ_out <- NULL
    numb_out <- NULL
    cat_out <- NULL
    numb_out_iDW <- NULL
    cat_out_iDW <- NULL
    extra_out <- NULL

    loi_rasts <- terra::rast(loi_rast_input, c(numb_rast, cat_rast))

    if (!is.null(iDW_rast_input)) {
      iDW_rasts <- terra::rast(iDW_rast_input, iDW_cols)
    }
    #if (!is.null(iDWSQ_rast_input)) {
    #  iDWSQ_rasts <- terra::rast(iDWSQ_rast_input, iDWSQ_cols)
    #}

    all_rasts <- list(loi_rasts, iDW_rasts) #iDWSQ_rasts
    all_rasts <- all_rasts[!sapply(all_rasts, is.null)]
    all_rasts <- terra::rast(all_rasts)

    fun <- c()
    if (median) {
      fun <- c(fun, "median")
    }
    if (!is.null(quantiles)) {
      fun <- c(fun, "quantile")
    }

    if (length(fun) > 0) {
      extra_out <- ihydro:::.extract_fun(
        all_rasts = all_rasts,
        subbasins = subbasins,
        x_cols = numb_rast,
        weight_cols = NULL,
        fun = fun,
        quantiles = quantiles,
        max_cells_in_memory = max_cells_in_memory,
        n_cores = n_cores,
        include_count = include_count
      )
    } else {
      iDW_out <- ihydro:::.extract_fun(
        all_rasts = all_rasts,
        subbasins = subbasins,
        x_cols = iDW_cols,
        weight_cols = NULL,
        default_value = 0,
        default_weight = NA_real_,
        fun = "sum",
        max_cells_in_memory = max_cells_in_memory,
        n_cores = n_cores,
        include_count = include_count
      )

      if (loi_summary) {
        if (length(numb_rast) > 0) {
          numb_out <- ihydro:::.extract_fun(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = c(numb_rast),
            weight_cols = NULL,
            default_value = NA_real_,
            default_weight = NA_real_,
            fun = c("sum", "mean", "stdev", "variance", "count", "min", "max"),
            max_cells_in_memory = max_cells_in_memory,
            n_cores = n_cores,
            include_count = FALSE
          )
        }

        if (length(cat_rast) > 0) {
          cat_out <- ihydro:::.extract_fun(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = c(cat_rast),
            weight_cols = NULL,
            default_value = 0,
            default_weight = NA_real_,
            fun = c("sum"),
            max_cells_in_memory = max_cells_in_memory,
            n_cores = n_cores,
            include_count = FALSE
          )
        }
      }

      if (length(iDW_cols) > 0) {
        if (length(numb_rast) > 0) {
          numb_out_iDW <- ihydro:::.extract_fun(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = c(numb_rast),
            weight_cols = iDW_cols,
            #sqweight_cols = iDWSQ_cols,
            default_value = NA_real_,
            default_weight = 0,
            fun = c("weighted_sum", "weighted_mean", "weighted_variance"),
            max_cells_in_memory = max_cells_in_memory,
            n_cores = n_cores,
            include_count = FALSE
          )
        }
        if (length(cat_rast) > 0) {
          cat_out_iDW <- ihydro:::.extract_fun(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = c(cat_rast),
            weight_cols = iDW_cols,
            #sqweight_cols = iDWSQ_cols,
            default_value = 0,
            default_weight = 0,
            fun = c("weighted_sum"),
            max_cells_in_memory = max_cells_in_memory,
            n_cores = n_cores,
            include_count = FALSE
          )
        }
      }
    }

    if (catch_source == "Catchment_poly") {
      link_id_otarget <- link_id
      link_id <- NA_character_
    }

    if (!is.null(progressor)) {
      progressor()
    }

    dplyr::bind_cols(
      tibble::tibble(
        link_id_otarget = link_id_otarget,
        subbasin_link_id = link_id
      ),
      iDW_out,
      #iDWSQ_out,
      numb_out,
      cat_out,
      numb_out_iDW,
      cat_out_iDW,
      extra_out
    )
  }
#)

#' Create subbasin extraction plan
#' @noRd
.extract_planner <- function(
    input,
    subb_ids,
    loi_rast_input,
    numb_rast,
    cat_rast,
    iDW_rast_input,
    iDW_cols,
    max_cells_in_memory = 3e+07,
    n_cores = 1L,
    chunks_per_worker = 1L,
    fun = c(
      "mean",
      "sd",
      "var",
      "cv",
      "min",
      "max",
      "sum",
      "count",
      "median",
      "quantile"
    ),
    quantiles = NULL,
    temp_dir = NULL,
    verbose = FALSE
) {

  fun <- match.arg(fun,several.ok = TRUE)

  safe_kmeans <- function(x, centers, ...) {
    x0 <- x

    x <- dplyr::distinct(x[, c("link_id", "cent_x", "cent_y")])
    if (nrow(x) <= centers) {
      return(
        list(
          cluster = seq_len(nrow(x))
        )
      )
    }
    x$clust <- unlist(kmeans(x[, c("cent_x", "cent_y")], centers, ...)$cluster)
    x <- dplyr::select(x, link_id, clust)
    out <- dplyr::left_join(
      x0,
      x,
      by = "link_id"
    )

    return(
      list(
        cluster = out$clust
      )
    )
  }

  max_cells_in_memory <- floor(max_cells_in_memory / n_cores)
  if (max_cells_in_memory >= .Machine$integer.max) {
    max_cells_in_memory <- .Machine$integer.max - 1
  }


  n_jobs <- n_cores * chunks_per_worker
  n_jobs <- min(n_jobs, length(unique(subb_ids$link_id)))

  ws_lumped <- "lumped" %in% iDW_cols
  ws_s <- iDW_cols[iDW_cols %in% c("iFLS", "HAiFLS")]
  ws_o <- iDW_cols[iDW_cols %in% c("iFLO", "HAiFLO")]

  subb_lookup <- ihydro::read_ihydro(
    input,
    "us_flowpaths"
  ) |>
    dplyr::select(
      final_link_id = pour_point_id,
      link_id = origin_link_id #subbasin ID
    ) |>
    dplyr::filter(
      final_link_id %in% subb_ids$link_id
    )

  subbasins <- sf::read_sf(
    input$outfile,
    query = build_sql_in(
      "Subbasins_poly",
      "link_id",
      as.numeric(unique(subb_lookup$link_id))
    )
  )

  subbasins$cent <- sf::st_centroid(subbasins$geom)
  subbasins$cent_x <- sf::st_coordinates(subbasins$cent)[, 1]
  subbasins$cent_y <- sf::st_coordinates(subbasins$cent)[, 2]

  tasks <- dplyr::select(
    tibble::as_tibble(subbasins),
    -geom,
    -cent
  )
  tasks$link_id <- as.character(tasks$link_id)

  tasks <- dplyr::left_join(
    tasks,
    dplyr::select(subb_ids, link_id), # unn_group
    by = "link_id"
  )

  tasks_s <- NULL
  tasks_o <- NULL
  tasks_lump <- NULL

  # develop ws_s tasks
  #count_with_o <- TRUE
  #if (length(ws_s) > 0L) {
  #count_with_o <- FALSE

  # group accordinging to location
  tasks_s_group <- safe_kmeans(
    tasks[, c("link_id", "cent_x", "cent_y")],
    n_jobs
  )$cluster

  tasks_s <- lapply(
    split(tasks$link_id, tasks_s_group),
    function(x) {
      list(
        input_file = input$outfile,
        link_id = x,
        loi_rast_input = loi_rast_input$outfile,
        loi_summary = TRUE,
        numb_rast = numb_rast,
        cat_rast = cat_rast,
        iDW_rast_input = iDW_rast_input$outfile,
        iDW_cols = c(ws_s, ws_o),
        max_cells_in_memory = max_cells_in_memory,
        n_cores = n_cores,
        include_count = TRUE #!count_with_o
      )
    }
  )
  #}

  # # develop ws_o tasks
  # if (length(ws_o) > 0L) {
  #   subb_lookup <- subb_lookup |>
  #     dplyr::left_join(
  #       tasks[, c("link_id", "cent_x", "cent_y")],
  #       by = c("link_id" = "link_id")
  #     ) |>
  #     dplyr::left_join(
  #       tasks[, c("link_id", "unn_group")],
  #       by = c("final_link_id" = "link_id")
  #     )
  #
  #   tasks_o <- split(subb_lookup, subb_lookup$unn_group)
  #
  #   # group accordinging to location by unn_group
  #   tasks_o <- lapply(tasks_o, function(x) {
  #     x$group <- safe_kmeans(
  #       x[, c("link_id", "cent_x", "cent_y")],
  #       n_jobs
  #     )$cluster
  #     x <- split(x, x$group)
  #     return(x)
  #   })
  #
  #   tasks_o <- unlist(tasks_o, recursive = F)
  #
  #   tasks_o <- lapply(
  #     tasks_o,
  #     function(x) {
  #       list(
  #         input_file = input$outfile,
  #         link_id = x$link_id,
  #         link_id_otarget = x$final_link_id,
  #         loi_rast_input = loi_rast_input$outfile,
  #         loi_summary = FALSE,
  #         numb_rast = numb_rast,
  #         cat_rast = cat_rast,
  #         iDW_rast_input = iDW_rast_input$outfile,
  #         iDW_cols = paste(ws_o, unique(x$unn_group), sep = "_unn_group"),
  #         max_cells_in_memory = max_cells_in_memory,
  #         n_cores = n_cores,
  #         include_count = count_with_o
  #       )
  #     }
  #   )
  # }

  # develop ws_lump tasks
  if ("median" %in% fun || !is.null(quantiles)) {
    tasks_lump <- dplyr::arrange(subb_ids, USChnLn_To)
    tasks_lump <- split(
      tasks_lump$link_id,
      rep(1:n_jobs, length.out = nrow(tasks_lump))
    )

    tasks_lump <- lapply(
      tasks_lump,
      function(x) {
        list(
          input_file = input$outfile,
          link_id = x,
          loi_rast_input = loi_rast_input$outfile,
          loi_summary = FALSE,
          numb_rast = numb_rast,
          cat_rast = NULL,
          iDW_rast_input = NULL,
          iDW_cols = NULL,
          catch_source = "Catchment_poly",
          median = "median" %in% fun,
          quantiles = quantiles,
          max_cells_in_memory = max_cells_in_memory,
          n_cores = n_cores,
          include_count = FALSE
        )
      }
    )
  }

  out <- list(ws_lump = tasks_lump, ws_s = tasks_s, ws_o = tasks_o)
  return(out)
}

#' Summarize subbasin extraction results at the catchment scale
#' @param extract_value List of tibbles from .extract_subbasins()
#' @param weighting_scheme Character vector, e.g. c("lumped", "iFLS", "HAiFLS", ...)
#' @param loi_numeric_stats Character vector, e.g. c("mean", "sd", "var", ...)
#' @param numeric_vars Character vector of numeric variable names (e.g. c("slope"))
#' @param cat_vars Character vector of categorical variable names (e.g. c("LC_1", ...))
#' @return Tibble with one row per catchment (link_id_otarget), columns for each summary
#' @noRd
.summarize_catchment <- function(
    extract_value,
    weighting_scheme,
    loi_numeric_stats,
    numeric_vars,
    cat_vars
) {
  # Helper: Pivot and clean
  clean_data <- function(df) {
    df[sapply(df, is.nan)] <- NA_real_
    df
  }

  # Summarize ws_lump (medians/quantiles only)
  summarize_ws_lump <- function(df) {
    res <- purrr::map(numeric_vars, function(var) {
      if ("median" %in% loi_numeric_stats) {
        setNames(
          list(df[[paste0("median.", var)]]),
          paste0(var, "_lumped_median")
        )
      } else {
        NULL
      }
    })
    tibble::as_tibble(purrr::compact(unlist(res, recursive = FALSE)))
  }

  # Summarize ws_s (S-targeted iDW and lumped stats)
  summarize_ws_s <- function(df) {
    res <- list()
    # Lumped stats
    if ("lumped" %in% weighting_scheme) {
      for (var in numeric_vars) {
        mean_col <- paste0("mean.", var)
        var_col <- paste0("variance.", var)
        n_col <- paste0("count.", var)
        cnt_col <- "count.count_internal"
        if (all(c(mean_col, var_col, n_col) %in% names(df))) {
          sub_mean <- df[[mean_col]]
          sub_var <- df[[var_col]]
          sub_n <- df[[n_col]]
          sub_cnt <- df[[cnt_col]]
          if ("mean" %in% loi_numeric_stats) {
            res[[paste0(var, "_lumped_mean")]] <- combine_lumped_mean(
              sub_mean,
              sub_n,
              sub_cnt
            )
          }
          if ("sd" %in% loi_numeric_stats) {
            res[[paste0(var, "_lumped_sd")]] <- combine_lumped_sample_sd(
              sub_mean,
              sub_var,
              sub_n,
              sub_cnt
            )
          }
          if ("var" %in% loi_numeric_stats) {
            res[[paste0(var, "_lumped_var")]] <- combine_lumped_sample_var(
              sub_mean,
              sub_var,
              sub_n,
              sub_cnt
            )
          }
          if ("min" %in% loi_numeric_stats) {
            res[[paste0(var, "_lumped_min")]] <- min(
              df[[paste0("min.", var)]],
              na.rm = TRUE
            )
          }
          if ("max" %in% loi_numeric_stats) {
            res[[paste0(var, "_lumped_max")]] <- max(
              df[[paste0("max.", var)]],
              na.rm = TRUE
            )
          }
          if ("sum" %in% loi_numeric_stats) {
            res[[paste0(var, "_lumped_sum")]] <- sum(
              df[[paste0("sum.", var)]],
              na.rm = TRUE
            )
          }
        }
      }
      for (cat in cat_vars) {
        sum_col <- paste0("sum.", cat)
        count_col <- "count.count_internal"
        if (sum_col %in% names(df) && count_col %in% names(df)) {
          res[[paste0(cat, "_lumped_prop")]] <- sum(
            df[[sum_col]],
            na.rm = TRUE
          ) /
            sum(df[[count_col]], na.rm = TRUE)
        }
      }
    }
    # S-targeted iDW
    for (ws in intersect(weighting_scheme, c("iFLS", "HAiFLS","iFLO", "HAiFLO"))) {
      for (var in numeric_vars) {
        mean_col <- paste0("weighted_mean.", var, ".", ws)
        var_col <- paste0("weighted_variance.", var, ".", ws)
        sum_col <- paste0("weighted_sum.", var, ".", ws)
        wt_col <- paste0("weighted_sum.", ws, ".", var)
        cnt_col <- "count.count_internal"
        if (all(c(mean_col, var_col, wt_col) %in% names(df))) {
          sub_mean <- df[[mean_col]]
          sub_var <- df[[var_col]]
          sub_wt <- df[[wt_col]]
          sub_cnt <- df[[cnt_col]]
          if ("mean" %in% loi_numeric_stats) {
            res[[paste0(var, "_", ws, "_mean")]] <- combine_weighted_mean(
              sub_mean,
              sub_wt,
              sub_cnt
            )
          }
          if ("sd" %in% loi_numeric_stats) {
            res[[paste0(var, "_", ws, "_sd")]] <- combine_weighted_sample_sd(
              sub_mean,
              sub_var,
              sub_wt,
              sub_cnt
            )
          }
          if ("var" %in% loi_numeric_stats) {
            res[[paste0(var, "_", ws, "_var")]] <- combine_weighted_sample_var(
              sub_mean,
              sub_var,
              sub_wt,
              sub_cnt
            )
          }
          if ("sum" %in% loi_numeric_stats && sum_col %in% names(df)) {
            res[[paste0(var, "_", ws, "_sum")]] <- sum(
              df[[sum_col]],
              na.rm = TRUE
            )
          }
        }
      }
      for (cat in cat_vars) {
        sum_col <- paste0("weighted_sum.", cat, ".", ws)
        wt_col <- paste0("sum.", ws)
        if (sum_col %in% names(df) && wt_col %in% names(df)) {
          res[[paste0(cat, "_", ws, "_prop")]] <- sum(
            df[[sum_col]],
            na.rm = TRUE
          ) /
            sum(df[[wt_col]], na.rm = TRUE)
        }
      }
    }
    tibble::as_tibble(res)
  }

  # # Summarize ws_o (O-targeted iDW only)
  # summarize_ws_o <- function(df) {
  #   gp <- unique(df$unn_group)
  #   gp <- gp[!is.na(gp)]
  #   gp <- paste0("_unn_group", gp)
  #   res <- list()
  #   for (ws in intersect(weighting_scheme, c("iFLO", "HAiFLO"))) {
  #     for (var in numeric_vars) {
  #       mean_col <- paste0("weighted_mean.", var, ".", ws, gp)
  #       var_col <- paste0("weighted_variance.", var, ".", ws, gp)
  #       sum_col <- paste0("weighted_sum.", var, ".", ws, gp)
  #       wt_col <- paste0("weighted_sum.", ws, gp, ".", var)
  #       cnt_col <- "count.count_internal"
  #       if (all(c(mean_col, var_col, wt_col) %in% names(df))) {
  #         sub_mean <- df[[mean_col]]
  #         sub_var <- df[[var_col]]
  #         sub_wt <- df[[wt_col]]
  #         sub_cnt <- df[[cnt_col]]
  #         if ("mean" %in% loi_numeric_stats) {
  #           res[[paste0(var, "_", ws, "_mean")]] <- combine_weighted_mean(
  #             sub_mean,
  #             sub_wt,
  #             sub_cnt
  #           )
  #         }
  #         if ("sd" %in% loi_numeric_stats) {
  #           res[[paste0(var, "_", ws, "_sd")]] <- combine_weighted_sample_sd(
  #             sub_mean,
  #             sub_var,
  #             sub_wt,
  #             sub_cnt
  #           )
  #         }
  #         if ("var" %in% loi_numeric_stats) {
  #           res[[paste0(var, "_", ws, "_var")]] <- combine_weighted_sample_var(
  #             sub_mean,
  #             sub_var,
  #             sub_wt,
  #             sub_cnt
  #           )
  #         }
  #         if ("sum" %in% loi_numeric_stats && sum_col %in% names(df)) {
  #           res[[paste0(var, "_", ws, "_sum")]] <- sum(
  #             df[[sum_col]],
  #             na.rm = TRUE
  #           )
  #         }
  #       }
  #     }
  #     for (cat in cat_vars) {
  #       sum_col <- paste0("weighted_sum.", cat, ".", ws, gp)
  #       wt_col <- paste0("sum.", ws, gp)
  #       if (sum_col %in% names(df) && wt_col %in% names(df)) {
  #         res[[paste0(cat, "_", ws, "_prop")]] <- sum(
  #           df[[sum_col]],
  #           na.rm = TRUE
  #         ) /
  #           sum(df[[wt_col]], na.rm = TRUE)
  #       }
  #     }
  #   }
  #   tibble::as_tibble(res)
  # }

  # Main logic
  result <- list()
  if (!is.null(extract_value$ws_lump)) {
    ws_lump <- clean_data(extract_value$ws_lump)
    result$ws_lump <- ws_lump |>
      dplyr::group_by(link_id_otarget) |>
      tidyr::nest() |>
      dplyr::mutate(summary = purrr::map(data, summarize_ws_lump))
  }
  if (!is.null(extract_value$ws_s)) {
    ws_s <- clean_data(extract_value$ws_s)
    result$ws_s <- ws_s |>
      dplyr::group_by(link_id_otarget) |>
      tidyr::nest() |>
      dplyr::mutate(summary = purrr::map(data, summarize_ws_s))
  }
  # if (!is.null(extract_value$ws_o)) {
  #   ws_o <- clean_data(extract_value$ws_o)
  #   result$ws_o <- ws_o |>
  #     dplyr::bind_rows(
  #       dplyr::select(ws_s, count.count_internal)
  #     ) |>
  #     dplyr::group_by(link_id_otarget) |>
  #     tidyr::nest() |>
  #     dplyr::mutate(summary = purrr::map(data, summarize_ws_o))
  # }

  # Unnest and join all summaries
  result <- purrr::map(
    result,
    ~ dplyr::select(.x, -data) |> tidyr::unnest(summary)
  )
  result <- result[sapply(result, nrow) != 0]
  out <- purrr::reduce(result, dplyr::left_join, by = "link_id_otarget")
  out <- dplyr::ungroup(out)
  dplyr::rename(out, link_id = link_id_otarget)
}
