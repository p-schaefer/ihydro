#' Attribute stream segments/sampling points with layers of interest
#'
#' Uses [hydroweight()] and [hydroweight_attributes()]
#' to attribute points with inverse-distance weighted landscape statistics.
#' Output columns are named to match [fasttrib_points()] conventions.
#'
#' @param input An `ihydro` object.
#' @param out_filename CSV output path.
#' @param loi_file Optional path to LOI GeoPackage or `ihydro` object.
#' @param loi_cols Character vector of LOI column names, or `NULL` for all.
#' @param sample_points Character vector of site IDs, or `NULL`.
#' @param link_id Character vector of link IDs, or `NULL`.
#' @param clip_region Optional region to exclude when summarising.
#' @param OS_combine Logical. Merge target_O and target_S?
#' @param target_o_type Character. Target geometry type.
#' @param weighting_scheme Character vector of weighting schemes.
#' @param loi_numeric_stats Character vector of summary statistics.
#' @param inv_function Inverse distance function.
#' @param temp_dir Temporary file directory.
#' @param return_products Logical. Return products?
#' @param verbose Logical.
#'
#' @return A data.frame of resulting attributes, also written to
#'   `out_filename`. Output column names match [fasttrib_points()]
#'   conventions: `pour_point_id` (lowercase), `status`,
#'   `{loi}_lumped_{stat}`, `{loi}_{scheme}_mean`,
#'   `{loi}_{scheme}_sd`, `{loi}_{scheme}_prop`.
#' @export
#'
#' @seealso [fasttrib_points()], [hydroweight()],
#'   [hydroweight_attributes()]
#'
#' @rdname attrib_points
attrib_points <- function(
  input,
  out_filename,
  loi_file = NULL,
  loi_cols = NULL,
  sample_points = NULL,
  link_id = NULL,
  clip_region = NULL,
  OS_combine = FALSE,
  target_o_type = c(
    "point",
    "segment_point",
    "segment_whole"
  ),
  weighting_scheme = c(
    "lumped",
    "iEucO",
    "iEucS",
    "iFLO",
    "iFLS",
    "HAiFLO",
    "HAiFLS"
  ),
  loi_numeric_stats = c(
    "distwtd_mean",
    "distwtd_sd",
    "mean",
    "sd",
    "median",
    "min",
    "max",
    "sum",
    "cell_count"
  ),
  inv_function = function(x) (x * 0.001 + 1)^-1,
  temp_dir = NULL,
  return_products = TRUE,
  verbose = FALSE
) {
  # -- Validate inputs ---
  check_ihydro(input)

  temp_dir <- ensure_temp_dir(temp_dir)
  target_o_type <- match.arg(target_o_type)
  weighting_scheme <- match.arg(
    weighting_scheme,
    several.ok = TRUE
  )
  loi_numeric_stats <- match.arg(
    loi_numeric_stats,
    several.ok = TRUE
  )
  loi_numeric_stats <- stats::setNames(
    loi_numeric_stats,
    loi_numeric_stats
  )

  db_fp <- input$outfile
  n_cores <- n_workers()

  # -- Configure external tools ---
  whitebox::wbt_options(
    exe_path = whitebox::wbt_exe_path(),
    verbose = verbose > 2,
    wd = temp_dir
  )

  old_terra_opts <- set_terra_options(
    n_cores = n_cores,
    temp_dir = temp_dir,
    verbose = verbose > 3
  )
  on.exit(restore_terra_options(old_terra_opts), add = TRUE)

  # -- Resolve LOI file ---
  loi_file <- .ap_resolve_loi(input, loi_file)
  loi_meta <- read_ihydro(loi_file, "loi_meta")

  if (is.null(loi_cols)) {
    loi_cols <- loi_meta$loi_var_nms
  }
  bad_cols <- loi_cols[!loi_cols %in% loi_meta$loi_var_nms]
  if (length(bad_cols) > 0L) {
    cli::cli_abort("LOI columns not found: {.val {bad_cols}}")
  }
  loi_meta <- dplyr::filter(
    loi_meta,
    loi_var_nms %in% loi_cols
  )

  # -- Resolve targets ---
  target_ids <- target_id_fun(
    db_fp = db_fp,
    sample_points = sample_points,
    link_id = link_id,
    segment_whole = target_o_type == "segment_whole",
    target_o_type = target_o_type
  )
  target_o <- target_o_fun(
    db_fp = db_fp,
    target_IDs = target_ids,
    target_o_type = target_o_type
  )

  # -- Clip region ---
  clip_region <- .ap_resolve_clip(clip_region, temp_dir)

  # -- Write DEM rasters to temp for worker access ---
  temp_dir_sub <- file.path(temp_dir, basename(tempfile()))
  dir.create(temp_dir_sub, recursive = TRUE)
  on.exit(
    suppressWarnings(
      unlink(temp_dir_sub, recursive = TRUE, force = TRUE)
    ),
    add = TRUE
  )

  dem_layers <- c(
    "dem_streams_d8_sub",
    "dem_final",
    "dem_accum_d8",
    "dem_d8"
  )
  for (lyr in dem_layers) {
    terra::writeRaster(
      read_ihydro(input, lyr),
      file.path(temp_dir_sub, paste0(lyr, ".tif")),
      datatype = "FLT4S"
    )
  }

  # -- Compute attributes in parallel ---
  if (verbose) {
    message("Calculating attributes")
  }

  raw_out <- .ap_dispatch(
    target_o = target_o,
    input = input,
    loi_file = loi_file,
    weighting_scheme = weighting_scheme,
    loi_cols = loi_cols,
    loi_numeric_stats = loi_numeric_stats,
    loi_meta = loi_meta,
    temp_dir_sub = temp_dir_sub,
    return_products = return_products,
    inv_function = inv_function,
    clip_region = clip_region,
    OS_combine = OS_combine,
    n_cores = n_cores,
    verbose = verbose
  )

  # -- Rename columns to fasttrib_points convention ---
  raw_out <- .ap_rename_columns(raw_out, weighting_scheme)

  # -- Join back to original target IDs ---
  target_ids_out <- target_id_fun(
    db_fp = db_fp,
    sample_points = sample_points,
    link_id = link_id,
    segment_whole = FALSE,
    target_o_type = target_o_type
  )

  final_out <- dplyr::left_join(
    target_ids_out,
    dplyr::mutate(
      raw_out,
      pour_point_id = as.character(pour_point_id)
    ),
    by = c("link_id" = "pour_point_id"),
    multiple = "all"
  )

  # -- Write CSV (exclude list-column products) ---
  data.table::fwrite(
    x = dplyr::select(
      final_out,
      -tidyselect::any_of("products")
    ),
    file = out_filename,
    buffMB = 128L,
    nThread = 1L,
    showProgress = FALSE
  )

  final_out
}


# ===========================================================
# Internal helpers -- NOT exported
# ===========================================================

# -- LOI / clip_region resolution ---

#' Resolve LOI file to an ihydro object
#' @noRd
.ap_resolve_loi <- function(input, loi_file) {
  if (is.null(loi_file)) {
    lyrs <- ihydro_layers(input)
    if (!"loi_meta" %in% lyrs$layer_name) {
      cli::cli_abort("No LOI data found in {.arg input}.")
    }
    return(as.ihydro(input$outfile))
  }
  if (inherits(loi_file, "ihydro")) {
    lf <- loi_file
  } else {
    lf <- as.ihydro(loi_file)
  }
  if (!"loi_meta" %in% ihydro_layers(lf)$layer_name) {
    cli::cli_abort(
      "No LOI data found in {.arg loi_file}."
    )
  }
  lf
}

#' Resolve clip_region to a file path (or NULL)
#' @noRd
.ap_resolve_clip <- function(clip_region, temp_dir) {
  clip_region <- process_input(
    clip_region,
    input_name = "clip_region"
  )
  if (is.null(clip_region)) {
    return(NULL)
  }
  if (inherits(clip_region, "SpatRaster")) {
    fp <- file.path(temp_dir, "clip_region.tif")
    terra::writeRaster(clip_region, fp, datatype = "FLT4S")
    return(fp)
  }
  if (inherits(clip_region, "SpatVector")) {
    fp <- file.path(temp_dir, "clip_region.shp")
    terra::writeVector(clip_region, fp)
    return(fp)
  }
  clip_region
}


# -- Parallel dispatch ---

#' Dispatch target_o rows to workers one task at a time
#'
#' Each row of target_o is one task. `scheduling = 1L` ensures
#' idle workers pick up the next task immediately (no batch
#' pre-assignment). Failed tasks are retried sequentially.
#'
#' @noRd
.ap_dispatch <- function(
  target_o,
  input,
  loi_file,
  weighting_scheme,
  loi_cols,
  loi_numeric_stats,
  loi_meta,
  temp_dir_sub,
  return_products,
  inv_function,
  clip_region,
  OS_combine,
  n_cores,
  verbose
) {
  # carrier::crate captures only explicitly-listed values,
  # preventing the entire parent environment from being
  # serialised to workers.
  worker <- carrier::crate(
    function(
      geom,
      link_id_val,
      input,
      loi_file,
      weighting_scheme,
      loi_cols,
      loi_numeric_stats,
      loi_meta,
      temp_dir_sub,
      return_products,
      inv_function,
      clip_region,
      OS_combine,
      n_cores
    ) {
      ihydro:::.ap_process_one(
        geom = geom,
        link_id_val = link_id_val,
        input = input,
        loi_file = loi_file,
        weighting_scheme = weighting_scheme,
        loi_cols = loi_cols,
        loi_numeric_stats = loi_numeric_stats,
        loi_meta = loi_meta,
        temp_dir_sub = temp_dir_sub,
        return_products = return_products,
        inv_function = inv_function,
        clip_region = clip_region,
        OS_combine = OS_combine,
        n_cores = n_cores
      )
    }
  )

  # Split target_o into individual rows
  geom_list <- split(
    target_o$geom,
    seq_len(nrow(target_o))
  )
  link_id_list <- as.list(target_o$link_id)

  args <- list(
    geom = geom_list,
    link_id_val = link_id_list,
    input = list(input),
    loi_file = list(loi_file),
    weighting_scheme = list(weighting_scheme),
    loi_cols = list(loi_cols),
    loi_numeric_stats = list(loi_numeric_stats),
    loi_meta = list(loi_meta),
    temp_dir_sub = list(temp_dir_sub),
    return_products = list(return_products),
    inv_function = list(inv_function),
    clip_region = list(clip_region),
    OS_combine = list(OS_combine),
    n_cores = list(n_cores)
  )

  # First pass: parallel, scheduling = 1L feeds one task
  # at a time to whichever worker is free
  results <- furrr::future_pmap(
    args,
    worker,
    .options = furrr::furrr_options(
      globals = FALSE,
      seed = NULL,
      scheduling = 4L
    )
  )
  out <- dplyr::bind_rows(results)

  # Retry Incomplete tasks sequentially
  incomplete_ids <- out$pour_point_id[
    out$status == "Incomplete"
  ]

  if (length(incomplete_ids) > 0L) {
    if (verbose) {
      message(
        "Retrying ",
        length(incomplete_ids),
        " point(s) sequentially."
      )
    }
    retry_idx <- which(
      target_o$link_id %in% incomplete_ids
    )
    retry_args <- lapply(args, function(a) {
      if (length(a) == length(geom_list)) {
        a[retry_idx]
      } else {
        a
      }
    })
    retry_results <- purrr::pmap(retry_args, worker)
    out_retry <- dplyr::bind_rows(retry_results)
    out <- dplyr::bind_rows(
      dplyr::filter(
        out,
        !pour_point_id %in% incomplete_ids
      ),
      out_retry
    )
  }

  out
}


# -- Single-point worker ---

#' Process one target point: hydroweight then
#' hydroweight_attributes
#'
#' Called inside a future worker (via carrier::crate).
#' Self-contained: sets per-worker terra options and
#' cleans up on exit.
#'
#' @noRd
.ap_process_one <- function(
  geom,
  link_id_val,
  input,
  loi_file,
  weighting_scheme,
  loi_cols,
  loi_numeric_stats,
  loi_meta,
  temp_dir_sub,
  return_products,
  inv_function,
  clip_region,
  OS_combine,
  n_cores
) {
  # Per-worker terra memory management
  old_terra_opts <- ihydro:::set_terra_options(
    n_cores = n_cores,
    temp_dir = temp_dir_sub,
    verbose = FALSE
  )
  on.exit(ihydro:::restore_terra_options(old_terra_opts), add = TRUE)

  # Build target_O sf from the geometry
  y <- sf::st_as_sf(geom) |>
    dplyr::mutate(link_id = link_id_val)

  # Per-task temp directory for hydroweight output
  hw_dir <- file.path(
    temp_dir_sub,
    basename(tempfile())
  )
  dir.create(hw_dir, recursive = TRUE)
  on.exit(
    suppressWarnings(
      unlink(hw_dir, recursive = TRUE, force = TRUE)
    ),
    add = TRUE
  )

  # Load DEM rasters (proxies, no data in RAM)
  target_S <- terra::rast(
    file.path(temp_dir_sub, "dem_streams_d8_sub.tif")
  )
  dem <- terra::rast(
    file.path(temp_dir_sub, "dem_final.tif")
  )
  flow_accum <- terra::rast(
    file.path(temp_dir_sub, "dem_accum_d8.tif")
  )

  # Get catchment
  catch <- tryCatch(
    ihydro::get_catchment(
      input,
      link_id = link_id_val
    ),
    error = function(cnd) NULL
  )
  if (is.null(catch)) {
    return(tibble::tibble(
      pour_point_id = link_id_val,
      status = "Incomplete"
    ))
  }

  # Compute inverse-distance weights
  hw <- tryCatch(
    suppressMessages(
      ihydro::hydroweight(
        hydroweight_dir = hw_dir,
        target_S = target_S,
        target_O = y,
        target_uid = link_id_val,
        OS_combine = OS_combine,
        dem = dem,
        flow_accum = flow_accum,
        clip_region = catch,
        weighting_scheme = weighting_scheme,
        inv_function = inv_function,
        return_products = return_products,
        wrap_return_products = TRUE,
        save_output = TRUE,
        clean_tempfiles = TRUE
      )
    ),
    error = function(cnd) NULL
  )
  if (is.null(hw)) {
    return(tibble::tibble(
      pour_point_id = link_id_val,
      status = "Incomplete"
    ))
  }

  # Compute attributes for numeric and categorical LOI
  num_vars <- loi_meta$loi_var_nms[
    loi_meta$loi_type == "num_rast"
  ]
  cat_vars <- loi_meta$loi_var_nms[
    loi_meta$loi_type == "cat_rast"
  ]

  hw_attr_num <- .ap_safe_hw_attr(
    loi_file,
    num_vars,
    loi_numeric = TRUE,
    loi_numeric_stats,
    catch,
    link_id_val,
    hw_dir,
    clip_region,
    return_products
  )
  hw_attr_cat <- .ap_safe_hw_attr(
    loi_file,
    cat_vars,
    loi_numeric = FALSE,
    loi_numeric_stats,
    catch,
    link_id_val,
    hw_dir,
    clip_region,
    return_products
  )

  # If both failed, mark incomplete
  if (is.null(hw_attr_num) && is.null(hw_attr_cat)) {
    return(tibble::tibble(
      pour_point_id = link_id_val,
      status = "Incomplete"
    ))
  }

  # Assemble products (if requested)
  p_out <- NULL
  if (return_products) {
    p0 <- hw
    p1 <- if (!is.null(hw_attr_num)) {
      p <- unlist(
        hw_attr_num$return_products,
        recursive = FALSE
      )
      if (!is.null(p)) {
        names(p) <- paste0(names(p), "_num")
      }
      p
    }
    p2 <- if (!is.null(hw_attr_cat)) {
      p <- unlist(
        hw_attr_cat$return_products,
        recursive = FALSE
      )
      if (!is.null(p)) {
        names(p) <- paste0(names(p), "_cat")
      }
      p
    }
    p_out <- c(p0, p1, p2)
    if (!is.null(p_out)) {
      p_out <- p_out[sort(names(p_out))]
    }
  }

  # Join numeric and categorical attribute tables
  attr_tables <- Filter(
    Negate(is.null),
    list(
      hw_attr_num$attribute_table,
      hw_attr_cat$attribute_table
    )
  )
  attr_joined <- purrr::reduce(
    attr_tables,
    dplyr::left_join,
    by = "pour_point_id"
  )

  dplyr::bind_cols(
    tibble::tibble(
      products = list(p_out)[1],
      status = "Complete"
    ),
    attr_joined
  ) |>
    dplyr::select(
      pour_point_id,
      status,
      tidyselect::everything()
    )
}


#' Safely call hydroweight_attributes
#'
#' Returns NULL on error or when loi_vars is empty.
#' @noRd
.ap_safe_hw_attr <- function(
  loi_file,
  loi_vars,
  loi_numeric,
  loi_numeric_stats,
  catch,
  link_id_val,
  hw_dir,
  clip_region,
  return_products
) {
  if (length(loi_vars) == 0L) {
    return(NULL)
  }
  tryCatch(
    suppressMessages(
      ihydro::hydroweight_attributes(
        loi = terra::rast(
          loi_file$outfile,
          lyrs = loi_vars
        ),
        loi_columns = loi_vars,
        loi_numeric = loi_numeric,
        loi_numeric_stats = loi_numeric_stats,
        roi = catch,
        roi_uid = link_id_val,
        roi_uid_col = "pour_point_id",
        distance_weights = file.path(
          hw_dir,
          paste0(link_id_val, "_inv_distances.zip")
        ),
        remove_region = clip_region,
        return_products = return_products
      )
    ),
    error = function(cnd) NULL
  )
}


# -- Column renaming ---

#' Rename hydroweight_attributes output columns to match
#' fasttrib_points conventions
#'
#' hydroweight_attributes produces:
#'   {loi}_mean, {loi}_sd, ... (lumped stats)
#'   {loi}_{scheme}_distwtd_mean (weighted mean)
#'   {loi}_{scheme}_distwtd_sd (weighted sd)
#'   {loi}_{scheme}_prop (categorical)
#'
#' fasttrib_points expects:
#'   {loi}_lumped_mean, {loi}_lumped_sd, ... (lumped)
#'   {loi}_{scheme}_mean (weighted, no "distwtd_")
#'   {loi}_{scheme}_sd (weighted, no "distwtd_")
#'   {loi}_lumped_prop (lumped categorical)
#'   {loi}_{scheme}_prop (already correct)
#'
#' @noRd
.ap_rename_columns <- function(df, weighting_scheme) {
  nms <- colnames(df)
  scheme_rx <- paste(weighting_scheme, collapse = "|")

  # 1. "{loi}_{scheme}_distwtd_mean" -> "{loi}_{scheme}_mean"
  nms <- gsub("_distwtd_mean$", "_mean", nms)

  # 2. "{loi}_{scheme}_distwtd_sd" -> "{loi}_{scheme}_sd"
  nms <- gsub("_distwtd_sd$", "_sd", nms)

  # 3. Bare lumped stats: "{loi}_mean" -> "{loi}_lumped_mean"
  #    Only if NOT preceded by a scheme name.
  lumped_stats <- c(
    "mean",
    "sd",
    "median",
    "min",
    "max",
    "sum",
    "cell_count",
    "NA_cell_count"
  )
  scheme_pat <- paste0("_(", scheme_rx, ")_")
  for (stat in lumped_stats) {
    sfx <- paste0("_", stat, "$")
    has_sfx <- grepl(sfx, nms)
    has_sch <- grepl(
      paste0(scheme_pat, stat, "$"),
      nms
    )
    is_lumped <- has_sfx & !has_sch
    nms[is_lumped] <- gsub(
      paste0("_", stat, "$"),
      paste0("_lumped_", stat),
      nms[is_lumped]
    )
  }

  # 4. Bare "{cat}_prop" -> "{cat}_lumped_prop"
  has_prop <- grepl("_prop$", nms)
  has_sch_prop <- grepl(
    paste0(scheme_pat, "prop$"),
    nms
  )
  is_bare_prop <- has_prop & !has_sch_prop & !grepl("_lumped_prop$", nms)
  nms[is_bare_prop] <- gsub(
    "_prop$",
    "_lumped_prop",
    nms[is_bare_prop]
  )

  colnames(df) <- nms
  df
}
