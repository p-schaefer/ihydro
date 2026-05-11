#' Attribute sampling points with distance-weighted landscape summaries (v2)
#'
#' Optimised version of fasttrib_points with two changes vs v1:
#'   1. Fewer exact_extract() calls in .extract_subbasins_v2(): numb + cat LOI
#'      rasters are stacked into a single call (unweighted and weighted), so
#'      coverage fractions are computed only twice instead of four times.
#'   2. Layer-duplication loops in .extract_fun_v2() are replaced with
#'      terra::subset() + vectorised names<- / terra::varnames<-.
#'
#' All other logic (including .summarize_catchment, .extract_planner) is
#' identical to the original.
#'
#' @inheritParams fasttrib_points
#' @inherit fasttrib_points return
#' @export
#' @rdname fasttrib_points_v2
fasttrib_points_v2 <- function(
  input,
  out_filename = NULL,
  loi_file = NULL,
  loi_cols = NULL,
  iDW_file = NULL,
  sample_points = NULL,
  link_id = NULL,
  weighting_scheme = c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO"),
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

  if (any(weighting_scheme != "lumped")) {
    if (is.null(iDW_file)) {
      temp_dir_sub <- file.path(temp_dir, basename(tempfile()))
      dir.create(temp_dir_sub, showWarnings = F, recursive = T)
      on.exit(unlink(temp_dir_sub, recursive = T, force = T), add = TRUE)
      args <- list(...)
      inv_function <- args$inv_function
      if (is.null(inv_function)) {
        inv_function <- function(x) (x * 0.001 + 1)^-1
      } else {
        inv_function <- eval(parse(
          text = paste0(deparse(inv_function), collapse = "")
        ))
      }

      fun_dep <- paste0(deparse(inv_function), collapse = "")

      if (verbose) {
        cli::cli_alert_info(c(
          "{.arg iDW_file} is Null. Calculating temporary iDW layers using: {.var {fun_dep}}.\n",
          "To use a different inverse distance function specify {.arg inv_function} in {.fun fasttrib_points_v2}\n",
          "or generate a permanent iDW file with {.fun prep_weights}"
        ))
      }

      iDW_file <- file.path(temp_dir, "temp_idw.gpkg")
      iDW_file <- prep_weights(
        input = input,
        output_filename = iDW_file,
        weighting_scheme = weighting_scheme[weighting_scheme != "lumped"],
        temp_dir = temp_dir_sub,
        verbose = verbose,
        inv_function = inv_function
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
  missing_layers <- missing_layers[!missing_layers %in% c("iFLO", "HAiFLO")]
  if (length(missing_layers) > 0L) {
    cli::cli_abort(c(
      "The following weighting layers are missing from iDW_file: ",
      "x" = "{.val {missing_layers}}"
    ))
  }

  # ─ Build subbasin lookup ────────────────────────

  unnest_catchment <- read_ihydro(input, "unnest_catchment")
  site_id_col <- read_ihydro(input, "site_id_col")[[1]]

  subb_ids <- read_ihydro(input, "stream_links_attr") |>
    dplyr::filter(link_lngth > 0) |>
    dplyr::filter(link_id %in% target_ids$link_id) |>
    dplyr::select(link_id, tidyselect::any_of(site_id_col), USChnLn_To) |>
    dplyr::left_join(
      dplyr::mutate(unnest_catchment, link_id = as.character(link_id)),
      by = c("link_id")
    )

  subb_ids <- subb_ids[!is.na(subb_ids$unn_group), ]

  subb_lookup <- ihydro::read_ihydro(
    input,
    "us_flowpaths"
  ) |>
    dplyr::select(
      final_link_id = pour_point_id,
      link_id = origin_link_id
    ) |>
    dplyr::filter(
      final_link_id %in% subb_ids$link_id
    )

  catch <- get_catchment(
    input = input,
    link_id = subb_ids$link_id,
    temp_dir = temp_dir,
    verbose = verbose,
    return = FALSE
  )

  numb_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "num_rast"]
  cat_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "cat_rast"]

  fun_sel <- unique(c("sum", loi_numeric_stats))

  if (is.null(n_batches)) {
    n_batches <- length(unique(subb_lookup$link_id))
  } else {
    stopifnot(is.numeric(n_batches))
    n_batches <- round(n_batches)
  }

  sub_summ <- .extract_planner(
    input = input,
    subb_ids = subb_ids,
    loi_rast_input = loi_file,
    numb_rast = numb_rast,
    cat_rast = cat_rast,
    iDW_rast_input = iDW_file,
    iDW_cols = weighting_scheme,
    mem_fraction = mem_fraction,
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

  progressr::handlers(progressr::handler_cli(
    format = "{cli::pb_bar} {cli::pb_percent} | {cli::pb_eta_str}"
  ))
  progressr::with_progress(enable = verbose, {
    total_tasks <- length(sub_summ)
    p <- progressr::progressor(steps = total_tasks)

    sub_summ <- lapply(sub_summ, function(args) {
      args$progressor <- p
      args
    })

    extract_execute <- lapply(sub_summ, function(args) {
      future::futureCall(
        .extract_subbasins_v2, # <-- v2 worker
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

  extract_value_comb <- list()
  extract_value_nms <- names(extract_value)
  ws_s_nms <- extract_value_nms[grepl("^ws_s\\.", extract_value_nms)]
  ws_lump_nms <- extract_value_nms[grepl("^ws_lump\\.", extract_value_nms)]

  extract_value_comb$ws_s <- tibble::tibble()
  extract_value_comb$ws_o <- tibble::tibble()
  extract_value_comb$ws_lump <- dplyr::bind_rows(extract_value[ws_lump_nms])

  if (length(ws_s_nms) > 0) {
    temp_comb <- list()
    ws_comb <- dplyr::bind_rows(extract_value[ws_s_nms])

    if (any(grepl("iFLO", weighting_scheme))) {
      o_targ <- dplyr::select(
        ws_comb,
        tidyselect::contains(c("iFLO", "HAiFLO")),
        subbasin_link_id,
        -link_id_otarget
      )

      wo_comb <- tidyr::pivot_longer(
        o_targ,
        tidyselect::contains(c("iFLO", "HAiFLO"))
      ) |>
        dplyr::mutate(
          type = dplyr::case_when(
            grepl("HAiFLO", name) ~ "HAiFLO",
            T ~ "iFLO"
          ),
          unn_group = strsplit(name, "iFLO_unn_group"),
          unn_group = purrr::map_chr(unn_group, ~ .x[[2]])
        ) |>
        dplyr::mutate(
          name2 = gsub("_unn_group.*\\.", ".", name),
          unn_group = gsub("\\..*$", "", unn_group),
          name = dplyr::case_when(
            name2 == name ~ gsub("_unn_group.*", "", name),
            T ~ name2
          )
        ) |>
        dplyr::select(unn_group, tidyselect::everything(), -name2, -type) |>
        tidyr::pivot_wider()

      unn_group_lookup <- dplyr::rename(
        subb_lookup,
        subbasin_link_id = link_id,
        link_id_otarget = final_link_id
      )

      unn_group_lookup <- dplyr::left_join(
        unn_group_lookup,
        dplyr::select(subb_ids, link_id_otarget = link_id, unn_group),
        by = "link_id_otarget"
      )

      temp_comb$ws_o <- dplyr::left_join(
        unn_group_lookup,
        wo_comb,
        by = c("subbasin_link_id", "unn_group")
      )
    }

    if (any(grepl("iFLS", weighting_scheme))) {
      temp_comb$ws_s <- dplyr::left_join(
        dplyr::rename(
          subb_lookup,
          subbasin_link_id = link_id,
          link_id_otarget = final_link_id
        ),
        dplyr::select(
          ws_comb,
          -tidyselect::contains(c(".iFLO", ".HAiFLO")),
          subbasin_link_id,
          -link_id_otarget
        ),
        by = "subbasin_link_id"
      )
    }

    extract_value_comb$ws_s <- purrr::reduce(
      temp_comb,
      dplyr::left_join,
      by = c("link_id_otarget", "subbasin_link_id")
    )
  }

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

  if (!is.null(out_filename)) {
    write.csv(result, out_filename)
  }

  return(result)
}


# ─────────────────────────────────────────────────────────────────────────────
# v2 internal helpers
# ─────────────────────────────────────────────────────────────────────────────

#' Extract raster attributes (v2) — vectorised layer duplication
#'
#' Identical contract to .extract_fun() but replaces the for-loops that build
#' `all_rast2` with terra::subset() + vectorised names<- / varnames<-.
#'
#' @noRd
.extract_fun_v2 <- function(
  all_rasts,
  subbasins,
  x_cols = NULL,
  weight_cols = NULL,
  default_value = NA_real_,
  default_weight = NA_real_,
  fun,
  quantiles = NULL,
  mem_fraction = 0.5,
  include_count = FALSE,
  n_cores = 1L,
  temp_dir_sub = NULL
) {
  if (is.null(x_cols) && is.null(weight_cols)) {
    return(NULL)
  }

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
  on.exit(
    terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE),
    add = TRUE
  )
  on.exit(ihydro:::restore_terra_options(old_terra_opts), add = TRUE)
  on.exit(
    suppressWarnings(unlink(temp_dir_sub, recursive = TRUE, force = TRUE)),
    add = TRUE
  )
  on.exit(gc(verbose = FALSE), add = TRUE)

  max_cells_in_memory <- ihydro:::.max_cells_in_memory_helper(
    mem_fraction = mem_fraction,
    n_cores = n_cores
  )

  if (include_count) {
    all_rasts[["count_internal"]] <- all_rasts[[1]]
    all_rasts[["count_internal"]][] <- 1
    x_cols <- c(x_cols, "count_internal")
    fun <- unique(c(fun, "count"))
  }

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
    # ── CHANGE 2a: replace for-loops with terra::subset() + vectorised rename ─

    # LOI layers duplicated: one copy per weight layer
    loi_idx <- rep(x_cols, each = length(weight_cols))
    loi_idxnm <- paste0(
      loi_idx,
      "_",
      rep(weight_cols, length.out = length(loi_idx))
    )

    all_rast2_loi <- terra::subset(all_rasts, loi_idx)
    names(all_rast2_loi) <- loi_idxnm

    # iDW layers duplicated: one copy per LOI layer
    iDW_idx <- rep(
      weight_cols,
      length.out = length(x_cols) * length(weight_cols)
    )
    iDW_idxnm <- paste0(iDW_idx, "_", rep(x_cols, each = length(weight_cols)))

    all_rast2_idw <- terra::subset(all_rasts, iDW_idx)
    names(all_rast2_idw) <- iDW_idxnm

    all_rast2 <- c(all_rast2_loi, all_rast2_idw)

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

    # rename columns back to canonical names
    final_names <- colnames(out[[1]])
    for (i in seq_along(loi_idxnm)) {
      final_names <- gsub(loi_idxnm[[i]], loi_idx[[i]], final_names)
    }
    for (i in seq_along(iDW_idxnm)) {
      final_names <- gsub(iDW_idxnm[[i]], iDW_idx[[i]], final_names)
    }
    colnames(out[[1]]) <- final_names

    # ── CHANGE 2b: weight-sum call — same refactor ────────────────────────────
    # Clamp LOI to 1 so weighted_sum gives sum-of-weights (accounting for NAs)
    all_rasts[[x_cols]] <- terra::clamp(all_rasts[[x_cols]], 1, 1)

    idw_comb <- c(weight_cols)
    iDW_idx2 <- rep(idw_comb, each = length(x_cols))
    iDW_idxnm2 <- paste0(
      iDW_idx2,
      "_",
      rep(x_cols, length.out = length(iDW_idx2))
    )

    loi_idx2 <- rep(x_cols, length.out = length(iDW_idx2))
    loi_idxnm2 <- paste0(loi_idx2, "_", rep(idw_comb, each = length(x_cols)))

    all_rast2b_idw <- terra::subset(all_rasts, iDW_idx2)
    names(all_rast2b_idw) <- iDW_idxnm2

    all_rast2b_loi <- terra::subset(all_rasts, loi_idx2)
    names(all_rast2b_loi) <- loi_idxnm2

    all_rast2b <- c(all_rast2b_idw, all_rast2b_loi)

    out[[2]] <- exactextractr::exact_extract(
      x = all_rast2b[[iDW_idxnm2]],
      y = subbasins,
      fun = "weighted_sum",
      weights = all_rast2b[[loi_idxnm2]],
      default_value = NA_real_,
      default_weight = 0,
      max_cells_in_memory = max_cells_in_memory,
      progress = FALSE,
      force_df = TRUE,
      full_colnames = TRUE
    )

    final_names <- colnames(out[[2]])
    for (i in seq_along(loi_idxnm2)) {
      final_names <- gsub(loi_idxnm2[[i]], loi_idx2[[i]], final_names)
    }
    for (i in seq_along(iDW_idxnm2)) {
      final_names <- gsub(iDW_idxnm2[[i]], iDW_idx2[[i]], final_names)
    }
    colnames(out[[2]]) <- final_names
  }

  dplyr::bind_cols(out)
}


#' Extract raster attributes for a set of subbasins (v2 — fewer exact_extract calls)
#'
#' Change 1: unweighted numb + cat LOI stacked into one exact_extract() call.
#'           weighted  numb + cat LOI stacked into one exact_extract() call.
#' Change 2: delegates to .extract_fun_v2() which uses terra::subset() instead
#'           of for-loops for layer duplication.
#'
#' @noRd
.extract_subbasins_v2 <- carrier::crate(
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
    catch_source = c("Subbasins_poly", "Catchment_poly"),
    median = FALSE,
    quantiles = NULL,
    mem_fraction = 0.5,
    n_cores = 1L,
    include_count = FALSE,
    progressor = NULL
  ) {
    catch_source <- match.arg(catch_source)

    subbasins <- sf::read_sf(
      input_file,
      query = ihydro:::build_sql_in(catch_source, "link_id", unique(link_id))
    )

    loi_rasts <- terra::rast(loi_rast_input, c(numb_rast, cat_rast))
    loi_rasts <- terra::crop(loi_rasts, terra::vect(subbasins))

    iDW_rasts <- NULL
    if (!is.null(iDW_rast_input)) {
      iDW_rasts <- lapply(iDW_cols, function(x) terra::rast(iDW_rast_input, x))
      iDW_rasts <- terra::rast(iDW_rasts)
      iDW_rasts <- terra::crop(iDW_rasts, terra::vect(subbasins))
    }

    all_rasts_list <- list(loi_rasts, iDW_rasts)
    all_rasts_list <- all_rasts_list[!sapply(all_rasts_list, is.null)]
    all_rasts <- terra::rast(all_rasts_list)

    iDW_out <- NULL
    loi_out <- NULL # replaces separate numb_out + cat_out
    loi_out_iDW <- NULL # replaces separate numb_out_iDW + cat_out_iDW
    extra_out <- NULL

    fun <- c()
    if (median) {
      fun <- c(fun, "median")
    }
    if (!is.null(quantiles)) {
      fun <- c(fun, "quantile")
    }

    if (length(fun) > 0) {
      extra_out <- ihydro:::.extract_fun_v2(
        all_rasts = all_rasts,
        subbasins = subbasins,
        x_cols = numb_rast,
        weight_cols = NULL,
        fun = fun,
        quantiles = quantiles,
        mem_fraction = mem_fraction,
        n_cores = n_cores,
        include_count = include_count
      )
    } else {
      # iDW weight sums (must stay separate — uses default_value = 0)
      iDW_out <- ihydro:::.extract_fun_v2(
        all_rasts = all_rasts,
        subbasins = subbasins,
        x_cols = iDW_cols,
        weight_cols = NULL,
        default_value = 0,
        default_weight = NA_real_,
        fun = "sum",
        mem_fraction = mem_fraction,
        n_cores = n_cores,
        include_count = include_count
      )

      # ── CHANGE 1a: single unweighted call for ALL LOI (numb + cat) ──────────
      all_loi <- c(numb_rast, cat_rast)
      if (loi_summary && length(all_loi) > 0) {
        # Use numb stats for all bands; cat columns will only use "sum" downstream
        # but computing the full set is cheap and keeps column names consistent.
        loi_out <- ihydro:::.extract_fun_v2(
          all_rasts = all_rasts,
          subbasins = subbasins,
          x_cols = all_loi,
          weight_cols = NULL,
          default_value = NA_real_, # numb default
          default_weight = NA_real_,
          fun = c("sum", "mean", "stdev", "variance", "count", "min", "max"),
          mem_fraction = mem_fraction,
          n_cores = n_cores,
          include_count = FALSE
        )
        # cat columns need default_value = 0; override the sum columns only
        if (length(cat_rast) > 0) {
          cat_out_only <- ihydro:::.extract_fun_v2(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = cat_rast,
            weight_cols = NULL,
            default_value = 0,
            default_weight = NA_real_,
            fun = "sum",
            mem_fraction = mem_fraction,
            n_cores = n_cores,
            include_count = FALSE
          )
          # Replace sum columns for cat variables with the correct default_value=0 version
          cat_sum_cols <- paste0("sum.", cat_rast)
          loi_out[cat_sum_cols] <- cat_out_only[cat_sum_cols]
        }
      }

      # ── CHANGE 1b: single weighted call for ALL LOI (numb + cat) ────────────
      if (length(iDW_cols) > 0 && length(all_loi) > 0) {
        loi_out_iDW <- ihydro:::.extract_fun_v2(
          all_rasts = all_rasts,
          subbasins = subbasins,
          x_cols = all_loi,
          weight_cols = iDW_cols,
          default_value = NA_real_,
          default_weight = 0,
          fun = c("weighted_sum", "weighted_mean", "weighted_variance"),
          mem_fraction = mem_fraction,
          n_cores = n_cores,
          include_count = FALSE
        )
        # cat weighted sums should use default_value = 0
        if (length(cat_rast) > 0) {
          cat_wt_only <- ihydro:::.extract_fun_v2(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = cat_rast,
            weight_cols = iDW_cols,
            default_value = 0,
            default_weight = 0,
            fun = "weighted_sum",
            mem_fraction = mem_fraction,
            n_cores = n_cores,
            include_count = FALSE
          )
          wt_sum_cols <- grep(
            paste0(
              "^weighted_sum\\.(",
              paste(cat_rast, collapse = "|"),
              ")\\."
            ),
            names(loi_out_iDW),
            value = TRUE
          )
          wt_sum_cols_new <- grep(
            paste0(
              "^weighted_sum\\.(",
              paste(cat_rast, collapse = "|"),
              ")\\."
            ),
            names(cat_wt_only),
            value = TRUE
          )
          loi_out_iDW[wt_sum_cols] <- cat_wt_only[wt_sum_cols_new]
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
      loi_out,
      loi_out_iDW,
      extra_out
    )
  }
)
