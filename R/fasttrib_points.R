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
#' @param loi_numeric_stats Character vector. Any of `"mean"`, `"sd"`, `"variance"`
#'   `"median"`, `"min"`, `"max"`, `"sum"`.
#' @param mem_fraction Numeric. Value between 0.1 and 0.9 (larger values give a warning).
#'   The fraction of RAM that may be used by the program.
#' @param quantile quantiles to be computed (`NULL` to skip quantile calcularion)
#' @param temp_dir Temporary directory for intermediate files.
#' @param verbose Logical.
#'
#' @return A data.frame of weighted attributes, also written to `out_filename`.
#' @export
#'
#' @seealso [process_hydrology()], [process_loi()], [prep_weights()],
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
    weighting_scheme = c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO"),
    loi_numeric_stats = c(
      "mean",
      "sd",
      "variance",
      "median",
      "min",
      "max",
      "sum"
    ),
    quantiles = NULL,
    mem_fraction = 0.5,
    n_batches = NULL,
    temp_dir = NULL,
    verbose = FALSE,
    ...
) {
  # ─ Validate inputs ───────────────────────────
  .check_ihydro(input)

  if (!is.logical(verbose)) {
    cli::cli_abort("{.arg verbose} must be {.code TRUE} or {.code FALSE}.")
  }
  if (!(mem_fraction > 0.1 & mem_fraction < 0.9)) {
    cli::cli_abort(
      "{.arg mem_fraction} must be between 0.1 and 0.9, not {mem_fraction}."
    )
  }
  if (!all(quantiles >= 0 & quantiles <= 1)) {
    cli::cli_abort("All values in {.arg quantiles} must be between 0 and 1.")
  }

  loi_numeric_stats <- match.arg(
    loi_numeric_stats,
    several.ok = TRUE,
    choices = c("mean", "sd", "variance", "median", "min", "max", "sum")
  )
  loi_numeric_stats <- stats::setNames(loi_numeric_stats, loi_numeric_stats)

  input_path <- input$outfile
  temp_dir <- .ensure_temp_dir(temp_dir)
  n_cores <- .n_workers()

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
          "To use a different inverse distance function specify {.arg inv_function} in {.fun fasttrib_points}\n",
          "or generate a permanent iDW file with {.fun prep_weights}"
        ))
      }

      iDW_file <- file.path(temp_dir, "temp_idw.gpkg")
      iDW_file <- prep_weights(
        input = input,
        output_filename = iDW_file,
        weighting_scheme = weighting_scheme[weighting_scheme != "lumped"],
        sample_points = sample_points,
        link_id = link_id,
        temp_dir = temp_dir_sub,
        verbose = verbose,
        inv_function = inv_function
      )
    }

    .check_ihydro(iDW_file)
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

  old_terra_opts <- .set_terra_options(
    n_cores = 1L,
    temp_dir = temp_dir,
    verbose = verbose > 3
  )
  on.exit(.restore_terra_options(old_terra_opts), add = TRUE)

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
  target_ids <- .target_id_fun(
    db_fp = input_path,
    sample_points = sample_points,
    link_id = link_id
  )

  # ─ Build subbasin lookup ────────────────────────

  unnest_catchment <- read_ihydro(input, "unnest_catchment")
  site_id_col <- read_ihydro(input, "site_id_col")[[1]]

  subb_ids <- read_ihydro(input, "stream_links_attr") |>
    dplyr::filter(link_lngth > 0) |> # small catchments without streams
    dplyr::filter(link_id %in% target_ids$link_id) |> # small catchments without streams
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
      link_id = origin_link_id #subbasin ID
    ) |>
    dplyr::filter(
      final_link_id %in% subb_ids$link_id
    )

  # ─ Resolve iDW file and prepare weights ─────────────────
  idw_layers <- ihydro_layers(iDW_file)$layer_name
  idw_layers <- idw_layers[!idw_layers %in% c("weight_meta", "DEM_Extent")]
  missing_layers <- setdiff(weighting_scheme, idw_layers)
  missing_layers <- missing_layers[missing_layers != "lumped"]
  if (any(grepl("iFLO", weighting_scheme))) {
    layers_needed <- lapply(
      weighting_scheme[grepl("iFLO", weighting_scheme)],
      function(x) paste0(x, "_unn_group", unique(subb_ids$unn_group))
    )
    layers_needed <- unlist(layers_needed)

    missing_layers <- missing_layers[!missing_layers %in% c("iFLO", "HAiFLO")]
    missing_layers <- c(missing_layers, setdiff(layers_needed, idw_layers))
  }

  if (length(missing_layers) > 0L) {
    cli::cli_abort(c(
      "The following weighting layers are missing from iDW_file: ",
      "x" = "{.val {missing_layers}}"
    ))
  }

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
  fun_sel <- unique(c("sum", loi_numeric_stats))

  if (is.null(n_batches)) {
    n_batches <- ceiling(length(unique(subb_lookup$link_id)) / 4)
  } else {
    if (!is.numeric(n_batches)) {
      cli::cli_abort("{.arg n_batches} must be numeric.")
    }
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
    quantiles = quantiles,
    temp_dir = temp_dir,
    verbose = verbose
  )
  sub_summ <- unlist(sub_summ, recursive = FALSE)

  if (verbose) {
    cli::cli_alert_info("Extracting attributes for each subbasin.")
  }

  subbs <- read_ihydro(input, "Subbasins_poly")
  subbs$link_id <- as.numeric(subbs$link_id) * 1000
  subbs <- terra::rasterize(
    terra::vect(subbs),
    read_ihydro(input, "dem_final"),
    field = "link_id"
  )
  subbs <- terra::wrap(subbs)

  loi_rasts <- terra::rast(loi_file$outfile)

  iDW_rasts <- NULL
  iDW_cols <- unique(unlist(lapply(sub_summ, function(x) x$iDW_cols)))
  if (!is.null(iDW_file) & length(iDW_cols) > 0) {
    iDW_rasts <- lapply(iDW_cols, function(x) terra::rast(iDW_file$outfile, x))
    iDW_rasts <- terra::rast(iDW_rasts)
  }

  all_rasts <- list(loi_rasts, iDW_rasts)
  all_rasts <- all_rasts[!sapply(all_rasts, is.null)]
  all_rasts <- terra::rast(all_rasts)
  all_rasts <- terra::wrap(all_rasts)

  gg <- gc(verbose = FALSE)

  progressr::handlers(progressr::handler_cli(
    format = "{cli::pb_bar} {cli::pb_percent} | {cli::pb_eta_str}"
  ))
  progressr::with_progress(enable = verbose, {
    total_tasks <- length(sub_summ)
    p <- progressr::progressor(steps = total_tasks)

    sub_summ <- lapply(sub_summ, function(args) {
      args$progressor <- p
      args$subbsn_rast <- subbs
      args$all_rasts <- all_rasts
      args
    })

    extract_execute <- lapply(sub_summ, function(args) {
      future::futureCall(
        .extract_subbasins,
        args = args,
        seed = NULL,
        globals = c("args"),
        packages = c(
          "stats",
          "sf",
          "terra",
          "exactextractr",
          "dplyr",
          "tidyr",
          "purrr",
          "Hmisc",
          "tibble",
          "ihydro"
        )
      )
    })
    extract_value <- lapply(extract_execute, future::value)
  })

  if (verbose) {
    cli::cli_alert_info("Combining subbasin attributes across catchments.")
  }

  # Combine extracted attributes into a single table
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
        batch,
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
      by = c("link_id_otarget", "subbasin_link_id","batch")
    )
  }

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

  if (!is.null(out_filename)) {
    write.csv(
      result,
      out_filename,
      row.names = FALSE
    )
  }

  return(result)
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

#' Extract raster attributes for a set of subbassins (for parallel)
#' @noRd
.extract_subbasins <- carrier::crate(
  function(
    input_file,
    link_id,
    all_rasts,
    subbsn_rast = NULL,
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

    extra_out <- NULL

    all_rasts <- terra::unwrap(all_rasts)
    all_rasts <- terra::subset(all_rasts, c(numb_rast, cat_rast, iDW_cols))

    fun <- c()
    if (median) {
      fun <- c(fun, "median")
    }
    if (!is.null(quantiles)) {
      fun <- c(fun, "quantile")
    }

    if (length(fun) > 0) {
      subbasins <- sf::read_sf(
        input_file,
        query = ihydro:::.build_sql_in(catch_source, "link_id", unique(link_id))
      )

      max_cells_in_memory <- ihydro:::.max_cells_in_memory_helper(
        mem_fraction = mem_fraction,
        n_cores = n_cores
      )

      extra_out <- exactextractr::exact_extract(
        x = all_rasts,
        y = subbasins,
        fun = fun,
        quantiles = quantiles,
        default_value = NA_real_,
        max_cells_in_memory = max_cells_in_memory,
        progress = FALSE,
        force_df = TRUE,
        full_colnames = TRUE
      )

      sub_tbl <- tibble::tibble(
        link_id_otarget = link_id,
        subbasin_link_id = NA_character_
      )
    } else {
      subbsn_rast <- terra::unwrap(subbsn_rast)
      pairs <- tibble::tibble(
        loi = c(numb_rast, cat_rast),
        iDW = list(iDW_cols)
      )
      pairs <- tidyr::unnest(pairs, iDW)
      pairs$new_cols <- paste(pairs[[1]], pairs[[2]], sep = ".")
      pairs$count_col <- paste(pairs[[2]], pairs[[1]], sep = ".")
      pairs <- tibble::as_tibble(pairs)

      max_cells_in_memory <- ihydro:::.max_cells_in_memory_helper(
        mem_fraction = mem_fraction,
        n_cores = n_cores,
        ncol = terra::nlyr(all_rasts)
      )

      stat_fns <- list(
        sum = ~ {
          if (sum(!is.na(.x)) == 0) {
            return(NA_real_)
          }
          sum(.x, na.rm = TRUE)
        },
        mean = ~ {
          if (sum(!is.na(.x)) == 0) {
            return(NA_real_)
          }
          mean(.x, na.rm = TRUE)
        },
        stdev = ~ {
          if (sum(!is.na(.x)) < 2) {
            return(NA_real_)
          }
          stats::sd(.x, na.rm = TRUE)
        },
        variance = ~ {
          if (sum(!is.na(.x)) < 2) {
            return(NA_real_)
          }
          stats::var(.x, na.rm = TRUE)
        },
        count = ~ sum(!is.na(.x)),
        min = ~ {
          if (sum(!is.na(.x)) == 0) {
            return(NA_real_)
          }
          min(.x, na.rm = TRUE)
        },
        max = ~ {
          if (sum(!is.na(.x)) == 0) {
            return(NA_real_)
          }
          max(.x, na.rm = TRUE)
        }
      )

      sub_sum <- ihydro:::.extract_subbasin_stats(
        link_id = link_id,
        subbsn_rast = subbsn_rast,
        all_rasts = all_rasts,
        numb_rast = numb_rast,
        cat_rast = cat_rast,
        iDW_cols = iDW_cols,
        pairs = pairs,
        stat_fns = stat_fns,
        max_cells_in_memory = max_cells_in_memory
      )

      extra_out <- sub_sum

      sub_tbl <- tibble::tibble(
        link_id_otarget = NA_character_,
        subbasin_link_id = extra_out$link_id
      )

      extra_out <- dplyr::select(extra_out, -link_id)
    }

    if (!is.null(progressor)) {
      progressor()
    }

    gg <- gc(verbose = FALSE)

    dplyr::bind_cols(
      sub_tbl,
      extra_out
    )
  }
)

#' Extract raster attributes for a single subbasin (with memory-safe batching)
#'
#' @param link_ids  Character vector of subbasin link IDs.
#' @param subbsn_rast  Unwrapped `SpatRaster` with subbasin IDs (link_id * 1000).
#' @param all_rasts  Unwrapped `SpatRaster` stack (LOI + iDW layers).
#' @param numb_rast  Character vector of numeric LOI layer names.
#' @param cat_rast   Character vector of categorical LOI layer names.
#' @param pairs      Tibble with columns `loi`, `iDW` (unnested), `new_cols`, `count_col`.
#' @param stat_fns   Named list of lambda functions passed to `dplyr::across`.
#' @param max_cells_in_memory  Integer upper bound on cells loaded at once.
#' @return A tibble with one row per `link_id`.
#' @noRd
.proc_batch <- function(
    link_id,
    sub_sel = NULL,
    subbsn_rast,
    all_rasts,
    numb_rast,
    cat_rast,
    iDW_cols,
    stat_fns,
    pairs
) {
  # Don't know why round is necessary here
  cells <- terra::cells(
    subbsn_rast,
    round(as.numeric(link_id) * 1000)
  )

  if (!is.null(sub_sel)) {
    cells$link_id <- cells$link_id[sub_sel]
  }

  vals <- terra::extract(
    all_rasts,
    cells$link_id
  )

  vals <- tibble::as_tibble(vals)

  out <- list()
  out[[length(out) + 1]] <- tibble::tibble(
    count.count_internal = length(cells$link_id)
  )

  if (length(numb_rast) > 0) {
    out[[length(out) + 1]] <- dplyr::summarise(
      tibble::as_tibble(vals),
      dplyr::across(
        dplyr::any_of(numb_rast),
        .fns = stat_fns,
        .names = "{.fn}.{.col}"
      )
    )
  }

  if (length(cat_rast) > 0) {
    out[[length(out) + 1]] <- dplyr::summarise(
      tibble::as_tibble(vals),
      dplyr::across(
        dplyr::any_of(cat_rast),
        .fns = stat_fns[names(stat_fns) == "sum"],
        .names = "{.fn}.{.col}"
      )
    )
  }

  if (length(iDW_cols) > 0) {
    out[[length(out) + 1]] <- dplyr::summarise(
      tibble::as_tibble(vals),
      dplyr::across(
        dplyr::any_of(iDW_cols),
        .fns = stat_fns[names(stat_fns) == "sum"],
        .names = "{.fn}.{.col}"
      )
    )
  }

  for (i in seq_len(nrow(pairs))) {
    v1 <- try(suppressWarnings(Hmisc::wtd.var(
      vals[[pairs[[1]][i]]],
      vals[[pairs[[2]][i]]],
      method = "ML" # This returns population level variance
    )),
    silent = TRUE)

    if (inherits(v1,"try-error")) {
      v1 <- NA_real_
    }

    out[[length(out) + 1]] <- tibble::tibble(
      weighted_mean = suppressWarnings(Hmisc::wtd.mean(
        vals[[pairs[[1]][i]]],
        vals[[pairs[[2]][i]]]
      )),
      weighted_variance = v1,
      weighted_sum = sum(
        vals[[pairs[[1]][i]]] * vals[[pairs[[2]][i]]],
        na.rm = T
      )
    )
    colnames(out[[length(out)]]) <- paste(
      colnames(out[[length(out)]]),
      pairs[[3]][[i]],
      sep = "."
    )

    out[[length(out) + 1]] <- tibble::tibble(
      weighted_sum = sum(
        !is.na(vals[[pairs[[1]][i]]]) * vals[[pairs[[2]][i]]],
        na.rm = T
      )
    )
    colnames(out[[length(out)]]) <- paste(
      colnames(out[[length(out)]]),
      pairs[[4]][[i]],
      sep = "."
    )
  }

  out <- purrr::list_cbind(out)
  out$link_id <- link_id

  return(out)
}


#' Extract and summarise per-subbasin raster statistics with memory-safe batching
#'
#' @param link_ids  Character vector of subbasin link IDs.
#' @param subbsn_rast  Unwrapped `SpatRaster` with subbasin IDs (link_id * 1000).
#' @param all_rasts  Unwrapped `SpatRaster` stack (LOI + iDW layers).
#' @param numb_rast  Character vector of numeric LOI layer names.
#' @param cat_rast   Character vector of categorical LOI layer names.
#' @param pairs      Tibble with columns `loi`, `iDW` (unnested), `new_cols`, `count_col`.
#' @param stat_fns   Named list of lambda functions passed to `dplyr::across`.
#' @param max_cells_in_memory  Integer upper bound on cells loaded at once.
#' @return A tibble with one row per `link_id`.
#' @noRd
.extract_subbasin_stats <- function(
    link_ids,
    subbsn_rast,
    all_rasts,
    numb_rast,
    cat_rast,
    iDW_cols,
    pairs,
    stat_fns,
    max_cells_in_memory
) {
  # Count cells per link_id
  cell_counts <- lapply(
    as.numeric(link_ids) * 1000,
    function(x) {
      # Don't know why round is necessary here
      length(terra::cells(subbsn_rast, round(x))$link_id)
    }
  )
  names(cell_counts) <- link_ids

  # Batch link_ids
  batches <- lapply(
    cell_counts,
    function(x) {
      split(
        seq_len(x),
        rep(
          1:ceiling(x / max_cells_in_memory),
          each = max_cells_in_memory,
          length.out = x
        )
      )
    }
  )

  # Process each batch
  batch_results <- purrr::map2(
    names(batches),
    batches,
    function(link_id, x) {
      lapply(x, function(y) {
        ihydro:::.proc_batch(
          link_id = link_id,
          sub_sel = y,
          subbsn_rast = subbsn_rast,
          all_rasts = all_rasts,
          numb_rast = numb_rast,
          cat_rast = cat_rast,
          iDW_cols = iDW_cols,
          stat_fns = stat_fns,
          pairs = pairs
        )
      })
    }
  )

  batch_results <- lapply(batch_results, purrr::list_rbind, names_to = "batch" )
  purrr::list_rbind(batch_results)
}

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
    mem_fraction = 0.5,
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
  fun <- match.arg(fun, several.ok = TRUE)

  safe_kmeans <- function(x, centers, ...) {
    x0 <- x
    x <- dplyr::group_by(x[, c("link_id", "cent_x", "cent_y")], link_id)
    x <- dplyr::summarise(
      x,
      cent_x = mean(cent_x),
      cent_y = mean(cent_y),
      .groups = "drop"
    )
    if (nrow(x) <= centers) {
      if (nrow(x) == nrow(x0)) {
        return(
          list(
            cluster = seq_len(nrow(x))
          )
        )
      } else {
        tmp <- x0 |>
          dplyr::group_by(link_id) |>
          dplyr::mutate(cluster = dplyr::cur_group_id())

        return(
          list(
            cluster = tmp$cluster
          )
        )
      }
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
    query = .build_sql_in(
      "Subbasins_poly",
      "link_id",
      as.numeric(unique(subb_lookup$link_id))
    )
  )

  n_jobs <- chunks_per_worker
  n_jobs <- min(n_jobs, nrow(subbasins))

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
    dplyr::select(subb_ids, link_id, unn_group),
    by = "link_id"
  )

  subb_unn_group <- subb_lookup |>
    dplyr::left_join(
      tasks[, c("link_id", "unn_group")],
      by = c("final_link_id" = "link_id")
    ) |>
    dplyr::select(-final_link_id)

  tasks_s <- NULL
  tasks_o <- NULL
  tasks_lump <- NULL

  if (length(ws_s) > 0L | length(ws_o) > 0L) {
    tasks_s_group <- safe_kmeans(
      tasks[, c("link_id", "cent_x", "cent_y")],
      n_jobs
    )$cluster

    tasks_s <- lapply(
      split(tasks$link_id, tasks_s_group),
      function(x) {
        unn_grp <- c()
        if (length(ws_o) > 0L) {
          unn_grp <- unique(subb_unn_group$unn_group[
            subb_unn_group$link_id %in% x
          ])
          unn_grp <- as.vector(sapply(ws_o, function(x) {
            paste0(x, "_unn_group", unn_grp)
          }))
        }

        list(
          input_file = input$outfile,
          link_id = x,
          loi_rast_input = loi_rast_input$outfile,
          loi_summary = TRUE,
          numb_rast = numb_rast,
          cat_rast = cat_rast,
          iDW_rast_input = iDW_rast_input$outfile,
          iDW_cols = c(ws_s, unn_grp),
          mem_fraction = mem_fraction,
          n_cores = n_cores,
          include_count = TRUE #!count_with_o #!count_with_o
        )
      }
    )
  }

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
          mem_fraction = mem_fraction,
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

  summarize_ws_so <- function(df) {
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
    # S/O-targeted iDW
    for (ws in intersect(
      weighting_scheme,
      c("iFLS", "HAiFLS", "iFLO", "HAiFLO")
    )) {
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
            res[[paste0(var, "_", ws, "_sd")]] <- combine_weighted_sample_sd( #combine_weighted_sample_sd
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

  # Main logic
  result <- list()
  if (!is.null(extract_value$ws_lump) && nrow(extract_value$ws_lump) > 0) {
    ws_lump <- clean_data(extract_value$ws_lump)
    result$ws_lump <- ws_lump |>
      dplyr::group_by(link_id_otarget) |>
      tidyr::nest() |>
      dplyr::mutate(summary = purrr::map(data, summarize_ws_lump))
  }
  if (!is.null(extract_value$ws_s) && nrow(extract_value$ws_s) > 0) {
    ws_s <- clean_data(extract_value$ws_s)
    result$ws_s <- ws_s |>
      dplyr::group_by(link_id_otarget) |>
      tidyr::nest() |>
      dplyr::mutate(summary = purrr::map(data, summarize_ws_so))
  }

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
