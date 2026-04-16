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
  out_filename,
  loi_file = NULL,
  loi_cols = NULL,
  iDW_file = NULL,
  store_iDW = FALSE,
  sample_points = NULL,
  link_id = NULL,
  target_o_type = c("point", "segment_point", "segment_whole"),
  weighting_scheme = NULL,
  loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
  inv_function = NULL,
  write_strategy = NULL,
  mem_fraction = 0.5,
  temp_dir = NULL,
  verbose = FALSE
) {
  # ─ Validate inputs ───────────────────────────
  check_ihydro(input)
  stopifnot(is.logical(store_iDW), length(store_iDW) == 1L)
  stopifnot(mem_fraction < 0.9 | mem_fraction > 0.1)

  if (!is.null(iDW_file)) {
    if (!is.null(store_iDW)) {
      cli::cli_warn(
        "Both {.arg iDW_file} and {.arg store_iDW} are set. Ignoring {.arg store_iDW} and using {.arg iDW_file} for weights storage."
      )
      store_iDW <- FALSE
    }
    if (!is.null(inv_function)) {
      cli::cli_warn(
        "Both {.arg iDW_file} and {.arg inv_function} are set. Ignoring {.arg inv_function} since weights will be read from {.arg iDW_file}."
      )
    }
    if (!is.null(weighting_scheme)) {
      cli::cli_warn(
        "Both {.arg iDW_file} and {.arg weighting_scheme} are set. Ignoring {.arg weighting_scheme} since weights will be read from {.arg iDW_file}."
      )
    }
    if (!is.null(write_strategy)) {
      cli::cli_warn(
        "Both {.arg iDW_file} and {.arg write_strategy} are set. Ignoring {.arg write_strategy} since weights will be read from {.arg iDW_file} and not written."
      )
    }
    if (!inherits(iDW_file, "ihydro") && !file.exists(iDW_file)) {
      cli::cli_abort(
        "{.arg iDW_file} does not exist at {.val {iDW_file}}."
      )
    }
  }

  if (is.null(inv_function)) {
    inv_function <- function(x) (x * 0.001 + 1)^-1
  }
  if (is.null(weighting_scheme)) {
    weighting_scheme <- c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO")
  }
  if (is.null(write_strategy)) {
    write_strategy <- c("sequential", "batched")
  }

  target_o_type <- match.arg(target_o_type)
  weighting_scheme <- match.arg(
    weighting_scheme,
    several.ok = TRUE,
    choices = c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO")
  )
  loi_numeric_stats <- match.arg(
    loi_numeric_stats,
    several.ok = TRUE,
    choices = c("mean", "sd", "median", "min", "max", "sum")
  )
  loi_numeric_stats <- stats::setNames(loi_numeric_stats, loi_numeric_stats)
  write_strategy <- match.arg(
    write_strategy,
    choices = c("sequential", "batched")
  )

  input_path <- input$outfile
  temp_dir <- ensure_temp_dir(temp_dir)
  n_cores <- n_workers()

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
  site_id_col <- read_ihydro(input, "site_id_col")[[1]]

  subb_ids <- read_ihydro(input, "stream_links_attr") |>
    dplyr::filter(link_lngth > 0) |> # small catchments without streams
    dplyr::filter(link_id %in% target_ids$link_id) |> # small catchments without streams
    dplyr::select(link_id, tidyselect::any_of(site_id_col), USChnLn_To) |>
    dplyr::left_join(
      dplyr::mutate(target_o_meta, link_id = as.character(link_id)),
      by = c("link_id")
    )

  subb_ids <- subb_ids[!is.na(subb_ids$unn_group), ]

  subb_lookup <- read_ihydro(
    input,
    "fcon_pwise_dist"
  ) |>
    dplyr::select(
      final_link_id = destination,
      link_id = origin #subbasin ID
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
  sub_summ <- .extract_planner(
    input = input,
    subb_ids = subb_ids,
    loi_rast_input = loi_file,
    numb_rast = numb_rast,
    cat_rast = cat_rast,
    iDW_rast_input = iDW_out,
    iDW_cols = weighting_scheme,
    max_cells_in_memory = max_cells_in_memory,
    n_cores = n_cores,
    chunks_per_worker = 4L,
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
    quantiles = NULL
  )
  sub_summ <- unlist(sub_summ, recursive = FALSE)

  extract_execute <- lapply(sub_summ, function(args) {
    future::futureCall(.extract_subbasins, args = args)
  })
  extract_value <- lapply(extract_execute, future::value)

  # Combine extracted attributes into a single table
  extract_value_comb <- list()
  extract_value_nms <- names(extract_value)
  ws_s_nms <- extract_value_nms[grepl("^ws_s\\.", extract_value_nms)]
  ws_o_nms <- extract_value_nms[grepl("^ws_o\\.", extract_value_nms)]
  ws_lump_nms <- extract_value_nms[grepl("^ws_lump\\.", extract_value_nms)]

  extract_value_comb$ws_s <- dplyr::bind_rows(extract_value[ws_s_nms])
  extract_value_comb$ws_0 <- dplyr::bind_rows(extract_value[ws_o_nms])
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

  if (length(ws_o_nms) > 0) {
    extract_value_comb$ws_0 <- extract_value_comb$ws_0 |>
      dplyr::left_join(
        dplyr::select(subb_ids, link_id, unn_group),
        by = c("link_id_otarget" = "link_id")
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
    dplyr::select(, link_id, tidyselect::any_of(site_id_col)) |>
    dplyr::left_join(
      result,
      by = "link_id"
    ) |>
    dplyr::mutate(
      status = "Complete",
      .after = 2
    )

  # # ─ Extract raster attributes ──────────────────────
  # temp_dir_sub <- file.path(temp_dir, basename(tempfile()))
  # dir.create(temp_dir_sub, recursive = TRUE)
  #
  # result <- list()
  #
  # # Identify numeric and categorical rasters from LOI metadata
  # numb_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "num_rast"]
  # cat_rast <- loi_meta$loi_var_nms[loi_meta$loi_type == "cat_rast"]
  #
  # # Format stats names for weighted extraction functions
  # weighted_stats <- loi_numeric_stats
  # weighted_stats_mod <- weighted_stats <- weighted_stats[
  #   weighted_stats %in% c("mean", "sd", "var")
  # ]
  # weighted_stats_mod[weighted_stats_mod == "mean"] <- "weighted_mean"
  # weighted_stats_mod[weighted_stats_mod == "sd"] <- "weighted_stdev"
  # weighted_stats_mod[weighted_stats_mod == "var"] <- "weighted_variance"
  #
  # lumped_stats <- loi_numeric_stats
  # lumped_stats[lumped_stats == "sd"] <- "stdev"
  # lumped_stats[lumped_stats == "var"] <- "variance"
  #
  # # Load LOI rasters
  # loi_rasts <- terra::rast(loi_path, loi_cols)
  #
  # # Unweighted LOI summaries (lumped) ---------------------------------------
  # lumped_out <- list()
  # if ("lumped" %in% weighting_scheme) {
  #   if (length(numb_rast) > 0) {
  #     if (verbose) {
  #       message("Calculating lumped numeric attributes")
  #     }
  #
  #     lumped_out$numb_rast <- exactextractr::exact_extract(
  #       loi_rasts[[numb_rast]],
  #       catch,
  #       fun = lumped_stats,
  #       progress = FALSE,
  #       full_colnames = TRUE,
  #       force_df = TRUE,
  #       max_cells_in_memory = max_cells_in_memory
  #     )
  #
  #     nms <- paste0(
  #       rep(
  #         numb_rast,
  #         length.out = length(loi_numeric_stats) * length(numb_rast)
  #       ),
  #       "_lumped_",
  #       rep(loi_numeric_stats, each = length(numb_rast))
  #     )
  #
  #     names(lumped_out$numb_rast) <- nms
  #
  #     lumped_out$numb_rast <- dplyr::mutate(
  #       lumped_out$numb_rast,
  #       link_id = catch$link_id,
  #       .before = 1
  #     )
  #
  #     # Convert population to sample sd
  #     if (any(loi_numeric_stats == "sd")) {
  #       catch_ncell <- read_ihydro(input, "stream_links_attr")
  #       catch_ncell <- dplyr::select(catch_ncell, link_id, sbbsn_area)
  #       cell_size <- terra::res(read_ihydro(input, "dem_final"))
  #       catch_ncell$sbbsn_area <- catch_ncell$sbbsn_area /
  #         (cell_size[1] * cell_size[2])
  #
  #       lumped_out$numb_rast <- lumped_out$numb_rast |>
  #         dplyr::left_join(
  #           catch_ncell,
  #           by = "link_id"
  #         ) |>
  #         dplyr::mutate(
  #           dplyr::across(
  #             tidyselect::ends_with("_lumped_sd"),
  #             ~ .x * sqrt((sbbsn_area - 1) / sbbsn_area)
  #           )
  #         ) |>
  #         dplyr::select(-sbbsn_area)
  #     }
  #   }
  #
  #   if (length(cat_rast) > 0) {
  #     if (verbose) {
  #       message("Calculating lumped categorical attributes")
  #     }
  #
  #     loi_rasts_cat <- loi_rasts[[cat_rast]]
  #     loi_rasts_cat[["catch_size"]] <- loi_rasts_cat[[1]]
  #     loi_rasts_cat[["catch_size"]][] <- 1
  #
  #     lumped_out$cat_rast <- exactextractr::exact_extract(
  #       loi_rasts_cat,
  #       catch,
  #       fun = c("sum"),
  #       progress = FALSE,
  #       full_colnames = TRUE,
  #       force_df = TRUE,
  #       max_cells_in_memory = max_cells_in_memory
  #     )
  #
  #     lumped_out$cat_rast <- lumped_out$cat_rast |>
  #       dplyr::mutate(
  #         dplyr::across(
  #           c(tidyselect::everything(), -sum.catch_size),
  #           ~ .x / sum.catch_size
  #         )
  #       ) |>
  #       dplyr::select(-sum.catch_size) |>
  #       dplyr::rename_with(
  #         ~ paste0(
  #           gsub("^sum\\.", "", .x),
  #           "_lumped_prop"
  #         )
  #       ) |>
  #       dplyr::mutate(
  #         link_id = catch$link_id,
  #         .before = 1
  #       )
  #   }
  #
  #   result$lumped <- purrr::reduce(lumped_out, dplyr::left_join, by = 'link_id')
  # }
  #
  # # Stream-targeted iDW rasters (iFLS / HAiFLS) -----------------------------
  # ws_s <- intersect(weighting_scheme, c("iFLS", "HAiFLS"))
  # iDWs_rasts <- NULL
  # iDWsLOI_rasts <- NULL
  # weighted_out <- list()
  #
  # if (length(ws_s) > 0L) {
  #   weighted_out$ws_s <- list()
  #   weighted_out$ws_sc <- list()
  #
  #   iDWs_rasts <- terra::rast(iDW_path, ws_s)
  #
  #   # -- Stream-targeted iDW numeric rasters (iFLS / HAiFLS) ----------------------
  #   if (length(numb_rast) > 0) {
  #     if (verbose) {
  #       message(
  #         "Calculating weighted numeric stream targeted (iFLS / HAiFLS) attributes"
  #       )
  #     }
  #
  #     weighted_out$ws_s <- .ft_numb_stats(
  #       ws = ws_s,
  #       loi_rasts = loi_rasts,
  #       iDW_rasts = iDWs_rasts,
  #       catch = catch,
  #       weighted_stats = weighted_stats_mod,
  #       numb_rast = numb_rast,
  #       progress = FALSE,
  #       full_colnames = TRUE,
  #       force_df = TRUE,
  #       max_cells_in_memory = max_cells_in_memory
  #     )
  #   }
  #
  #   # -- Stream-targeted iDW categorical rasters (iFLS / HAiFLS) -----------------
  #   if (length(cat_rast) > 0) {
  #     if (verbose) {
  #       message(
  #         "Calculating weighted categorical stream targeted (iFLS / HAiFLS) attributes"
  #       )
  #     }
  #
  #     weighted_out$ws_sc <- .ft_cat_stats(
  #       ws = ws_s,
  #       loi_rasts = loi_rasts,
  #       iDW_rasts = iDWs_rasts,
  #       catch = catch,
  #       cat_rast = cat_rast,
  #       progress = FALSE,
  #       full_colnames = TRUE,
  #       force_df = TRUE,
  #       max_cells_in_memory = max_cells_in_memory
  #     )
  #   }
  #
  #   result$weighted_s <- purrr::reduce(
  #     weighted_out,
  #     dplyr::left_join,
  #     by = 'link_id'
  #   )
  # }
  #
  # # Site-targeted iDW rasters (iFLO/ HAiFLO) -----------------------------
  # ws_o <- intersect(weighting_scheme, c("iFLO", "HAiFLO"))
  #
  # if (length(ws_o) > 0) {
  #   result_int <- list()
  #
  #   cnt <- 0
  #   tot <- length(unique(subb_ids$unn_group))
  #   for (i in unique(subb_ids$unn_group)) {
  #     cnt <- cnt + 1
  #     weighted_out$ws_s <- list()
  #     weighted_out$ws_sc <- list()
  #
  #     unn_group <- i
  #
  #     catch_sub <- subb_ids |>
  #       dplyr::filter(unn_group == i) |>
  #       dplyr::pull(link_id)
  #
  #     catch_sub <- dplyr::filter(
  #       catch,
  #       link_id %in% catch_sub
  #     )
  #
  #     o_lyrs <- paste0(
  #       rep(ws_o, each = length(unique(unn_group))),
  #       "_unn_group",
  #       rep(unique(unn_group), times = length(ws_o))
  #     )
  #     iDWo_rasts <- terra::rast(iDW_path, o_lyrs)
  #     names(iDWo_rasts) <- gsub(
  #       paste0("_unn_group", unn_group),
  #       "",
  #       names(iDWo_rasts)
  #     )
  #
  #     # -- Outlet-targeted iDW numeric rasters (iFLO / HAiFLO) ----------------------
  #     if (length(numb_rast) > 0) {
  #       if (verbose) {
  #         message(
  #           paste0(
  #             "Calculating weighted numeric outlet targeted (iFLO / HAiFLO) attributes: ",
  #             cnt,
  #             "/",
  #             tot
  #           )
  #         )
  #       }
  #
  #       weighted_out$ws_s <- .ft_numb_stats(
  #         ws = ws_o,
  #         loi_rasts = loi_rasts,
  #         iDW_rasts = iDWo_rasts,
  #         catch = catch_sub,
  #         weighted_stats = weighted_stats_mod,
  #         numb_rast = numb_rast,
  #         progress = FALSE,
  #         full_colnames = TRUE,
  #         force_df = TRUE,
  #         max_cells_in_memory = max_cells_in_memory
  #       )
  #     }
  #
  #     # -- Outlet-targeted iDW categorical rasters (iFLO / HAiFLO) -----------------
  #     if (length(cat_rast) > 0) {
  #       if (verbose) {
  #         message(
  #           paste0(
  #             "Calculating weighted categorical outlet targeted (iFLO / HAiFLO) attributes: ",
  #             cnt,
  #             "/",
  #             tot
  #           )
  #         )
  #       }
  #
  #       weighted_out$ws_sc <- .ft_cat_stats(
  #         ws = ws_o,
  #         loi_rasts = loi_rasts,
  #         iDW_rasts = iDWo_rasts,
  #         catch = catch_sub,
  #         cat_rast = cat_rast,
  #         progress = FALSE,
  #         full_colnames = TRUE,
  #         force_df = TRUE,
  #         max_cells_in_memory = max_cells_in_memory
  #       )
  #     }
  #
  #     result_int[[i]] <- purrr::reduce(
  #       weighted_out,
  #       dplyr::left_join,
  #       by = 'link_id'
  #     )
  #   }
  #
  #   result$weighted_o <- dplyr::bind_rows(
  #     result_int
  #   )
  # }
  #
  # result <- purrr::reduce(
  #   result,
  #   dplyr::left_join,
  #   by = 'link_id'
  # )
  #
  # result <- tibble::as_tibble(result)
  #
  # result <- dplyr::left_join(
  #   result,
  #   dplyr::select(
  #     subb_ids,
  #     link_id,
  #     tidyselect::any_of(site_id_col)
  #   ),
  #   by = "link_id"
  # )
  #
  # result <- dplyr::mutate(
  #   result,
  #   dplyr::across(
  #     c(tidyselect::everything(), -link_id, -tidyselect::any_of(site_id_col)),
  #     ~ case_when(
  #       is.nan(.x) ~ 0,
  #       is.na(.x) ~ 0,
  #       T ~ .x
  #     )
  #   )
  # )
  #
  # result <- dplyr::mutate(
  #   result,
  #   status = "Complete"
  # )
  #
  # result <- dplyr::select(
  #   result,
  #   link_id,
  #   tidyselect::any_of(site_id_col),
  #   status,
  #   tidyselect::everything()
  # )

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


#' # ────────────────────────────────────────
#' # Raster extraction pipeline
#' # ────────────────────────────────────────
#'
#' #' Weighted stats with IDW and numeric rasters
#' #' @noRd
#' .ft_numb_stats <- function(
#'   ws,
#'   loi_rasts,
#'   iDW_rasts,
#'   catch,
#'   weighted_stats,
#'   numb_rast,
#'   ...
#' ) {
#'   weighted_out <- list()
#'   for (i in ws) {
#'     weighted_out[[i]] <- exactextractr::exact_extract(
#'       loi_rasts[[numb_rast]],
#'       catch,
#'       fun = weighted_stats,
#'       weights = iDW_rasts[[i]],
#'       default_weight = 0,
#'       ...
#'     )
#'
#'     nms <- paste0(
#'       rep(
#'         numb_rast,
#'         length.out = length(weighted_stats) * length(i) * length(numb_rast)
#'       ),
#'       "_",
#'       rep(i, length.out = length(numb_rast) * length(weighted_stats)),
#'       "_",
#'       rep(weighted_stats, each = length(numb_rast) * length(i))
#'     )
#'
#'     names(weighted_out[[i]]) <- nms
#'
#'     weighted_out[[i]] <- dplyr::mutate(
#'       weighted_out[[i]],
#'       link_id = catch$link_id,
#'       .before = 1
#'     )
#'   }
#'   purrr::reduce(
#'     weighted_out,
#'     dplyr::left_join,
#'     by = 'link_id'
#'   )
#' }
#'
#' #' Weighted stats with IDW and categorical rasters
#' #' @noRd
#' .ft_cat_stats <- function(
#'   ws,
#'   loi_rasts,
#'   iDW_rasts,
#'   catch,
#'   cat_rast,
#'   ...
#' ) {
#'   w_sum <- exactextractr::exact_extract(
#'     iDW_rasts,
#'     catch,
#'     fun = "sum",
#'     ...
#'   )
#'   w_sum <- dplyr::mutate(
#'     w_sum,
#'     link_id = catch$link_id,
#'     .before = 1
#'   ) |>
#'     dplyr::rename_with(~ gsub("^sum\\.", "", .x))
#'
#'   weighted_out <- list()
#'   for (i in ws) {
#'     weighted_out[[i]] <- exactextractr::exact_extract(
#'       loi_rasts[[cat_rast]],
#'       catch,
#'       fun = "weighted_sum",
#'       weights = iDW_rasts[[i]],
#'       default_weight = 0,
#'       ...
#'     )
#'
#'     weighted_out[[i]] <- weighted_out[[i]] / w_sum[[i]]
#'
#'     nms <- paste0(
#'       rep(cat_rast, length.out = length(i) * length(cat_rast)),
#'       "_",
#'       rep(i, length.out = length(cat_rast)),
#'       "_prop"
#'     )
#'
#'     names(weighted_out[[i]]) <- nms
#'
#'     weighted_out[[i]] <- dplyr::mutate(
#'       weighted_out[[i]],
#'       link_id = catch$link_id,
#'       .before = 1
#'     )
#'   }
#'   purrr::reduce(
#'     weighted_out,
#'     dplyr::left_join,
#'     by = 'link_id'
#'   )
#' }

#' Extract raster attributes
#' @noRd
.extract_fun <- function(
  all_rasts,
  subbasins,
  x_cols = NULL,
  weight_cols = NULL,
  default_value = NA_real_,
  default_weight = NA_real_,
  fun,
  quantiles = NULL,
  max_cells_in_memory = 3e+07,
  include_count = FALSE
) {
  if (is.null(x_cols) && is.null(weight_cols)) {
    return(NULL)
  }

  weights_resolve <- function(all_rasts, weight_cols = NULL) {
    if (is.null(weight_cols)) {
      return(NULL)
    }
    all_rasts[[weight_cols]]
  }

  if (include_count) {
    all_rasts[["count_internal"]] <- all_rasts[[1]]
    all_rasts[["count_internal"]][] <- 1
    x_cols <- c(x_cols, "count_internal")
    fun <- unique(c(fun, "count"))
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
    for (i in weight_cols) {
      out[[i]] <- exactextractr::exact_extract(
        x = all_rasts[[x_cols]],
        y = subbasins,
        fun = fun,
        quantiles = quantiles,
        weights = weights_resolve(all_rasts, i),
        default_value = default_value,
        default_weight = default_weight,
        max_cells_in_memory = max_cells_in_memory,
        progress = FALSE,
        force_df = TRUE,
        full_colnames = TRUE
      )
    }
  }

  dplyr::bind_cols(out)
}


#' Extract raster attributes for a set of subbassins (for parallel)
#' @noRd
.extract_subbasins <- carrier::crate(
  function(
    input_file,
    link_id,
    link_id_otarget = NA_character_,
    loi_rast_input,
    loi_summary = TRUE,
    numb_rast,
    cat_rast,
    iDW_rast_input,
    iDW_cols,
    catch_source = c("Subbasins_poly", "Catchment_poly"),
    median = FALSE,
    quantiles = NULL,
    max_cells_in_memory = 3e+07,
    include_count = FALSE
  ) {
    catch_source <- match.arg(catch_source)

    subbasins <- sf::read_sf(
      input_file,
      query = ihydro:::build_sql_in(catch_source, "link_id", unique(link_id))
    )

    iDW_rasts <- NULL
    loi_rasts <- NULL
    iDW_out <- NULL
    numb_out <- NULL
    cat_out <- NULL
    numb_out_iDW <- NULL
    cat_out_iDW <- NULL
    extra_out <- NULL

    if (!is.null(loi_rast_input)) {
      loi_rasts <- terra::rast(loi_rast_input, c(numb_rast, cat_rast))
    }
    if (!is.null(iDW_rast_input)) {
      iDW_rasts <- terra::rast(iDW_rast_input, iDW_cols)
    }
    all_rasts <- list(loi_rasts, iDW_rasts)
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
        include_count = include_count
      )
    } else {
      iDW_out <- ihydro:::.extract_fun(
        all_rasts = all_rasts,
        subbasins = subbasins,
        x_cols = iDW_cols,
        weight_cols = NULL,
        fun = "sum",
        max_cells_in_memory = max_cells_in_memory,
        include_count = include_count
      )

      if (loi_summary) {
        if (length(numb_rast) > 0) {
          numb_out <- ihydro:::.extract_fun(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = c(numb_rast),
            weight_cols = NULL,
            fun = c("sum", "mean", "stdev", "variance", "count", "min", "max"),
            max_cells_in_memory = max_cells_in_memory,
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
            fun = c("sum", "count"),
            max_cells_in_memory = max_cells_in_memory,
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
            default_value = 0,
            default_weight = 0,
            fun = c("weighted_sum", "weighted_mean", "weighted_variance"),
            max_cells_in_memory = max_cells_in_memory,
            include_count = FALSE
          )
        }
        if (length(cat_rast) > 0) {
          cat_out_iDW <- ihydro:::.extract_fun(
            all_rasts = all_rasts,
            subbasins = subbasins,
            x_cols = c(cat_rast),
            weight_cols = iDW_cols,
            default_value = 0,
            default_weight = 0,
            fun = c("weighted_sum"),
            max_cells_in_memory = max_cells_in_memory,
            include_count = FALSE
          )
        }
      }
    }

    if (catch_source == "Catchment_poly") {
      link_id_otarget <- link_id
      link_id <- NA_character_
    }

    dplyr::bind_cols(
      tibble::tibble(
        link_id_otarget = link_id_otarget,
        subbasin_link_id = link_id
      ),
      iDW_out,
      numb_out,
      cat_out,
      numb_out_iDW,
      cat_out_iDW,
      extra_out
    )
  }
)

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
  quantiles = NULL
) {
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

  n_jobs <- n_cores * chunks_per_worker
  n_jobs <- min(n_jobs, length(unique(subb_ids$link_id)))

  ws_lumped <- "lumped" %in% iDW_cols
  ws_s <- iDW_cols[iDW_cols %in% c("iFLS", "HAiFLS")]
  ws_o <- iDW_cols[iDW_cols %in% c("iFLO", "HAiFLO")]

  subb_lookup <- ihydro::read_ihydro(
    input,
    "fcon_pwise_dist"
  ) |>
    dplyr::select(
      final_link_id = destination,
      link_id = origin #subbasin ID
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
    dplyr::select(subb_ids, link_id, unn_group),
    by = "link_id"
  )

  tasks_s <- NULL
  tasks_o <- NULL
  tasks_lump <- NULL

  # develop ws_s tasks
  count_with_o <- TRUE
  if (length(ws_s) > 0L) {
    count_with_o <- FALSE

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
          iDW_cols = ws_s,
          max_cells_in_memory = max_cells_in_memory,
          include_count = !count_with_o
        )
      }
    )
  }

  # develop ws_o tasks
  if (length(ws_o) > 0L) {
    subb_lookup <- subb_lookup |>
      dplyr::left_join(
        tasks[, c("link_id", "cent_x", "cent_y")],
        by = c("link_id" = "link_id")
      ) |>
      dplyr::left_join(
        tasks[, c("link_id", "unn_group")],
        by = c("final_link_id" = "link_id")
      )

    tasks_o <- split(subb_lookup, subb_lookup$unn_group)

    # group accordinging to location by unn_group
    tasks_o <- lapply(tasks_o, function(x) {
      x$group <- safe_kmeans(
        x[, c("link_id", "cent_x", "cent_y")],
        n_jobs
      )$cluster
      x <- split(x, x$group)
      return(x)
    })

    tasks_o <- unlist(tasks_o, recursive = F)

    tasks_o <- lapply(
      tasks_o,
      function(x) {
        list(
          input_file = input$outfile,
          link_id = x$link_id,
          link_id_otarget = x$final_link_id,
          loi_rast_input = loi_rast_input$outfile,
          loi_summary = FALSE,
          numb_rast = numb_rast,
          cat_rast = cat_rast,
          iDW_rast_input = iDW_rast_input$outfile,
          iDW_cols = paste(ws_o, unique(x$unn_group), sep = "_unn_group"),
          max_cells_in_memory = max_cells_in_memory,
          include_count = count_with_o
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
          max_cells_in_memory = max_cells_in_memory,
          include_count = FALSE
        )
      }
    )
  }

  out <- list(ws_s = tasks_s, ws_o = tasks_o, ws_lump = tasks_lump)
  return(out)
}

.catchment_assemble <- function(
  input,
  subb_ids
) {
  # subbasin_funs <- fun[
  #   fun %in% c("mean", "sd", "var", "cv", "min", "max", "sum", "count")
  # ]
  # names(subbasin_funs) <- subbasin_funs
  # subbasin_funs <- dplyr::case_when(
  #   subbasin_funs == "sd" ~ "stdev",
  #   subbasin_funs == "var" ~ "variance",
  #   subbasin_funs == "cv" ~ "coefficient_of_variation",
  #   T ~ subbasin_funs
  # )
  # subbasin_wfuns <- fun[fun %in% c("mean", "sd", "var", "sum")]
  # names(subbasin_wfuns) <- subbasin_wfuns
  # subbasin_wfuns <- dplyr::case_when(
  #   subbasin_wfuns == "mean" ~ "weighted_mean",
  #   subbasin_wfuns == "sd" ~ "weighted_stdev",
  #   subbasin_wfuns == "var" ~ "weighted_variance",
  #   subbasin_wfuns == "sum" ~ "weighted_sum",
  #   T ~ subbasin_wfuns
  # )
  # cathcment_funs <- fun[fun %in% c("median", "quantile")]
  # names(cathcment_funs) <- cathcment_funs
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
  extract_value <- lapply(extract_value, function(x) {
    x[sapply(x, is.nan)] <- 0
    x <- tidyr::pivot_longer(
      x,
      c(-link_id_otarget, -subbasin_link_id, -tidyselect::any_of("unn_group"))
    )

    if ("unn_group" %in% colnames(x)) {
      x$name <- gsub(
        paste0(
          paste0("_unn_group", unique(x$unn_group)),
          collapse = "|"
        ),
        "",
        x$name
      )
    }
    return(x)
  })

  result <- list()

  # Whole catchment summaries (median and quantiles) ------------------------
  if (!is.null(extract_value$ws_lump)) {
    result$ws_lump <- extract_value$ws_lump |>
      dplyr::group_by(link_id_otarget) |>
      tidyr::nest()

    result$ws_lump$summary <- lapply(
      result$ws_lump$data,
      function(x) {
        res <- list()

        if ("lumped" %in% weighting_scheme) {
          x <- x[!is.na(x$value), ]
          x <- tidyr::pivot_wider(
            x,
            values_fill = 0
          )

          for (var in numeric_vars) {
            if ("median" %in% loi_numeric_stats) {
              # pull from median subtable
              res[[paste0(var, "_lumped_median")]] <- x[[paste0(
                "median.",
                var
              )]]
            }
            # TODO: add quantiles
          }
        }
        tibble::as_tibble(res)
      }
    )
  }

  # Lumped and S-targeted iDW summaries --------------------------------------------------------
  if (!is.null(extract_value$ws_s)) {
    result$ws_s_lump <- extract_value$ws_s |>
      dplyr::group_by(link_id_otarget) |>
      tidyr::nest()

    result$ws_s_lump$summary <- lapply(
      result$ws_s_lump$data,
      function(x) {
        res <- list()

        if ("lumped" %in% weighting_scheme) {
          x <- x[!is.na(x$value), ]
          x <- tidyr::pivot_wider(
            x,
            values_fill = 0
          )

          for (var in numeric_vars) {
            if ("mean" %in% loi_numeric_stats) {
              sum_x <- sum(x[[paste0("sum.", var)]], na.rm = TRUE)
              sum_n <- sum(x[[paste0("count.", var)]], na.rm = TRUE)
              res[[paste0(var, "_lumped_mean")]] <- sum_x / sum_n
            }

            if ("var" %in% loi_numeric_stats || "sd" %in% loi_numeric_stats) {
              means <- x[[paste0("mean.", var)]]
              vars <- x[[paste0("variance.", var)]] # chunk population variances
              ns <- x[[paste0("count.", var)]]
              overall_mean <- sum(means * ns, na.rm = TRUE) /
                sum(ns, na.rm = TRUE)
              numerator <- sum((ns - 1) * vars, na.rm = TRUE) +
                sum(ns * (means - overall_mean)^2, na.rm = TRUE)
              denominator <- sum(ns, na.rm = TRUE) - 1
              sample_var <- numerator / denominator
              if ("var" %in% loi_numeric_stats) {
                # sample variance (unbiased)
                res[[paste0(var, "_lumped_var")]] <- sample_var
              }
              if ("sd" %in% loi_numeric_stats) {
                # sample standard deviation (unbiased)
                res[[paste0(var, "_lumped_sd")]] <- sqrt(sample_var)
              }
            }

            if ("min" %in% loi_numeric_stats) {
              res[[paste0(var, "_lumped_min")]] <- min(
                x[[paste0("min.", var)]],
                na.rm = TRUE
              )
            }
            if ("max" %in% loi_numeric_stats) {
              res[[paste0(var, "_lumped_max")]] <- max(
                x[[paste0("max.", var)]],
                na.rm = TRUE
              )
            }
            if ("sum" %in% loi_numeric_stats) {
              res[[paste0(var, "_lumped_sum")]] <- sum(
                x[[paste0("sum.", var)]],
                na.rm = TRUE
              )
            }
          }

          # Categorical proportions
          for (cat in cat_vars) {
            sum_col <- paste0("sum.", cat)
            count_col <- paste0("count.", "count_internal")
            if (sum_col %in% names(x) && count_col %in% names(x)) {
              res[[paste0(cat, "_lumped_prop")]] <- sum(
                x[[sum_col]],
                na.rm = TRUE
              ) /
                sum(x[[count_col]], na.rm = TRUE)
            }
          }
        }

        # Weighted summaries
        ws <- setdiff(weighting_scheme, c("lumped", "iFLO", "HAiFLO"))
        for (ws in setdiff(weighting_scheme, "lumped")) {
          ws_data_col <- ws
          for (var in numeric_vars) {
            mean_col <- paste0("weighted_mean.", var, ".", ws_data_col)
            var_col <- paste0("weighted_variance.", var, ".", ws_data_col)
            sum_col <- paste0("weighted_sum.", var, ".", ws_data_col)
            wt_col <- paste0("sum.", ws_data_col)
            if (all(c(mean_col, var_col, wt_col) %in% names(x))) {
              sub_mean <- x[[mean_col]]
              sub_mean <- sub_mean[!is.na(sub_mean)]
              sub_var <- x[[var_col]]
              sub_var <- sub_var[!is.na(sub_var)]
              sub_wt <- x[[wt_col]]
              sub_wt <- sub_wt[!is.na(sub_wt)]
              if ("mean" %in% loi_numeric_stats) {
                res[[paste0(var, "_", ws, "_mean")]] <- combine_weighted_mean(
                  sub_mean,
                  sub_wt
                )
              }
              if ("sd" %in% loi_numeric_stats) {
                res[[paste0(var, "_", ws, "_sd")]] <- combine_weighted_sd(
                  sub_mean,
                  sub_var,
                  sub_wt
                )
              }
              if ("var" %in% loi_numeric_stats) {
                res[[paste0(var, "_", ws, "_var")]] <- combine_weighted_var(
                  sub_mean,
                  sub_var,
                  sub_wt
                )
              }
              if ("sum" %in% loi_numeric_stats && sum_col %in% names(x)) {
                res[[paste0(var, "_", ws, "_sum")]] <- sum(
                  x[[sum_col]],
                  na.rm = TRUE
                )
              }
            }
          }
          # Categorical proportions (weighted)
          for (cat in cat_vars) {
            sum_col <- paste0("weighted_sum.", cat, ".", ws)
            wt_col <- paste0("sum.", ws)
            if (sum_col %in% names(x) && wt_col %in% names(x)) {
              res[[paste0(cat, "_", ws, "_prop")]] <- sum(
                x[[sum_col]],
                na.rm = TRUE
              ) /
                sum(x[[wt_col]], na.rm = TRUE)
            }
          }
        }
        tibble::as_tibble(res)
      }
    )
  }

  # O-targeted iDW summaries ------------------------
  if (!is.null(extract_value$ws_0)) {
    result$ws_0 <- extract_value$ws_0 |>
      dplyr::group_by(link_id_otarget) |>
      tidyr::nest()

    result$ws_0$summary <- lapply(
      result$ws_0$data,
      function(x) {
        res <- list()

        x <- x[!is.na(x$value), ]
        x <- tidyr::pivot_wider(
          x,
          values_fill = 0
        )

        # Weighted summaries
        for (ws in setdiff(weighting_scheme, c("lumped", "iFLS", "HAiFLS"))) {
          ws_data_col <- ws
          for (var in numeric_vars) {
            mean_col <- paste0("weighted_mean.", var, ".", ws_data_col)
            var_col <- paste0("weighted_variance.", var, ".", ws_data_col)
            sum_col <- paste0("weighted_sum.", var, ".", ws_data_col)
            wt_col <- paste0("sum.", ws_data_col)
            if (all(c(mean_col, var_col, wt_col) %in% names(x))) {
              sub_mean <- x[[mean_col]]
              sub_mean <- sub_mean[!is.na(sub_mean)]
              sub_var <- x[[var_col]]
              sub_var <- sub_var[!is.na(sub_var)]
              sub_wt <- x[[wt_col]]
              sub_wt <- sub_wt[!is.na(sub_wt)]
              if ("mean" %in% loi_numeric_stats) {
                res[[paste0(var, "_", ws, "_mean")]] <- combine_weighted_mean(
                  sub_mean,
                  sub_wt
                )
              }
              if ("sd" %in% loi_numeric_stats) {
                res[[paste0(var, "_", ws, "_sd")]] <- combine_weighted_sd(
                  sub_mean,
                  sub_var,
                  sub_wt
                )
              }
              if ("var" %in% loi_numeric_stats) {
                res[[paste0(var, "_", ws, "_var")]] <- combine_weighted_var(
                  sub_mean,
                  sub_var,
                  sub_wt
                )
              }
              if ("sum" %in% loi_numeric_stats && sum_col %in% names(x)) {
                res[[paste0(var, "_", ws, "_sum")]] <- sum(
                  x[[sum_col]],
                  na.rm = TRUE
                )
              }
            }
          }
          # Categorical proportions (weighted)
          for (cat in cat_vars) {
            sum_col <- paste0("weighted_sum.", cat, ".", ws)
            wt_col <- paste0("sum.", ws)
            if (sum_col %in% names(x) && wt_col %in% names(x)) {
              res[[paste0(cat, "_", ws, "_prop")]] <- sum(
                x[[sum_col]],
                na.rm = TRUE
              ) /
                sum(x[[wt_col]], na.rm = TRUE)
            }
          }
        }
        tibble::as_tibble(res)
      }
    )
  }

  result <- lapply(
    result,
    function(x) {
      x |>
        dplyr::select(-data) |>
        tidyr::unnest(summary)
    }
  )

  result <- purrr::reduce(result, dplyr::left_join, by = "link_id_otarget")
  result <- dplyr::ungroup(result)
  result <- dplyr::rename(result, link_id = link_id_otarget)
  return(result)
}
