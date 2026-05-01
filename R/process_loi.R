#' Process layers of interest (LOI) for stream network attribution
#'
#' Standardises numeric and categorical landscape layers to match the DEM's
#' resolution, extent, and CRS. Categorical layers are one-hot encoded.
#'
#' @details
#' This function is designed to process any landscape layers of interest (LOI) that
#' users want to attribute to the stream network. It can handle both numeric and
#' categorical raster inputs, as well as vector inputs that will be rasterized.
#' The processed layers are aligned to the DEM used for hydrological processing,
#' ensuring consistency in resolution, extent, and CRS. Categorical layers are one-hot encoded
#' to create separate binary rasters for each category. The processed layers are saved to a
#' GeoPackage for efficient storage and retrieval.
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
#' @param input `NULL` or an `ihydro` object. If `NULL`, `dem` is required.
#' @param dem `NULL` or a DEM (path, `SpatRaster`, etc.).
#' @param clip_region Optional clip polygon/raster.
#' @param num_inputs Named list of numeric raster layers.
#' @param cat_inputs Named list of categorical raster layers.
#' @param output_filename `NULL` or path ending in `.gpkg`. If `NULL`, must be appended to `input`.
#' @param variable_names Named list of variable names to keep per input.
#' @param return_products Logical. Return raster products?
#' @param temp_dir Temporary file directory.
#' @param write_strategy Character. How processed rasters are written to the
#'   GeoPackage. `"sequential"` (default) waits for all parallel workers to
#'   finish, then writes every raster to the GeoPackage in one pass — simplest
#'   and avoids nested-future connection warnings. `"batched"` processes layers
#'   in chunks of size `n_cores`, writing completed rasters to the GeoPackage
#'   between chunks to reduce peak temporary disk usage.
#' @param verbose Logical.
#' @param overwrite Logical. Overwrite existing layers?
#'
#' @return An `ihydro` object.
#' The contained rasters are the processed LOI layers, aligned to the DEM and saved in the specified GeoPackage. If `return_products` is `TRUE`, the rasters are also returned as `SpatRaster` objects within the `ihydro` object for immediate use in R.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
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
#'
#' # Numeric Raster
#'
#' whitebox::wbt_slope(
#'  dem = file.path(ex_loc, "toy_dem.tif"),
#'  output = file.path(ex_loc, "slope.tif")
#' )
#'
#' # Combine loi layers
#' output_filename_loi <- file.path(ex_loc, "Processed_loi.gpkg")
#'
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
#' ihydro::ihydro_layers(loi_combined, n = Inf)
#'
#' }
#'

process_loi <- function(
  input = NULL,
  dem = NULL,
  clip_region = NULL,
  num_inputs = list(),
  cat_inputs = list(),
  output_filename = NULL,
  variable_names = NULL,
  return_products = FALSE,
  write_strategy = c("sequential", "batched"),
  temp_dir = NULL,
  verbose = FALSE,
  overwrite = TRUE
) {
  # ── Validation ──────────────────────────────────────────────────────────
  write_strategy <- match.arg(write_strategy)

  if (!is.null(input)) {
    check_ihydro(input)
  }
  stopifnot(
    is.logical(return_products),
    is.logical(verbose),
    is.logical(overwrite)
  )

  temp_dir <- ensure_temp_dir(temp_dir)

  if (!inherits(num_inputs, "list")) {
    cli::cli_abort("{.arg num_inputs} must be a named list.")
  }
  if (length(num_inputs) > 0 && is.null(names(num_inputs))) {
    cli::cli_abort("Objects in {.arg num_inputs} must be named.")
  }
  if (!inherits(cat_inputs, "list")) {
    cli::cli_abort("{.arg cat_inputs} must be a named list.")
  }
  if (length(cat_inputs) > 0 && is.null(names(cat_inputs))) {
    cli::cli_abort("Objects in {.arg cat_inputs} must be named.")
  }
  if (any(names(num_inputs) %in% names(cat_inputs))) {
    cli::cli_abort(
      "{.arg num_inputs} and {.arg cat_inputs} cannot share names."
    )
  }
  if (is.null(variable_names)) {
    if (verbose) {
      message(
        "No variables specified in 'variable_names'; all variables will be used."
      )
    }
  }
  if (is.null(input) && is.null(dem)) {
    cli::cli_abort("Either {.arg input} or {.arg dem} must be specified.")
  }
  if (is.null(input) && is.null(output_filename)) {
    cli::cli_abort(
      "Either {.arg input} or {.arg output_filename} must be specified."
    )
  }

  if (!is.null(input) && is.null(output_filename)) {
    output_filename <- input$outfile
  }
  output_filename <- normalizePath(output_filename, mustWork = FALSE)
  if (!grepl("\\.gpkg$", output_filename)) {
    cli::cli_abort("{.arg output_filename} must end in {.val .gpkg}.")
  }
  if (dirname(output_filename) == temp_dir) {
    cli::cli_abort("{.arg output_filename} must not be inside {.arg temp_dir}.")
  }

  dir.create(dirname(output_filename), showWarnings = FALSE, recursive = TRUE)

  whitebox::wbt_options(
    exe_path = whitebox::wbt_exe_path(),
    verbose = verbose > 2,
    wd = temp_dir
  )
  terra::terraOptions(verbose = verbose > 3, tempdir = temp_dir)

  # ── Prepare DEM ─────────────────────────────────────────────────────────
  if (verbose) {
    message("Preparing DEM")
  }
  if (!is.null(input)) {
    dem <- read_ihydro(input, "dem_final")
  } else {
    dem <- process_input(dem, working_dir = temp_dir)
    if (!inherits(dem, "SpatRaster")) {
      cli::cli_abort("{.arg dem} must resolve to a {.cls SpatRaster}.")
    }
  }
  terra::writeRaster(
    dem,
    file.path(temp_dir, "dem_final.tif"),
    overwrite = TRUE,
    gdal = "COMPRESS=NONE",
    datatype = "FLT4S"
  )

  # Clip region
  if (is.null(clip_region)) {
    clip_region <- terra::as.polygons(
      terra::rast(terra::ext(dem), crs = terra::crs(dem))
    )
  } else {
    clip_region <- process_input(
      clip_region,
      working_dir = temp_dir,
      align_to = terra::vect(
        "POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",
        crs = terra::crs(dem)
      )
    )
  }
  terra::writeVector(
    clip_region,
    file.path(temp_dir, "clip_region.shp"),
    overwrite = TRUE
  )

  # ── Materialise inputs to disk ──────────────────────────────────────────
  num_inputs <- purrr::map(num_inputs, materialise_input, temp_dir = temp_dir)
  cat_inputs <- purrr::map(cat_inputs, materialise_input, temp_dir = temp_dir)

  all_names <- c(names(num_inputs), names(cat_inputs))
  all_inputs <- c(num_inputs, cat_inputs)
  all_types <- c(
    rep("num_rast", length(num_inputs)),
    rep("cat_rast", length(cat_inputs))
  )

  if (!is.null(variable_names)) {
    lyr_variables <- variable_names[all_names]

    for (nm in names(variable_names)) {
      ip <- process_input(all_inputs[[nm]])
      if (!variable_names[[nm]] %in% names(ip)) {
        cli::cli_abort(
          "Variable {.val {variable_names[[nm]]}} not found in input {.val {nm}}."
        )
      }
    }
  } else {
    lyr_variables <- rep(list(NA_character_), length(all_names))
  }
  lyr_variables[sapply(lyr_variables, is.null)] <- list(NA_character_)

  # ── Create gpkg if needed ───────────────────────────────────────────────
  if (!file.exists(output_filename)) {
    dem |>
      terra::ext() |>
      terra::as.polygons(crs = terra::crs(dem)) |>
      sf::st_as_sf() |>
      sf::write_sf(
        output_filename,
        layer = "DEM_Extent",
        append = TRUE,
        delete_layer = FALSE,
        delete_dsn = FALSE
      )
  }

  # ── Process LOI (parallel) ─────────────────────────────────────────────
  n_cores <- n_workers()
  max_cores_opt <- getOption("parallelly.maxWorkers.localhost")
  on.exit(options(parallelly.maxWorkers.localhost = max_cores_opt), add = TRUE)
  options(parallelly.maxWorkers.localhost = n_cores)

  if (n_cores > 1) {
    oplan <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(oplan), add = TRUE)
  }

  temp_dir_save <- file.path(temp_dir, basename(tempfile()))
  dir.create(temp_dir_save)

  if (verbose) {
    message("Processing layers of interest (LOI) across ", n_cores, " cores")
  }

  ip <- list(
    lyr_nms = as.list(all_names),
    lyr = all_inputs,
    lyr_variables = lyr_variables,
    rln = as.list(all_types),
    temp_dir = rep(list(temp_dir), length(all_names)),
    temp_dir_save = rep(list(temp_dir_save), length(all_names)),
    output_filename = rep(list(output_filename), length(all_names)),
    overwrite = rep(list(overwrite), length(all_names)),
    p = rep(list(NULL), length(all_names))
  )

  if (write_strategy == "batched") {
    # ── Batched: process in chunks, drain to gpkg between chunks ────────
    chunks <- split(
      seq_along(all_names),
      ceiling(seq_along(all_names) / max(n_cores, 1L))
    )
    future_proc <- vector("list", length(all_names))

    progressr::with_progress(enable = verbose, {
      p <- progressr::progressor(steps = length(all_names))
      ip$p <- rep(list(p), length(all_names))

      for (chunk in chunks) {
        ip_chunk <- purrr::map(ip, `[`, chunk)

        results_chunk <- furrr::future_pmap(
          ip_chunk,
          .options = furrr::furrr_options(
            globals = FALSE,
            seed = NULL,
            scheduling = 4L
          ),
          process_single_loi_worker
        )
        future_proc[chunk] <- results_chunk

        # Drain this batch's rasters to gpkg immediately
        fl <- list.files(temp_dir_save, "\\.tif$", full.names = TRUE)
        for (f in fl) {
          r <- terra::rast(f)
          if (verbose) {
            message("Writing: ", names(r))
          }
          res <- try(write_raster_gpkg(r, output_filename), silent = TRUE)
          if (inherits(res, "try-error")) {
            msg <- conditionMessage(attr(res, "condition"))
            if (!msg %in% c("stoi", "stol")) stop(msg)
          }
          file.remove(f)
        }
      }
    })
  } else {
    # ── Sequential: process all, then write ─────────────────────────────
    progressr::with_progress(enable = verbose, {
      p <- progressr::progressor(steps = length(all_names))
      ip$p <- rep(list(p), length(all_names))

      future_proc <- furrr::future_pmap(
        ip,
        .options = furrr::furrr_options(
          globals = FALSE,
          seed = NULL,
          scheduling = 4L
        ),
        process_single_loi_worker
      )
    })

    if (verbose) {
      message("Writing outputs")
    }
    fl <- list.files(temp_dir_save, "\\.tif$", full.names = TRUE)
    for (f in fl) {
      r <- terra::rast(f)
      if (verbose) {
        message("Writing: ", names(r))
      }
      res <- try(write_raster_gpkg(r, output_filename), silent = TRUE)
      if (inherits(res, "try-error")) {
        msg <- conditionMessage(attr(res, "condition"))
        if (!msg %in% c("stoi", "stol")) stop(msg)
      }
      file.remove(f)
    }
  }

  ot <- future_proc
  names(ot) <- all_names

  meta <- purrr::map_dfr(
    ot,
    ~ tibble::tibble(
      loi_lyr_nms = .x[["lyr_nms"]],
      loi_var_nms = .x[["lyr_variables"]],
      loi_type = .x[["rln"]]
    )
  )

  if (verbose) {
    message("Saving metadata")
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), output_filename)
  dplyr::copy_to(
    con,
    meta,
    "loi_meta",
    overwrite = TRUE,
    temporary = FALSE,
    analyze = TRUE,
    in_transaction = TRUE
  )
  DBI::dbDisconnect(con)

  # ── Return ──────────────────────────────────────────────────────────────
  output <- if (is.null(input)) list(outfile = output_filename) else input

  if (return_products) {
    ott <- purrr::map(
      meta |> split(meta$loi_type),
      ~ purrr::map(.$loi_var_nms, ~ terra::rast(output_filename, lyrs = .))
    )
    rast_lists <- list(
      num_inputs = try(
        terra::rast(
          unlist(
            ott["num_rast"],
            use.names = FALSE
          )
        ),
        silent = TRUE
      ),
      cat_inputs = try(
        terra::rast(
          unlist(
            ott["cat_rast"],
            use.names = FALSE
          )
        ),
        silent = TRUE
      )
    )
    rast_lists <- rast_lists[!sapply(rast_lists, inherits, "try-error")]
    output <- c(rast_lists, output)
  }

  clean_temp_files(temp_dir)
  structure(output, class = "ihydro")
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Write in-memory R objects to disk so workers can read them
#' @noRd
materialise_input <- function(x, temp_dir) {
  if (inherits(x, "character")) {
    return(x)
  }
  if (inherits(x, c("sf", "sfc", "sfg", "SpatVector"))) {
    fp <- file.path(temp_dir, paste0(basename(tempfile()), ".shp"))
    if (inherits(x, "SpatVector")) {
      x <- sf::st_as_sf(x)
    }
    sf::write_sf(x, fp)
    return(fp)
  }
  if (inherits(x, "SpatRaster")) {
    fp <- file.path(temp_dir, paste0(basename(tempfile()), ".tif"))
    terra::writeRaster(x, fp, datatype = "FLT4S")
    return(fp)
  }
  x
}

#' Worker function for processing a single LOI
#' @noRd
process_single_loi_worker <- carrier::crate(
  function(
    lyr_nms,
    lyr,
    lyr_variables,
    rln,
    temp_dir,
    temp_dir_save,
    output_filename,
    overwrite,
    p
  ) {
    #suppressPackageStartupMessages(library(sf)) # Not sure why, but this is necessary
    #options(dplyr.summarise.inform = FALSE, scipen = 999)
    `%>%` <- magrittr::`%>%`

    temp_sub <- file.path(temp_dir, basename(tempfile()))
    dir.create(temp_sub)

    resample_method <- ifelse(grepl("num_rast", rln), "bilinear", "near")
    if (all(is.na(lyr_variables))) {
      lyr_variables <- NULL
    }

    suppressMessages(
      output <- ihydro:::process_input(
        input = unlist(lyr),
        input_variable_names = if (!is.null(lyr_variables)) {
          unlist(lyr_variables)
        } else {
          NULL
        },
        align_to = file.path(temp_dir, "dem_final.tif"),
        clip_region = file.path(temp_dir, "clip_region.shp"),
        resample_type = resample_method,
        working_dir = temp_sub
      )
    )

    output <- terra::split(output, names(output))
    names(output) <- sapply(output, names)

    for (nm in names(output)) {
      terra::writeRaster(
        output[[nm]],
        file.path(temp_dir_save, paste0(names(output[[nm]]), ".tif")),
        datatype = "FLT4S",
        todisk = TRUE,
        overwrite = TRUE,
        gdal = "COMPRESS=NONE"
      )
    }

    unlink(temp_sub, recursive = TRUE, force = TRUE)
    p()

    list(
      lyr_nms = lyr_nms,
      lyr_variables = sapply(output, names),
      rln = rln
    )
  }
)

#' Drain temp raster files into gpkg while future is running
#' @noRd
drain_temp_rasters <- function(
  temp_dir_save,
  output_filename,
  future_proc,
  verbose
) {
  future_status <- future::futureOf(future_proc)

  write_available <- function() {
    fl <- list.files(temp_dir_save, "\\.tif$", full.names = TRUE)
    fl <- fl[file.mtime(fl) < Sys.time() - 60]
    for (f in fl) {
      r <- try(terra::rast(f), silent = TRUE)
      if (inherits(r, "try-error")) {
        next
      }
      if (verbose) {
        message("Writing: ", names(r))
      }
      res <- try(write_raster_gpkg(r, output_filename), silent = TRUE)
      if (inherits(res, "try-error")) {
        msg <- conditionMessage(attr(res, "condition"))
        if (!msg %in% c("stoi", "stol")) stop(msg)
      }
      file.remove(f)
    }
  }

  while (!future::resolved(future_status)) {
    Sys.sleep(0.5)
    write_available()
  }

  Sys.sleep(5)

  # Check for errors
  stop_on_future_errors(future_proc)
  #if (length(future_proc$result$conditions) > 0) {
  #  errs <- purrr::keep(
  #    purrr::map(future_proc$result$conditions, "condition"),
  #    ~ inherits(., "error")
  #  )
  #  if (length(errs) > 0) stop(paste(errs, collapse = "\n"))
  #}

  # Final drain
  fl <- list.files(temp_dir_save, "\\.tif$", full.names = TRUE)
  for (f in fl) {
    r <- terra::rast(f)
    if (verbose) {
      message("Writing: ", names(r))
    }
    res <- try(write_raster_gpkg(r, output_filename), silent = TRUE)
    if (inherits(res, "try-error")) {
      msg <- conditionMessage(attr(res, "condition"))
      if (!msg %in% c("stoi", "stol")) stop(msg)
    }
    file.remove(f)
  }
}
