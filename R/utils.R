#' Extract a table layer from an ihydro GeoPackage, returning a tibble
#' @noRd
.ihydro_extract_tbl <- function(x, layer, collect) {
  .check_ihydro(x)
  .check_ihydro_layer(x, layer)
  con <- DBI::dbConnect(RSQLite::SQLite(), x$outfile)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  out <- dplyr::tbl(con, layer)

  if (collect) {
    out <- dplyr::collect(out)
  }

  return(out)
}

#' Extract a raster layer from an ihydro GeoPackage, returning a SpatRaster
#' @noRd
.ihydro_extract_rast <- function(x, layer) {
  .check_ihydro(x)
  .check_ihydro_layer(x, layer)

  out <- lapply(layer, function(lyr) terra::rast(x$outfile, lyr))
  out <- terra::rast(out)

  return(out)
}

#' Validate that an object is class ihydro
#' @noRd
.check_ihydro <- function(x, arg = deparse(substitute(x))) {
  if (!inherits(x, "ihydro")) {
    cli::cli_abort("{.arg {arg}} must be of class {.cls ihydro}.")
  }
  if (!file.exists(x$outfile)) {
    cli::cli_abort("GeoPackage file not found at {.path {x$outfile}}.")
  }
  invisible(x)
}

#' Validate that ihydro layers exist
#' @noRd
.check_ihydro_layer <- function(x, layer) {
  .check_ihydro(x)
  layer_list <- ihydro_layers(x)
  missing <- layer[!layer %in% layer_list$layer_name]
  if (length(missing) > 0) {
    if (all(missing == "site_id_col")) {
      return(invisible(x))
    } else {
      cli::cli_abort(
        "{.qty {length(missing)}} Layer{?s} {.val {missing}} not found in ihydro object"
      )
    }
  }
  invisible(x)
}

#' Read the site_id_col metadata from the database
#' @noRd
.read_site_id_col <- function(db_fp) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  dplyr::collect(dplyr::tbl(con, "site_id_col"))$site_id_col
}

#' Resolve target link IDs from sample points and/or link IDs
#' @noRd
.target_id_fun <- function(
  db_fp,
  sample_points = NULL,
  link_id = NULL
) {
  if (is.character(db_fp)) {
    db_fp <- as_ihydro(db_fp)
  }
  site_id_col <- read_ihydro(db_fp, "site_id_col")$site_id_col

  all_points <- read_ihydro(db_fp, "stream_links_attr") |>
    dplyr::mutate(
      dplyr::across(
        c(link_id, tidyselect::any_of(site_id_col)),
        as.character
      )
    ) |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(site_id_col),
        ~ dplyr::na_if(., "")
      )
    )

  sample_points <- as.character(sample_points)
  link_id <- as.character(link_id)

  # When neither filter is provided, return all link_ids

  if (length(sample_points) == 0 && length(link_id) == 0) {
    target_ids <- all_points |>
      tibble::as_tibble() |>
      dplyr::select(link_id, tidyselect::any_of(site_id_col))
  } else {
    parts <- list()

    if (site_id_col != "link_id" && length(sample_points) > 0) {
      parts[[1]] <- all_points |>
        tibble::as_tibble() |>
        dplyr::select(link_id, tidyselect::any_of(site_id_col)) |>
        dplyr::filter(!!rlang::sym(site_id_col) %in% sample_points)
    }

    if (length(link_id) > 0) {
      link_id_sel <- link_id
      parts[[length(parts) + 1]] <- all_points |>
        tibble::as_tibble() |>
        dplyr::select(link_id, tidyselect::any_of(site_id_col)) |>
        dplyr::filter(link_id %in% link_id_sel)
    }

    target_ids <- dplyr::bind_rows(parts)
  }

  target_ids <- dplyr::distinct(target_ids)

  dplyr::mutate(target_ids, link_id = as.character(link_id))
}

#' Build a SQL IN-clause query for reading from gpkg
#' @noRd
.build_sql_in <- function(table, column, values) {
  vals <- paste0("'", values, "'", collapse = ", ")
  paste0(
    "SELECT `link_id`, `geom` FROM `",
    table,
    "` WHERE (`",
    column,
    "` IN (",
    vals,
    "))"
  )
}

# ── GPKG raster writer ────────────────────────────────────────────────────────

#' Write a SpatRaster layer into a GeoPackage
#' @noRd
.write_raster_gpkg <- function(
  x,
  filename,
  layer_name = NULL,
  na_val = -9999.9999
) {
  if (is.null(layer_name)) {
    layer_name <- names(x)
  }

  # Ensure the sentinel value isn't already present

  x <- terra::classify(
    x,
    rbind(
      c(na_val, na_val + .Machine$double.eps),
      c(NA_real_, na_val)
    )
  )

  # x[x == na_val] <- na_val + .Machine$double.eps
  # x[is.na(x)] <- na_val
  terra::NAflag(x) <- na_val

  terra::writeRaster(
    x = x,
    filename = filename,
    NAflag = na_val,
    filetype = "GPKG",
    datatype = "FLT4S",
    gdal = c(
      "APPEND_SUBDATASET=YES",
      paste0("RASTER_TABLE=", layer_name),
      paste0("RASTER_IDENTIFIER=", layer_name),
      paste0("RASTER_DESCRIPTION=", layer_name),
      paste0("FIELD_NAME=", layer_name),
      paste0("UOM=", layer_name),
      paste0("QUANTITY_DEFINITION=", layer_name),
      "TILE_FORMAT=TIFF",
      "METADATA_TABLES=YES",
      "COMPRESS=LZW",
      "TILED=YES",
      "BAND_COUNT=1"
    )
  )

  # Update null sentinel in coverage metadata

  con <- DBI::dbConnect(RSQLite::SQLite(), filename)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(
    con,
    paste0(
      "UPDATE gpkg_2d_gridded_coverage_ancillary SET data_null = ",
      na_val,
      " WHERE tile_matrix_set_name = '",
      layer_name,
      "'"
    )
  )

  invisible(NULL)
}


# ── Temp-directory helper ─────────────────────────────────────────────────────

#' Ensure a temporary directory exists, return its normalised path
#' @noRd
.ensure_temp_dir <- function(temp_dir = NULL) {
  if (is.null(temp_dir)) {
    temp_dir <- tempfile()
  }
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  normalizePath(temp_dir)
}

# ── Worker-count helper ──────────────────────────────────────────────────────

#' Determine the number of available workers from the current future plan
#' @noRd
.n_workers <- function() {
  n <- future::nbrOfWorkers()
  if (is.infinite(n)) {
    n <- future::availableCores(logical = FALSE)
  }
  max(n, 1L)
}


# ── Terra memory management helpers ──────────────────────────────────────────

#' Save current terraOptions and set per-worker memory limits
#'
#' Computes a safe per-worker memory cap based on total system RAM divided by
#' `n_cores`, then calls `terra::terraOptions()` with `memmax` and `memfrac`.
#' Returns the previous options so they can be restored later.
#'
#' @param n_cores Integer. Number of parallel workers sharing the machine.
#' @param temp_dir Character. Directory for terra temporary files.
#' @param verbose Logical. Forwarded to `terra::terraOptions(verbose = .)`.
#' @param mem_fraction Numeric (0, 1]. Fraction of total RAM available to all
#'   workers combined.  Defaults to 0.8 (80 %).
#'
#' @return A list of previous terraOptions, suitable for `.restore_terra_options()`.
#' @noRd
.set_terra_options <- function(
  n_cores = 1L,
  temp_dir = tempdir(),
  verbose = FALSE,
  mem_fraction = 0.5
) {
  old_opts <- terra::terraOptions(print = FALSE)

  total_target_frac <- mem_fraction
  per_worker_frac <- total_target_frac / max(1L, as.integer(n_cores))
  per_worker_frac <- max(0.05, min(total_target_frac, per_worker_frac))

  max_mem <- .avail_mem(
    n_cores = n_cores,
    mem_fraction = per_worker_frac,
    min_per_worker_gb = 2.5
  )

  progress <- 0

  terra::terraOptions(
    memmax = max_mem / 1024^3,
    memfrac = per_worker_frac,
    tempdir = temp_dir,
    verbose = verbose,
    progress = progress,
    parallel = n_cores == 1
  )

  old_opts
}

#' Restore terraOptions from a snapshot
#'
#' Safely restores options previously captured by `.set_terra_options()`,
#' filtering out read-only fields that `terraOptions()` does not accept as
#' inputs.
#'
#' @param old_opts List returned by `terra::terraOptions(print = FALSE)`.
#' @noRd
.restore_terra_options <- function(old_opts) {
  opts <- old_opts
  # Remove read-only / non-settable fields
  opts$metadata <- NULL
  opts$names <- NULL
  if (identical(opts$filetype, "")) {
    opts$filetype <- NULL
  }
  do.call(terra::terraOptions, opts)
}

# Calculate an appropriate per-worker memory limit based on total system RAM and number of workers
#' @noRd
.avail_mem <- function(
  n_cores = 1L,
  mem_fraction = 0.5,
  min_per_worker_gb = 2.5
) {
  max_mem <- as.numeric(memuse::Sys.meminfo()$totalram) * mem_fraction
  max_mem <- max_mem / max(1L, as.integer(n_cores))
  max_mem <- max(min_per_worker_gb * 1024^3, max_mem) # enforce a minimum per-worker memory
  return(max_mem)
}


#' @noRd
.max_cells_in_memory_helper <- function(
  n_cores = 1L,
  mem_fraction = 0.5,
  ncol = NULL
) {
  max_mem <- memuse::Sys.meminfo()$totalram * mem_fraction
  if (is.null(ncol)) {
    max_square <- memuse::howmany(max_mem)
  } else {
    max_square <- memuse::howmany(max_mem, ncol = ncol + 3)
    max_square[[2]] <- 1
  }

  max_cells_in_memory <- max_square[[1]] * max_square[[2]]
  max_cells_in_memory <- floor(max_cells_in_memory / n_cores)

  if (max_cells_in_memory >= .Machine$integer.max) {
    max_cells_in_memory <- .Machine$integer.max - 1
  }
  return(max_cells_in_memory)
}

