# ── ihydro class ──────────────────────────────────────────────────────────────

#' Create an ihydro object from a GeoPackage file
#'
#' Loads an existing ihydro GeoPackage (`.gpkg`) file and wraps it in an
#' `ihydro` S3 object. The resulting object stores the file path and can be
#' passed to any other function in the package (e.g., [generate_vectors()],
#' [trace_flowpaths()], [process_loi()]).
#'
#' `as_ihydro()` is provided as a snake_case alias for convenience.
#'
#' @param file Character. File path to a GeoPackage file (`.gpkg`). The file
#'   must already exist on disk.
#'
#' @return An object of class `ihydro`, which is a list with element `outfile`
#'   containing the normalised path to the GeoPackage.
#'
#' @seealso [ihydro_layers()] to inspect the contents of an `ihydro` object,
#'   [process_flowdir()] or [process_hydrology()] to create a new GeoPackage
#'   from a DEM.
#'
#' @examples
#' \dontrun{
#' # Open an existing ihydro GeoPackage
#' my_hydro <- as.ihydro("path/to/hydrology.gpkg")
#'
#' # Inspect the layers stored inside
#' ihydro_layers(my_hydro)
#'
#' # The snake_case alias works identically
#' my_hydro <- as_ihydro("path/to/hydrology.gpkg")
#' }
#'
#' @export
as.ihydro <- function(file) {
  file <- normalizePath(file, mustWork = FALSE)
  if (!file.exists(file)) {
    cli::cli_abort("{.arg file} does not exist: {.path {file}}")
  }
  if (!grepl("\\.gpkg$", file)) {
    cli::cli_abort("{.arg file} must end in {.val .gpkg}.")
  }
  structure(list(outfile = file), class = "ihydro")
}

#' @export
#' @rdname as.ihydro
as_ihydro <- as.ihydro

# ── Layer listing ─────────────────────────────────────────────────────────────

#' List layers stored in an ihydro GeoPackage
#'
#' Queries the GeoPackage database to enumerate every raster, vector, and
#' plain-table layer it contains, along with a logical grouping label
#' (`data_group`). This is useful for inspecting what has been computed so
#' far and for verifying that all required layers are present before running
#' downstream functions.
#'
#' @details
#' Layers are classified into the following groups:
#' \describe{
#'   \item{`hydro`}{DEM-derived rasters and stream network vectors
#'     (e.g., `dem_final`, `stream_links`, `Subbasins_poly`).}
#'   \item{`flow_path`}{Upstream / downstream lookup tables (`us_flowpaths`,
#'     `ds_flowpaths`).}
#'   \item{`pwise_dist`}{Pairwise distance matrices (`fcon_pwise_dist`,
#'     `funcon_pwise_dist`).}
#'   \item{`sample_points`}{Original and snapped sampling-point layers.}
#'   \item{`loi`}{Layers of interest (landscape rasters processed by
#'     [process_loi()]).}
#'   \item{`iDW`}{Inverse distance weight rasters (`iFLS`, `HAiFLO`, etc.).}
#'   \item{`meta`}{Metadata tables (`site_id_col`, `loi_meta`, etc.).}
#' }
#'
#' @param input An object of class `ihydro`, as returned by [as.ihydro()],
#'   [process_flowdir()], or any other ihydro function.
#'
#' @return A [tibble][tibble::tibble] with columns:
#' \describe{
#'   \item{`layer_name`}{Name of the layer in the GeoPackage.}
#'   \item{`data_type`}{One of `"Raster"`, `"Vector"`, or `"Table"`.}
#'   \item{`data_group`}{Logical grouping label (see **Details**).}
#' }
#'
#' @seealso [as.ihydro()] to create an `ihydro` object.
#'
#' @examples
#' \dontrun{
#' hydro <- as.ihydro("my_watershed.gpkg")
#'
#' # List everything
#' ihydro_layers(hydro)
#'
#' # Filter to just the hydrology layers
#' ihydro_layers(hydro) |>
#'   dplyr::filter(data_group == "hydro")
#' }
#'
#' @export
ihydro_layers <- function(input) {
  check_ihydro(input)

  con <- DBI::dbConnect(RSQLite::SQLite(), input$outfile)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Raster layers
  rast_lyrs <- dplyr::tbl(con, "gpkg_contents") |>
    dplyr::collect() |>
    dplyr::filter(grepl("gridded", data_type)) |>
    dplyr::transmute(layer_name = table_name, data_type = "Raster")

  # Vector layers (identified by rtree spatial index)
  all_tables <- DBI::dbListTables(con)

  vect_lyrs <- tibble::tibble(
    layer_name = gsub(
      "^rtree_|_geom$",
      "",
      all_tables[
        grepl("^rtree_", all_tables) &
          grepl("_geom$", all_tables)
      ]
    ),
    data_type = "Vector"
  )

  # Plain tables (everything else)
  exclude <- grepl("^rtree_|^gpkg_|^sqlite_|_geom$", all_tables) |
    all_tables %in% rast_lyrs$layer_name
  tbl_lyrs <- tibble::tibble(
    layer_name = all_tables[!exclude],
    data_type = "Table"
  ) |>
    dplyr::filter(
      !layer_name %in% vect_lyrs$layer_name,
      !layer_name %in% rast_lyrs$layer_name,
    )

  out <- dplyr::bind_rows(rast_lyrs, vect_lyrs, tbl_lyrs) |>
    dplyr::mutate(
      data_group = dplyr::case_when(
        grepl("^dem_", layer_name) ~ "hydro",
        grepl("iFLS|HAiFLS|iFLO|HAiFLO", layer_name) ~ "iDW",
        layer_name %in%
          c(
            "site_id_col",
            "DEM_Extent",
            "loi_meta",
            "target_o_meta"
          ) ~ "meta",
        layer_name %in%
          c("snapped_points", "original_points") ~ "sample_points",
        layer_name %in%
          c(
            "stream_links",
            "stream_lines",
            "stream_points",
            "Subbasins_poly",
            "Catchment_poly",
            "stream_points_attr",
            "stream_links_attr"
          ) ~ "hydro",
        layer_name %in% c("ds_flowpaths", "us_flowpaths") ~ "flow_path",
        layer_name %in%
          c("funcon_pwise_dist", "fcon_pwise_dist") ~ "pwise_dist",
        TRUE ~ "loi"
      )
    ) |>
    dplyr::arrange(data_type, data_group)

  out <- as.data.frame(out)

  return(out)
}

# ── Extract Layer ─────────────────────────────────────────────────────────────

#' Extract layers stored in an ihydro GeoPackage
#'
#' Reads a specified layer from the ihydro GeoPackage and returns it as an R object (see return value for details).
#'
#' @param ihydro An object of class `ihydro`, as returned by [as.ihydro()],
#'   [process_flowdir()], or any other ihydro function.
#' @param layer Character. Name of the layer to read from the GeoPackage. Use [ihydro_layers()] to see what layers are available.
#' @param collect Logical. Whether to use [dplyr::collect()] to retieve data into a local tibble. Only relevant when `layer` is a table.
#'
#' @return The requested layer, read from the GeoPackage. The type of the returned object depends on the layer:
#' - Raster layers are returned as `SpatRaster` objects.
#' - Vector layers are returned as `SpatVector` objects.
#' - Plain tables are returned as `tibble` data frames.
#'
#'
#' @seealso [ihydro_layers()] to create an `ihydro` object.
#'
#' @examples
#' \dontrun{
#' hydro <- as.ihydro("my_watershed.gpkg")
#'
#' # List everything
#' ihydro_layers(hydro)
#'
#' # Return the stream lines
#' read_ihydro(hydro, "stream_lines")
#' }
#'
#' @export
read_ihydro <- function(ihydro, layer, collect = TRUE) {
  stopifnot(is.character(layer), length(layer) == 1L)
  check_ihydro(ihydro)
  check_ihydro_layer(ihydro, layer)

  layer_list <- ihydro_layers(ihydro)

  layer_sel <- dplyr::filter(layer_list, layer_name == layer)

  switch(
    layer_sel$data_type,
    "Raster" = terra::rast(ihydro$outfile, lyrs = layer),
    "Vector" = sf::read_sf(ihydro$outfile, layer = layer),
    "Table" = ihydro_extract_tbl(ihydro, layer, collect),
    cli::cli_abort("Unknown data type for layer {.val {layer}}")
  )
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Validate that an object is class ihydro
#' @noRd
ihydro_extract_tbl <- function(x, layer, collect) {
  check_ihydro(x)
  check_ihydro_layer(x, layer)
  con <- DBI::dbConnect(RSQLite::SQLite(), x$outfile)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  out <- dplyr::tbl(con, layer)

  if (collect) {
    out <- dplyr::collect(out)
  }

  return(out)
}

#' Validate that an object is class ihydro
#' @noRd
check_ihydro <- function(x, arg = deparse(substitute(x))) {
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
check_ihydro_layer <- function(x, layer) {
  check_ihydro(x)
  layer_list <- ihydro_layers(x)
  missing <- layer[!layer %in% layer_list$layer_name]
  if (length(missing) > 0) {
    cli::cli_abort(
      "{cli::qty(length(missing))} Layer{?s} {.val {missing}} not found in ihydro object"
    )
  }
  invisible(x)
}

#' Read the site_id_col metadata from the database
#' @noRd
read_site_id_col <- function(db_fp) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  dplyr::collect(dplyr::tbl(con, "site_id_col"))$site_id_col
}

#' Resolve target link IDs from sample points and/or link IDs
#' @noRd
target_id_fun <- function(
  db_fp,
  sample_points = NULL,
  link_id = NULL,
  segment_whole = FALSE,
  target_o_type = c(
    "point",
    "segment_point",
    "segment_whole"
  )
) {
  target_o_type <- match.arg(target_o_type)

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  site_id_col <- dplyr::collect(dplyr::tbl(con, "site_id_col"))$site_id_col

  all_points <- dplyr::collect(dplyr::tbl(con, "stream_links_attr")) |>
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

  if (segment_whole) {
    target_ids <- target_ids |>
      dplyr::select(link_id) |>
      dplyr::mutate(link_id = as.character(floor(as.numeric(link_id))))
  }

  if (target_o_type != "point") {
    available <- sf::read_sf(
      db_fp,
      query = build_sql_in("stream_lines", "link_id", target_ids$link_id)
    )
    target_ids <- dplyr::filter(target_ids, link_id %in% available$link_id)
  }

  dplyr::mutate(target_ids, link_id = as.character(link_id))
}

#' Read the target geometry (point, segment, or whole-segment)
#' @noRd
target_o_fun <- function(
  db_fp,
  target_IDs,
  target_o_type = c(
    "point",
    "segment_point",
    "segment_whole"
  )
) {
  target_o_type <- match.arg(target_o_type)
  ids <- target_IDs$link_id

  if (target_o_type == "point") {
    layer <- "stream_links"
  } else {
    layer <- "stream_lines"
  }

  target_o <- sf::read_sf(db_fp, query = build_sql_in(layer, "link_id", ids))

  if (target_o_type == "segment_whole") {
    target_o <- target_o |>
      dplyr::select(link_id) |>
      dplyr::mutate(link_id = as.character(floor(as.numeric(link_id)))) |>
      dplyr::filter(link_id %in% ids) |>
      dplyr::summarise(geom = sf::st_union(geom), .by = link_id)
  }

  target_o
}

#' Build a SQL IN-clause query for reading from gpkg
#' @noRd
build_sql_in <- function(table, column, values) {
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
write_raster_gpkg <- function(
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
      c(NA_real_,na_val)
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
ensure_temp_dir <- function(temp_dir = NULL) {
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
n_workers <- function() {
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
#' @return A list of previous terraOptions, suitable for `restore_terra_options()`.
#' @noRd
set_terra_options <- function(
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
    parallel = FALSE
  )

  old_opts
}

#' Restore terraOptions from a snapshot
#'
#' Safely restores options previously captured by `set_terra_options()`,
#' filtering out read-only fields that `terraOptions()` does not accept as
#' inputs.
#'
#' @param old_opts List returned by `terra::terraOptions(print = FALSE)`.
#' @noRd
restore_terra_options <- function(old_opts) {
  opts <- old_opts
  # Remove read-only / non-settable fields
  opts$metadata <- NULL
  opts$names <- NULL
  if (identical(opts$filetype, "")) {
    opts$filetype <- NULL
  }
  do.call(terra::terraOptions, opts)
}


# ── Access Example Files ────────────────────────────────────────────────────

#' Access example data layers included in the package
#'
#' @param layer Character. Name of the example layer to load. Must be one of:
#'   "elev_ned_30m.tif", "geology.shp", "lakes.shp", "landuse_r.tif",
#'   "pointsources.shp", "sites_nc.shp", "streams.shp"
#'
#' @return A `SpatRaster` or `SpatVector` object, depending on the layer.
#'
#' @details Original source of data: https://github.com/MiKatt/openSTARS
#'
#' @examples
#' \dontrun{
#' ex_data("elev_ned_30m.tif") # example DEM raster
#' ex_data("geology.shp") # example geology vector
#' ex_data("lakes.shp") # example lakes vector
#' ex_data("landuse_r.tif") # example land use raster
#' ex_data("pointsources.shp") # example point sources vector
#' ex_data("sites_nc.shp") # example sampling sites vector
#' ex_data("streams.shp") # example stream network vector
#' }
#'
#' @export
ex_data <- function(
  layer = c(
    "elev_ned_30m.tif",
    "geology.shp",
    "lakes.shp",
    "landuse_r.tif",
    "pointsources.shp",
    "sites_nc.shp",
    "streams.shp"
  )
) {
  layer <- match.arg(layer)

  ex_files <- system.file(
    "extdata",
    "openSTARS_data.zip",
    package = "ihydro"
  ) |>
    normalizePath() |>
    utils::unzip(list = TRUE)

  ex_files <- ex_files$Name[grepl("\\.shp|\\.tif", ex_files$Name)]
  stopifnot(layer %in% basename(ex_files))

  final_file <- file.path(
    "/vsizip",
    system.file("extdata", "openSTARS_data.zip", package = "ihydro"),
    layer
  )

  file_ext <- function(x) {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }

  switch(
    file_ext(layer),
    "shp" = terra::vect(final_file),
    "tif" = terra::rast(final_file)
  )
}

# Parse conditions (now captures call and any saved trace) and stop helper
#' @noRd
parse_future_conditions <- function(fut) {
  conds <- tryCatch(fut$result$conditions, error = function(e) NULL)
  if (is.null(conds) || length(conds) == 0) {
    return(tibble::tibble())
  }

  get_cnd <- function(x) if (!is.null(x$condition)) x$condition else x

  extract_trace <- function(cnd) {
    # try common places future/rlang store traces
    tr <- NULL
    if (!is.null(attr(cnd, "trace"))) {
      tr <- attr(cnd, "trace")
    }
    if (is.null(tr) && !is.null(cnd$trace)) {
      tr <- cnd$trace
    }
    if (is.null(tr) && !is.null(attr(cnd, "traceback"))) {
      tr <- attr(cnd, "traceback")
    }

    if (!is.null(tr)) {
      if (is.character(tr)) {
        return(paste0(tr, collapse = "\n"))
      }
      # best-effort stringification for arbitrary trace objects
      return(paste(capture.output(str(tr, max.level = 2)), collapse = "\n"))
    }

    # fallback to a deparsed call if present
    cl <- tryCatch(conditionCall(cnd), error = function(e) NULL)
    if (!is.null(cl)) {
      return(paste(deparse(cl), collapse = ""))
    }
    NA_character_
  }

  tibble::tibble(
    index = seq_along(conds),
    class = vapply(
      conds,
      function(x) paste(class(get_cnd(x)), collapse = "|"),
      character(1)
    ),
    message = vapply(
      conds,
      function(x) {
        cnd <- get_cnd(x)
        msg <- tryCatch(conditionMessage(cnd), error = function(e) NULL)
        if (is.null(msg) || length(msg) == 0) {
          NA_character_
        } else {
          paste(msg, collapse = "\n")
        }
      },
      character(1)
    ),
    call = vapply(
      conds,
      function(x) {
        cl <- tryCatch(conditionCall(get_cnd(x)), error = function(e) NULL)
        if (is.null(cl)) NA_character_ else paste(deparse(cl), collapse = "")
      },
      character(1)
    ),
    trace = vapply(conds, function(x) extract_trace(get_cnd(x)), character(1)),
    is_error = vapply(
      conds,
      function(x) inherits(get_cnd(x), "error"),
      logical(1)
    ),
    is_warning = vapply(
      conds,
      function(x) inherits(get_cnd(x), "warning"),
      logical(1)
    ),
    is_message = vapply(
      conds,
      function(x) inherits(get_cnd(x), "message"),
      logical(1)
    ),
    raw = conds
  )
}

# Stop if any errors present; include call/trace in message
#' @noRd
stop_on_future_errors <- function(
  fut,
  collapse = "\n\n---\n\n",
  show_trace = TRUE
) {
  df <- parse_future_conditions(fut)
  if (nrow(df) == 0) {
    return(invisible(TRUE))
  }
  errs <- df[df$is_error, , drop = FALSE]
  if (nrow(errs) == 0) {
    return(invisible(TRUE))
  }

  msgs <- vapply(
    seq_len(nrow(errs)),
    function(i) {
      parts <- c(paste0(
        "Future error [",
        errs$index[i],
        "]: ",
        errs$message[i]
      ))
      if (!is.na(errs$call[i]) && nzchar(errs$call[i])) {
        parts <- c(parts, paste0("  call: ", errs$call[i]))
      }
      if (show_trace && !is.na(errs$trace[i]) && nzchar(errs$trace[i])) {
        parts <- c(parts, "  trace:", errs$trace[i])
      }
      paste(parts, collapse = "\n")
    },
    character(1)
  )

  stop(paste(msgs, collapse = collapse), call. = FALSE)
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

.rem_na_helper <- function(...) {
  args <- list(...)
  nms <- match.call(expand.dots = FALSE)$`...`
  names(args) <- nms
  out <- tibble::as_tibble(args)
  out <- out[!apply(out, 1, function(row) any(is.na(row))), , drop = FALSE]
  out <- out[out$cnt != 1, ]
  as.list(out)
}

#' Combine unweighted means across chunks (variance input)
#' @noRd
combine_lumped_mean <- function(sub_mean, sub_n, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_n, cnt)
  sub_mean <- cln$sub_mean
  sub_n <- cln$sub_n
  sum(sub_mean * sub_n) / sum(sub_n)
}

#' Combine unweighted population variances across chunks (denominator N, variance input)
#' @noRd
combine_lumped_var <- function(sub_mean, sub_var, sub_n, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_n, cnt)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_n <- cln$sub_n

  mu_pooled <- sum(sub_mean * sub_n) / sum(sub_n)
  sum(sub_n * (sub_var + (sub_mean - mu_pooled)^2)) / sum(sub_n)
}

#' Combine unweighted population SD across chunks (denominator N, variance input)
#' @noRd
combine_lumped_sd <- function(sub_mean, sub_var, sub_n, cnt) {
  sqrt(combine_lumped_var(sub_mean, sub_var, sub_n, cnt))
}

#' Combine unweighted sample variances across chunks (denominator N-1, variance input)
#' @noRd
combine_lumped_sample_var <- function(sub_mean, sub_var, sub_n, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_n, cnt)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_n <- cln$sub_n

  mu_pooled <- sum(sub_mean * sub_n) / sum(sub_n)
  # unbiased sample variance for each chunk
  sub_sample_var <- sub_var * sub_n / (sub_n - 1)
  numerator <- sum((sub_n - 1) * sub_sample_var) +
    sum(sub_n * (sub_mean - mu_pooled)^2)
  denominator <- sum(sub_n) - 1
  numerator / denominator
}

#' Combine unweighted sample SD across chunks (denominator N-1, variance input)
#' @noRd
combine_lumped_sample_sd <- function(sub_mean, sub_var, sub_n, cnt) {
  sqrt(combine_lumped_sample_var(sub_mean, sub_var, sub_n, cnt))
}

#' Combine weighted means across chunks (variance input)
#' @noRd
combine_weighted_mean <- function(sub_mean, sub_wt, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_wt, cnt)
  sub_mean <- cln$sub_mean
  sub_wt <- cln$sub_wt
  sum(sub_mean * sub_wt) / sum(sub_wt)
}

#' Combine weighted population variances across chunks (denominator W, variance input)
#' @noRd
combine_weighted_var <- function(sub_mean, sub_var, sub_wt, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_wt, cnt)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_wt <- cln$sub_wt

  mu_pooled <- sum(sub_mean * sub_wt) / sum(sub_wt)
  sum(sub_wt * (sub_var + (sub_mean - mu_pooled)^2)) / sum(sub_wt)
}

#' Combine weighted population SD across chunks (denominator W, variance input)
#' @noRd
combine_weighted_sd <- function(sub_mean, sub_var, sub_wt, cnt) {
  sqrt(combine_weighted_var(sub_mean, sub_var, sub_wt, cnt))
}

#' Combine weighted sample variances across chunks (unbiased, denominator W-1, variance input)
#' @noRd
combine_weighted_sample_var <- function(sub_mean, sub_var, sub_wt, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_wt, cnt)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_wt <- cln$sub_wt

  mu_pooled <- sum(sub_mean * sub_wt) / sum(sub_wt)
  # unbiased sample variance for each chunk
  sub_sample_var <- sub_var * sub_wt / (sub_wt - 1)
  W <- sum(sub_wt)
  numerator <- sum((sub_wt - 1) * sub_sample_var) +
    sum(sub_wt * (sub_mean - mu_pooled)^2)
  denominator <- W - 1
  numerator / denominator
}

#' Combine weighted sample SD across chunks (unbiased, denominator W-1, variance input)
#' @noRd
combine_weighted_sample_sd <- function(sub_mean, sub_var, sub_wt, cnt) {
  sqrt(combine_weighted_sample_var(sub_mean, sub_var, sub_wt, cnt))
}


