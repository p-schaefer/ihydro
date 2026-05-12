#' Generate all hydrology products from a DEM
#'
#' Convenience wrapper that calls [process_flowdir()], [generate_vectors()],
#' and [trace_flowpaths()] in sequence.
#'
#' @inheritParams process_flowdir
#' @inheritParams generate_vectors
#' @inheritParams trace_flowpaths
#'
#' @return An `ihydro` object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
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
#' }
#'

process_hydrology <- function(
  dem,
  output_filename,
  threshold,
  burn_streams = NULL,
  burn_depth = NULL,
  burn_slope_dist = NULL,
  burn_slope_depth = NULL,
  min_length = NULL,
  depression_corr = c(NULL, "fill", "breach"),
  extra_attr = c(
    "link_slope",
    "cont_slope",
    "USChnLn_To",
    "Elevation",
    "StOrd_Hack",
    "StOrd_Str",
    "StOrd_Hort",
    "StOrd_Shr"
  ),
  points = NULL,
  snap_distance = NULL,
  break_on_noSnap = TRUE,
  site_id_col = NULL,
  pwise_dist = FALSE,
  pwise_all_links = FALSE,
  return_products = FALSE,
  temp_dir = NULL,
  compress = FALSE,
  verbose = FALSE
) {
  # ── Validation ──────────────────────────────────────────────────────────
  if (!is.integer(threshold)) {
    if (as.integer(threshold) != threshold) {
      cli::cli_alert_info("{.arg threshold} was converted to an integer.")
    }
    threshold <- as.integer(threshold)
  }
  if (
    !is.null(burn_streams) && (is.null(burn_depth) && is.null(burn_slope_depth))
  ) {
    cli::cli_abort(
      "{.arg burn_depth} and/or {.arg burn_slope_depth} must be provided when {.arg burn_streams} is present"
    )
  }
  if (!is.null(burn_streams)) {
    if (!is.null(burn_depth) && !is.numeric(burn_depth)) {
      cli::cli_abort("{.arg burn_depth} must be numeric.")
    }
    if (!is.null(burn_slope_dist) && !is.numeric(burn_slope_dist)) {
      cli::cli_abort("{.arg burn_slope_dist} must be numeric.")
    }
    if (!is.null(burn_slope_depth) && !is.numeric(burn_slope_depth)) {
      cli::cli_abort("{.arg burn_slope_depth} must be numeric.")
    }
  }

  if (!is.null(min_length) && !is.numeric(min_length)) {
    cli::cli_abort("{.arg min_length} must be numeric.")
  }
  stopifnot(
    is.logical(return_products),
    is.logical(compress),
    is.logical(verbose)
  )
  burn_depth <- as.integer(burn_depth)

  stopifnot(
    is.logical(pwise_dist),
    is.logical(return_products),
    is.logical(compress),
    is.logical(break_on_noSnap),
    is.logical(verbose)
  )

  if (!is.null(points)) {
    points <- sf::st_as_sf(process_input(points))
    validate_site_id_col(site_id_col)
    if (!site_id_col %in% names(points)) {
      cli::cli_abort("{.arg site_id_col} must be a column in {.arg points}.")
    }
  }

  temp_dir <- .ensure_temp_dir(temp_dir)
  output_filename <- normalizePath(output_filename, mustWork = FALSE)
  if (!grepl("\\.gpkg$", output_filename)) {
    cli::cli_abort("{.arg output_filename} must end in {.val .gpkg}.")
  }
  extra_attr <- match.arg(extra_attr, several.ok = TRUE)

  # ── Pipeline ────────────────────────────────────────────────────────────
  if (verbose) {
    message("Processing flow direction")
  }
  hydro_out <- process_flowdir(
    dem = dem,
    threshold = threshold,
    burn_streams = burn_streams,
    burn_depth = burn_depth,
    min_length = min_length,
    depression_corr = depression_corr,
    return_products = return_products,
    output_filename = output_filename,
    temp_dir = temp_dir,
    compress = compress,
    verbose = verbose
  )

  if (verbose) {
    message("Generating vectors")
  }
  hydro_out <- generate_vectors(
    input = hydro_out,
    extra_attr = extra_attr,
    points = points,
    site_id_col = site_id_col,
    snap_distance = snap_distance,
    break_on_noSnap = break_on_noSnap,
    return_products = return_products,
    temp_dir = temp_dir,
    compress = compress,
    verbose = verbose
  )

  if (verbose) {
    message("Tracing flowpaths")
  }
  hydro_out <- trace_flowpaths(
    input = hydro_out,
    return_products = return_products,
    temp_dir = temp_dir,
    pwise_dist = pwise_dist,
    pwise_all_links = pwise_all_links,
    verbose = verbose
  )

  hydro_out
}
