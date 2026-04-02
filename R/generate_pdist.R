#' Generate pairwise distances
#'
#' Computes flow-connected and flow-unconnected pairwise distances between
#' stream links, including proportion of shared catchment area.
#'
#' @details
#' This function computes pairwise distances between stream links based on the flow paths traced by [trace_flowpaths()]. It generates two tables in the GeoPackage database: `fcon_pwise_dist` for flow-connected pairs of links, and `funcon_pwise_dist` for flow-unconnected pairs. For flow-connected pairs, it calculates the directed path length along the stream network, as well as the proportion of shared catchment area between the origin and destination links. For flow-unconnected pairs, it calculates the undirected path length as the sum of the directed path lengths from each link to their nearest common downstream link. The function is designed to be efficient and can handle large networks, but the runtime may increase with the size and complexity of the stream network. It is recommended to run this function on a machine with sufficient resources and to use the `verbose` option to monitor progress.
#'
#' @param input An `ihydro` object from [trace_flowpaths()].
#' @param pwise_all_links Logical. Calculate all flow un-connected pairwise distances?
#' @param verbose Logical.
#'
#' @return The input `ihydro` object with distance tables added.
#' The contained tables in the output GeoPackage are:
#' - `funcon_pwise_dist`: Pairwise flow-unconnected distances between links (if `pwise_dist = TRUE`).
#' - `fcon_pwise_dist`: Pairwise flow-connected distances between links (if `pwise_dist = TRUE`).
#' The function returns an `ihydro` object containing references to the flowpath tables, which can be accessed using [read_ihydro()] or listed with [ihydro_layers()].
#'
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
#' ihydro_obj <- ihydro::process_flowdir(
#'    dem = ex_dem,
#'    threshold = 1000,
#'    burn_streams = ex_streams,
#'    burn_depth = 5,
#'    burn_slope_dist = 250,
#'    burn_slope_depth = 5,
#'    min_length = 3,
#'    depression_corr = "breach",
#'    output_filename = output_gpkg,
#'    return_products = TRUE,
#'    verbose = TRUE
#' )
#'
#' ihydro::generate_vectors(
#'    input = ihydro_obj,
#'    points = ex_points,
#'    site_id_col = "site_id",
#'    snap_distance = 150,
#'    break_on_noSnap = FALSE,
#'    return_products = TRUE,
#'    verbose = TRUE
#' )
#'
#' ihydro::trace_flowpaths(
#'   input = ihydro_obj,
#'   pwise_dist = TRUE,
#'   pwise_all_links = FALSE,
#'   return_products = TRUE,
#'   verbose = TRUE
#' )
#'
#' ihydro::generate_pdist(
#'   input = ihydro_obj,
#'   pwise_all_links = TRUE,
#'   return_products = TRUE,
#'   verbose = TRUE
#' )
#'
#' ihydro::ihydro_layers(ihydro_obj, n = Inf)
#'
#' pwise_dist <- dplyr::bind_rows(
#'   ihydro::read_ihydro(hydro_out,"fcon_pwise_dist") |>
#'     dplyr::mutate(dist_type = "Flow Connected"),
#'   ihydro::read_ihydro(hydro_out,"funcon_pwise_dist") |>
#'     dplyr::mutate(dist_type = "Flow Unconnected")
#' )
#'
#' lines_out <- ihydro::read_ihydro(hydro_out,"stream_lines") |>
#'   dplyr::mutate(link_id = as.character(link_id)) |>
#'   dplyr::left_join(
#'     pwise_dist |>
#'       dplyr::filter(dist_type == "Flow Connected") |>
#'       dplyr::filter(origin == "200"),
#'     by = c("link_id" = "destination")
#' )
#'
#' mapview::mapview(lines_out, zcol = "directed_path_length", layer.name = "Flow-connected distance to link 200")
#'
#' }
#'

generate_pdist <- function(
  input,
  verbose = FALSE,
  return_products = FALSE,
  pwise_all_links = FALSE
) {
  check_ihydro(input)

  stopifnot(is.logical(verbose), is.logical(pwise_all_links))

  if (pwise_all_links) {
    message("pwise_all_links = TRUE can be very slow for large datasets.")
  }

  db_fp <- input$outfile

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  site_id_col <- read_site_id_col(db_fp)
  if (pwise_all_links) {
    site_id_col <- "link_id"
  }

  # ── Flow-connected distances ────────────────────────────────────────────
  DBI::dbExecute(con, "DROP TABLE IF EXISTS fcon_pwise_dist")

  if (verbose) {
    message("Calculating flow-connected distances")
  }

  # Pre-compute catchment areas per pour point
  us_catchment <- dplyr::tbl(con, "us_flowpaths") |>
    dplyr::left_join(
      dplyr::tbl(con, "stream_links_attr") |>
        dplyr::select(link_id, sbbsn_area),
      by = c("origin_link_id" = "link_id")
    ) |>
    dplyr::group_by(pour_point_id) |>
    dplyr::summarize(catchment_area = sum(sbbsn_area, na.rm = TRUE)) |>
    dplyr::ungroup()

  dplyr::tbl(con, "stream_links_attr") |>
    dplyr::filter(!is.na(dbplyr::sql(site_id_col))) |>
    dplyr::select(link_id) |>
    dplyr::rename(origin_link_id = link_id) |>
    dplyr::left_join(
      dplyr::tbl(con, "ds_flowpaths"),
      by = "origin_link_id"
    ) |>
    dplyr::left_join(
      dplyr::tbl(con, "stream_links_attr") |>
        dplyr::select(link_id, link_lngth, USChnLn_Fr),
      by = c("destination_link_id" = "link_id")
    ) |>
    dplyr::rename(
      origin = origin_link_id,
      destination = destination_link_id
    ) |>
    dplyr::group_by(origin) |>
    dbplyr::window_order(USChnLn_Fr) |>
    dplyr::mutate(directed_path_length = cumsum(link_lngth)) |>
    dplyr::ungroup() |>
    dplyr::select(origin, destination, directed_path_length) |>
    dplyr::distinct() |>
    dplyr::left_join(
      us_catchment |>
        dplyr::rename(
          origin = pour_point_id,
          origin_catchment = catchment_area
        ),
      by = "origin"
    ) |>
    dplyr::left_join(
      us_catchment |>
        dplyr::rename(
          destination = pour_point_id,
          destination_catchment = catchment_area
        ),
      by = "destination"
    ) |>
    dplyr::mutate(
      prop_shared_catchment = dplyr::case_when(
        directed_path_length > 0 ~
          as.numeric(origin_catchment) / as.numeric(destination_catchment),
        TRUE ~ 0
      ),
      prop_shared_logcatchment = dplyr::case_when(
        directed_path_length > 0 ~
          log(as.numeric(origin_catchment)) /
          log(as.numeric(destination_catchment)),
        TRUE ~ 0
      ),
      undirected_path_length = directed_path_length
    ) |>
    dplyr::select(-origin_catchment, -destination_catchment) |>
    dplyr::compute(
      name = "fcon_pwise_dist",
      temporary = FALSE,
      overwrite = TRUE,
      indexes = c("origin", "destination")
    )

  # ── Flow-unconnected distances ──────────────────────────────────────────
  DBI::dbExecute(con, "DROP TABLE IF EXISTS funcon_pwise_dist")

  if (verbose) {
    message("Calculating flow-unconnected distances")
  }

  fcon <- dplyr::tbl(con, "fcon_pwise_dist")
  ds <- dplyr::tbl(con, "ds_flowpaths")

  dplyr::full_join(
    fcon |>
      dplyr::select(
        midpoint = destination,
        origin_p1 = origin,
        directed_path_length_p1 = directed_path_length
      ),
    fcon |>
      dplyr::select(
        midpoint = destination,
        origin_p2 = origin,
        directed_path_length_p2 = directed_path_length
      ),
    by = "midpoint"
  ) |>
    dplyr::filter(origin_p1 != origin_p2) |>
    dplyr::anti_join(
      ds |>
        dplyr::select(
          origin_p2 = destination_link_id,
          origin_p1 = origin_link_id
        ),
      by = c("origin_p1", "origin_p2")
    ) |>
    dplyr::anti_join(
      ds |>
        dplyr::select(
          origin_p1 = destination_link_id,
          origin_p2 = origin_link_id
        ),
      by = c("origin_p1", "origin_p2")
    ) |>
    dplyr::mutate(
      undirected_path_length = as.numeric(directed_path_length_p1) +
        as.numeric(directed_path_length_p2)
    ) |>
    dplyr::select(
      origin = origin_p1,
      destination = origin_p2,
      undirected_path_length
    ) |>
    dplyr::mutate(
      directed_path_length = 0,
      prop_shared_catchment = 0,
      prop_shared_logcatchment = 0
    ) |>
    dplyr::group_by(origin, destination) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ min(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    dplyr::compute(
      name = "funcon_pwise_dist",
      temporary = FALSE,
      overwrite = TRUE,
      indexes = c("origin", "destination")
    )

  if (return_products) {
    input$funcon_pwise_dist <- read_ihydro(input, "funcon_pwise_dist")
    input$fcon_pwise_dist <- read_ihydro(input, "fcon_pwise_dist")
  }

  structure(input, class = "ihydro")
}
