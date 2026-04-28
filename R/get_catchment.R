#' Generate upstream catchment areas
#'
#' Retrieves complete upstream catchment polygons for specified reaches or
#' sampling points. Results are cached in the `Catchment_poly` layer.
#'
#' @param input An `ihydro` object.
#' @param sample_points Character vector of site IDs (must exist in the `site_id_col` column provided in [generate_vectors()]), or `NULL`.
#' @param link_id Character vector of link IDs, or `NULL`.
#' @param temp_dir Temporary file directory.
#' @param verbose Logical.
#'
#' @return An `sf` polygon object of upstream catchments.
#'
#' @export
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
#' catch_poly <- ihydro::get_catchment(
#'   input = ihydro_obj,
#'   sample_points = c("101.1"),
#'   link_id = c("800"),
#' )
#'
#' mapview::mapview(catch_poly) +
#'   mapview::mapview(ihydro::read_ihydro(ihydro_obj, "streams"))
#'
#' }
#'

get_catchment <- function(
  input,
  sample_points = NULL,
  link_id = NULL,
  temp_dir = NULL,
  verbose = FALSE,
  return = TRUE
) {
  check_ihydro(input)
  whitebox::wbt_options(exe_path = whitebox::wbt_exe_path(), verbose = FALSE)

  n_cores <- n_workers()
  max_cores_opt <- getOption("parallelly.maxWorkers.localhost")
  on.exit(options(parallelly.maxWorkers.localhost = max_cores_opt), add = TRUE)
  options(parallelly.maxWorkers.localhost = n_cores)

  db_fp <- input$outfile
  temp_dir <- ensure_temp_dir(temp_dir)

  ihydro_tbl <- ihydro_layers(input)

  # Check for existing catchments
  existing_catch <- tibble::tibble(link_id = character(0))
  site_id_col <- read_site_id_col(db_fp)

  if ("Catchment_poly" %in% ihydro_tbl$layer_name) {
    existing_catch <- read_ihydro(input, "Catchment_poly") |>
      tibble::as_tibble() |>
      dplyr::select(link_id)
  }

  # Resolve targets
  target_ids <- target_id_fun(
    db_fp = db_fp,
    sample_points = sample_points,
    link_id = link_id
  )

  already_done <- dplyr::filter(
    target_ids,
    link_id %in% existing_catch$link_id
  )
  still_needed <- dplyr::filter(
    target_ids,
    !link_id %in% existing_catch$link_id
  )

  # If all already cached, return from gpkg
  if (nrow(still_needed) == 0) {
    if (!return) {
      return(invisible(NULL))
    }

    return(sf::read_sf(
      db_fp,
      query = build_sql_in("Catchment_poly", "link_id", already_done$link_id)
    ))
  }

  if (verbose) {
    message(
      "Generating catchments for ",
      nrow(still_needed),
      " targets (",
      nrow(already_done),
      " already cached)"
    )
  }
  # ── Compute missing catchments in parallel ──────────────────────────────
  progressr::with_progress(enable = verbose, {
    p <- progressr::progressor(steps = nrow(still_needed))

    out <- still_needed |>
      dplyr::mutate(
        db_loc = list(db_fp),
        p = list(p)
      ) |>
      dplyr::mutate(
        catch = furrr::future_pmap(
          list(link_id = link_id, db_loc = db_loc, p = p),
          .options = furrr::furrr_options(
            globals = FALSE,
            seed = NULL,
            scheduling = 4L
          ),
          compute_single_catchment
        )
      )

    # out <- out |>
    #   dplyr::filter(!is.na(link_id)) |>
    #   dplyr::select(catch) |>
    #   tidyr::unnest(catch) |>
    #   dplyr::filter(!is.na(link_id))
    #
    # out <- terra::vect(out)
    #
    # out <- lapply(1:nrow(out), function(x) terra::rasterize(out[x,],r,field = "link_id"))
    # out <- lapply(out,function(x) terra::as.polygons(x))
    #
    # out <- sf::st_as_sf(out)
    #
  })

  out <- out |>
    dplyr::filter(!is.na(link_id)) |>
    dplyr::select(catch) |>
    tidyr::unnest(catch) |>
    dplyr::filter(!is.na(link_id))

  # Cache results
  sf::write_sf(
    out,
    db_fp,
    layer = "Catchment_poly",
    delete_dsn = FALSE,
    delete_layer = FALSE,
    append = TRUE
  )

  # Combine with previously cached

  if (nrow(already_done) > 0) {
    out <- dplyr::bind_rows(
      out,
      sf::read_sf(
        db_fp,
        query = build_sql_in(
          "Catchment_poly",
          "link_id",
          already_done$link_id
        )
      )
    )
  }

  if (!return) {
    return(invisible(NULL))
  }
  sf::st_as_sf(out)
}


# ── Worker for single-catchment computation ───────────────────────────────────

#' @noRd
compute_single_catchment <- carrier::crate(
  function(link_id, db_loc, p) {
    #suppressPackageStartupMessages(library(sf))
    #options(dplyr.summarise.inform = FALSE, scipen = 999)
    `%>%` <- magrittr::`%>%`

    con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)

    query <- dplyr::tbl(con, "us_flowpaths") %>%
      dplyr::filter(pour_point_id %in% local(link_id)) %>%
      dplyr::rename(link_id = origin_link_id) %>%
      dplyr::left_join(
        dplyr::tbl(con, "Subbasins_poly") %>% dplyr::select(link_id, geom),
        by = "link_id"
      ) %>%
      dplyr::show_query() %>%
      utils::capture.output() %>%
      utils::tail(-1) %>%
      paste(collapse = "")

    DBI::dbDisconnect(con)

    s1 <- suppressWarnings(sf::read_sf(db_loc, query = query))

    out_geom <- s1 %>%
      dplyr::select(-link_id, link_id = pour_point_id) %>%
      #sf::st_buffer(units::as_units(0.01, sf::st_crs(s1)$units)) %>%
      sf::st_union() %>%
      sfheaders::sf_remove_holes() %>%
      sf::st_cast("POLYGON")

    out_sf <- sf::st_sf(
      link_id = link_id,
      geom = out_geom,
      crs = sf::st_crs(s1)
    )

    p()

    return(out_sf[1, , drop = FALSE])
  }
)
