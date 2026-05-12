#' Convert iHydro stream data to SSN format
#'
#' This function converts stream network data from an iHydro object into a format
#' compatible with the Spatial Stream Network (SSN) modeling framework. It processes
#' the stream lines, prepares site data for observed and predicted locations, and assembles
#' the SSN object with the appropriate attributes for edges and sites. The function also
#' handles topology errors that may arise during the conversion process.
#'
#' @param ihydro_obj An iHydro object
#' @param obs A data frame of observed site data with a link_id column
#' @param pred A data frame of predicted site data with a link_id column
#' @param lsn_path Pathname to a directory in character format specifying where to store the outputs. The directory will be created if it does not already exist..
#' @param verbose Logical.
#'
#' @return An SSN object containing the stream network and site data
#' @export
#'

ihydro_to_ssn <- function(
    ihydro_obj,
    obs = NULL,
    pred = NULL,
    lsn_path = NULL,
    verbose = FALSE
) {
  .check_ihydro(ihydro_obj)

  n_cores <- .n_workers()
  site_list <- .prep_sites(ihydro_obj, obs, pred)

  if (is.null(lsn_path)) {
    lsn_path <- tempdir()
  }

  lsn.path <- file.path(lsn_path, "lsn_files")
  dir.create(lsn.path, showWarnings = FALSE)

  edges <- .prep_edges(
    ihydro_obj,
    n_cores = n_cores,
    lsn_path = lsn_path,
    verbose = verbose
  )

  if (!is.null(site_list$pred)) {
    site_list[["pred"]] <- SSNbler::sites_to_lsn(
      sites = site_list[["pred"]],
      edges = edges,
      lsn_path = lsn_path,
      file_name = "pred",
      snap_tolerance = Inf,
      save_local = TRUE,
      overwrite = TRUE,
      verbose = verbose
    )
  }

  if (!is.null(site_list$obs)) {
    site_list[["obs"]] <- SSNbler::sites_to_lsn(
      sites = site_list[["obs"]],
      edges = edges,
      lsn_path = lsn_path,
      file_name = "obs",
      snap_tolerance = Inf,
      save_local = TRUE,
      overwrite = TRUE,
      verbose = verbose
    )
  }

  snn_site_list <- SSNbler::updist_sites(
    sites = site_list,
    edges = edges,
    length_col = "Length",
    save_local = TRUE,
    lsn_path = lsn_path
  )

  snn_site_list <- SSNbler::afv_sites(
    sites = snn_site_list,
    edges = edges,
    afv_col = "afvArea",
    save_local = TRUE,
    lsn_path = lsn_path
  )

  # Assemble SSN object
  ssn_obj <- SSNbler::ssn_assemble(
    edges = edges,
    lsn_path = lsn_path,
    obs_sites = snn_site_list$obs,
    preds_list = snn_site_list[names(snn_site_list) != "obs"],
    ssn_path = paste0(lsn_path, "/ihydro.ssn"),
    import = TRUE,
    check = TRUE,
    afv_col = "afvArea",
    overwrite = TRUE,
    verbose = verbose
  )

  return(ssn_obj)
}


#' Convert iHydro stream data to lsn format
#' @param ihydro_obj An iHydro object
#' @param n_cores Number of cores to use for parallel processing
#' @param lsn_path Directory to save intermediate files
#' @param verbose Whether to print progress messages
#' @return A sf object of edges in lsn format
#' @noRd
.prep_edges <- function(
    ihydro_obj,
    n_cores = 1L,
    lsn_path = NULL,
    verbose = FALSE
) {
  .check_ihydro(ihydro_obj)

  stream_points <- .prep_ihydro_points(ihydro_obj)

  lns <- dplyr::mutate(
    read_ihydro(ihydro_obj, "stream_lines"),
    link_id = as.character(link_id)
  )
  lns0 <- lns <- dplyr::left_join(
    lns,
    stream_points,
    by = "link_id"
  )

  res <- terra::res(
    read_ihydro(ihydro_obj, "dem_final")
  )

  edges <- SSNbler::lines_to_lsn(
    streams = lns,
    lsn_path = lsn_path,
    check_topology = TRUE,
    snap_tolerance = sqrt(res[[1]])/8,
    topo_tolerance = sqrt(res[[1]])/8,
    overwrite = TRUE,
    use_parallel = n_cores > 1,
    no_cores = n_cores,
    verbose = verbose
  )

  err_fl <- file.path(lsn_path, "node_errors.gpkg")
  if (file.exists(err_fl)) {
    if (verbose) {
      message("Attempting to fix topology errors")
    }

    ers <- sf::read_sf(err_fl)

    ft <- file.remove(err_fl)

    ers_buff <- sf::st_buffer(ers, res[[1]] / 4)
    ers_lns <- suppressWarnings(sf::st_intersection(lns, ers_buff))
    ers_buff <- suppressWarnings(sf::st_cast(ers_buff, "LINESTRING"))
    ers_pnt <- suppressWarnings(sf::st_intersection(lns, ers_buff))

    parent_segs <- ers_pnt |>
      dplyr::group_by(pointid) |>
      tidyr::nest() |>
      dplyr::mutate(
        min_dist = purrr::map(
          data,
          ~ {
            dmat <- sf::st_distance(.x)
            dmat <- as.data.frame(as.table(dmat))
            dmat$Freq <- as.numeric(dmat$Freq)
            dmat <- dmat[dmat$Freq != 0, ]
            dmat <- dmat[which.min(dmat$Freq), ]
            out <- .x[as.numeric(as.character(c(dmat$Var1, dmat$Var2))), ]
            link_to_move <- out$link_id[rank(out$USChnLn_Fr)][1]
            coord_to_move_to <- sf::st_coordinates(
              out$geom[rank(out$USChnLn_Fr)][2]
            )
            list(
              link_to_move = link_to_move,
              coord_to_move_to = coord_to_move_to,
              link_to_recieve = out$link_id[rank(out$USChnLn_Fr)][2]
            )
          }
        )
      )

    for (i in parent_segs$min_dist) {
      trib_to_move <- lns[lns$link_id == i[[1]], ]
      coords <- sf::st_coordinates(trib_to_move)
      coords[nrow(coords), 1:2] <- unlist(i[[2]])
      coords <- sf::st_linestring(coords[, 1:2])
      lns$geom[lns$link_id == i[[1]]] <- coords

      parent_split <- tibble::as_tibble(lns[lns$link_id %in% i[[3]], ])
      parent_line <- sf::st_coordinates(parent_split$geom)[, 1:2]
      parent_line <- sf::st_linestring(parent_line)
      snap_point <- sf::st_point(i[[2]])
      parent_line <- sf::st_snap(parent_line, snap_point, sqrt(res[[1]])/16)
      parent_line <- sf::st_coordinates(parent_line)[, 1:2]

      parent_index <- apply(parent_line, 1, function(x) all(x == i[[2]]))
      parent_index <- which(parent_index)
      parent_line <- list(
        parent_line[1:parent_index, 1:2],
        parent_line[parent_index:nrow(parent_line), 1:2]
      )

      parent_line <- lapply(parent_line, sf::st_linestring)
      parent_line <- tibble::tibble(geom = sf::st_as_sfc(parent_line))

      parent_split <- dplyr::select(parent_split, -geom)
      parent_split <- dplyr::bind_rows(parent_split, parent_split)
      parent_split <- dplyr::bind_cols(
        parent_split,
        parent_line
      ) |>
        sf::st_as_sf()

      sf::st_crs(parent_split) <- sf::st_crs(lns)
      lns <- lns[!lns$link_id %in% i[[3]], ]
      lns <- dplyr::bind_rows(lns, parent_split)
    }

    #lns_inter <- sf::st_snap(lns, sf::st_union(lns), sqrt(res[[1]])/16)

    # Edges
    edges <- SSNbler::lines_to_lsn(
      streams = lns,
      lsn_path = lsn_path,
      check_topology = TRUE,
      snap_tolerance = sqrt(res[[1]])/8,
      topo_tolerance = sqrt(res[[1]])/8,
      overwrite = TRUE,
      use_parallel = n_cores > 1,
      no_cores = n_cores,
      verbose = verbose
    )

    if (file.exists(err_fl)) {
      ers2 <- sf::read_sf(err_fl)
      ers2 <- tibble::as_tibble(ers2)

      if (nrow(ers2) > 0) {
        rlang::abort(
          message = "The following topology errors could not be resolved:",
          body = c(ers2 |> capture.output())
        )
      }
    }
  }

  edges <- SSNbler::updist_edges(
    edges = edges,
    save_local = TRUE,
    lsn_path = lsn_path,
    calc_length = TRUE,
    verbose = verbose
  )

  edges <- SSNbler::afv_edges(
    edges = edges,
    infl_col = "area",
    segpi_col = "areaPI",
    afv_col = "afvArea",
    lsn_path = lsn_path
  )

  return(edges)

}

#' Prepare site data for conversion to lsn format
#' @param ihydro_obj An iHydro object
#' @param obs A data frame of observed site data with a link_id column
#' @param pred A data frame of predicted site data with a link_id column
#' @return A list of data frames for observed and predicted sites, with stream link attributes joined
#' @noRd
.prep_sites <- function(
    ihydro_obj,
    obs = NULL,
    pred = NULL
) {
  site_list <- list(
    obs = obs,
    pred = pred
  )
  if (is.null(obs) & is.null(pred)) {
    return(site_list)
  }

  sites <- read_ihydro(ihydro_obj, "stream_links")

  if (!is.null(obs)) {
    if (inherits(obs, "sf")) {
      cli::cli_abort("{.arg obs} must not be a class {.cls sf}")
    }
    if (!inherits(obs, "data.frame")) {
      cli::cli_abort("{.arg obs} must be a {.cls data.frame}")
    }
    if (!"link_id" %in% colnames(obs)) {
      cli::cli_abort("{.arg obs} must contain a {.code link_id} column")
    }
    if (any(!obs$link_id %in% sites$link_id)) {
      missing_links <- obs$link_id[!obs$link_id %in% sites$link_id]
      cli::cli_abort(c(
        "{.arg obs} contains {.code link_id}s that are not in the iHydro object:",
        "x" = "{.code link_id}s: {.val {missing_links}}"
      ))
    }

    obs$link_id <- as.character(obs$link_id)

    site_list$obs <- dplyr::right_join(
      sites,
      obs,
      by = "link_id"
    )
  }
  if (!is.null(pred)) {
    if (inherits(pred, "sf")) {
      cli::cli_abort("{.arg obs} must not be a class {.cls sf}")
    }
    if (!inherits(pred, "data.frame")) {
      cli::cli_abort("{.arg pred} must be a {.cls data.frame}")
    }
    if (!"link_id" %in% colnames(pred)) {
      cli::cli_abort("{.arg pred} must contain a {.code link_id} column")
    }
    if (any(!pred$link_id %in% sites$link_id)) {
      missing_links <- pred$link_id[!pred$link_id %in% sites$link_id]
      cli::cli_abort(c(
        "{.arg pred} contains {.code link_id}s that are not in the iHydro object:",
        "x" = "{.code link_id}s: {.val {missing_links}}"
      ))
    }

    pred$link_id <- as.character(pred$link_id)

    site_list$pred <- dplyr::right_join(
      sites,
      pred,
      by = "link_id"
    )
  }

  return(site_list)
}

#' Prepare stream points for conversion to lsn format
#' @param ihydro_obj An iHydro object
#' @param keep_site_id_col Whether to keep the site ID column in the output
#' @param keep_dir_cols Whether to keep the directional columns (e.g. ustrib_, uslink_, dstrib_, dslink_) in the output
#' @return A data frame of stream points with attributes joined from the iHydro object
#' @noRd
.prep_ihydro_points <- function(
    ihydro_obj,
    keep_site_id_col = TRUE,
    keep_dir_cols = TRUE
) {
  site_id_col <- read_ihydro(ihydro_obj, "site_id_col")[[1]]

  sites <- read_ihydro(ihydro_obj, "stream_links_attr") |>
    dplyr::filter(link_lngth > 0)

  subb_lookup <- ihydro::read_ihydro(
    ihydro_obj,
    "us_flowpaths"
  )

  stream_points <- dplyr::left_join(
    subb_lookup,
    sites,
    by = c("origin_link_id" = "link_id")
  ) |>
    dplyr::group_by(pour_point_id) |>
    dplyr::summarise(
      area = sum(sbbsn_area, na.rm = TRUE)
    )

  stream_points <- dplyr::rename(
    stream_points,
    link_id = pour_point_id
  )

  stream_points <- dplyr::left_join(
    stream_points,
    sites,
    by = "link_id"
  )

  stream_points <- dplyr::select(
    stream_points,
    link_id,
    tidyselect::any_of(site_id_col),
    tidyselect::everything(),
    -ID,
    -link_type,
    -link_class,
  )

  if (!keep_site_id_col) {
    stream_points <- dplyr::select(
      stream_points,
      -tidyselect::any_of(site_id_col)
    )
  }

  if (!keep_dir_cols) {
    stream_points <- dplyr::select(
      stream_points,
      -tidyselect::starts_with(c("ustrib_", "uslink_", "dstrib_", "dslink_"))
    )
  }

  return(stream_points)
}
