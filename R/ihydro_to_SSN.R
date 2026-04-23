.ihydro_to_lsn <- function(
    ihydro_obj,
    pred = NULL,
    obs = NULL,
    temp_dir = NULL,
    verbose = FALSE
) {
  check_ihydro(ihydro_obj)
  if (!is.null(pred)) {
    if (!inherits(pred, "data.frame") | "link_id" %in% colnames(pred)) {
      cli::cli_abort("{.arg pred} must be a dataframe with a `link_id` column")
    }
  }
  if (!is.null(obs)) {
    if (!inherits(obs, "data.frame") | "link_id" %in% colnames(obs)) {
      cli::cli_abort("{.arg obs} must be a dataframe with a `link_id` column")
    }
  }

  n_cores <- n_workers()
  on.exit(options(parallelly.maxWorkers.localhost = max_cores_opt), add = TRUE)

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
    #-tidyselect::starts_with(c("ustrib_", "uslink_", "dstrib_", "dslink_"))
  )

  lns <- dplyr::mutate(
    read_ihydro(ihydro_obj, "stream_lines"),
    link_id = as.character(link_id)
  )
  lns0 <- lns <- dplyr::left_join(
    lns,
    stream_points,
    by = "link_id"
  )

  if (is.null(temp_dir)) {
    temp_dir <- tempdir()
  }

  lsn.path <- file.path(temp_dir, "lsn_files")
  dir.create(lsn.path, showWarnings = FALSE)

  res <- terra::res(
    read_ihydro(ihydro_obj, "dem_final")
  )

  edges <- SSNbler::lines_to_lsn(
    streams = lns,
    lsn_path = lsn.path,
    check_topology = TRUE,
    snap_tolerance = sqrt(res[[1]]),
    topo_tolerance = sqrt(res[[1]]),
    overwrite = TRUE,
    use_parallel = n_cores > 1,
    no_cores = n_cores,
    verbose = verbose
  )

  err_fl <- file.path(lsn.path,"node_errors.gpkg")
  while (file.exists(err_fl)) {
    if (verbose){
      message("Attempting to fix topology errors")
    }

    ers <- sf::read_sf(err_fl)

    ft <- file.remove(err_fl)

    ers_buff <- sf::st_buffer(ers,res[[1]])
    ers_lns <- suppressWarnings(sf::st_intersection(lns,ers_buff))
    ers_buff <- suppressWarnings(sf::st_cast(ers_buff,"LINESTRING"))
    ers_pnt <- suppressWarnings(sf::st_intersection(lns,ers_buff))

    parent_segs <- ers_pnt |>
      dplyr::group_by(pointid) |>
      tidyr::nest() |>
      dplyr::mutate(
        min_dist = purrr::map(data,~{
          dmat <- sf::st_distance(.x)
          dmat <- as.data.frame(as.table(dmat))
          dmat$Freq <- as.numeric(dmat$Freq)
          dmat <- dmat[dmat$Freq != 0,]
          dmat <- dmat[which.min(dmat$Freq),]
          out <- .x[as.numeric(as.character(c(dmat$Var1,dmat$Var2))),]
          link_to_move <- out$link_id[rank(out$USChnLn_Fr)][1]
          coord_to_move_to <- sf::st_coordinates(out$geom[rank(out$USChnLn_Fr)][2])
          list(
            link_to_move = link_to_move,
            coord_to_move_to = coord_to_move_to,
            link_to_recieve = out$link_id[rank(out$USChnLn_Fr)][2]
          )
        })
      )

    for (i in parent_segs$min_dist) {
      trib_to_move <- lns[lns$link_id == i[[1]],]
      coords <- sf::st_coordinates(trib_to_move)
      coords[nrow(coords),1:2] <- unlist(i[[2]])
      coords <- sf::st_linestring(coords[,1:2])
      lns$geom[lns$link_id == i[[1]]] <- coords

      parent_split <- tibble::as_tibble(lns[lns$link_id %in% i[[3]],])
      parent_line <- sf::st_coordinates(parent_split$geom)[,1:2]
      parent_line <-  sf::st_linestring(parent_line)
      snap_point <-  sf::st_point(i[[2]])
      parent_line <- sf::st_snap(parent_line,snap_point,sqrt(res[[1]]))

      parent_index <- apply(parent_line,1,function(x) all(x==i[[2]]))
      parent_index <- which(parent_index)
      parent_line <- list(
        parent_line[1:parent_index,1:2],
        parent_line[parent_index:nrow(parent_line),1:2]
      )

      parent_line <- lapply(parent_line, sf::st_linestring)
      parent_line <- tibble::tibble(geom = sf::st_as_sfc(parent_line))

      parent_split <- dplyr::select(parent_split, - geom)
      parent_split <- dplyr::bind_rows(parent_split,parent_split)
      parent_split <- dplyr::bind_cols(
        parent_split,
        parent_line
      ) |>
        sf::st_as_sf()

      sf::st_crs(parent_split) <- sf::st_crs(lns)
      lns <- lns[!lns$link_id %in% i[[3]],]
      lns <- dplyr::bind_rows(lns,parent_split)
    }

    lns_inter <- sf::st_snap(lns, sf::st_union(lns), sqrt(res))

    # Edges
    edges <- SSNbler::lines_to_lsn(
      streams = lns_inter,
      lsn_path = lsn.path,
      check_topology = TRUE,
      snap_tolerance = sqrt(res[[1]]),
      topo_tolerance = sqrt(res[[1]]),
      overwrite = TRUE,
      use_parallel = n_cores > 1,
      no_cores = n_cores,
      verbose = verbose
    )
  }

  edges <- SSNbler::updist_edges(
    edges = edges,
    save_local = TRUE,
    lsn_path = lsn.path,
    calc_length = TRUE,
    verbose = verbose
  )

  edges <- SSNbler::afv_edges(
    edges = edges,
    infl_col = "area",
    segpi_col = "areaPI",
    afv_col = "afvArea",
    lsn_path = lsn.path
  )

  return(edges)

  #
  #
  #
  #
  # # Points
  # site_list <- list(
  #   pred = NULL,
  #   obs = NULL
  # )
  # site_list[["pred"]] <- SSNbler::sites_to_lsn(
  #   sites = pts,
  #   edges = edges,
  #   lsn_path = lsn.path,
  #   file_name = "pred",
  #   snap_tolerance = Inf,
  #   save_local = TRUE,
  #   overwrite = TRUE,
  #   verbose = verbose
  # )
  #
  # if (!is.null(obs_obj)) {
  #   obs_obj <- dplyr::left_join(
  #     obs_obj,
  #     dplyr::select(
  #       stream_points,
  #       link_id,
  #       tidyselect::any_of(
  #         colnames(stream_points)[
  #           !colnames(stream_points) %in% colnames(obs_obj)
  #         ]
  #       )
  #     ),
  #     by = "link_id"
  #   )
  #
  #   site_list[["obs"]] <- SSNbler::sites_to_lsn(
  #     sites = obs_obj,
  #     edges = edges,
  #     lsn_path = lsn.path,
  #     file_name = "pred",
  #     snap_tolerance = Inf,
  #     save_local = TRUE,
  #     overwrite = TRUE,
  #     verbose = verbose
  #   )
  # }
  #
  # site.list <- SSNbler::updist_sites(
  #   sites = site_list,
  #   edges = edges,
  #   length_col = "Length",
  #   save_local = TRUE,
  #   lsn_path = lsn.path
  # )
  #
  #
  # site.list <- SSNbler::afv_sites(
  #   sites = site.list,
  #   edges = edges,
  #   afv_col = "afvArea",
  #   save_local = TRUE,
  #   lsn_path = lsn.path
  # )
  #
  # # Assemble SSN object
  # ssn_obj <- SSNbler::ssn_assemble(
  #   edges = edges,
  #   lsn_path = lsn.path,
  #   obs_sites = site.list$obs,
  #   preds_list = site.list[!"obs" %in% names(site.list)],
  #   ssn_path = paste0(lsn.path, "/ihydro.ssn"),
  #   import = TRUE,
  #   check = TRUE,
  #   afv_col = "afvArea",
  #   overwrite = TRUE,
  #   verbose = verbose
  # )
  #
  # ssn_obj
}
