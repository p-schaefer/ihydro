
#' Generate upstream catchment areas for points
#'
#' Catchments are automatically added to Catchment_poly layer in the geopackage for easier future retrieval
#'
#' @param input output of `process_hydrology()` or one containing `generate_subbasins()` and `trace_flowpaths()`
#' @param sample_points character or NULL. IDs of unique station identifiers provided in 'site_id_col' of `generate_vectors()` to calculate catchments for.
#' @param link_id character or NULL. 'link_id's of reaches to calculate catchments for.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#'
#' @return sf polygon of upstream catchments
#'
#' @export


get_catchment<-function(
    input,
    sample_points=NULL,
    link_id=NULL,
    temp_dir=NULL
) {
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  options(scipen = 999)
  options(dplyr.summarise.inform = FALSE,future.rng.onMisuse="ignore")
  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),verbose=F)

  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1
  max_cores_opt<-getOption("parallelly.maxWorkers.localhost")
  options(parallelly.maxWorkers.localhost = n_cores)

  ihydro_tbl<-ihydro_layers(input)

  existing_catch<-tibble::tibble(link_id=NA_character_)[F,]

  db_loc<-db_fp<-zip_loc<-input$outfile

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)

  site_id_col<-dplyr::collect(dplyr::tbl(con,"site_id_col"))$site_id_col

  all_points<-dplyr::collect(dplyr::tbl(con,"stream_links_attr")) %>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),~dplyr::na_if(.,"")))

  if ("Catchment_poly" %in% ihydro_tbl$layer_name) {
    existing_catch<-dplyr::collect(dplyr::tbl(con,"Catchment_poly") %>% select(link_id))
  }

  DBI::dbDisconnect(con)

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  tdir<-file.path(temp_dir,basename(tempfile()))
  dir.create(tdir)

  target_IDs<-target_id_fun(
    db_fp=db_fp,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=F
  )


  existing_IDs<-target_IDs %>%
    filter(link_id %in% existing_catch$link_id)

  target_IDs<-target_IDs%>%
    filter(!link_id %in% existing_catch$link_id)

  if (nrow(target_IDs)==0) {
    return(sf::read_sf(
      db_loc,
      query=paste0("SELECT `link_id`, `geom` FROM `Catchment_poly` WHERE (`link_id` IN (",paste0("'",existing_IDs$link_id,collapse = "',"),"'))")
    ))
  }


  progressr::with_progress(enable=T,{
    p <- progressr::progressor(steps = nrow(target_IDs))

    out<-target_IDs %>%
      dplyr::mutate(db_loc=list(db_loc),
                    p=list(p)) %>%
      dplyr::mutate(catch=furrr::future_pmap(
        list(link_id=link_id,db_loc=db_loc,p=p),
        .options = furrr::furrr_options(globals = F),
        carrier::crate(
          function(link_id,db_loc,p) {
            options(dplyr.summarise.inform = FALSE)
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)

            out<-dplyr::tbl(con,"us_flowpaths") %>%
              dplyr::filter(pour_point_id %in% local(link_id)) %>%
              dplyr::rename(link_id=origin_link_id) %>%
              dplyr::left_join(dplyr::tbl(con,"Subbasins_poly") %>%
                                 dplyr::select(link_id,geom),
                               by="link_id")  %>%
              dplyr::show_query() %>%
              utils::capture.output() %>%
              utils::tail(-1) %>%
              paste(collapse = "")

            DBI::dbDisconnect(con)

            p()
            #browser()

            s1<-suppressWarnings(sf::read_sf(db_loc,
                                         query=out))

            s1 %>%
              dplyr::select(-link_id,link_id=pour_point_id) %>%
              sf::st_buffer(units::as_units(0.01,sf::st_crs(s1)$units)) %>%
              sf::st_union() %>%
              #sfheaders::sf_remove_holes() %>%
              nngeo::st_remove_holes() %>%
              sf::st_cast("POLYGON") %>%
              .[1]

          }
        )
      ))

  })

  out<-out %>%
    dplyr::select(tidyselect::any_of("link_id"),tidyselect::any_of("catch")) %>%
    tidyr::unnest(catch) %>%
    dplyr::rename(geom=catch)

  sf::write_sf(
    out,
    db_loc,
    layer="Catchment_poly",
    delete_dsn =F,
    delete_layer=F,
    append = T
  )

  if (nrow(existing_IDs)>0) {
    out<-out %>%
      dplyr::bind_rows(
        sf::read_sf(
          db_loc,
          query=paste0("SELECT `link_id`, `geom` FROM `Catchment_poly` WHERE (`link_id` IN (",paste0("'",existing_IDs$link_id,collapse = "',"),"'))")
        )
      )
  }

  options(parallelly.maxWorkers.localhost=max_cores_opt)

  return(sf::st_as_sf(out))

}
