# necessary temporarily: https://github.com/tidyverse/dtplyr/issues/398
.datatable.aware <- TRUE

# ihydro_layers<-function(x) UseMethod("ihydro_layers")
# target_id_fun<-function(x) UseMethod("target_id_fun")
# target_o_fun<-function(x) UseMethod("target_o_fun")

#' ihydro Class
#'
#' Creates an object of class "ihydro".
#'
#' @param file filepath to a geopackage file (.gpkg)
#'
#' @return an object of class ihydro
#' @export
#'

as.ihydro<-function(file){
  file<-normalizePath(file,mustWork =F)
  if (!file.exists(file)) stop("'file' does not exist")
  if (!grepl("\\.gpkg$",file)) stop("'file' must be a file path ending in '.gpkg'")
  out<-list(outfile=file)
  class(out)<-"ihydro"
  return(out)
}

#' @export
as_ihydro<-function(...) as.ihydro(...)

# print.ihydro<-function(x){
#   ihydro_layers(x)
# }

#' Describe layers objects of class ihydro
#'
#' @param input an object of class ihydro
#'
#' @return table of objects in the ihydro geopackage
#' @export
#'

ihydro_layers<-function(input) {
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")

  con<-DBI::dbConnect(RSQLite::SQLite(),input$outfile)

  rast_lyrs<-dplyr::tbl(con,"gpkg_contents") %>%
    dplyr::collect() %>%
    dplyr::select(layer_name=table_name,
                  data_type) %>%
    dplyr::filter(grepl("gridded",data_type)) %>%
    dplyr::mutate(data_type="Raster")

  lyr<-DBI::dbListTables(con)

  DBI::dbDisconnect(con)

  vect_lyrs<-tibble::tibble(
    layer_name=lyr[grepl("^rtree_",lyr) & grepl("_geom$",lyr)],
    data_type="Vector"
  ) %>%
    mutate(layer_name=gsub("^rtree_|_geom$","",layer_name))

  lyr_data<-lyr[!grepl("^rtree_|^gpkg_|^sqlite_|_geom$",lyr) & !lyr %in% rast_lyrs$layer_name]

  lyr_data<-tibble::tibble(
    layer_name=lyr_data,
    data_type="Table"
  )

  out<-dplyr::bind_rows(rast_lyrs,vect_lyrs) %>%
    dplyr::bind_rows(lyr_data)  %>%
    dplyr::mutate(data_group=dplyr::case_when(
      grepl("^dem_",layer_name) ~ "hydro",
      grepl("iFLS|HAiFLS|iFLO|HAiFLO",layer_name) ~ "iDW",
      layer_name %in% c("site_id_col","DEM_Extent","loi_meta","target_o_meta") ~ "meta",
      layer_name %in% c("snapped_points","original_points") ~ "sample_points",
      layer_name %in% c("stream_links","stream_lines","stream_points",
                        "Subbasins_poly","Catchment_poly",
                        "stream_points_attr","stream_links_attr") ~ "hydro",
      layer_name %in% c("ds_flowpaths","us_flowpaths") ~ "flow_path",
      layer_name %in% c("funcon_pwise_dist","fcon_pwise_dist") ~ "pwise_dist",
      T ~ "loi"
    )) %>%
    dplyr::arrange(data_type,data_group)

  # lyr<-suppressWarnings(sf::st_layers(input$outfile))
  #
  # vect_lyrs<-tibble::tibble(
  #   layer_name=unlist(lyr[[1]]),
  #   data_type=unlist(lyr[[2]])
  # ) %>%
  #   dplyr::mutate(data_type=dplyr::case_when(
  #     grepl("Polygon",data_type) ~ "Polygon",
  #     grepl("Line",data_type) ~ "Line",
  #     is.na(data_type) ~ "Table",
  #     T ~ data_type
  #   )) %>%
  #   dplyr::mutate(data_group=dplyr::case_when(
  #     layer_name %in% c("site_id_col","DEM_Extent","loi_meta","target_o_meta") ~ "meta",
  #     layer_name %in% c("snapped_points","original_points") ~ "sample_points",
  #     layer_name %in% c("stream_links","stream_lines","stream_points",
  #                       "Subbasins_poly","Catchment_poly",
  #                       "stream_points_attr","stream_links_attr") ~ "hydro",
  #     layer_name %in% c("ds_flowpaths","us_flowpaths") ~ "flow_path",
  #     layer_name %in% c("funcon_pwise_dist","fcon_pwise_dist") ~ "pwise_dist"
  #   ))



  return(out)
}


target_id_fun<-function(
    db_fp,
    sample_points=NULL,
    link_id=NULL,
    segment_whole=F,
    target_o_type=c("point","segment_point","segment_whole")
) {

  target_o_type<-match.arg(target_o_type,several.ok = F)

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)

  site_id_col<-dplyr::collect(dplyr::tbl(con,"site_id_col"))$site_id_col

  all_points<-dplyr::collect(dplyr::tbl(con,"stream_links_attr")) %>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),~dplyr::na_if(.,"")))

  DBI::dbDisconnect(con)

  # Get target link_id ------------------------------------------------------
  sample_points<-as.character(sample_points)
  link_id<-as.character(link_id)
  if (length(sample_points)==0 & length(link_id)==0) {
    #message("`sample_points` and `link_id` are NULL, all `link_id`s will evaluated")
    target_IDs<-all_points %>%
      tibble::as_tibble() %>%
      dplyr::select(link_id,tidyselect::any_of(site_id_col))
  } else {
    if (site_id_col!="link_id" & length(sample_points)>0){
      target_IDs<-all_points %>%
        tibble::as_tibble() %>%
        dplyr::select(link_id,tidyselect::any_of(site_id_col)) %>%
        dplyr::filter(!!rlang::sym(site_id_col) %in% sample_points)
    } else {
      target_IDs<-NULL
    }

    if (length(link_id)>0){
      link_id_sel<-link_id
      target_IDs<-dplyr::bind_rows(
        target_IDs,
        all_points %>%
          tibble::as_tibble() %>%
          dplyr::select(link_id,tidyselect::any_of(site_id_col)) %>%
          dplyr::filter(link_id %in% link_id_sel)
      )
    }
  }

  target_IDs<-target_IDs %>%
    dplyr::distinct()

  if (segment_whole) {
    target_IDs<-target_IDs %>%
      dplyr::select(link_id) %>%
      dplyr::mutate(link_id=as.character(floor(as.numeric(link_id))))
  }

  if (target_o_type!="point") {
    target_O<-sf::read_sf(db_fp,
                          query=paste0("SELECT `link_id`, `geom` FROM `stream_lines` WHERE (`link_id` IN (",paste0("'",target_IDs$link_id,collapse = "',"),"'))")
    )

    target_IDs<-target_IDs %>%
      dplyr::filter(link_id %in% target_O$link_id)

  }

  target_IDs<-target_IDs %>%
    dplyr::mutate(
      link_id=as.character(link_id)
    )

  return(target_IDs)
}


target_o_fun<-function(
    db_fp,
    target_IDs,
    target_o_type=c("point","segment_point","segment_whole")
) {

  # Select correct target for O -------------------------------------
  if (target_o_type=="point"){
    target_O<-sf::read_sf(db_fp,
                          query=paste0("SELECT `link_id`, `geom` FROM `stream_links` WHERE (`link_id` IN (",paste0("'",target_IDs$link_id,collapse = "',"),"'))")
    )
  } else {
    if (target_o_type=="segment_point"){
      target_O<-sf::read_sf(db_fp,
                            query=paste0("SELECT `link_id`, `geom` FROM `stream_lines` WHERE (`link_id` IN (",paste0("'",target_IDs$link_id,collapse = "',"),"'))")
      )
    } else {
      #if (verbose) message("Merging stream segments")

      target_O<-sf::read_sf(db_fp,
                            query=paste0("SELECT `link_id`, `geom` FROM `stream_lines` WHERE (`link_id` IN (",paste0("'",target_IDs$link_id,collapse = "',"),"'))")
      ) %>%
        dplyr::select(link_id) %>%
        dplyr::mutate(link_id=as.character(floor(as.numeric(link_id)))) %>%
        dplyr::filter(link_id %in% target_IDs$link_id) %>%
        dplyr::group_by(link_id) %>%
        dplyr::summarize(geom=sf::st_union(geom)) %>%
        dplyr::ungroup()

    }
  }

  return(target_O)

}


if (F) {
  lines_out<-hydro_out$stream_lines %>%
    dplyr::filter(link_id %in% us_647$origin_link_id #|
                  #link_id %in% ds_1071.1$destination_link_id
    )
  sub_out<-hydro_out$subbasins %>%
    dplyr::filter(link_id %in% us_647$origin_link_id #|
                  #link_id %in% ds_1071.1$destination_link_id
    )

  t1<-tm_shape(sub_out) + tm_polygons(col="white",alpha =0.8,legend.show=F,lwd =2) +
    tm_shape(lines_out) + tm_lines(col="blue",alpha =1,legend.show=F,lwd =4) +
    tm_shape(hydro_out$links %>% filter(link_id %in% c("647"))) + #,"1071.1"
    tm_dots(legend.show=F,size=1.5,border.col="black",border.alpha=1,border.lwd=1) +
    tm_layout(frame = FALSE,bg.color="transparent")

  tmap_save (t1, filename = "man/figures/logo_img.png", device=png, bg="transparent")

  hexSticker::sticker("man/figures/logo_img.png",
                      package="ihydro",
                      p_color = "#080808",
                      p_family = "Courgette",
                      p_fontface = "bold",
                      p_x = 1,
                      p_y = 1.6,
                      h_size = 1.7,
                      h_fill = "#51996d",
                      h_color = "#0f1f15",
                      p_size=20,
                      s_x=1,
                      s_y=0.9,
                      s_width=.8,
                      filename="man/figures/logo.png")
}


# GPKG raster writer helper -----------------------------------------------

writeRaster_fun <- function(x, filename,layer_name=NULL,NA_val = -9999.9999){
  if (is.null(layer_name)) layer_name <- names(x)
  x[x==NA_val] <- NA_val + .Machine$double.eps
  x[is.na(x)] <- NA_val
  terra::NAflag(x) <- NA_val

  out <- terra::writeRaster(
    x = x,
    filename = filename,
    NAflag = NA_val,
    filetype = "GPKG",
    datatype = "FLT4S",
    gdal = c("APPEND_SUBDATASET=YES",
             paste0("RASTER_TABLE=",layer_name,""),
             paste0("RASTER_IDENTIFIER=",layer_name,""),
             paste0("RASTER_DESCRIPTION=",layer_name,""),
             paste0("FIELD_NAME=",layer_name,""),
             paste0("UOM=",layer_name,""),
             paste0("QUANTITY_DEFINITION=",layer_name,""),
             "TILE_FORMAT=TIFF",
             "METADATA_TABLES=YES",
             "COMPRESS=LZW",
             "TILED=YES",
             "BAND_COUNT=1"
    )
  )

  con<-DBI::dbConnect(RSQLite::SQLite(),filename)
  DBI::dbExecute(con, paste0("UPDATE gpkg_2d_gridded_coverage_ancillary SET data_null = ",NA_val," WHERE tile_matrix_set_name = '",layer_name,"'"))
  DBI::dbDisconnect(con)

  return(out)
}
