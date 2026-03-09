
#' Generate and attribute stream line and points along the lines
#'
#' @param input resulting object from `generate_subbasins()`
#' @param extra_attr character. One or more of c("link_slope", "cont_slope", "USChnLn_To", "Elevation", "StOrd_Hack", "StOrd_Str", "StOrd_Hort", "StOrd_Shr"). Optional attributes to add to stream vector outputs.
#' @param points character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/points.shp"), or any GIS data object that will be converted to spatial points. Points representing sampling locations.
#' @param site_id_col character. Variable name in `points` that corresponds to unique site identifiers. This column will be included in all vector geospatial analysis products. Note, if multiple points have the same `site_id_col`, their centroid will be used and returned; if multiple points overlap after snapping, only the first is used.
#' @param snap_distance integer. Maximum distance which points will be snapped to stream lines in map units
#' @param break_on_noSnap logical. Should the function stop if any points are not snapped to a stream segment (i.e., are beyon snap_distance)
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param compress logical. Should output rasters be compressed, slower but more space efficient.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#'
#' @export

attrib_streamline<-function(
    input,
    extra_attr=c(
      "link_slope",
      "cont_slope",
      "USChnLn_To",
      "Elevation",
      "StOrd_Hack",
      "StOrd_Str",
      "StOrd_Hort",
      "StOrd_Shr"
    ),
    points=NULL,
    site_id_col=NULL,
    snap_distance=10L,
    break_on_noSnap=T,
    return_products=F,
    temp_dir=NULL,
    compress=F,
    verbose=F
) {
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

  if (!is.null(points) & is.null(site_id_col)) stop("`site_id_col` can not be NULL if `points` is provided")

  extra_attr<-match.arg(extra_attr,several.ok = T)

  if (!is.null(points) & !is.numeric(snap_distance)) stop("'snap_distance' must be an integer value")
  if (!is.null(points) & !is.logical(break_on_noSnap)) stop("'break_on_noSnap' must be logical")

  if (!is.null(points)) {
    if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
    if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")
    if (site_id_col=="link_id") stop("'site_id_col' cannot be 'link_id'")
  } else {
    site_id_col<-"link_id"
  }

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.logical(compress)) stop("'compress' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),
                        verbose=verbose>2,
                        wd=temp_dir,
                        compress_rasters =compress)

  terra::terraOptions(verbose = verbose>3,
                      tempdir = temp_dir
  )

  gdal_arg<-NULL
  if (compress){
    gdal_arg<-"COMPRESS=NONE"
  }

  options(scipen = 999)
  options(dplyr.summarise.inform = FALSE)

  ds_fp<-db_loc<-zip_loc<-input$outfile
  #fl<-utils::unzip(list=T,zip_loc)
  #db_loc<-file.path(gsub(basename(zip_loc),"",zip_loc),gsub(".zip",".db",basename(zip_loc)))
  # db_loc<-zip_loc
  # ds_fp<-db_loc
  # if (file.exists(ds_fp)) {
  #   stop(paste0("sqlite database: ",ds_fp," Already Exists, please delete the file before proceeding."))
  #   #warning(paste0("sqlite database: ",ds_fp," Already Exists, it was deleted and replaced."))
  #   #file.remove(ds_fp)
  # }

  # utils::unzip(zip_loc,
  #       c("dem_d8.tif","dem_streams_d8.tif","dem_final.tif"),
  #       exdir=temp_dir,
  #       overwrite=T,
  #       junkpaths=T)

  for (i in c("dem_d8","dem_streams_d8","dem_final")){
    terra::writeRaster(
      terra::rast(zip_loc,lyrs=i),
      file.path(temp_dir,paste0(i,".tif")),overwrite=T,gdal=gdal_arg
    )
  }

  dem_final<-terra::rast(file.path(temp_dir,"dem_final.tif"))
  names(dem_final)<-"Elevation"
  terra::writeRaster(dem_final,file.path(temp_dir,"Elevation.tif"),overwrite=T,gdal=gdal_arg)

  #subb<-st_as_sf(vect(file.path("/vsizip",zip_loc,"Subbasins_poly.shp")))

  # Attribute Stream network ------------------------------------------------

  # # Essential attributes
  if (verbose) message("Calculating stream link attributes")

  whitebox::wbt_stream_link_identifier(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="link_id.tif"
  )

  whitebox::wbt_stream_link_class(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="link_class.tif"
  )

  whitebox::wbt_stream_link_length(
    d8_pntr= "dem_d8.tif",
    linkid="link_id.tif",
    output="link_lngth.tif"
  )

  whitebox::wbt_tributary_identifier(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="trib_id.tif"
  )

  whitebox::wbt_farthest_channel_head(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="USChnLn_Fr.tif"
  )

  # Extra Attributes --------------------------------------------------------
  if ("link_slope" %in% extra_attr)
    whitebox::wbt_stream_link_slope( # slope gradient in degrees
      d8_pntr= "dem_d8.tif",
      linkid="link_id.tif",
      dem= "Elevation.tif",
      output="link_slope.tif"
    )
  # # convert above to % change, %slope= tan(Angle in degrees*pi/180)*100
  if ("cont_slope" %in% extra_attr)
    whitebox::wbt_stream_slope_continuous(
      d8_pntr= "dem_d8.tif",
      streams= "dem_streams_d8.tif",
      dem= "Elevation.tif",
      output="cont_slope.tif"
    )
  # # convert above to % change, %slope= tan(Angle in degrees*pi/180)*100

  if ("StOrd_Hack" %in% extra_attr)
    whitebox::wbt_hack_stream_order(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Hack.tif"
    )
  if ("StOrd_Str" %in% extra_attr)
    whitebox::wbt_strahler_stream_order(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Str.tif"
    )

  if ("StOrd_Hort" %in% extra_attr)
    whitebox::wbt_horton_stream_order(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Hort.tif"
    )

  if ("StOrd_Shr" %in% extra_attr)
    whitebox::wbt_shreve_stream_magnitude(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Shr.tif"
    )

  if ("USChnLn_To" %in% extra_attr)
    whitebox::wbt_length_of_upstream_channels(
      d8_pntr= "dem_d8.tif",
      streams= "dem_streams_d8.tif",
      output="USChnLn_To.tif"
    )



  ## Generate stream lines from raster ---------------
  whitebox::wbt_raster_streams_to_vector(
    streams = "link_id.tif",
    d8_pntr = "dem_d8.tif",
    output = "strm_link_id.shp"
  )
  strm<-sf::read_sf(file.path(temp_dir, "strm_link_id.shp")) %>%
    dplyr::select(STRM_VAL)
  sf::st_crs(strm)<-terra::crs(dem_final)
  colnames(strm)[1]<-"link_id"

  #saveRDS(strm,file.path(temp_dir, "strm_link_id.rds"))

  # Generate point attributions ---------------------------------------------
  if (verbose) message("Extracting stream link attributes")

  id1<-terra::rast(file.path(temp_dir,"link_id.tif"))
  id2<-terra::as.points(id1)
  id21<-sf::st_as_sf(id2)

  attr_main<-c(
    file.path(temp_dir,"link_class.tif"),
    file.path(temp_dir,"link_lngth.tif"),
    file.path(temp_dir,"trib_id.tif"),
    file.path(temp_dir,"USChnLn_Fr.tif")
  )

  extra_attr<-file.path(temp_dir,paste0(extra_attr,".tif"))
  extra_attr<-extra_attr[file.exists(extra_attr)]

  extra_attr<-c(attr_main,extra_attr)
  names(extra_attr)<-gsub("\\.tif","",basename(extra_attr))

  attr<-lapply(extra_attr,function(x) terra::rast(x))

  id3<-suppressWarnings(suppressMessages(terra::extract(do.call(c,stats::setNames(attr,NULL)),id2)))

  id4<-id3 %>%
    dplyr::mutate(link_type=dplyr::case_when(
      link_class==1~"Exterior Link",
      link_class==2~"Interior Link",
      link_class==3~"Source Node (head water)",
      link_class==4~"Link Node",
      link_class==5~"Sink Node",
    )) %>%
    tibble::tibble() %>%
    dplyr::select(ID,link_class,link_type,tidyselect::everything())

  names(id4)<-abbreviate(names(id4),10)

  final_points<-dplyr::bind_cols(id21,id4)


  # Instert Points if present -----------------------------------------------

  if (!is.null(points)){
    #browser()
    points<-hydroweight::process_input(points,
                                       align_to  = terra::vect(utils::head(final_points[,1])),
                                       input_name="points")

    sf::write_sf(sf::st_as_sf(points) %>%
                   dplyr::select(tidyselect::any_of(site_id_col),tidyselect::everything()),
                 file.path(temp_dir,"original_points.shp"))

    points<-sf::st_as_sf(points) %>%
      dplyr::group_by(!!rlang::sym(site_id_col)) %>%
      dplyr::summarize(dplyr::across(tidyselect::contains("geom"),sf::st_union)) %>%
      sf::st_centroid() %>%
      dplyr::ungroup()

    if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")

    message("Snapping Points")
    snapped_points<-points %>%
      sf::st_join(final_points %>%
                    dplyr::filter(!is.na(link_type)) %>%
                    dplyr::filter(!link_type %in% c("Link Node","Sink Node","Source Node (head water)")) %>%
                    dplyr::group_by(link_id) %>%
                    dplyr::filter(USChnLn_Fr!=max(USChnLn_Fr)) %>%
                    dplyr::select(ID,link_id),
                  join=nngeo::st_nn,
                  k=1,
                  maxdist = snap_distance,
                  progress =T
      ) %>%
      tibble::as_tibble() %>%
      dplyr::select(-geometry) %>%
      dplyr::left_join(final_points %>%
                         #filter(!link_type %in% c("Link Node","Sink Node","Source Node (head water)")) %>%
                         dplyr::select(ID,link_id),
                       by = c("ID", "link_id")) %>%
      sf::st_as_sf()


    if (break_on_noSnap){
      if (any(is.na(snapped_points$link_id))) stop(paste0("The following points could not be snapped: ", paste0(snapped_points[[site_id_col]][is.na(snapped_points$link_id)],collapse = ", ") ))
    }

    if (any(is.na(snapped_points$link_id))) warning(paste0("The following points could not be snapped and were not included: ", paste0(snapped_points[[site_id_col]][is.na(snapped_points$link_id)],collapse = ", ") ))

    snapped_points<-snapped_points %>%
      dplyr::filter(!is.na(link_id)) %>%
      dplyr::select(tidyselect::any_of(site_id_col),tidyselect::everything()) %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(dplyr::across(tidyselect::everything(),~utils::head(.,1))) %>%
      dplyr::ungroup()

    new_final_points<-final_points
    new_final_points$link_class[new_final_points$ID %in% snapped_points$ID]<-6
    new_final_points$link_type[new_final_points$ID %in% snapped_points$ID]<-"Sample Point"
    new_final_points<-new_final_points %>%
      dplyr::left_join(
        snapped_points %>% tibble::as_tibble() %>% dplyr::select(ID,tidyselect::any_of(site_id_col)),
        by="ID"
      ) %>%
      dplyr::filter(!is.na(!!rlang::sym(site_id_col))) %>%
      dplyr::group_by(link_id) %>%
      dplyr::arrange(link_id,dplyr::desc(USChnLn_Fr)) %>%
      dplyr::mutate(link_id_new=dplyr::row_number()) %>%
      dplyr::mutate(link_id_new=formatC(link_id_new,width=nchar(max(link_id_new)),format="d",flag=0)) %>%
      dplyr::mutate(link_id_new=as.numeric(paste0(paste0(link_id),".",link_id_new)))

    final_points<-final_points %>%
      dplyr::filter(!ID %in% new_final_points$ID) %>%
      dplyr::bind_rows(new_final_points) %>%
      dplyr::group_by(link_id) %>%
      dplyr::arrange(link_id,dplyr::desc(USChnLn_Fr)) %>%
      dplyr::mutate(link_id_new=dplyr::case_when(
        dplyr::row_number()==1 & is.na(link_id_new) ~ link_id,
        T ~ link_id_new
      )) %>%
      tidyr::fill(link_id_new,.direction = "down") %>%
      dplyr::select(-link_id) %>%
      dplyr::rename(link_id=link_id_new) %>%
      dplyr::select(link_id,tidyselect::everything()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(ID)

    snapped_points<-final_points %>%
      dplyr::select(link_id,tidyselect::any_of(site_id_col)) %>%
      dplyr::filter(!is.na(link_id)) %>%
      dplyr::filter(!dplyr::if_any(tidyselect::any_of(site_id_col),is.na))

    sf::write_sf(snapped_points,file.path(temp_dir,"snapped_points.shp"))


    # Make new stream raster --------------------------------------------------

    strm<-final_points %>%
      dplyr::select(link_id) %>%
      dplyr::arrange(link_id) %>%
      terra::vect() %>%
      terra::rasterize(y=dem_final,field="link_id") %>%
      terra::writeRaster(file.path(temp_dir,"new_stream_layer.tif"),overwrite=T,gdal="COMPRESS=NONE")

    whitebox::wbt_raster_streams_to_vector(
      streams = file.path(temp_dir,"new_stream_layer.tif"),
      d8_pntr= file.path(temp_dir,"dem_d8.tif"),
      output = file.path(temp_dir,"new_stream_layer.shp")
    )

    # for some reason wbt_raster_streams_to_vector() rounds numbers weirdly
    un_ID<-unique(final_points$link_id)
    strm<-sf::read_sf(file.path(temp_dir,"new_stream_layer.shp")) %>%
      dplyr::select(STRM_VAL) %>%
      dplyr::rename(link_id=STRM_VAL) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(link_id=un_ID[which.min(abs(link_id-un_ID))]) %>%
      dplyr::ungroup()
    sf::st_crs(strm)<-terra::crs(dem_final)
  }

  sf::write_sf(strm,file.path(temp_dir,"stream_lines.shp"))

  # Add columns: for us and ds link_id and trib_id ---------------------------------

  st_r<-terra::rast(file.path(temp_dir, "dem_streams_d8.tif"))
  d8_pntr<-terra::rast(file.path(temp_dir, "dem_d8.tif"))

  # Upstream ----------------------------------------------------------
  if (verbose) message("Identifying Upstream Links")
  #browser()

  nodes<-final_points %>%
    dplyr::group_by(link_id) %>%
    dplyr::filter(USChnLn_Fr==min(USChnLn_Fr)) %>%
    dplyr::ungroup() %>%
    # filter(link_class %in% c(4,5)) %>%
    dplyr::select(ID) %>%
    dplyr::arrange(ID) %>%
    terra::vect()

  cv <- terra::cells(d8_pntr,nodes)
  all_cell <- terra::cells(d8_pntr,terra::vect(final_points))
  all_cell <- tibble::as_tibble(dplyr::mutate(final_points,cell=all_cell[,2]))

  cv1 <- tibble::tibble(target=cv[,2]) %>%
    dplyr::bind_cols(terra::adjacent(d8_pntr, cells=.$target, directions="queen") %>%
                       data.frame() %>%
                       stats::setNames(c("TL","TM","TR","ML","MR","BL","BM","BR"))) %>%
    tidyr::gather("direction","cell_num",-target) %>%
    dplyr::arrange(target) %>%
    dplyr::filter(!is.na(cell_num),!is.nan(cell_num)) %>%
    dplyr::mutate(on_stream=data.frame(terra::extract(st_r,terra::xyFromCell(st_r,.$cell_num)))$dem_streams_d8) %>%
    dplyr::filter(on_stream==1) %>%
    dplyr::left_join(
      all_cell %>% dplyr::select(link_id,cell),
      by=c("cell_num"="cell")
    )%>%
    dplyr::left_join(
      all_cell %>% dplyr::select(trib_id,cell),
      by=c("cell_num"="cell")
    )%>%
    dplyr::left_join(
      all_cell %>% dplyr::select(target_link_id=link_id,cell),
      by=c("target"="cell")
    )%>%
    dplyr::left_join(
      all_cell %>% dplyr::select(target_trib_id=trib_id,cell),
      by=c("target"="cell")
    )

  cv2<-cv1  %>%
    dplyr::mutate(flow_dir=terra::extract(d8_pntr,terra::xyFromCell(d8_pntr,.$cell_num))$dem_d8) %>%
    dplyr::mutate(flow_in=dplyr::case_when(
      direction=="TL" & flow_dir == 4 ~ T,
      direction=="TM" & flow_dir == 8 ~ T,
      direction=="TR" & flow_dir == 16 ~ T,
      direction=="ML" & flow_dir == 2 ~ T,
      direction=="MR" & flow_dir == 32 ~ T,
      direction=="BL" & flow_dir == 1 ~ T,
      direction=="BM" & flow_dir == 128 ~ T,
      direction=="BR" & flow_dir == 64 ~ T,
      T ~ F
    )) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(n_inflows=sum(flow_in)) %>%
    dplyr::ungroup()

  final_us<-cv2 %>%
    dplyr::filter(flow_in) %>%
    dplyr::select(n_inflows,link_id,trib_id,target_link_id,target_trib_id) %>%
    dplyr::group_by(target_link_id,target_trib_id) %>%
    dplyr::mutate(nm=dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::rename(ustrib_id = trib_id,
                  uslink_id=link_id) %>%
    tidyr::pivot_wider(id_cols = c(target_link_id,target_trib_id),
                       names_from = c(nm),
                       names_sep = "",
                       values_from = c(ustrib_id,uslink_id)) %>%
    dplyr::rename(link_id =target_link_id,
                  trib_id=target_trib_id)

  names(final_us)<-abbreviate(names(final_us),10)


  # Downstream -----------------------------------------------------------
  if (verbose) message("Identifying Downstream Links")

  nodes<-final_points %>%
    dplyr::group_by(link_id) %>%
    dplyr::filter(USChnLn_Fr==max(USChnLn_Fr)) %>%
    dplyr::ungroup() %>%
    dplyr::select(ID) %>%
    dplyr::arrange(ID) %>%
    terra::vect()

  cv <- terra::cells(d8_pntr,nodes)

  cv1 <- tibble::tibble(target=cv[,2]) %>%
    dplyr::bind_cols(terra::adjacent(d8_pntr, cells=.$target, directions="queen")%>%
                       data.frame() %>%
                       stats::setNames(c("TL","TM","TR","ML","MR","BL","BM","BR"))) %>%
    tidyr::gather("direction","cell_num",-target) %>%
    dplyr::arrange(target) %>%
    dplyr::mutate(on_stream=terra::extract(st_r,terra::xyFromCell(st_r,.$cell_num))$dem_streams_d8) %>%
    dplyr::filter(on_stream==1) %>%
    dplyr::left_join(
      all_cell %>% dplyr::select(link_id,cell),
      by=c("cell_num"="cell")
    )%>%
    dplyr::left_join(
      all_cell %>% dplyr::select(trib_id,cell),
      by=c("cell_num"="cell")
    )%>%
    dplyr::left_join(
      all_cell %>% dplyr::select(target_link_id=link_id,cell),
      by=c("target"="cell")
    )%>%
    dplyr::left_join(
      all_cell %>% dplyr::select(target_trib_id=trib_id,cell),
      by=c("target"="cell")
    )

  cv2<-cv1  %>%
    dplyr::mutate(flow_dir=terra::extract(d8_pntr,terra::xyFromCell(d8_pntr,.$target))$dem_d8) %>%
    dplyr::mutate(flow_in=dplyr::case_when(
      direction=="TL" & on_stream & flow_dir == 64 ~ T,
      direction=="TM" & on_stream & flow_dir == 128 ~ T,
      direction=="TR" & on_stream & flow_dir == 1 ~ T,
      direction=="ML" & on_stream & flow_dir == 32 ~ T,
      direction=="MR" & on_stream & flow_dir == 2 ~ T,
      direction=="BL" & on_stream & flow_dir == 16 ~ T,
      direction=="BM" & on_stream & flow_dir == 8 ~ T,
      direction=="BR" & on_stream & flow_dir == 4 ~ T,
      T ~ F
    )) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(n_inflows=sum(flow_in)) %>%
    dplyr::ungroup()

  final_ds<-cv2 %>%
    dplyr::filter(flow_in) %>%
    dplyr::select(n_inflows,link_id,trib_id,target_link_id,target_trib_id) %>%
    dplyr::group_by(target_link_id,target_trib_id) %>%
    dplyr::mutate(nm=dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::rename(dstrib_id = trib_id,
                  dslink_id=link_id) %>%
    tidyr::pivot_wider(id_cols = c(target_link_id,target_trib_id),
                       names_from = c(nm),
                       names_sep = "",
                       values_from = c(dstrib_id,dslink_id)) %>%
    dplyr::rename(link_id =target_link_id,
                  trib_id=target_trib_id)

  names(final_ds)<-abbreviate(names(final_ds),10)

  # Putting it all together -------------------------------------------------
  links<-final_points %>%
    dplyr::group_by(link_id) %>%
    dplyr::filter(USChnLn_Fr==max(USChnLn_Fr)) %>%
    dplyr::ungroup()

  final_links<-links %>%
    dplyr::full_join(final_us, by = c("link_id", "trib_id")) %>%
    dplyr::full_join(final_ds, by = c("link_id", "trib_id")) #

  # check both ID and link_id are unique
  check_link_id<-any(duplicated(final_links$link_id))
  check_id<-any(duplicated(final_links$ID))
  if (check_link_id | check_id) warning("Some link_id's and/or ID's are duplicated, this may indicate an issue with the upstream/downstream IDs")


  #browser()
  sf::write_sf(final_links %>% dplyr::select(link_id),file.path(temp_dir,"stream_links.shp"))
  sf::write_sf(final_points %>% dplyr::select(ID),file.path(temp_dir,"stream_points.shp"))

  con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)

  # Only keep data with an associated stream line
  # SOmetimes you can get a sink at the edge of a DEM with no assicaited stream

  #lns<-sf::read_sf(file.path(temp_dir,paste0("stream_lines",".shp")))

  ot<-final_links %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry) %>%
    #dplyr::filter(link_id %in% lns$link_id) %>%
    dplyr::copy_to(df=.,
                   con,
                   "stream_links_attr",
                   overwrite =T,
                   temporary =F,
                   indexes=c("link_id","trib_id"),
                   analyze=T,
                   in_transaction=T)


  ot<-final_points %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry) %>%
    #dplyr::filter(link_id %in% lns$link_id) %>%
    dplyr::copy_to(df=.,
                   con,
                   "stream_points_attr",
                   overwrite =T,
                   temporary =F,
                   unique_indexes=c("ID"),
                   indexes=c("link_id","trib_id"),
                   analyze=T,
                   in_transaction=T)




  # data.table::fwrite(tibble::tibble(site_id_col=site_id_col),file.path(temp_dir,"site_id_col.csv"))

  # Generate Output ---------------------------------------------------------
  if (verbose) message("Generating Output")

  # dist_list_out<-c(
  #   list.files(temp_dir,"site_id_col"),
  #   list.files(temp_dir,"snapped_points"),
  #   list.files(temp_dir,"original_points"),
  #   list.files(temp_dir,"stream_links"),
  #   list.files(temp_dir,"stream_lines"),
  #   list.files(temp_dir,"stream_points")
  # )

  # dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))
  #
  # out_file<-zip_loc

  t1<-dplyr::copy_to(df=tibble::tibble(site_id_col=site_id_col),
                     con,
                     "site_id_col",
                     overwrite =T,
                     temporary =F,
                     analyze=T,
                     in_transaction=T)

  for (i in c("snapped_points",
              "original_points",
              "stream_links",
              "stream_lines",
              "stream_points")){
    if (file.exists(file.path(temp_dir,paste0(i,".shp")))){
      t1<-sf::write_sf(
        sf::read_sf(file.path(temp_dir,paste0(i,".shp"))), #%>%
          #dplyr::filter(dplyr::if_any(.cols=tidyselect::any_of("link_id"),~.x %in% lns$link_id)),
        db_loc,
        layer=i,
        delete_layer=T,
        #layer_options = "OVERWRITE=true",
        append = T
      )
    }
  }

  DBI::dbDisconnect(con)

  # utils::zip(out_file,
  #     unlist(dist_list_out),
  #     flags = '-r9Xjq'
  # )

  output<-input[!names(input) %in% c("stream_lines",
                                     "links",
                                     "points",
                                     "snapped_points",
                                     "original_points"
  )]

  if (return_products){
    output<-c(
      list(
        stream_lines=strm,
        links=final_links,
        points=final_points,
        snapped_points=snapped_points,
        original_points=points
      ),
      output
    )
  }

  #output$db_loc<-db_loc

  suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive = T)))

  class(output)<-"ihydro"
  return(output)
}
