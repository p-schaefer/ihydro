#' Quickly attribute stream segments/sampling points with layers of interest (loi)
#'
#' @param input output from `process_hydrology()` (if `process_loi()` was not run on `process_hydrology()`, `loi_file` must be specified)
#' @param out_filename Output file name, ending in '.csv'.
#' @param loi_file filepath of `process_loi()` output (optional, will search for loi layers in 'input' if not specified).
#' @param loi_cols character or NULL. Names of loi layers to include in summary. If NULL, all layers used.
#' @param iDW_file filepath of `prep_weights()` output (optional, will search for weight layers in 'input' if not specified).
#' @param store_iDW logical, if 'iDW_file' is NULL, weights will be added to 'input'
#' @param sample_points character or NULL. IDs of unique station identifiers provided in 'site_id_col' of `generate_vectors()` to calculate attributes for.
#' @param link_id character or NULL. 'link_id's of reaches to calculate attributes for.
#' @param target_o_type character. One of: c("point","segment_point","segment_whole"). Target for"iFLO", and "HAiFLO" weighting schemes. "Point" represents the sampling point on the stream, "segment_point" represents the upstream segment of the sampling points, and "segment_whole" will target the entire reach, regardless of where sampling occurred.
#' @param weighting_scheme character. One or more weighting schemes: c("lumped","iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param loi_numeric_stats character. One or more of c("mean", "sd", "median", "min", "max", "sum"). Only distance-weighted versions of mean and SD are returned for all weighting schemes except lumped.
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#' @param backend character. One of "data.table" or "SQLite"
#' @param SQLdb_path character. If backend = "SQLite" directory to store temporary SQL tables, or ":memory:" for an in-memory SQL table
#'
#' @return A data.frame of weighted attributes for the requested areas
#' @export
#'

fasttrib_points<-function(
    input,
    out_filename,
    loi_file=NULL,
    loi_cols=NULL,
    iDW_file=NULL,
    store_iDW=F,
    sample_points=NULL,
    link_id=NULL,
    target_o_type=c("point","segment_point","segment_whole"),
    weighting_scheme =  c("lumped",  "iFLS", "HAiFLS","iFLO",  "HAiFLO"),
    loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir=NULL,
    verbose=F,
    backend=c("tibble","SQLite","data.table"),
    SQLdb_path=":memory:"
){
  backend<-match.arg(backend)

  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  if (inherits(loi_file,"ihydro")) loi_file<-loi_file$outfile
  if (inherits(iDW_file,"ihydro")) iDW_file<-iDW_file$outfile

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir,recursive = T)
  temp_dir<-normalizePath(temp_dir)

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),
                        verbose=verbose>2,
                        wd=temp_dir)

  terra::terraOptions(verbose = verbose>3,
                      tempdir = temp_dir
  )

  target_o_type<-match.arg(target_o_type,several.ok = F)
  weighting_scheme<-match.arg(weighting_scheme,several.ok = T)
  loi_numeric_stats<-match.arg(loi_numeric_stats,several.ok = T)

  weighting_scheme_s<-weighting_scheme[grepl("FLS|iEucS",weighting_scheme)]
  weighting_scheme_o<-weighting_scheme[!grepl("lumped|FLS|iEucS",weighting_scheme)]
  lumped_scheme<-"lumped" %in% weighting_scheme
  #if (length(weighting_scheme_o)>0) message("Calculation for iFLO, and HAiFLO are slow")

  loi_numeric_stats<-stats::setNames(loi_numeric_stats,loi_numeric_stats)

  attr_db_loc<-db_loc<-dw_dir<-zip_loc<-input$outfile

  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  if (!is.logical(store_iDW)) stop("'store_iDW' must be logical")

  # Get Target IDs ----------------------------------------------------------
  target_IDs<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=target_o_type=="segment_whole",
    target_o_type=target_o_type
  )

  lyrs<-ihydro_layers(input)

  if (is.null(loi_file)) {
    if (!"loi_meta" %in% lyrs$layer_name) stop("no 'loi' data found in 'input'")
    loi_file<-input$outfile
    loi_file<-as.ihydro(loi_file)
  } else {
    loi_file<-as.ihydro(loi_file)
    loi_lyrs<-ihydro_layers(loi_file)
    if (!"loi_meta" %in% loi_lyrs$layer_name) stop("no 'loi' data found in 'loi_lyrs'")
  }
  # TODO: Add checks to make sure requested loi_cols are in input
  loi_meta<-sf::read_sf(loi_file,"loi_meta")
  if (is.null(loi_cols)) loi_cols<-loi_meta$loi_var_nms
  if (any(!loi_cols %in% loi_meta$loi_var_nms)) stop("some `loi_cols` not present in `loi_file`")
  loi_meta<-loi_meta %>%
    dplyr::filter(loi_var_nms %in% loi_cols)


  # Not 100% sure this is going to work
  if (is.null(iDW_file)) {
    if (store_iDW){
      iDW_file<-input$outfile

    } else {
      iDW_file<-tempfile()
      dir.create(iDW_file,recursive = T)
      iDW_file<-file.path(iDW_file,"temp_iDW.gpkg")
    }
  }

  if (verbose) message("Preparing inverse distance weights")
  target_o_meta<-target_IDs
  target_o_meta$unn_group<-"1"

  idw_lyrs<-NULL
  if (!all(weighting_scheme=="lumped")) {
    iDW_out<-prep_weights(
      input=input,
      output_filename=iDW_file,
      sample_points=sample_points,
      link_id=link_id,
      target_o_type=target_o_type,
      weighting_scheme =   weighting_scheme[!grepl("lumped",weighting_scheme)],
      inv_function = inv_function,
      temp_dir=temp_dir,
      verbose=verbose
    )

    if (!inherits(iDW_out,"ihydro")) stop(iDW_out)

    idw_lyrs<-ihydro_layers(iDW_out)

    # TODO: Add checks to make sure requested weighting_scheme are in input
    if (any(!sapply(weighting_scheme[!grepl("lumped",weighting_scheme)],function(x) any(grepl(x,idw_lyrs$layer_name))))) stop("no 'iDW' data found in 'input'")

    target_o_meta<-sf::read_sf(iDW_file,"target_o_meta")
  }





  #browser()

  # Get subbasin IDs --------------------------------------------------------
  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)

  subb_IDs<-dplyr::tbl(con,"us_flowpaths") %>%
    dplyr::filter(pour_point_id %in% local(target_IDs$link_id)) %>%
    dplyr::select(pour_point_id,origin_link_id) %>%
    dplyr::distinct() %>%
    dplyr::collect() %>%
    dplyr::left_join(
      target_o_meta %>%
        dplyr::mutate(link_id=as.character(link_id)),
      by=c("pour_point_id"="link_id")
    ) %>%
    dplyr::arrange(origin_link_id)

  DBI::dbDisconnect(con)

  # Setup output ------------------------------------------------------------
  temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
  dir.create(temp_dir_sub,recursive = T)

  o1<-extract_raster_attributes(
    input=input,
    iDW_file=iDW_file,
    loi_file=loi_file,
    weighting_scheme=weighting_scheme,
    loi_numeric_stats=loi_numeric_stats,
    loi_cols=loi_cols,
    subb_IDs=subb_IDs,
    temp_dir_sub=temp_dir_sub,
    verbose=verbose,
    backend=backend,
    SQLdb_path=SQLdb_path
  )

  target_IDs_out<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=F,
    target_o_type=target_o_type
  )

  #browser()
  final_out<-target_IDs_out %>%
    dplyr::left_join(
      o1 %>% dplyr::mutate(pour_point_id=as.character(pour_point_id)),
      by=c("link_id"="pour_point_id"),
      multiple = "all"
    )

  data.table::fwrite(x=final_out,
                     file=out_filename,
                     buffMB = 128L,
                     nThread = 1,
                     showProgress = F)

  return(final_out)


}




extract_raster_attributes<-function(
    input,
    iDW_file,
    loi_file,
    weighting_scheme,
    loi_cols,
    subb_IDs,
    loi_numeric_stats,
    temp_dir_sub,
    verbose,
    backend=c("data.table","SQLite"),
    SQLdb_path=":memory:"
){
  if (inherits(loi_file,"ihydro")) loi_file<-loi_file$outfile
  if (inherits(iDW_file,"ihydro")) iDW_file<-iDW_file$outfile

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)
  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  loi_meta<-sf::read_sf(loi_file,"loi_meta")
  loi_meta<-loi_meta %>% dplyr::filter(loi_var_nms %in% loi_cols)

  if (verbose) message("Calculating Catchments")

  input_poly<-ihydro::get_catchment(input,link_id=unique(subb_IDs$pour_point_id))

  if (verbose) message("Calculating Weighted Attributes")


  # Calculate in parallel ---------------------------------------------------

  progressr::with_progress(enable=T,{
    p <- progressr::progressor(steps = (length(unique(subb_IDs$unn_group))))

    ot<-subb_IDs %>%
      dplyr::group_by(unn_group) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(core=rep(1:n_cores,length.out=dplyr::n())) %>%
      dplyr::group_by(core) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        out=furrr::future_pmap(
          #out=purrr::pmap(
          list(x=data,
               input=list(as.ihydro(input$outfile)),
               iDW_file=list(iDW_file),
               loi_file=list(loi_file),
               loi_meta=list(loi_meta),
               weighting_scheme=list(weighting_scheme),
               loi_numeric_stats=list(loi_numeric_stats),
               loi_cols=list(loi_cols),
               p=list(p),
               temp_dir_sub=list(temp_dir_sub),
               n_cores=list(n_cores),
               backend=list(backend),
               SQLdb_path=list(SQLdb_path)
          ),
          .options = furrr::furrr_options(globals = F),
          carrier::crate(
            function(x,
                     input,
                     iDW_file,
                     loi_file,
                     loi_meta,
                     weighting_scheme,
                     loi_numeric_stats,
                     loi_cols,
                     p,
                     temp_dir_sub,
                     n_cores,
                     backend,
                     SQLdb_path) {
              #browser()
              ihydro::.summ_fn(
                x=x,
                input=input,
                iDW_file=iDW_file,
                loi_file=loi_file,
                loi_meta=loi_meta,
                weighting_scheme=weighting_scheme,
                loi_numeric_stats=loi_numeric_stats,
                loi_cols=loi_cols,
                p=p,
                temp_dir_sub=temp_dir_sub,
                n_cores=n_cores,
                backend=backend,
                SQLdb_path=SQLdb_path
              )

            }
          )
        ))

  })


  out<-ot %>%
    dplyr::select(out) %>%
    tidyr::unnest(cols=out)



  # Any that couldn't be calculated in parallel -----------------------------
  incomplete_proc<-out$pour_point_id[out$status=="Incomplete"]

  if (length(incomplete_proc)>0){
    if (n_cores>1) {
      message(
        paste0(
          "The following link_ids could not be processed in parallel due to insufficient memory: ",
          paste(out$pour_point_id[out$status=="Incomplete"],collapse = ","),
          ". Attempting to recalculate sequentially."
        )
      )

      progressr::with_progress(enable=T,{
        p <- progressr::progressor(steps = length(incomplete_proc))

        ot2<-subb_IDs %>%
          dplyr::filter(pour_point_id %in% incomplete_proc) %>%
          dplyr::group_by(unn_group) %>%
          tidyr::nest() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(core=rep(1,length.out=dplyr::n())) %>%
          dplyr::group_by(core) %>%
          tidyr::nest() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            out=purrr::pmap(
              list(x=data,
                   input=list(as.ihydro(input$outfile)),
                   iDW_file=list(iDW_file),
                   loi_file=list(loi_file),
                   loi_meta=list(loi_meta),
                   weighting_scheme=list(weighting_scheme),
                   loi_numeric_stats=list(loi_numeric_stats),
                   loi_cols=list(loi_cols),
                   p=list(p),
                   temp_dir_sub=list(temp_dir_sub),
                   n_cores=list(1),
                   backend=list(backend),
                   SQLdb_path=list(SQLdb_path)
              ),
              carrier::crate(
                function(x,
                         input,
                         iDW_file,
                         loi_file,
                         loi_meta,
                         weighting_scheme,
                         loi_numeric_stats,
                         loi_cols,
                         p,
                         temp_dir_sub,
                         n_cores,
                         backend,
                         SQLdb_path) {

                  ihydro::.summ_fn(
                    x=x,
                    input=input,
                    iDW_file=iDW_file,
                    loi_file=loi_file,
                    loi_meta=loi_meta,
                    weighting_scheme=weighting_scheme,
                    loi_numeric_stats=loi_numeric_stats,
                    loi_cols=loi_cols,
                    p=p,
                    temp_dir_sub=temp_dir_sub,
                    n_cores=n_cores,
                    backend=backend,
                    SQLdb_path=SQLdb_path
                  )
                }
              )
            ))
      })

      out2<-ot2 %>%
        dplyr::select(out) %>%
        tidyr::unnest(cols=out)

      out<-bind_rows(out,out2)

      incomplete_proc<-out$pour_point_id[out$status=="Incomplete"]
    }
  }

  if (length(incomplete_proc)>0) {
    warning(
      paste0(
        "The following link_ids could not be processed due to insufficient memory: ",
        paste(out$pour_point_id[out$status=="Incomplete"],collapse = ",")
      )
    )

    # TODO: If all else fails, extract the data using subcatchments, write to
    #        sqlite database and do calculations there...
  }


  return(out)

}

#' Process gpkg rasters for attribute calculation
#'
#' @keywords internal
#' @return a data.frame of attributes
#' @export
#'
.summ_fn<-function(x,
                   input,
                   iDW_file,
                   loi_file,
                   loi_meta,
                   weighting_scheme,
                   loi_numeric_stats,
                   loi_cols,
                   p,
                   temp_dir_sub,
                   n_cores,
                   backend=c("tibble","data.table","SQLite"),
                   SQLdb_path=":memory:"
){
  backend<-match.arg(backend)

  #sys.mem<-(memuse::Sys.meminfo()$freeram/n_cores)*0.9
  sys.mem<-memuse::Sys.meminfo()$freeram*0.9
  max_pixels_mem<-as.numeric(sys.mem)/12
  max_pixels_mem<-min(max_pixels_mem,3e+06) # this isn't quite working right

  if (F) {
    # The below section is for using exact_extract
    # loi_numeric_stats[loi_numeric_stats=="sd"]<-"stdev"
    # names(loi_numeric_stats)[names(loi_numeric_stats)=="sd"]<-"stdev"
    #
    # if (length(weighting_scheme[weighting_scheme!="lumped"])>0) {
    #   #names(loi_numeric_stats)<-paste0("lumped_",names(loi_numeric_stats))
    #   if ("mean" %in% loi_numeric_stats) loi_numeric_stats<-c(loi_numeric_stats,weighted_mean="weighted_mean")
    #   if ("stdev" %in% loi_numeric_stats) loi_numeric_stats<-c(loi_numeric_stats,weighted_stdev="weighted_stdev")
    # }
    #
    # nm_fn<-function(c_nm,weight_nm,loi_numeric_stats,cat=F,append_cols="link_id"){
    #   loi_numeric_stats_loc<-lapply(loi_numeric_stats,function(x) grepl(paste0("^",x,"\\."),c_nm) )
    #
    #   c_nm_new<-gsub(paste0(paste0(loi_numeric_stats,"."),collapse = "|"),"",c_nm)
    #
    #   w_loc<-loi_numeric_stats_loc[grepl("weighted",names(loi_numeric_stats_loc))]
    #
    #   if (length(w_loc)>0){
    #     w_loc<-purrr::reduce(w_loc,function(x,y) x|y)
    #   } else {
    #     w_loc<-rep(F,length(c_nm_new))
    #   }
    #
    #   c_nm_new[w_loc]<-paste0(c_nm_new[w_loc],paste0("_",weight_nm))
    #   c_nm_new[!w_loc & !c_nm_new %in% append_cols]<-paste0(c_nm_new[!w_loc & !c_nm_new %in% append_cols],paste0("_","lumped"))
    #
    #   names(loi_numeric_stats_loc)<-gsub("stdev","sd",names(loi_numeric_stats_loc))
    #   if (cat) {
    #     names(loi_numeric_stats_loc)<-gsub("mean","prop",names(loi_numeric_stats_loc))
    #   }
    #
    #   for (i in names(loi_numeric_stats_loc)){
    #     c_nm_new[loi_numeric_stats_loc[[i]]]<-paste0(c_nm_new[loi_numeric_stats_loc[[i]]],
    #                                                  "_",
    #                                                  gsub("weighted_","",i))
    #   }
    #
    #   return(c_nm_new)
    #
    # }
  }


  sys.disk<-NULL

  if (F & "SQLite" %in% backend){ #not using this currently
    sys.disk <- system(paste0("df -BM ",temp_dir_sub),intern=T)
    sys.disk <- strsplit(sys.disk[length(sys.disk)], "[ ]+")[[1]]
    sys.disk <- as.numeric(gsub("\\D", "", sys.disk[4]))/1024
  }

  options(dplyr.summarise.inform = FALSE)
  options(scipen = 999)
  `%>%` <- magrittr::`%>%`

  iDWs_rasts<-NULL
  iDWo_rasts<-NULL

  loi_rasts<-terra::rast(loi_file,loi_cols)

  if (length(weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])>0){
    iDWs_rasts<-terra::rast(iDW_file,weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])
  }

  purrr::pmap_dfr(
    list(
      y=x$data,
      yy=x$unn_group,
      loi_rasts=list(loi_rasts),
      iDW_file=list(iDW_file),
      iDWs_rasts=list(iDWs_rasts),
      iDWo_rasts=list(iDWo_rasts),
      loi_cols=list(loi_cols),
      input=list(input),
      weighting_scheme=list(weighting_scheme),
      loi_meta=list(loi_meta),
      loi_numeric_stats=list(loi_numeric_stats),
      backend=list(backend),
      SQLdb_path=list(SQLdb_path),
      sys.mem=list(sys.mem),
      temp_dir_sub=list(temp_dir_sub),
      p=list(p)
    ),
    carrier::crate(
      function(y,
               yy,
               loi_rasts,
               iDW_file,
               iDWs_rasts,
               iDWo_rasts,
               loi_cols,
               input,
               weighting_scheme,
               loi_meta,
               loi_numeric_stats,
               backend,
               SQLdb_path,
               sys.mem,
               temp_dir_sub,
               p){
        #browser()

        input_poly<-ihydro::get_catchment(input,link_id=unique(y$pour_point_id))

        if (length(weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")])>0){
          iDWo_rasts<-terra::rast(iDW_file,
                                  unlist(purrr:::map(unique(yy),
                                                     ~paste0(
                                                       paste0(
                                                         weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")],
                                                         "_unn_group"
                                                       ),
                                                       .
                                                     )))
          )

          names(iDWo_rasts)<-gsub(paste0("_unn_group",yy),"",names(iDWo_rasts))

        }

        if (F){
          # The below section is for using exact_extract
          # iDW_rasts<-c(iDWs_rasts,iDWo_rasts)
          #
          # num_out<-NULL
          # cat_out<-NULL
          #
          # # Numeric Rasters
          # if (length(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"])>0){
          #   num_out<-purrr::map(names(iDW_rasts),function(idw){
          #     loi_numeric_stats_sub<-loi_numeric_stats
          #
          #     if (!("lumped" %in% weighting_scheme & idw == weighting_scheme[!weighting_scheme %in% "lumped"][[1]])) {
          #       loi_numeric_stats_sub<-loi_numeric_stats_sub[grepl("^weighted_",loi_numeric_stats_sub)]
          #     }
          #
          #     ot<-exactextractr::exact_extract(
          #       x=terra::subset(loi_rasts,loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
          #       y=input_poly,
          #       weights=terra::subset(iDW_rasts,idw),
          #       append_cols="link_id",
          #       fun=loi_numeric_stats_sub,
          #       default_value=NA_real_,
          #       default_weight=0,
          #       max_cells_in_memory=max_pixels_mem,
          #       progress=F
          #     )
          #
          #     colnames(ot)<-nm_fn(c_nm=colnames(ot),
          #                         weight_nm=idw,
          #                         loi_numeric_stats=loi_numeric_stats_sub,
          #                         cat=F
          #     )
          #
          #     return(ot)
          #
          #   }) %>%
          #     purrr::reduce(dplyr::left_join,by="link_id")
          # }
          #
          # # Categoric Rasters
          # if (length(loi_meta$loi_var_nms[loi_meta$loi_type=="cat_rast"])>0){
          #   cat_out<-purrr::map(names(iDW_rasts),function(idw){
          #     loi_numeric_stats_sub<-loi_numeric_stats[grepl("mean",loi_numeric_stats)]
          #
          #     if (!("lumped" %in% weighting_scheme & idw == weighting_scheme[!weighting_scheme %in% "lumped"][[1]])) {
          #       loi_numeric_stats_sub<-loi_numeric_stats_sub[grepl("^weighted_",loi_numeric_stats_sub)]
          #     }
          #
          #     ot<-exactextractr::exact_extract(
          #       x=terra::subset(loi_rasts,loi_meta$loi_var_nms[loi_meta$loi_type=="cat_rast"]),
          #       y=input_poly,
          #       weights=terra::subset(iDW_rasts,idw),
          #       append_cols="link_id",
          #       fun=loi_numeric_stats_sub,
          #       default_value=0,
          #       default_weight=0,
          #       max_cells_in_memory=max_pixels_mem,
          #       progress=F
          #     )
          #
          #     colnames(ot)<-nm_fn(c_nm=colnames(ot),
          #                         weight_nm=idw,
          #                         loi_numeric_stats=loi_numeric_stats_sub,
          #                         cat=T
          #     )
          #
          #     return(ot)
          #
          #   }) %>%
          #     purrr::reduce(dplyr::left_join,by="link_id")
          # }
          #
          # out<-list(num_out,cat_out)
          # out<-out[!sapply(out,is.null)]
          # out<-purrr::reduce(out,dplyr::full_join,by="link_id") %>%
          #   dplyr::mutate(status="Complete",
          #                 pour_point_id=link_id) %>%
          #   dplyr::select(
          #     tidyselect::any_of("pour_point_id"),
          #     tidyselect::any_of("status"),
          #     tidyselect::contains(loi_meta$loi_var_nms)
          #   )
          #
          # p()
          #
          # return(out)
        }

        input_rasts<-c(loi_rasts,iDWs_rasts,iDWo_rasts)

        loi_types<-table(loi_meta$loi_type[loi_meta$loi_var_nms %in% loi_cols])
        num_rast_analysis_cols<-loi_types["num_rast"]*length(weighting_scheme[weighting_scheme!="lumped"])*3
        cat_rast_analysis_cols<-loi_types["cat_rast"]*length(weighting_scheme[weighting_scheme!="lumped"])*2
        analysis_cols<-sum(c(num_rast_analysis_cols,cat_rast_analysis_cols),na.rm=T)

        #Estimates the number of rows that could fit into memory without analysis
        max.obj.fulldata<-memuse::howmany(sys.mem,
                                          ncol=length(loi_cols)+length(weighting_scheme[weighting_scheme!="lumped"]))


        # Summary Function --------------------------------------------------------
        ot<-purrr::pmap_dfr(
          list(
            sub_poly=split(input_poly,seq_along(input_poly$link_id)),
            input_rasts=list(input_rasts),
            loi_rasts=list(loi_rasts),
            iDWs_rasts=list(iDWs_rasts),
            iDWo_rasts=list(iDWo_rasts),
            loi_cols=list(loi_cols),
            input=list(input),
            weighting_scheme=list(weighting_scheme),
            loi_meta=list(loi_meta),
            loi_numeric_stats=list(loi_numeric_stats),
            backend=list(backend),
            SQLdb_path=list(SQLdb_path),
            sys.mem=list(sys.mem),
            temp_dir_sub=list(temp_dir_sub),
            loi_types=list(loi_types),
            num_rast_analysis_cols=list(num_rast_analysis_cols),
            cat_rast_analysis_cols=list(cat_rast_analysis_cols),
            analysis_cols=list(analysis_cols),
            max.obj.fulldata=list(max.obj.fulldata),
            p=list(p)
          ),
          carrier::crate(
            function(sub_poly,
                     input_rasts,
                     loi_rasts,
                     iDWs_rasts,
                     iDWo_rasts,
                     loi_cols,
                     input,
                     weighting_scheme,
                     loi_meta,
                     loi_numeric_stats,
                     backend,
                     SQLdb_path,
                     sys.mem,
                     temp_dir_sub,
                     loi_types,
                     num_rast_analysis_cols,
                     cat_rast_analysis_cols,
                     analysis_cols,
                     max.obj.fulldata,
                     p) {
              temp_fl<-tempfile(pattern ="ihdyro",tmpdir =temp_dir_sub, fileext = ".tif")
              #browser()
              sub_id<-sub_poly$link_id
              sub_poly_rast<-terra::rasterize(terra::vect(sub_poly),input_rasts,fun=sum,field=1,filename=temp_fl)
              sub_poly_rast<-terra::cells(sub_poly_rast)

              #Estimates the number of rows that could fit into memory and do analysis
              max.obj.fullanalysis_row<-memuse::howmany(sys.mem,
                                                        ncol=analysis_cols)
              max.obj.fullanalysis_col<-memuse::howmany(sys.mem,
                                                        nrow=length(sub_poly_rast))

              # Estimate table size if reading entire dataframe into memory
              obj.size.fulldata<-memuse::howbig(nrow=length(sub_poly_rast),
                                                ncol=length(loi_cols)+length(weighting_scheme[weighting_scheme!="lumped"]),
                                                unit="GiB")

              # Estimate table size if reading entire dataframe into memory
              obj.size.fullanalysis<-memuse::howbig(nrow=length(sub_poly_rast),
                                                    ncol=analysis_cols,
                                                    unit="GiB")

              # Estimate table size if single loi and IDW into memory
              obj.size.single<-memuse::howbig(nrow=length(sub_poly_rast),
                                              ncol=5,
                                              unit="GiB")

              if (obj.size.single>sys.mem) { # object won't fit into memory
                # TODO: process rasters in chunks and write to sqlite, then do calculations there
                out<-tibble::tibble(pour_point_id=sub_id,status="Incomplete")

              } else {
                if (obj.size.fullanalysis<sys.mem) { # object will fit entirely into memory
                  # browser()

                  ot<-terra::extract(
                    input_rasts,
                    sub_poly_rast
                  )

                  out<-ihydro::.attr_fn(ot,
                                        point_id=sub_id,
                                        weighting_scheme2=weighting_scheme,
                                        loi_meta2=loi_meta,
                                        loi_cols2=loi_cols,
                                        loi_numeric_stats2=loi_numeric_stats,
                                        backend=backend,
                                        SQLdb_path=SQLdb_path)

                } else {
                  if (max.obj.fullanalysis_col[2]>(length(weighting_scheme[weighting_scheme!="lumped"])*3)) { # all IDW plus at least some loi will fit into memory

                    weighting_scheme2<-weighting_scheme
                    if (all(weighting_scheme2=="lumped")) {
                      weighting_scheme2<-"lumped"
                    } else {
                      weighting_scheme2<-weighting_scheme[weighting_scheme!="lumped"]
                    }

                    if (all(weighting_scheme2=="lumped")){
                      ot_idw<-data.frame(lumped=1)
                    } else {
                      ot_idw<-terra::extract(
                        terra::subset(input_rasts,c(weighting_scheme2)[c(weighting_scheme2) %in% names(input_rasts)]),
                        sub_poly_rast
                      )


                      if (ncol(ot_idw)==1) colnames(ot_idw)<-weighting_scheme2
                    }

                    remaining_cols<-max.obj.fullanalysis_col[2]-(length(weighting_scheme[weighting_scheme!="lumped"])*3)
                    remaining_cols_ratio<-remaining_cols/analysis_cols

                    if (remaining_cols_ratio>1){
                      remaining_cols<-length(loi_cols)
                    } else {
                      remaining_cols<-floor(length(loi_cols)*remaining_cols_ratio*0.5)
                    }

                    loi_cols_split<-split(loi_cols,
                                          rep(1:length(loi_cols),
                                              length.out=length(loi_cols),
                                              each =remaining_cols)
                    )

                    out<-purrr::pmap(
                      list(
                        loi_sub=loi_cols_split,
                        input_rasts=list(input_rasts),
                        sub_poly_rast=list(sub_poly_rast),
                        sub_id=list(sub_id),
                        weighting_scheme=list(weighting_scheme),
                        loi_meta=list(loi_meta),
                        loi_numeric_stats=list(loi_numeric_stats),
                        backend=list(backend),
                        SQLdb_path=list(SQLdb_path)
                      ),
                      carrier::crate(
                        function(loi_sub,
                                 input_rasts,
                                 sub_poly_rast,
                                 sub_id,
                                 weighting_scheme,
                                 loi_meta,
                                 loi_numeric_stats,
                                 backend,
                                 SQLdb_path){

                          ot<-terra::extract(
                            terra::subset(input_rasts,c(loi_sub)[c(loi_sub) %in% names(input_rasts)]),
                            sub_poly_rast
                          )

                          if (ncol(ot)==1) colnames(ot)<-loi_sub

                          out<-ihydro::.attr_fn(dplyr::bind_cols(ot,ot_idw),
                                                point_id=sub_id,
                                                weighting_scheme2=weighting_scheme,
                                                loi_meta2=loi_meta,
                                                loi_cols2=loi_sub,
                                                loi_numeric_stats2=loi_numeric_stats,
                                                backend=backend,
                                                SQLdb_path=SQLdb_path)

                          return(out)
                        })) %>%
                      purrr::reduce(dplyr::left_join,by=c("pour_point_id","status"))

                  } else {
                    # one IDW plus at one loi will fit into memory

                    weighting_scheme2<-weighting_scheme
                    if (all(weighting_scheme2=="lumped")) {
                      weighting_scheme2<-"lumped"
                    } else {
                      weighting_scheme2<-weighting_scheme[weighting_scheme!="lumped"]
                    }

                    out<- purrr::pmap(
                      list(
                        sub_weighting_scheme=weighting_scheme2,
                        loi_cols=list(loi_cols),
                        input_rasts=list(input_rasts),
                        sub_poly_rast=list(sub_poly_rast),
                        sub_id=list(sub_id),
                        weighting_scheme=list(weighting_scheme),
                        loi_meta=list(loi_meta),
                        loi_numeric_stats=list(loi_numeric_stats),
                        backend=list(backend),
                        SQLdb_path=list(SQLdb_path)
                      ),
                      carrier::crate(
                        function(sub_weighting_scheme,
                                 loi_cols,
                                 input_rasts,
                                 sub_poly_rast,
                                 sub_id,
                                 weighting_scheme,
                                 loi_meta,
                                 loi_numeric_stats,
                                 backend,
                                 SQLdb_path) {

                          if (sub_weighting_scheme=="lumped"){
                            ot_idw<-data.frame(lumped=1)
                          } else {
                            ot_idw<-terra::extract(
                              terra::subset(input_rasts,c(sub_weighting_scheme)[c(sub_weighting_scheme) %in% names(input_rasts)]),
                              sub_poly_rast
                            )

                            if (ncol(ot_idw)==1) colnames(ot_idw)<-sub_weighting_scheme
                          }


                          purrr::pmap(
                            list(
                              loi_sub=loi_cols,
                              ot_idw=list(ot_idw),
                              sub_weighting_scheme=list(sub_weighting_scheme),
                              input_rasts=list(input_rasts),
                              sub_poly_rast=list(sub_poly_rast),
                              sub_id=list(sub_id),
                              weighting_scheme=list(weighting_scheme),
                              loi_meta=list(loi_meta),
                              loi_numeric_stats=list(loi_numeric_stats),
                              backend=list(backend),
                              SQLdb_path=list(SQLdb_path)
                            ),
                            carrier::crate(
                              function(loi_sub,
                                       ot_idw,
                                       sub_weighting_scheme,
                                       input_rasts,
                                       sub_poly_rast,
                                       sub_id,
                                       weighting_scheme,
                                       loi_meta,
                                       loi_numeric_stats,
                                       backend,
                                       SQLdb_path){

                                ot<-terra::extract(
                                  terra::subset(input_rasts,c(loi_sub)[c(loi_sub) %in% names(input_rasts)]),
                                  sub_poly_rast
                                )

                                if (ncol(ot)==1) colnames(ot)<-loi_sub

                                with_lumped<-"lumped" %in% weighting_scheme

                                if (with_lumped & sub_weighting_scheme==weighting_scheme2[1]){
                                  sub_weighting_scheme<-unique(c("lumped",sub_weighting_scheme))
                                }

                                out<-ihydro::.attr_fn(dplyr::bind_cols(ot,ot_idw),
                                                      point_id=sub_id,
                                                      weighting_scheme2=sub_weighting_scheme,
                                                      loi_meta2=loi_meta,
                                                      loi_cols2=loi_sub,
                                                      loi_numeric_stats2=loi_numeric_stats,
                                                      backend=backend,
                                                      SQLdb_path=SQLdb_path)
                                return(out)
                              })) %>%
                            purrr::reduce(dplyr::left_join,by=c("pour_point_id","status"))
                        })) %>%
                      purrr::reduce(dplyr::left_join,by=c("pour_point_id","status"))

                  }
                }
              }

              rm(sub_poly_rast)
              file.remove(temp_fl)
              #r<-terra::tmpFiles(current=FALSE, orphan=TRUE, old=FALSE, remove=TRUE)

              return(out)

            }))


        p()
        return(ot)
      }))
}

#' Calculate loi attributes from input rasters
#'
#' @keywords internal
#' @return a data.frame of attributes
#' @export
#'
.attr_fn<-function(df=NULL,
                   SQLdb_path=":memory:",
                   tbl_name=basename(tempfile(pattern = "ihydro", fileext = ".sql")),
                   point_id,
                   weighting_scheme2,
                   loi_meta2,
                   loi_cols2,
                   loi_numeric_stats2,
                   backend=c("tibble","data.table","SQLite")
){
  stopifnot(!is.null(df))
  options(dplyr.summarise.inform = FALSE)
  options(scipen = 999)
  `%>%` <- magrittr::`%>%`
  `:=` <- data.table::`:=`

  backend<-match.arg(backend)

  df<-df %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),~dplyr::if_else(is.nan(.),NA_real_,.)))

  loi_meta2<-dplyr::filter(loi_meta2,
                           loi_var_nms %in% loi_cols2,
                           loi_var_nms %in% colnames(df)
  )

  weighting_scheme2<-weighting_scheme2[weighting_scheme2 %in% c("lumped",colnames(df))]

  if (backend=="data.table") {

    df<-df %>%
      dtplyr::lazy_dt()
  }
  if (backend=="SQLite") {
    stopifnot(!is.null(SQLdb_path))
    if (SQLdb_path!=":memory:"){
      if (!dir.exists(SQLdb_path)) dir.create(SQLdb_path)
      SQLdb_path<-tempfile(tmpdir = SQLdb_path,pattern="ihydro",fileext = ".sql")
    }

    con<-DBI::dbConnect(RSQLite::SQLite(),SQLdb_path)

    if (SQLdb_path!=":memory:"){
      res <- DBI::dbSendQuery(con, paste0("PRAGMA temp_store = 1;"))
      res <- DBI::dbSendQuery(con, paste0("PRAGMA temp_store_directory='",gsub(basename(SQLdb_path),"",SQLdb_path),"';"))
      res <- DBI::dbSendQuery(con, paste0("PRAGMA SQLITE_TMPDIR='",gsub(basename(SQLdb_path),"",SQLdb_path),"';"))
    }

    df1<-dplyr::copy_to(con,df,name = tbl_name)

    df<-dplyr::tbl(con,tbl_name)

  }


  if ("iFLS" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2), # tidyselect::any_of
                                  ~.*(!!rlang::sym("iFLS")),.names=paste0("{.col}_","iFLS")))
  }
  if ("iFLO" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2), #tidyselect::any_of
                                  ~.*(!!rlang::sym("iFLO")),.names=paste0("{.col}_","iFLO")))
  }
  if ("HAiFLS" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2), #tidyselect::any_of
                                  ~.*(!!rlang::sym("HAiFLS")),.names=paste0("{.col}_","HAiFLS")))
  }
  if ("HAiFLO" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2), #tidyselect::any_of
                                  ~.*(!!rlang::sym("HAiFLO")),.names=paste0("{.col}_","HAiFLO")))
  }

  df<-df %>%
    dplyr::compute()


  weighted_mean_out<-NULL
  lumped_mean_out<-NULL
  weighted_sd_out<-NULL
  lumped_sd_out<-NULL
  min_out<-NULL
  max_out<-NULL
  count_out<-NULL
  median_out<-NULL
  sum_out<-NULL

  numb_rast<-loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]
  cat_rast<-loi_meta2$loi_var_nms[loi_meta2$loi_type=="cat_rast"]

  #browser()


  # Lumped Summaries --------------------------------------------------------

  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="mean")) {
    lumped_mean_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(
        dplyr::across(!!(numb_rast), #tidyselect::any_of
                      ~sum(.,na.rm=T)/sum(!is.na(.),na.rm=T)
        ),
        dplyr::across(!!(cat_rast), #tidyselect::any_of
                      ~sum(.,na.rm=T)/dplyr::n()
        )
      ) %>%
      dplyr::collect()

    if (length(numb_rast)>0) {
      lumped_mean_out<-lumped_mean_out %>%
        dplyr::rename_with(.cols=tidyselect::any_of(numb_rast),~paste0(.,"_lumped_mean"))
    }

    if (length(cat_rast)>0) {
      lumped_mean_out<-lumped_mean_out %>%
        dplyr::rename_with(.cols=tidyselect::any_of(cat_rast),~paste0(.,"_lumped_prop")) %>%
        dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.)))
    }

  }

  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="min") & length(numb_rast)>0){
    min_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~min(.,na.rm=T))) %>% #tidyselect::any_of
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_min"))

  }
  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="max" & length(numb_rast)>0)){
    max_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~max(.,na.rm=T))) %>% #tidyselect::any_of
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_max"))

  }
  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="count")){
    count_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::everything(),~sum(!is.na(.),na.rm=T))) %>% #tidyselect::everything
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_count"))

  }

  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="sum" & length(numb_rast)>0)){
    sum_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~sum(.,na.rm=T))) %>% #tidyselect::any_of
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_sum"))
  }

  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="median" & length(numb_rast)>0)){
    median_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::median(.,na.rm=T))) %>% #tidyselect::any_of
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_median"))
  }
  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2 %in% c("sd","stdev") & length(numb_rast)>0)){
    lumped_sd_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::sd(.,na.rm=T))) %>% #tidyselect::any_of
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_sd"))
  }

  # Weighted summaries -----------------------------------------------------------

  if (length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 & (any(loi_numeric_stats2 %in% c("mean")) | length(cat_rast)>0)) {
    weighted_mean_out<-df %>%
      dplyr::summarize(
        dplyr::across(tidyselect::ends_with(paste0("_iFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)), #tidyselect::any_of
        dplyr::across(tidyselect::ends_with(paste0("_HAiFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)), #tidyselect::any_of
        dplyr::across(tidyselect::ends_with(paste0("_iFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)), #tidyselect::any_of
        dplyr::across(tidyselect::ends_with(paste0("_HAiFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T)) #tidyselect::any_of
      ) %>%
      dplyr::collect()

    if (length(numb_rast)>0) {
      weighted_mean_out<-weighted_mean_out %>%
        dplyr::rename_with(.cols=tidyselect::starts_with(paste0(numb_rast,"_")),~paste0(.,"_mean"))
    }

    if (length(cat_rast)>0) {
      weighted_mean_out<-weighted_mean_out %>%
        dplyr::rename_with(.cols=tidyselect::starts_with(paste0(cat_rast,"_")),~paste0(.,"_prop")) %>%
        dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.)))
    }

  }


  if (
    length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 &
    any(loi_numeric_stats2 %in% c("sd","stdev")) &
    length(numb_rast)>0
  ) {

    weighted_sd_out<-df %>%
      dplyr::select(
        tidyselect::starts_with(numb_rast),
        tidyselect::any_of(weighting_scheme2)
      )

    if ("iFLS" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~(!!rlang::sym("iFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)))^2)),
                                    .names="{.col}_iFLS_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~ ((sum(!!rlang::sym("iFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLS")!=0,na.rm=T)) * sum(!!rlang::sym("iFLS"),na.rm=T),
                                    .names="{.col}_iFLS_term2"
                      ))
    }
    if ("HAiFLS" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~(!!rlang::sym("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)))^2)),
                                    .names="{.col}_HAiFLS_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~ ((sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLS"),na.rm=T),
                                    .names="{.col}_HAiFLS_term2"
                      ))
    }
    if ("iFLO" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~(!!rlang::sym("iFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)))^2)),
                                    .names="{.col}_iFLO_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~ ((sum(!!rlang::sym("iFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLO")!=0,na.rm=T)) * sum(!!rlang::sym("iFLO"),na.rm=T),
                                    .names="{.col}_iFLO_term2"
                      ))
    }
    if ("HAiFLO" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~(!!rlang::sym("HAiFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T)))^2)),
                                    .names="{.col}_HAiFLO_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),#tidyselect::any_of
                                    ~ ((sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLO"),na.rm=T),
                                    .names="{.col}_HAiFLO_term2"
                      ))
    }

    weighted_sd_out<- weighted_sd_out%>%
      dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),~sum(.,na.rm=T)),
                       dplyr::across(tidyselect::ends_with("_term2"),~.[1])
      ) %>%
      dplyr::collect() %>%
      # The below is some rearranging
      tidyr::pivot_longer(cols=c(tidyselect::everything())) %>%
      dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",2)[,1],
                    term=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",2)[,2]) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
      dplyr::group_by(name) %>%
      dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = name,values_from=sd)

  }


  final_list<-list(
    lumped_mean_out,
    lumped_sd_out,
    min_out,
    max_out,
    count_out,
    median_out,
    sum_out,
    weighted_mean_out,
    weighted_sd_out
  )
  final_list<-final_list[!sapply(final_list,is.null)]
  final_list<-final_list[sapply(final_list,nrow)>0]

  final_out<-dplyr::bind_cols(tibble::tibble(pour_point_id=point_id,status="Complete"),final_list)

  final_out<-final_out %>%
    dplyr::select(
      tidyselect::any_of("pour_point_id"),
      tidyselect::any_of("status"),
      tidyselect::contains(loi_meta2$loi_var_nms)
    )

  if (backend=="SQLite") {
    DBI::dbDisconnect(con)

    if (SQLdb_path!=":memory:"){
      file.remove(SQLdb_path)
    }
  }

  return(final_out)
}



