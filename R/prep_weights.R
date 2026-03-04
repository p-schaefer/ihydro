
#' Prepares hydroweights and adds them to existing zip file
#'
#' @param input an object of class ihydro.
#' @param sample_points character or NULL. IDs of unique station identifiers priveded in 'site_id_col' of `generate_vectors()`
#' @param link_id character or NULL. 'link_id's of reaches to calculate attributes for.
#' @param target_o_type character. One of: c("point","segment_point","segment_whole"). Target for iEucO" "iFLO", and "HAiFLO" weighting schemes. "Point" represents the sampling point on the stream, "segment_point" represents the upstream segment of the sampling points, and "segment_whole" will target the entire reach, regardless of where sampling occurred.
#' @param weighting_scheme character. One or more weighting schemes: c( "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return input
#' @export
#'

prep_weights<-function(
    input,
    output_filename=NULL,
    sample_points=NULL,
    link_id=NULL,
    target_o_type=c("point","segment_point","segment_whole"),
    weighting_scheme =  c("iFLS", "HAiFLS","iFLO",  "HAiFLO"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir=NULL,
    verbose=F
){
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  if (inherits(output_filename,"ihydro")) output_filename<-output_filename$outfile

  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1
  max_cores_opt<-getOption("parallelly.maxWorkers.localhost")
  options(parallelly.maxWorkers.localhost = n_cores)

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  if (is.null(target_o_type)) target_o_type<-"point"
  if (length(target_o_type)>1) target_o_type<-target_o_type[[1]]
  target_o_type<-match.arg(target_o_type,several.ok = F)

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),
                        verbose=verbose>2,
                        wd=temp_dir)

  terra::terraOptions(verbose = verbose>3,
                      tempdir = temp_dir
  )

  # Determine output --------------------------------------------------------

  db_loc<-zip_loc<-input$outfile

  if (is.null(output_filename) | output_filename==db_loc) {
    output_filename<-db_loc

    if (!file.exists(output_filename)){
      t1<-try(dir.create(gsub(basename(output_filename),"",output_filename)),silent=T)
    }

  } else {
    if (!grepl("\\.gpkg$",output_filename)) stop("output_filename must be a character ending in '.gpkg'")

    if (!file.exists(output_filename)){
      t1<-try(dir.create(gsub(basename(output_filename),"",output_filename)),silent=T)

      dem<-terra::rast(db_loc,"dem_final")
      dem_ext<-dem %>%
        terra::ext() %>%
        terra::as.polygons(crs=terra::crs(dem)) %>%
        sf::st_as_sf() %>%
        sf::write_sf(output_filename,
                     layer="DEM_Extent",
                     append = T,
                     delete_layer = F,
                     delete_dsn = F)
    }
  }

  output_filename<-as.ihydro(output_filename)
  lyrs<-ihydro_layers(output_filename)

  # Organize weighting schemes ----------------------------------------------
  weighting_scheme<-match.arg(weighting_scheme,several.ok = T)
  weighting_scheme<-weighting_scheme[weighting_scheme!="lumped"]

  weighting_scheme_s<-weighting_scheme[grepl("FLS",weighting_scheme)]
  if (all(weighting_scheme_s %in% lyrs$layer_name)) message("Stream Targeted Weights already exist in output and won't be recalculated")
  weighting_scheme_s<-weighting_scheme_s[!weighting_scheme_s %in% lyrs$layer_name]

  weighting_scheme_o<-weighting_scheme[!grepl("lumped|FLS",weighting_scheme)]

  # Determine what sites/layers needs to be calculated ----------------------
  target_IDs<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=target_o_type=="segment_whole",
    target_o_type=target_o_type
  )

  if (any(lyrs=="target_o_meta")) {
    avail_weights<-sf::read_sf(output_filename,"target_o_meta") %>%
      dplyr::mutate(link_id=as.character(link_id)) %>%
      dplyr::full_join(target_IDs %>% dplyr::select(link_id),by="link_id")

    avail_weights<-lapply(setNames(weighting_scheme_o,weighting_scheme_o),function(x){
      avail_weights %>%
        dplyr::mutate(weight=x) %>%
        dplyr::mutate(layer_name=paste0(x,"_unn_group",unn_group)) %>%
        dplyr::mutate(layer_name_exists = layer_name %in% lyrs$layer_name)
    }) %>%
      dplyr::bind_rows()
  } else {
    avail_weights<-lapply(setNames(weighting_scheme_o,weighting_scheme_o),function(x){
      tibble::tibble(link_id=target_IDs$link_id) %>%
        dplyr::mutate(weight=x) %>%
        dplyr::mutate(layer_name=NA_character_) %>%
        dplyr::mutate(layer_name_exists = F)
    }) %>%
      dplyr::bind_rows()
  }


  # Check if remaining function needs to be run: ----------------------------
  if (length(weighting_scheme_s)==0 & all(avail_weights$layer_name_exists==T)){
    if (!inherits(output_filename,"ihydro")) output_filename<-ihydro::as_ihydro(output_filename)
    return(output_filename)
  }

  temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
  dir.create(temp_dir_sub)


  # Write rasters to temporary locatio --------------------------------------
  for (i in c("dem_streams_d8",
              "dem_final",
              "dem_accum_d8",
              "dem_d8")) {
    terra::writeRaster(
      terra::rast(zip_loc,i),
      file.path(temp_dir_sub,paste0(i,".tif")),
      overwrite=T
    )
  }

  # Calculate weighted S-target distances -------------------------------------

  if (length(weighting_scheme_s) > 0) {
    temp_dir_sub2<-file.path(temp_dir_sub,basename(tempfile()))
    dir.create(temp_dir_sub2)

    if (verbose) message("Generating Stream Targeted Weights")
    suppressMessages(
      hw_streams<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub2,
                                           target_O = NULL,
                                           target_S = file.path(temp_dir_sub,paste0("dem_streams_d8",".tif")),
                                           target_uid = 'ALL',
                                           OS_combine = FALSE,
                                           dem=file.path(temp_dir_sub,paste0("dem_final",".tif")),
                                           flow_accum = file.path(temp_dir_sub,paste0("dem_accum_d8",".tif")),
                                           weighting_scheme = weighting_scheme_s[!weighting_scheme_s %in% lyrs$layer_name],
                                           inv_function = inv_function,
                                           clean_tempfiles=F,
                                           return_products = T,
                                           wrap_return_products=F,
                                           save_output=F#,
                                           #temp_dir=temp_dir_sub2
      )
    )

    for (i in hw_streams) {
      o <- writeRaster_fun(i,output_filename$outfile)
    }

    rm(hw_streams)
    t1<-try((unlink(temp_dir_sub2,force = T,recursive = T)),silent=T)
  }

  # Calculate weighted O-target distances -------------------------------------

  rand_id <- function(n = 1,d=12) {
    paste0(replicate(d, sample(c(LETTERS,letters,0:9),n,TRUE)),collapse = "")
  }

  if (length(weighting_scheme_o) > 0) {
    #browser()

    temp_dir_sub2<-file.path(temp_dir_sub,basename(tempfile()))
    dir.create(temp_dir_sub2)

    target_o_meta<-NULL
    if ("target_o_meta" %in% lyrs$layer_name){
      target_o_meta<-sf::read_sf(output_filename,"target_o_meta") %>%
        dplyr::filter(paste0(weighting_scheme_o,"_unn_group",unn_group) %in% lyrs$layer_name)
    }

    target_O<-target_o_fun(
      db_fp=db_loc,
      target_IDs=target_IDs,
      target_o_type=target_o_type
    ) %>%
      dplyr::filter(!link_id %in% target_o_meta$link_id)

    all_points<-target_o_fun(
      db_fp=db_loc,
      target_IDs=target_IDs,
      target_o_type="point"
    )%>%
      dplyr::filter(!link_id %in% target_o_meta$link_id)

    sf::write_sf(all_points %>% dplyr::select(link_id),
                 file.path(temp_dir_sub2,"pour_points.shp"),
                 overwrite=T)

    if (nrow(target_O)>0){
      if (verbose) message("Generating Site Targeted Weights")
      if (verbose) message("Calculation for iFLO, and HAiFLO can be slow")

      # Unnest basins -----------------------------------------------------------

      # Use unnest basins to find catchments that don't overlap
      # Use asynchronous evaluation (if parallel backend registered)
      # to clear up rasters as they are made
      future_unnest<-future::future({
        whitebox::wbt_unnest_basins(
          d8_pntr=file.path(temp_dir_sub,"dem_d8.tif"),
          pour_pts=file.path(temp_dir_sub2,"pour_points.shp"),
          output=file.path(temp_dir_sub2,"unnest.tif")
        )
      })

      future_unnest_status <- future::futureOf(future_unnest)

      rast_out<-list()
      while(!future::resolved(future_unnest_status)){
        Sys.sleep(0.2)
        fl_un<-list.files(temp_dir_sub2,"unnest_",full.names = T)

        if (length(fl_un)==0) next

        rast_all<-purrr::map(fl_un,function(x) try(terra::rast(x),silent=T))
        rast_all<-rast_all[!sapply(rast_all,function(x) inherits(x,"try-error"))]

        if (length(rast_all)>0){
          rast_out<-c(rast_out,purrr::map(rast_all,terra::unique))
          suppressWarnings(file.remove(unlist(purrr::map(rast_all,terra::sources))))
        }
      }

      if (length(future_unnest$result$conditions)>0){
        err<-future_unnest$result$conditions[[1]]$condition
        if (inherits(err,"error")){
          stop(err)
        }
      }

      fl_un<-list.files(temp_dir_sub2,"unnest_",full.names = T)
      fl_un<-fl_un[grepl(".tif",fl_un)]
      rast_all<-purrr::map(fl_un,function(x) try(terra::rast(x),silent=T))
      rast_all<-rast_all[!sapply(rast_all,function(x) inherits(x,"try-error"))]
      if (length(rast_all)>0){
        rast_out<-c(rast_out,purrr::map(rast_all,terra::unique))
        suppressWarnings(file.remove(unlist(purrr::map(rast_all,terra::sources))))
      }

      link_id_split<-purrr::map(rast_out,~all_points$link_id[.[[1]]])

      # Calculating O target weights --------------------------------------------

      n_cores_2<-n_cores

      if (n_cores>1) {
        n_cores_2<-n_cores_2-1
        oplan <- future::plan(list(future::tweak(future::multisession, workers = 2),
                                   future::tweak(future::multisession, workers = n_cores_2)))
        on.exit(future::plan(oplan), add = TRUE)
      }

      target_O_sub<-purrr::map(link_id_split,~dplyr::filter(target_O,link_id %in% .) %>%
                                 dplyr::select(link_id) %>%
                                 dplyr::mutate(unn_group=rand_id(1)))
      splt_lst<-suppressWarnings(terra::split(target_O_sub,1:n_cores_2))
      splt_lst<-splt_lst[!sapply(splt_lst,is.null)]


      future_proc<-future::future({
        options(parallelly.maxWorkers.localhost = n_cores_2+1)
        hw_o_targ<-furrr::future_pmap(
          #hw_o_targ<-purrr::pmap(
          list(x=splt_lst,
               output_filename=list(output_filename),
               temp_dir_sub=list(temp_dir_sub),
               temp_dir_sub2=list(temp_dir_sub2),
               weighting_scheme_o=list(weighting_scheme_o),
               inv_function=list(inv_function),
               verbose=list(verbose)
          ),
          .options = furrr::furrr_options(globals = F),
          carrier::crate(
            function(x,
                     output_filename,
                     temp_dir_sub,
                     temp_dir_sub2,
                     weighting_scheme_o,
                     inv_function,
                     verbose
            ){
              options(dplyr.summarise.inform = FALSE)
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              target_S <- terra::rast(file.path(temp_dir_sub,"dem_streams_d8.tif"))
              dem <- terra::rast(file.path(temp_dir_sub,"dem_final.tif"))
              flow_accum <- terra::rast(file.path(temp_dir_sub,"dem_accum_d8.tif"))

              o_out<-purrr::map(x,function(y){
                if(verbose) {
                  message(paste("Calculating weights for:",y$unn_group[[1]])) #

                }
                temp_dir_sub_sub<-file.path(temp_dir_sub2,basename(tempfile()))
                dir.create(temp_dir_sub_sub)
                suppressMessages(
                  hw_o<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub_sub,
                                                 target_O = y,
                                                 target_S = target_S,
                                                 target_uid = paste0("unnest_group_",y$unn_group[[1]]),
                                                 OS_combine = FALSE,
                                                 dem=dem,
                                                 flow_accum = flow_accum,
                                                 weighting_scheme = weighting_scheme_o,
                                                 inv_function = inv_function,
                                                 clean_tempfiles=F,
                                                 return_products = T,
                                                 wrap_return_products=F,
                                                 save_output=F#,
                                                 #temp_dir=temp_dir_sub_sub
                  )
                )

                for (i in hw_o) {
                  #i[is.na(i)]<-(-9999)

                  t1<-terra::writeRaster(
                    i,
                    file.path(temp_dir_sub2,paste0(names(i),"_unn_group",y$unn_group[[1]],".tif")),
                    datatype="FLT4S",
                    todisk = T,
                    overwrite=T,
                    gdal="COMPRESS=NONE"
                  )


                }

                rm(hw_o)
                file.remove(list.files(temp_dir_sub_sub,recursive=T,full.names = T))
                t1<-try((unlink(temp_dir_sub_sub,force = T,recursive = T)),silent=T)
                #p()
                return(NULL)

                #return(unlist(rout))
              })

            }))
      })

      progressr::with_progress(enable=T,{
        p <- progressr::progressor(steps = (length(target_O_sub)*length(weighting_scheme_o)))
        future_proc_status <- future::futureOf(future_proc)

        while(!future::resolved(future_proc_status)){
          Sys.sleep(0.2)
          fl_un<-list.files(temp_dir_sub2,full.names = T)
          fl_un<-fl_un[grepl(paste0(weighting_scheme_o,collapse = "|"),fl_un)]

          fl_un_time<-file.mtime(fl_un)
          fl_un<-fl_un[fl_un_time<(Sys.time()-60)]

          if (length(fl_un)>0) {
            rast_all<-purrr::map(fl_un,function(y) {
              Sys.sleep(0.2)
              x<-try(terra::rast(y),
                     silent = T)

              if (inherits(x,"try-error")) return(NULL)


              ot <- try(
                writeRaster_fun(x,
                                output_filename$outfile,
                                gsub("\\.tif$","",basename(y))
                                ),
                silent=T)

              if (inherits(ot,"try-error")) {
                if (!attr(ot,"condition")$message %in% c("stoi","stol")){
                  stop(attr(ot,"condition")$message)
                }
              }

              p()
              return(y)
            })

            rast_all<-rast_all[!sapply(rast_all,is.null)]

            try(suppressWarnings(file.remove(unlist(rast_all))),silent=T)

          }

        }

        Sys.sleep(60)

        if (length(future_proc$result$conditions)>0){
          err<-lapply(future_proc$result$conditions,function(x) x$condition)
          err<-err[sapply(err,function(x) inherits(x,"error"))]
          if (length(err)>0){
            return(err)
          }
        }

        fl_un<-list.files(temp_dir_sub2,full.names = T)
        fl_un<-fl_un[grepl(paste0(weighting_scheme_o,collapse = "|"),fl_un)]

        if (length(fl_un)>0) {
          rast_all<-purrr::map(fl_un,function(y) {

            x<-terra::rast(y)

            ot <- try(
              writeRaster_fun(x,
                              output_filename$outfile,
                              gsub("\\.tif$","",basename(y))
                              ),
              silent=T)

            if (inherits(ot,"try-error")) {
              if (!attr(ot,"condition")$message %in% c("stoi","stol")){
                stop(attr(ot,"condition")$message)
              }
            }

            p()
            return(y)
          })
          rast_all<-rast_all[!sapply(rast_all,is.null)]

          try(suppressWarnings(file.remove(unlist(rast_all))),silent=T)
        }


        t1<-try((unlink(temp_dir_sub2,force = T,recursive = T)),silent=T)
      })

      if (verbose) message("Saving meta data")

      dplyr::bind_rows(target_O_sub) %>%
        tibble::as_tibble() %>%
        dplyr::select(link_id,unn_group) %>%
        sf::write_sf(
          output_filename$outfile,
          layer="target_o_meta",
          append=T,
          delete_layer = F,
          delete_dsn=F
        )

    }


  }

  unlink(temp_dir,recursive = T,force = T)

  output<-input
  if (output$outfile != output_filename$outfile) output<-list(outfile=output_filename$outfile)

  #class(output)<-"ihydro"

  if (!inherits(output,"ihydro")) output<-as.ihydro(output$outfile)

  options(parallelly.maxWorkers.localhost=max_cores_opt)

  return(output)


}

