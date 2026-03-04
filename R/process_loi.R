
#' Processes layers of interest (loi) for attribution of stream network
#'
#' @param input \code{NULL} or output from `process_hydrology()` or `process_flowdir()`. If \code{NULL}, `dem` must be specified.
#' @param dem \code{NULL} or character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster. If \code{NULL}, input must be specified.
#' @param clip_region character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of ESRI Shapefile type or GeoTiFF type only. Region over which loi are returned. defaults to `dem` extent.
#' @param num_inputs named list containing file paths (with extensions), and/or R GIS objects to be coerced into \code{SpatRaster}. Layers of interest to be summarized numerically (e.g., with mean, SD, etc.)
#' @param cat_inputs named list containing file paths (with extensions), and/or R GIS objects to be coerced into \code{SpatRaster}. Layers of interest to be summarized categorical (e.g., with with proportions or weighted proportions in the catchment)
#' @param output_filename \code{NULL} or character (full file path with extensions). If `input` is provided, loi are saved in the same .zip output folder regadless of `output_filename`. If `input` not provided, file path to write resulting .zip file.
#' @param variable_names named list containing variables from `num_inputs` and `cat_inputs` to include in output. If not specied, all variables are used.
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export

process_loi<-function(
    input=NULL,
    dem=NULL,
    clip_region=NULL,
    num_inputs=list(),
    cat_inputs=list(),
    output_filename=NULL,
    variable_names=NULL,
    return_products=F,
    temp_dir=NULL,
    verbose=F,
    overwrite=T
) {
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

  if (!is.null(input)) {
    if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  }

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.logical(overwrite)) stop("'overwrite' must be logical")
  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  if (!inherits(num_inputs,"list")) stop("'num_inputs' should be a named list")
  if (length(num_inputs)>0 & is.null(names(num_inputs))) stop("objects in 'num_inputs' must be named")
  if (!inherits(cat_inputs,"list")) stop("'cat_inputs' should be a named list")
  if (length(cat_inputs)>0 & is.null(names(cat_inputs))) stop("objects in 'cat_inputs' must be named")

  if (any(names(num_inputs) %in% names(cat_inputs)) |
      any(names(cat_inputs) %in% names(num_inputs))) stop("'num_inputs' and 'cat_inputs' cannot share names")

  if (!is.null(variable_names)) {
    # Add a check to make sure names match
  } else {
    message("No variables specified in 'variable_names', all variables in 'num_inputs' and 'cat_inputs' will be used")
  }

  if (is.null(input) & is.null(dem)) stop("Either 'input' or 'dem' must be specifed")
  if (is.null(input) & is.null(output_filename)) stop("Either 'input' or 'output_filename' must be specifed")

  if (!is.null(input) & is.null(output_filename)) {
    output_filename<-input$outfile
  }
  output_filename<-normalizePath(output_filename,mustWork =F)
  #if (!grepl("\\.zip$",output_filename)) stop("output_filename must be a character ending in '.zip'")
  if (!grepl("\\.gpkg$",output_filename)) stop("output_filename must be a character ending in '.gpkg'")

  if (gsub(basename(output_filename),"",output_filename) == temp_dir) stop("'output_filename' should not be in the same directory as 'temp_dir'")

  tt<-suppressWarnings(dir.create(gsub(basename(output_filename),"",output_filename)))

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),
                        verbose=verbose>2,
                        wd=temp_dir)

  terra::terraOptions(verbose = verbose>3,
                      tempdir = temp_dir
  )

  # Prepare DEM -------------------------------------------------------------
  if (verbose) message("Preparing DEM")
  if (!is.null(input)){
    zip_loc<-input$outfile

    dem<-terra::rast(input$outfile,lyrs="dem_final")
    terra::writeRaster(dem,file.path(temp_dir,"dem_final.tif"),overwrite=T,gdal="COMPRESS=NONE")

  } else {
    dem<-hydroweight::process_input(dem,input_name="dem",working_dir=temp_dir)
    if (!inherits(dem,"SpatRaster")) stop("dem must be a class 'SpatRaster'")
    terra::writeRaster(dem,file.path(temp_dir,"dem_final.tif"),overwrite=T,gdal="COMPRESS=NONE")
  }

  if (is.null(clip_region)) {
    clip_region<-terra::as.polygons(terra::rast(terra::ext(dem),crs=terra::crs(dem)))
    terra::writeVector(clip_region,file.path(temp_dir,"clip_region.shp"),overwrite=T)
  } else {
    clip_region<-hydroweight::process_input(clip_region,
                                            input_name="clip_region",
                                            working_dir=temp_dir,
                                            align_to=terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",crs=terra::crs(dem)))

    terra::writeVector(clip_region,file.path(temp_dir,"clip_region.shp"),overwrite=T)
  }

  target_crs<-terra::crs(dem)

  # Prepare loi -------------------------------------------------------------

  num_inputs<-purrr::map(num_inputs,function(x){
    if (inherits(x,"character")) return(x)
    if (inherits(x,c("sf","sfc","sfg","SpatVector"))) {
      fp<-file.path(temp_dir,paste0(basename(tempfile()),".shp"))
      if (inherits(x,"SpatVector")) x<-sf::st_as_sf(x)
      sf::write_sf(x,fp)
      x<-fp
      #return(terra::wrap(terra::vect(x)))
    }
    if (inherits(x,c("SpatRaster"))) {
      fp<-file.path(temp_dir,paste0(basename(tempfile()),".tif"))
      terra::writeRaster(x,fp)
      x<-fp
      #return(terra::wrap(x))
    }
    return(x)
  })

  cat_inputs<-purrr::map(cat_inputs,function(x){
    if (inherits(x,"character")) return(x)
    if (inherits(x,c("sf","sfc","sfg","SpatVector"))) {
      fp<-file.path(temp_dir,paste0(basename(tempfile()),".shp"))
      if (inherits(x,"SpatVector")) x<-sf::st_as_sf(x)
      sf::write_sf(x,fp)
      x<-fp
      #return(terra::wrap(terra::vect(x)))
    }
    if (inherits(x,c("SpatRaster"))) {
      fp<-file.path(temp_dir,paste0(basename(tempfile()),".tif"))
      terra::writeRaster(x,fp)
      x<-fp
      #return(terra::wrap(x))
    }
    return(x)
  })

  inputs_list<-list(
    lyr_nms = as.list(c(names(num_inputs),names(cat_inputs))),
    lyr = c(num_inputs,cat_inputs),
    lyr_variables = c(variable_names[names(num_inputs)],variable_names[names(cat_inputs)]),
    rln = as.list(c(rep("num_rast",length(num_inputs)), rep("cat_rast",length(cat_inputs)))),
    temp_dir = as.list(rep(temp_dir,length(c(num_inputs,cat_inputs)))),
    overwrite=as.list(rep(overwrite,length(c(num_inputs,cat_inputs)))),
    output_filename = as.list(rep(output_filename,length(c(num_inputs,cat_inputs))))
  )

  if (!is.null(inputs_list$lyr_variables[sapply(inputs_list$lyr_variables,is.null)])){
    inputs_list$lyr_variables[sapply(inputs_list$lyr_variables,is.null)]<-NA_character_
  } else {
    inputs_list$lyr_variables<-rep(list(NA_character_),length(inputs_list$lyr_nms))
  }


  if (!file.exists(output_filename)){
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

  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1
  max_cores_opt<-getOption("parallelly.maxWorkers.localhost")
  options(parallelly.maxWorkers.localhost = n_cores)

  n_cores_2<-n_cores

  if (n_cores>1) {
    n_cores_2<-n_cores_2-1
    oplan <- future::plan(list(future::tweak(future::multisession, workers = 2),
                               future::tweak(future::multisession, workers = n_cores_2)))
    on.exit(future::plan(oplan), add = TRUE)
  }

  temp_dir_save<-file.path(temp_dir,basename(tempfile()))
  dir.create(temp_dir_save)


  # Check if loi exists -----------------------------------------------------

  # if (file.exists(output_filename)){
  #   con<-DBI::dbConnect(RSQLite::SQLite(),output_filename)
  #   tl<-DBI::dbListTables(con)
  #
  #   tbl_lst<-dplyr::tbl(con,"gpkg_tile_matrix_set") %>%
  #     dplyr::collect()
  #
  #   if (i %in% tbl_lst$table_name) {
  #     if (overwrite) {
  #       DBI::dbSendStatement(con, paste0("DELETE FROM gpkg_tile_matrix_set where table_name=","'",i,"'"))
  #       DBI::dbSendStatement(con, paste0("DELETE FROM gpkg_tile_matrix where table_name=","'",i,"'"))
  #       DBI::dbSendStatement(con, paste0("DELETE FROM gpkg_2d_gridded_tile_ancillary where tpudt_name=","'",i,"'"))
  #       DBI::dbSendStatement(con, paste0("DELETE FROM gpkg_2d_gridded_coverage_ancillary where tile_matrix_set_name=","'",i,"'"))
  #       DBI::dbSendStatement(con, paste0("DELETE FROM gpkg_contents where table_name=","'",i,"'"))
  #       DBI::dbSendStatement(con, paste0("DROP TABLE ",i))
  #
  #     } else {
  #       warning(paste0("Layer ",i," could not be writted to gpkg because a table with that name already exists."))
  #       p()
  #       return(NULL)
  #     }
  #   }
  #
  #   DBI::dbDisconnect(con)
  # }


  progressr::with_progress(enable=T,{
    message("Processing loi")
    p <- progressr::progressor(steps = (length(num_inputs) + length(cat_inputs)))

    ip<-c(inputs_list,
          list(p=rep(list(p),length(inputs_list$output_filename)),
               temp_dir_save=rep(list(temp_dir_save),length(inputs_list$output_filename))
          ))

    future_proc<-future::future({
      options(parallelly.maxWorkers.localhost = n_cores_2+1)
      ot<-furrr::future_pmap(ip,
                             #purrr::pmap(ip,
                             .options=furrr::furrr_options(globals = F),
                             carrier::crate(
                               function(lyr_nms=lyr_nms,
                                        lyr=lyr,
                                        lyr_variables=lyr_variables,
                                        gdal_arg=gdal_arg,
                                        rln=rln,
                                        temp_dir=temp_dir,
                                        temp_dir_save=temp_dir_save,
                                        output_filename=output_filename,
                                        overwrite=overwrite,
                                        p=p){
                                 #browser()

                                 purrr::pmap(list(lyr_nms=lyr_nms,
                                                  lyr=lyr,
                                                  lyr_variables=lyr_variables,
                                                  rln=rln,
                                                  temp_dir=temp_dir,
                                                  temp_dir_save=temp_dir_save,
                                                  output_filename=output_filename),
                                             function(lyr_nms,
                                                      lyr,
                                                      lyr_variables,
                                                      rln,
                                                      temp_dir,
                                                      temp_dir_save,
                                                      output_filename){
                                               options(scipen = 999)
                                               `%>%` <- magrittr::`%>%`

                                               temp_temp_dir<-file.path(temp_dir,basename(tempfile()))
                                               dir.create(temp_temp_dir)

                                               resaml<-ifelse(grepl("num_rast",rln),"bilinear","near")

                                               if (all(is.na(lyr_variables))) lyr_variables<-NULL

                                               #browser()

                                               suppressMessages(
                                                 output<-hydroweight::process_input(
                                                   input=unlist(lyr),
                                                   input_name = unlist(lyr_nms),
                                                   #variable_names=unlist(lyr_variables),
                                                   #target=file.path(temp_dir,"dem_final.tif"),
                                                   input_variable_names = unlist(lyr_variables),
                                                   align_to = file.path(temp_dir,"dem_final.tif"),
                                                   clip_region = file.path(temp_dir,"clip_region.shp"),
                                                   resample_type = resaml,
                                                   working_dir=temp_temp_dir
                                                 )
                                               )

                                               out_files<-file.path(temp_dir,paste0(rln,"_",names(output),".tif"))
                                               names(out_files)<-names(output)

                                               output<-terra::split(output,names(output))
                                               names(output)<-sapply(output,names)

                                               for (i in names(output)){

                                                 t1<-terra::writeRaster(
                                                   output[[i]],
                                                   file.path(temp_dir_save,paste0(names(output[[i]]),".tif")),
                                                   datatype="FLT4S",
                                                   todisk = T,
                                                   overwrite=T,
                                                   gdal="COMPRESS=NONE"
                                                 )


                                               }

                                               file.remove(list.files(temp_temp_dir,recursive=T,full.names = T))
                                               unlink(temp_temp_dir,recursive=T,force=T)

                                               p()

                                               out<-list(
                                                 lyr_nms=lyr_nms,
                                                 lyr_variables=sapply(output,names),
                                                 rln=rln
                                               )

                                               return(out)
                                             })



                               }))
    })

  })

  if (verbose) message("Generating Outputs")

  future_proc_status <- future::futureOf(future_proc)
  while(!future::resolved(future_proc_status)){
    Sys.sleep(0.5)
    fl<-list.files(temp_dir_save,".tif",full.names = T)
    fl<-fl[grepl(".tif$",fl)]

    fl_un_time<-file.mtime(fl)
    fl<-fl[fl_un_time<Sys.time()-60]


    for (x in fl) {
      tot<-try(terra::rast(x),silent=T)
      if (inherits(tot,"try-error")) next()
      if (verbose) message(paste0("Writing: ",names(tot)))

      tott<-try(writeRaster_fun(tot,output_filename),silent=T)

      if (inherits(tott,"try-error")) {
        if (!attr(tott,"condition")$message %in% c("stoi","stol")){
          stop(attr(tott,"condition")$message)
        }
      }

      file.remove(x)
    }
  }

  Sys.sleep(60)

  if (length(future_proc$result$conditions)>0){
    err<-lapply(future_proc$result$conditions,function(x) x$condition)
    err<-err[sapply(err,function(x) inherits(x,"error"))]
    if (length(err)>0){
      stop(paste0(unlist(err)))
    }
  }

  fl<-list.files(temp_dir_save,".tif",full.names = T)
  fl<-fl[grepl(".tif$",fl)]

  for (x in fl) {
    tot<-terra::rast(x)
    if (verbose) message(paste0("Writing: ",names(tot)))
    tott<-try(writeRaster_fun(tot,output_filename),silent=T)

    if (inherits(tott,"try-error")) {
      if (!attr(tott,"condition")$message %in% c("stoi","stol")){
        stop(attr(tott,"condition")$message)
      }
    }

    file.remove(x)
  }


  ot<-future::value(future_proc)
  names(ot)<-unlist(inputs_list$lyr_nms)

  ot<-lapply(ot,function(x) unlist(x,recursive=F))

  meta<-map_dfr(ot,~tibble::tibble(
    loi_lyr_nms=.$lyr_nms,
    loi_var_nms=.$lyr_variables,
    loi_type=.$rln
  ))

  if (verbose) message("Saving meta data")

  con<-DBI::dbConnect(RSQLite::SQLite(),output_filename)

  tmp<-copy_to(dest=con,
               df=meta,
               name="loi_meta",
               overwrite =T,
               temporary =F,
               analyze=T,
               in_transaction=T)

  DBI::dbDisconnect(con)

  # Generate Output ---------------------------------------------------------

  output<-input

  if (is.null(output)){
    output<-list(outfile=output_filename)
  }

  if (return_products){

    ott<-purrr::map(meta %>% split(.,.$loi_type),~purrr::map(.$loi_var_nms,~terra::rast(output_filename,lyrs=.)))
    ott2<-list(
      num_inputs=try(terra::rast(unlist(ott["num_rast"],use.names = F)),silent=T),
      cat_inputs=try(terra::rast(unlist(ott["cat_rast"],use.names = F)),silent=T)
    )

    ott2<-ott2[!sapply(ott2,function(x) inherits(x,"try-error"))]

    output<-c(
      ott2,
      output
    )
  }

  suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive = T)))

  class(output)<-"ihydro"

  options(parallelly.maxWorkers.localhost=max_cores_opt)

  return(output)
}
