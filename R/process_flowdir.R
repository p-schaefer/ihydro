
#' Process flow direction/accumulation and extracts streams from DEM
#'
#' If `return_products` = T, Packed-SpatRaster objects are returned. If `burn_streams` is specified, [whitebox::wbt_fill_burn] is used to burn the stream layer into the DEM. This is followed by [whitebox::wbt_feature_preserving_smoothing] to smooth final DEM. The `depression_corr` argument can be used to apply filling or breaching approaches to enforce flow accumulation from internal pit cells.
#'
#' @param dem character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster.
#' @param threshold integer. Flow accumulation threshold for stream initiation.
#' @param burn_streams character. (full file path with extension, e.g., "C:/Users/Administrator/Desktop/input.shp"), sf, SpatVector, PackedSpatVector, RasterLayer, SpatRaster, or PackedSpatRaster. Stream vector to burn into DEM.
#' @param burn_depth numeric. Depth (in meters) to burn stream into the DEM.
#' @param min_length numeric. Minimum tributary length, shorter 1st order tributaries are removed.
#' @param depression_corr NULL or character. One of c("fill","breach"), specifying whether depressions should be filled or breached. NULL will perform neither, if DEM is already corrected.
#' @param output_filename character. Full file path (with extension, e.g., "C:/Users/Administrator/Desktop/out.zip") to write resulting .zip file.
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param compress logical. Should output rasters be compressed, slower but more space efficient.
#' @param verbose logical.
#'
#' @seealso [whitebox::wbt_d8_pointer], [whitebox::wbt_d8_flow_accumulation], [whitebox::wbt_extract_streams]
#'
#' @return An object of class ihydro. If \code{return_products = TRUE}, all geospatial analysis products are returned as well.
#' @export
#'

process_flowdir<-function(
    dem,
    threshold,
    burn_streams=NULL,
    burn_depth=NULL,
    min_length=NULL,
    depression_corr=c(NULL,"fill","breach"),
    output_filename,
    return_products=F,
    temp_dir=NULL,
    compress=F,
    verbose=F
) {
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  if (!is.integer(threshold)) stop("'threshold' must be an integer value")
  if (!is.null(burn_depth) & !is.numeric(burn_depth)) stop("'burn_depth' must be an numeric value")
  if (!is.null(min_length) & !is.numeric(min_length)) stop("'min_length' must be an numeric value")

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(compress)) stop("'compress' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)
  output_filename<-normalizePath(output_filename,mustWork =F)
  #if (!grepl("\\.zip$",output_filename)) stop("output_filename must be a character ending in '.zip'")
  if (!grepl("\\.gpkg$",output_filename)) stop("'output_filename' must be a file path ending in '.gpkg'")
  if (file.exists(output_filename)) stop("'output_filename' already exists")

  out_file<-output_filename
  out_dir<-gsub(basename(out_file),"",out_file)
  if (!dir.exists(out_dir)) dir.create(out_dir)

  depression_corr<-match.arg(depression_corr)

  if (gsub(basename(output_filename),"",output_filename) == temp_dir) stop("'output_filename' should not be in the same directory as 'temp_dir'")

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

  dem<-hydroweight::process_input(dem,input_name="dem",working_dir=temp_dir)
  if (!inherits(dem,"SpatRaster")) stop("dem must be a class 'SpatRaster'")
  target_crs<-terra::crs(dem)

  terra::writeRaster(dem,file.path(temp_dir,"dem_final.tif"),overwrite=T,gdal=gdal_arg)

  if (!is.null(burn_streams)){

    resol<-terra::res(dem)[1]

    burn_streams<-hydroweight::process_input(burn_streams,
                                             align_to = terra::as.lines(terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",
                                                                         crs = target_crs)),
                                             clip_region = terra::as.polygons(terra::ext(dem),crs = target_crs), #+c(-resol*1.5,-resol*1.5,-resol*1.5,-resol*1.5)
                                             input_name="burn_streams",
                                             working_dir=temp_dir)

    terra::writeVector(burn_streams,file.path(temp_dir,"stream_final.shp"),overwrite=T)
    # Fillburn is a bit too aggressive in its burning (burns too deep)
    # wbt_fill_burn(
    #   dem="dem_final.tif",
    #   streams="stream_final.shp",
    #   output="dem_final.tif"
    # )

    r1<-dem

    sr1<-terra::mask(r1,burn_streams %>% sf::st_as_sf() %>% sf::st_buffer(resol*3) %>% terra::vect())
    r1[!is.na(sr1)]<-sr1-ceiling(burn_depth/3)
    sr1<-terra::mask(r1,burn_streams %>% sf::st_as_sf() %>% sf::st_buffer(resol*2) %>% terra::vect())
    r1[!is.na(sr1)]<-sr1-ceiling(burn_depth/3)
    sr1<-terra::mask(r1,burn_streams %>% sf::st_as_sf() %>% sf::st_buffer(resol*1) %>% terra::vect())
    r1[!is.na(sr1)]<-sr1-ceiling(burn_depth/3)
    terra::writeRaster(r1,file.path(temp_dir,"dem_final.tif"),overwrite=T,gdal=gdal_arg)

    whitebox::wbt_feature_preserving_smoothing(
      dem="dem_final.tif",
      output="dem_final.tif"
    )
  }

  if (!is.null(depression_corr)){
    if (depression_corr=="fill") {
      whitebox::wbt_fill_depressions(
        dem="dem_final.tif",
        output="dem_final.tif"
      )
    }

    if (depression_corr=="breach") {
      whitebox::wbt_breach_depressions(
        dem="dem_final.tif",
        output="dem_final.tif",
        fill_pits=T
      )
    }
  }


  # Process flow dirrection -------------------------------------------------
  ## Generate d8 pointer
  if (verbose) message("Generating d8 pointer")
  whitebox::wbt_d8_pointer(
    dem="dem_final.tif",
    output="dem_d8.tif"
  )

  ## Generate d8 flow accumulation in units of cells
  if (verbose) message("Generating d8 flow accumulation")
  whitebox::wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
    output = "dem_accum_d8.tif",
    out_type = "cells",
    pntr=T
  )
  whitebox::wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
    output = "dem_accum_d8_sca.tif",
    out_type = "sca",
    pntr=T
  )


  ## Generate streams with a stream initiation threshold
  if (verbose) message("Extracting Streams")
  whitebox::wbt_extract_streams(
    flow_accum = "dem_accum_d8.tif",
    output = "dem_streams_d8.tif",
    threshold = threshold
  )

  if (!is.null(min_length)){
    if (verbose) message("Trimming Streams")
    whitebox::wbt_remove_short_streams(
      d8_pntr="dem_d8.tif",
      streams="dem_streams_d8.tif",
      output="dem_streams_d8.tif",
      min_length=min_length
    )
  }

  # Prepare Output ----------------------------------------------------------
  #browser()

  dist_list_out<-list(
    "dem_final.tif",
    "dem_d8.tif",
    "dem_accum_d8.tif",
    "dem_accum_d8_sca.tif",
    "dem_streams_d8.tif"
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-output_filename
  out_dir<-gsub(basename(out_file),"",out_file)
  if (!dir.exists(out_dir)) dir.create(out_dir)

  if (file.exists(out_file)) {
    out_file<-file.path(out_dir,paste0(gsub("\\.gpkg","",basename(out_file)),"_",paste0(basename(tempfile()), ".gpkg")))
    warning(paste0("Target .gpkg file already exists. Saving as: ",out_file))

    # out_file<-file.path(out_dir,paste0(gsub("\\.zip","",basename(out_file)),"_",paste0(basename(tempfile()), ".zip")))
    # warning(paste0("Target .zip file already exists. Saving as: ",out_file))
  }

  if (verbose) message("Generating Output")

  #browser()

  dem_ext<-terra::rast(dist_list_out[[1]])

  dem_ext<-dem_ext %>%
    terra::ext() %>%
    terra::as.polygons(crs=terra::crs(dem_ext)) %>%
    sf::st_as_sf() %>%
    sf::write_sf(out_file,
                 layer="DEM_Extent",
                 append = T,
                 delete_layer = F,
                 delete_dsn = F)

  for (i in dist_list_out){
    tmp<-terra::rast(unlist(i))

    if (verbose) message(paste0("Writing: ",names(tmp)))

    ot <- writeRaster_fun(tmp,out_file)
  }

  # utils::zip(out_file,
  #     unlist(dist_list_out),
  #     flags = '-r9Xjq'
  # )

  output<-list(
    outfile=out_file
  )

  if (return_products){
    output<-c(
      lapply(stats::setNames(dist_list_out,basename(unlist(dist_list_out))),function(x) terra::wrap(terra::rast(x))),
      output
    )
  }
  suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive=T)))

  class(output)<-"ihydro"
  return(output)
}
