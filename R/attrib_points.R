
#' Attributes stream segments/sampling points with layers of interest (loi)
#'
#' @param input an object of class ihydro
#' @param out_filename filepath of csv to save resulting attributes
#' @param loi_file filepath of `process_loi()` output
#' @param loi_cols character or NULL. Names of loi layers to include in summary. If NULL, all layers used.
#' @param sample_points character or NULL. IDs of unique station identifiers provided in 'site_id_col' of `generate_vectors()` to calculate attributes for.
#' @param link_id character or NULL. 'link_id's of reaches to calculate attributes for.
#' @param target_o_type character. One of: c("point","segment_point","segment_whole"). Target for "iEucO", "iFLO", and "HAiFLO" weighting schemes. "Point" represents the sampling point on the stream, "segment_point" represents the upstream segment of the sampling points, and "segment_whole" will target the entire reach, regardless of where sampling occurred.
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param loi_numeric_stats character. One or more of c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"). Those without distwtd_ are simple "lumped" statistics.
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param clip_region character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Regions to remove when summarizing the attributes (e.g., remove lake from catchment)
#' @param OS_combine logical. Should target_O and target_S be merged as targets for iEucS, iFLS, and/or HAiFLS? Use \code{TRUE} or \code{FALSE}. This allows cells surrounding \code{target_O} to flow directly into \code{target_O} rather than be forced through \code{target_S}.
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return A tibble of resulting attributes. If \code{return_products = TRUE}, all geospatial analysis products are returned as well.
#'
#' @export


attrib_points<-function(
    input,
    out_filename,
    loi_file=NULL,
    loi_cols=NULL,
    sample_points=NULL,
    link_id=NULL,
    clip_region=NULL,
    OS_combine=F,
    target_o_type=c("point","segment_point","segment_whole"),
    weighting_scheme = c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS"),
    loi_numeric_stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir=NULL,
    return_products=T,
    verbose=F
){
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  if (inherits(loi_file,"ihydro")) loi_file<-loi_file$outfile

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
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

  loi_numeric_stats<-stats::setNames(loi_numeric_stats,loi_numeric_stats)

  attr_db_loc<-db_loc<-dw_dir<-zip_loc<-input$outfile

  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

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

  # Get Target IDs ----------------------------------------------------------
  target_IDs<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=target_o_type=="segment_whole",
    target_o_type=target_o_type
  )

  target_O<-target_o_fun(
    db_fp=db_loc,
    target_IDs=target_IDs,
    target_o_type=target_o_type
  )

  # Setup remove_region -----------------------------------------------------
  clip_region<-hydroweight::process_input(clip_region,input_name="clip_region")
  if (!is.null(clip_region)){
    if (inherits(clip_region,"SpatRaster")) {
      terra::writeRaster(clip_region,file.path(temp_dir,"clip_region.tif"))
      clip_region<-file.path(temp_dir,"clip_region.tif")
    }
    if (inherits(clip_region,"SpatVector")) {
      terra::writeVector(clip_region,file.path(temp_dir,"clip_region.shp"))
      clip_region<-file.path(temp_dir,"clip_region.shp")
    }
  }

  temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
  dir.create(temp_dir_sub)

  for (i in c("dem_streams_d8",
              "dem_final",
              "dem_accum_d8",
              "dem_d8")) {
    terra::writeRaster(
      terra::rast(zip_loc,i),
      file.path(temp_dir_sub,paste0(i,".tif"))
    )
  }

  if (verbose) message("Calculating Attributes")
  progressr::with_progress(enable=verbose,{
    p <- progressr::progressor(steps = nrow(target_O))

    attr_out<-target_O %>%
      dplyr::mutate(core=rep(1:n_cores,length.out=dplyr::n())) %>%
      dplyr::group_by(core) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        attr=furrr::future_pmap(
          #attr = purrr::pmap(
          list(
            x=data,
            p=list(p),
            input=list(input),
            loi_file=list(loi_file),
            weighting_scheme=list(weighting_scheme),
            loi_cols=list(loi_cols),
            loi_numeric_stats=list(loi_numeric_stats),
            temp_dir_sub=list(temp_dir_sub),
            return_products=list(return_products),
            inv_function=list(inv_function),
            clip_region=list(clip_region),
            OS_combine=list(OS_combine),
            loi_meta=list(loi_meta)
          ),
          .options = furrr_options(globals = F),
          #carrier::crate(
          function(x,
                   p,
                   input,
                   loi_file,
                   weighting_scheme,
                   loi_cols,
                   loi_numeric_stats,
                   temp_dir_sub,
                   return_products,
                   inv_function,
                   clip_region,
                   OS_combine,
                   loi_meta
          ){
            suppressPackageStartupMessages(library(sf))
            options(dplyr.summarise.inform = FALSE)
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            target_S <- terra::rast(file.path(temp_dir_sub,"dem_streams_d8.tif"))
            dem <- terra::rast(file.path(temp_dir_sub,"dem_final.tif"))
            flow_accum <- terra::rast(file.path(temp_dir_sub,"dem_accum_d8.tif"))

            purrr::pmap_dfr(
              list(
                y=split(x$geom,x$link_id),
                yy=split(x$link_id,x$link_id)
              ),
              function(y,yy){

                y<-sf::st_as_sf(y) %>%
                  dplyr::mutate(link_id=yy)

                temp_dir_sub_sub<-file.path(temp_dir_sub,basename(tempfile()))
                dir.create(temp_dir_sub_sub)

                catch<-ihydro::get_catchment(input,link_id =yy)

                suppressMessages(
                  hw<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub_sub,
                                               target_S = target_S,
                                               target_O = y,
                                               target_uid = yy,
                                               OS_combine = OS_combine,
                                               dem=dem,
                                               flow_accum = flow_accum,
                                               clip_region=catch,
                                               weighting_scheme = weighting_scheme,
                                               inv_function = inv_function,
                                               return_products=return_products,
                                               wrap_return_products = T,
                                               save_output = T,
                                               clean_tempfiles = T)
                )

                hw_attr_num<-NULL
                hw_attr_cat<-NULL
                if (any(loi_meta$loi_type=="num_rast")) {
                  suppressMessages(
                    hw_attr_num<-hydroweight::hydroweight_attributes(loi=terra::rast(loi_file$outfile,lyrs=loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     loi_columns = loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"],
                                                                     loi_numeric=T,
                                                                     loi_numeric_stats = loi_numeric_stats,
                                                                     roi = catch,
                                                                     roi_uid=yy,
                                                                     roi_uid_col = "pour_point_ID",
                                                                     distance_weights=file.path(temp_dir_sub_sub,paste0(yy,"_inv_distances.zip")),
                                                                     remove_region = clip_region,
                                                                     return_products = return_products)
                  )
                }

                if (any(loi_meta$loi_type=="cat_rast")) {
                  suppressMessages(
                    hw_attr_cat<-hydroweight::hydroweight_attributes(loi=terra::rast(loi_file$outfile,lyrs=loi_meta$loi_var_nms[loi_meta$loi_type=="cat_rast"]),
                                                                     loi_columns = loi_meta$loi_var_nms[loi_meta$loi_type=="cat_rast"],
                                                                     loi_numeric=F,
                                                                     loi_numeric_stats = loi_numeric_stats,
                                                                     roi = catch,
                                                                     roi_uid=yy,
                                                                     roi_uid_col = "pour_point_ID",
                                                                     distance_weights=file.path(temp_dir_sub_sub,paste0(yy,"_inv_distances.zip")),
                                                                     remove_region = clip_region,
                                                                     return_products = return_products)
                  )
                }

                #p()

                p_out <- NULL
                if (return_products){
                  p0<-hw
                  p1<-unlist(hw_attr_num$return_products,recursive=F)
                  if (!is.null(p1)) names(p1)<-paste0(names(hw_attr_num$return_products),"_num")
                  p2<-unlist(hw_attr_cat$return_products,recursive=F)
                  if (!is.null(p2)) names(p2)<-paste0(names(hw_attr_cat$return_products),"_cat")
                  p_out<-c(p0,p1,p2)
                  p_out<-p_out[sort(names(p_out))]
                }

                return(
                  dplyr::bind_cols(
                    tibble::tibble(products=list(p_out)[1]),
                    purrr::reduce(
                      list(hw_attr_num$attribute_table,
                           hw_attr_cat$attribute_table),
                      dplyr::left_join,
                      by="pour_point_ID")
                  ) %>%
                    dplyr::select(pour_point_ID,tidyselect::everything())
                )

              })
          }
          #)
        ))
  })

  final_out<-attr_out %>%
    dplyr::select(attr) %>%
    tidyr::unnest(cols=c("attr"))

  target_IDs_out<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=F,
    target_o_type=target_o_type
  )

  final_out<-target_IDs_out %>%
    dplyr::left_join(
      final_out %>%
        dplyr::mutate(pour_point_ID=as.character(pour_point_ID)),
      by=c("link_id"="pour_point_ID")
    )

  data.table::fwrite(x=final_out %>% dplyr::select(-tidyselect::any_of("products")),
                     file=out_filename,
                     buffMB = 128L,
                     nThread = 1,
                     showProgress = F)

  unlink(temp_dir_sub,recursive = T,force=T)

  return(final_out)

}


