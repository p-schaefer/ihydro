
library(whitebox)
## Load libraries
# remotes::install_github("p-schaefer/ihydro)
library(ihydro)
library(tmap)
library(furrr)
library(whitebox)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)

# Many function in 'ihydro' can be run in parallel internally.

# Parallelization is done through the future package, so all parallel backends
# should be supported i.e.,:
plan(multisession(workers=4))

## Generate save_dir as a temporary directory
save_dir <- tempdir()


## Import toy_dem from openSTARS package
# devtools::install_github("MiKatt/openSTARS", ref = "dev")

ot<-system.file("extdata", "nc", "elev_ned_30m.tif", package = "openSTARS") %>%
  rast()

crs(ot)<-crs(rast(system.file("extdata", "nc", "landuse_r.tif", package = "openSTARS")))
writeRaster(ot,file.path(save_dir, "toy_dem.tif"),overwrite=T)

toy_dem<-rast(file.path(save_dir, "toy_dem.tif"))

## Identify some sampling points
system.file("extdata", "nc", "sites_nc.shp", package = "openSTARS") %>%
  vect() %>%
  st_as_sf() %>%
  st_transform(st_crs(toy_dem)) %>%
  vect() %>%
  writeVector(file.path(save_dir, "sites.shp"),overwrite=T)


output_filename_hydro<-file.path(save_dir,"Processed_Hydrology.gpkg")

# Generates d8 flow direction and accumulation, extracts streams at a specified
# flow accumulation threshold
hydro_out<-process_flowdir(
  dem=toy_dem,
  burn_streams=system.file("extdata", "nc", "streams.shp", package = "openSTARS"),
  burn_depth=5,
  min_length=3,
  depression_corr="breach",
  threshold=100L,
  return_products=F,
  output_filename=output_filename_hydro,
  temp_dir=NULL,
  verbose=T
)

hydro_out<-generate_vectors(
  input=hydro_out,
  points=file.path(save_dir, "sites.shp"), # These are optional sampling locations
  site_id_col="site_id", # Column name in points layer that corresponds to
  #                      # unique IDs that will be available in data products
  snap_distance=100L, # points that are more than 100m from closest stream are excluded
  break_on_noSnap =F, # default is to stop when any points don't snap, this will ignore that
  return_products=F,
  temp_dir=NULL,
  verbose=T
)

hydro_out<-trace_flowpaths(
  input=hydro_out,
  return_products=F,
  pwise_dist=F, # This will calculate all downstream flow-connected pairwise distances
  pwise_all_links=F, # This will calculate all flow un-connected pairwise distances (can be very time consuming)
  temp_dir=NULL,
  verbose=T
)

## Predictors from openSTARS
system.file("extdata", "nc", "landuse_r.tif", package = "openSTARS") %>%
  rast() %>%
  setNames("LC") %>%
  writeRaster(file.path(save_dir, "LC.tif"),overwrite=T)

landuse_r_path <-file.path(save_dir, "LC.tif")
geology_path<-system.file("extdata", "nc", "geology.shp", package = "openSTARS")
pointsources_path<-system.file("extdata", "nc", "pointsources.shp", package = "openSTARS")

read_sf(pointsources_path) %>%
  mutate(pointsource="pontsrc") %>%
  st_buffer(60) %>%
  write_sf(file.path(save_dir, "pointsources.shp"),overwrite=T)

pointsources_path<-file.path(save_dir, "pointsources.shp")


# Numeric Raster

wbt_slope(
  dem = file.path(save_dir, "toy_dem.tif"),
  output = file.path(save_dir, "slope.tif")
)

# Combine loi layers
output_filename_loi<-file.path(save_dir,"Processed_loi.gpkg")

# This function standardizes numeric and categorical loi layers.

loi_combined<-process_loi(
  dem=toy_dem,
  num_inputs=list(# Can be given as a mixture of input types (file paths, or any sf or terra format)
    slope=file.path(save_dir, "slope.tif")
  ),
  cat_inputs=list(# Can be given as a mixture of input types (file paths, or any sf or terra format)
    landcover=landuse_r_path,
    geology=geology_path,
    pointsources=pointsources_path
  ),
  variable_names=list( # any unlisted inputs will be used in their entirety
    geology="GEO_NAME", # names listed here will subset those attributes or layers from the inputs
    pointsources="pontsrc"
  ),
  output_filename=output_filename_loi,
  return_products=F,
  temp_dir=NULL,
  verbose=T
)

pw<-prep_weights(
    input=hydro_out,
    output_filename=file.path(tempdir(),"prepped_weights.gpkg"),
    sample_points=c("1","25","4"),
    target_o_type=c("point"),
    weighting_scheme =  c("iFLS", "HAiFLS","iFLO",  "HAiFLO"),
    verbose=T
)

final_attributes_sub<-fasttrib_points(
  input=hydro_out,
  out_filename="sample_points_wgtattr.csv",
  loi_file=output_filename_loi, # specify loi file, if NULL, function will look for loi in 'input'
  loi_cols=NULL,                # Specify loi columns to use, if NULL, all will be used
  iDW_file=file.path(tempdir(),"prepped_weights.gpkg"), # Leaving this as NULL will look for iDW in 'input' and calculate any not available
  store_iDW=T,   # This will save the distance weights to the input or iDW_file if it is specified
  sample_points=c("1","25","4"),
  link_id=NULL,
  target_o_type=c("point"),
  weighting_scheme =  c("lumped",  "iFLS", "HAiFLS","iFLO",  "HAiFLO"),#
  loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
  inv_function = function(x) {
    (x * 0.001 + 1)^-1
  },
  temp_dir=NULL,
  verbose=T,
  backend="SQLite"
)

final_attributes_sub<-fasttrib_points(
  input=hydro_out,
  out_filename="sample_points_wgtattr.csv",
  loi_file=output_filename_loi, # specify loi file, if NULL, function will look for loi in 'input'
  loi_cols=NULL,                # Specify loi columns to use, if NULL, all will be used
  iDW_file=file.path(tempdir(),"prepped_weights.gpkg"), # Leaving this as NULL will look for iDW in 'input' and calculate any not available
  store_iDW=T,   # This will save the distance weights to the input or iDW_file if it is specified
  sample_points=c("1","25","80"),
  link_id=NULL,
  target_o_type=c("point"),
  weighting_scheme =  c("lumped",  "iFLS", "HAiFLS","iFLO",  "HAiFLO"),#
  loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
  inv_function = function(x) {
    (x * 0.001 + 1)^-1
  },
  temp_dir=NULL,
  verbose=T,
  backend="SQLite"
)

final_attributes_sub<-fasttrib_points(
  input=hydro_out,
  out_filename="sample_points_wgtattr.csv",
  loi_file=output_filename_loi, # specify loi file, if NULL, function will look for loi in 'input'
  loi_cols=NULL,                # Specify loi columns to use, if NULL, all will be used
  iDW_file=NULL, # Leaving this as NULL will look for iDW in 'input' and calculate any not available
  store_iDW=T,   # This will save the distance weights to the input or iDW_file if it is specified
  sample_points=c("1","25","4"),
  link_id=NULL,
  target_o_type=c("point"),
  weighting_scheme =  c("lumped",  "iFLS", "HAiFLS","iFLO",  "HAiFLO"),#
  loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
  inv_function = function(x) {
    (x * 0.001 + 1)^-1
  },
  temp_dir=NULL,
  verbose=T,
  backend="SQLite"
)

# mtcars2<-mtcars
#
#
# a1<-carrier::crate(function(x=mtcars2,col1=c("mpg","cyl")) {
#   `%>%` <- magrittr::`%>%`
#   x %>%
#     dtplyr::lazy_dt() %>%
#     dplyr::mutate(dplyr::across(tidyselect::any_of(col1), ~.*!!rlang::sym("gear"),.names=paste0("{.col}_","gear"))) %>%
#     dplyr::as_tibble()
# })
#
# a1(mtcars2,col1=c("mpg","cyl"))
