## code to prepare `dem`, `streams`, `landuse`, `geology`, and `pointsources` datasets goes here

## Generate save_dir as a temporary directory
save_dir <- tempdir()
final_dir <- file.path("inst", "extdata")
if (!dir.exists(final_dir)) {
  dir.create(final_dir, recursive = TRUE)
}

## Import toy_dem from openSTARS package
download.file(
  "https://github.com/MiKatt/openSTARS/archive/refs/heads/master.zip",
  file.path(save_dir, "openSTARS.zip")
)
unzip(file.path(save_dir, "openSTARS.zip"), exdir = save_dir)

fl <- list.files(
  file.path(save_dir, "openSTARS-master", "inst", "extdata", "nc"),
  full.names = TRUE
)

fl <- fl[grepl("\\.tif$|\\.shp$", fl)]

targ_crs <- terra::crs(terra::rast(fl[grepl("elev",fl)]))

for (i in fl) {
  if (grepl("\\.tif$", i)) {
    r <- terra::rast(i)
    r <- terra::project(r, targ_crs)
    terra::writeRaster(r,
                       i,
                       overwrite = TRUE,
                       datatype = "FLT4S")
  }
  if (grepl("\\.shp$", i)) {
    v <- terra::vect(i)
    v <- terra::project(v, targ_crs)
    terra::writeVector(v, i, overwrite = TRUE)
  }
}

zip(
  zipfile = file.path(
    file.path(save_dir, "openSTARS-master", "inst", "extdata", "nc"),
    "openSTARS_data.zip"
  ),
  files = list.files(
    file.path(save_dir, "openSTARS-master", "inst", "extdata", "nc"),
    full.names = TRUE
  ),
  flags = '-r9Xj'
)

file.copy(
  file.path(
    file.path(save_dir, "openSTARS-master", "inst", "extdata", "nc"),
    "openSTARS_data.zip"
  ),
  file.path(final_dir, "openSTARS_data.zip"),
  overwrite = TRUE
)
