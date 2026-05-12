library(testthat)

skip_on_cran()

# Tests for process_loi, prep_weights, fasttrib_points / attrib_points

test_that("process_loi and fasttrib_points produce attributes for sample points", {
  skip_if_not_installed("whitebox")

  dem <- ihydro::ex_data("elev_ned_30m.tif")
  streams <- ihydro::ex_data("streams.shp")
  points <- ihydro::ex_data("sites_nc.shp")
  landuse <- ihydro::ex_data("landuse_r.tif")
  geology <- ihydro::ex_data("geology.shp")
  pointsrc <- ihydro::ex_data("pointsources.shp")

  ex_loc <- tempdir()
  terra::writeRaster(landuse, file.path(ex_loc, "LC.tif"), overwrite = TRUE)
  terra::writeVector(geology, file.path(ex_loc, "geology.shp"), overwrite = TRUE)
  terra::writeVector(pointsrc, file.path(ex_loc, "pointsources.shp"), overwrite = TRUE)

  out_gpkg <- file.path(tempdir(), "test_loi_hydro.gpkg")
  if (file.exists(out_gpkg)) file.remove(out_gpkg)

  hydro_out <- ihydro::process_hydrology(
    dem = dem,
    threshold = 200L,
    burn_streams = streams,
    burn_depth = 5,
    depression_corr = "breach",
    return_products = TRUE,
    output_filename = out_gpkg,
    verbose = FALSE
  )

  output_filename_loi <- file.path(ex_loc, "Processed_loi.gpkg")
  if (file.exists(output_filename_loi)) file.remove(output_filename_loi)

  loi_combined <- ihydro::process_loi(
    dem = dem,
    num_inputs = list(slope = file.path(ex_loc, "LC.tif")),
    cat_inputs = list(landcover = file.path(ex_loc, "LC.tif"),
                      geology = file.path(ex_loc, "geology.shp"),
                      pointsources = file.path(ex_loc, "pointsources.shp")),
    variable_names = list(geology = "GEO_NAME", pointsources = "psource"),
    output_filename = output_filename_loi,
    return_products = TRUE,
    verbose = FALSE
  )

  expect_s3_class(loi_combined, "ihydro")
  layers_loi <- ihydro::ihydro_layers(loi_combined)
  expect_true("loi_meta" %in% layers_loi$layer_name)

  # Now compute weights and attributes for a small set of sample points
  fast_out <- ihydro::fasttrib_points(
    input = hydro_out,
    out_filename = file.path(tempdir(), "attr.csv"),
    loi_file = loi_combined,
    link_id = c("1", "25", "80"),
    weighting_scheme = c("lumped", "iFLS"),
    loi_numeric_stats = c("mean", "sd"),
    inv_function = function(x) (x * 0.001 + 1)^-1,
    verbose = FALSE
  )

  expect_true(is.data.frame(fast_out))
  expect_true(all(c("link_id") %in% colnames(fast_out)))

})
