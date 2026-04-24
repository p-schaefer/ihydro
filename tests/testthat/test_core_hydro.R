library(testthat)

skip_on_cran()

test_that("process_flowdir produces expected ihydro layers and files", {
  skip_if_not_installed("whitebox")
  dem <- ihydro::ex_data("elev_ned_30m.tif")
  streams <- ihydro::ex_data("streams.shp")

  out_gpkg <- file.path(tempdir(), "test_proc_flow.gpkg")
  if (file.exists(out_gpkg)) file.remove(out_gpkg)

  res <- ihydro::process_flowdir(
    dem = dem,
    threshold = 100L,
    burn_streams = streams,
    burn_depth = 5,
    depression_corr = "breach",
    return_products = TRUE,
    output_filename = out_gpkg,
    verbose = FALSE
  )

  expect_s3_class(res, "ihydro")
  expect_true(file.exists(res$outfile))

  layers <- ihydro::ihydro_layers(res)
  expect_true("dem_d8" %in% layers$layer_name)
  expect_true("dem_streams_d8" %in% layers$layer_name)

  # reading a raster product works
  acc <- ihydro::read_ihydro(res, "dem_accum_d8")
  expect_true(inherits(acc, "SpatRaster") || inherits(acc, "RasterLayer"))
})
