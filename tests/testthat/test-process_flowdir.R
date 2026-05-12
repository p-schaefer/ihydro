testthat::test_that("process_flowdir runs and produces expected layers", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("whitebox")
  testthat::skip_if_not_installed("terra")

  dem <- ihydro::ex_data("elev_ned_30m.tif")
  streams <- ihydro::ex_data("streams.shp")
  out_gpkg <- file.path(tempdir(), "test_hydro.gpkg")

  res <- ihydro::process_flowdir(
    dem = dem,
    burn_streams = streams,
    burn_depth = 5,
    depression_corr = "breach",
    threshold = 100L,
    return_products = FALSE,
    output_filename = out_gpkg,
    temp_dir = tempfile(),
    verbose = FALSE
  )

  expect_s3_class(res, "ihydro")
  layers <- ihydro::ihydro_layers(res)
  expect_true(any(grepl("dem_d8|dem_streams_d8|dem_accum_d8", layers$layer_name)))
})
