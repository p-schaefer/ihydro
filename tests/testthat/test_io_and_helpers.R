library(testthat)
skip_on_cran()

test_that("ex_data returns example files and process_input handles inputs", {
  dem <- ihydro::ex_data("elev_ned_30m.tif")
  expect_true(!is.null(dem))

  shp <- ihydro::ex_data("sites_nc.shp")
  expect_true(!is.null(shp))

  # process_input should accept a path, an sf object, and return an sf when appropriate
  sf_in <- ihydro:::process_input(shp)
  expect_true(inherits(sf_in, "SpatVector"))

  # process_input on dem returns a terra SpatRaster / RasterLayer or path
  dem_in <- ihydro:::process_input(dem)
  expect_true(inherits(dem_in, "SpatRaster") || inherits(dem_in, "RasterLayer") || is.character(dem_in))
})


test_that("as_ihydro/as.ihydro and ihydro_layers/read_ihydro basic behaviour", {
  skip_if_not_installed("whitebox")

  dem <- ihydro::ex_data("elev_ned_30m.tif")
  streams <- ihydro::ex_data("streams.shp")
  out_gpkg <- file.path(tempdir(), "test_io_hydro.gpkg")
  if (file.exists(out_gpkg)) file.remove(out_gpkg)

  hydro_out <- ihydro::process_flowdir(
    dem = dem,
    threshold = 100L,
    burn_streams = streams,
    burn_depth = 5,
    depression_corr = "breach",
    return_products = TRUE,
    output_filename = out_gpkg,
    verbose = FALSE
  )

  expect_s3_class(hydro_out, "ihydro")
  expect_true(file.exists(hydro_out$outfile))

  ih <- ihydro::as.ihydro(hydro_out$outfile)
  expect_s3_class(ih, "ihydro")

  layers <- ihydro::ihydro_layers(ih)
  expect_true(is.data.frame(layers))
  expect_true(nrow(layers) > 0)

  # read a raster layer that should exist
  if ("dem_accum_d8" %in% layers$layer_name) {
    r <- ihydro::read_ihydro(ih, "dem_accum_d8")
    expect_true(inherits(r, "SpatRaster") || inherits(r, "RasterLayer"))
  }
})
