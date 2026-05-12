testthat::skip_on_cran()

testthat::test_that("ex_data provides example data and as.ihydro works", {
  # ex_data should return paths/objects
  dem <- ihydro::ex_data("elev_ned_30m.tif")
  pts <- ihydro::ex_data("sites_nc.shp")
  expect_true(inherits(dem, "SpatRaster") || inherits(dem, "character"))
  expect_true(inherits(pts, "SpatVector") || inherits(pts, "character"))

  # create a minimal gpkg using process_loi or process_flowdir later; instead
  # verify as.ihydro can wrap an existing gpkg path (create temporary gpkg)
  tmpgpkg <- tempfile(fileext = ".gpkg")
  dem |>
    terra::ext() |>
    terra::as.polygons(crs = terra::crs(dem)) |>
    sf::st_as_sf() |>
    sf::write_sf(tmpgpkg, layer = "DEM_Extent")

  ih <- ihydro::as.ihydro(tmpgpkg)
  expect_s3_class(ih, "ihydro")
  layers <- ihydro::ihydro_layers(ih)
  expect_true("DEM_Extent" %in% layers$layer_name)
})
