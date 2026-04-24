# tests/testthat/test-process_input.R

test_that("process_input handles SpatRaster input", {
  r <- terra::rast(nrows = 10, ncols = 10, vals = 1:100)
  terra::crs(r) <- "EPSG:4326"
  out <- process_input(r)
  expect_s4_class(out, "SpatRaster")
  expect_equal(terra::crs(out), terra::crs(r))
})

test_that("process_input handles SpatVector input", {
  v <- terra::vect(
    matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE),
    type = "points",
    crs = "EPSG:4326"
  )
  out <- process_input(v)
  expect_s4_class(out, "SpatVector")
  expect_equal(terra::crs(out), terra::crs(v))
})

test_that("process_input subsets variables", {
  r <- terra::rast(nrows = 5, ncols = 5, nlyrs = 2, vals = 1:50)
  names(r) <- c("a", "b")
  terra::crs(r) <- "EPSG:4326"
  out <- process_input(r, input_variable_names = "a")
  expect_equal(names(out), "a")
})

test_that("process_input aligns raster to target", {
  r1 <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:4326")
  r2 <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:3857")
  out <- process_input(r1, align_to = r2)
  expect_equal(terra::crs(out), terra::crs(r2))
})

test_that("process_input errors on missing CRS", {
  r <- terra::rast(nrows = 5, ncols = 5, vals = 1:25)
  terra::crs(r) <- NULL
  expect_error(process_input(r), "crs")
})

test_that("process_input errors on missing variable", {
  r <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:4326")
  expect_error(process_input(r, variable_names = "not_a_var"), "not in input")
})

test_that("process_input rasterizes SpatVector with resample_type = 'bilinear'", {
  # Create a SpatVector (points with numeric attribute)
  v <- terra::vect(
    data.frame(x = c(0, 1), y = c(0, 1), val = c(10, 20)),
    geom = c("x", "y"),
    crs = "EPSG:4326"
  )
  r <- terra::rast(
    nrows = 10,
    ncols = 10,
    xmin = -1,
    xmax = 2,
    ymin = -1,
    ymax = 2,
    crs = "EPSG:4326"
  )
  out <- process_input(
    v,
    variable_names = "val",
    target = r,
    resample_type = "bilinear"
  )
  expect_s4_class(out, "SpatRaster")
  expect_equal(names(out), "val")
  expect_true(all(!is.na(terra::values(out))))
})

test_that("process_input rasterizes SpatVector with resample_type = 'near'", {
  # Create a SpatVector (points with categorical attribute)
  v <- terra::vect(
    data.frame(x = c(0, 1), y = c(0, 1), cat = c("a", "b")),
    geom = c("x", "y"),
    crs = "EPSG:4326"
  )
  r <- terra::rast(
    nrows = 10,
    ncols = 10,
    xmin = -1,
    xmax = 2,
    ymin = -1,
    ymax = 2,
    crs = "EPSG:4326"
  )
  out <- process_input(
    v,
    variable_names = "cat",
    target = r,
    resample_type = "near"
  )
  expect_s4_class(out, "SpatRaster")
  # Should have one layer per unique category
  expect_true(all(c("cat_a", "cat_b") %in% names(out)))
})

test_that("process_input projects raster with resample_type = 'bilinear'", {
  r1 <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:4326")
  r2 <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:3857")
  out <- process_input(r1, target = r2, resample_type = "bilinear")
  expect_equal(terra::crs(out), terra::crs(r2))
})

test_that("process_input projects raster with resample_type = 'near'", {
  r1 <- terra::rast(
    nrows = 5,
    ncols = 5,
    vals = rep(1:5, each = 5),
    crs = "EPSG:4326"
  )
  r2 <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:3857")
  out <- process_input(r1, target = r2, resample_type = "near")
  expect_equal(terra::crs(out), terra::crs(r2))
})
