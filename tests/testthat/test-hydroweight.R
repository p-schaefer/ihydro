test_that("hydroweight returns expected structure for lumped", {
  dir <- tempdir()
  dem <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:4326")
  flow_accum <- terra::rast(
    nrows = 5,
    ncols = 5,
    vals = 1:25,
    crs = "EPSG:4326"
  )
  target_O <- terra::vect(
    matrix(c(0, 0), ncol = 2),
    type = "points",
    crs = "EPSG:4326"
  )
  uid <- "testsite"
  result <- hydroweight(
    hydroweight_dir = dir,
    dem = dem,
    flow_accum = flow_accum,
    target_O = target_O,
    target_uid = uid,
    weighting_scheme = "lumped",
    save_output = FALSE,
    return_products = TRUE,
    clean_tempfiles = TRUE
  )
  expect_type(result, "list")
  expect_true("lumped" %in% names(result))
  expect_s4_class(result$lumped, "SpatRaster")
})

test_that("hydroweight returns all requested schemes", {
  dir <- tempdir()
  dem <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:4326")
  flow_accum <- terra::rast(
    nrows = 5,
    ncols = 5,
    vals = 1:25,
    crs = "EPSG:4326"
  )
  target_O <- terra::vect(
    matrix(c(0, 0), ncol = 2),
    type = "points",
    crs = "EPSG:4326"
  )
  target_S <- terra::vect(
    matrix(c(1, 1), ncol = 2),
    type = "points",
    crs = "EPSG:4326"
  )
  uid <- "testsite"
  schemes <- c("lumped", "iEucO", "iEucS")
  result <- hydroweight(
    hydroweight_dir = dir,
    dem = dem,
    flow_accum = flow_accum,
    target_O = target_O,
    target_S = target_S,
    target_uid = uid,
    weighting_scheme = schemes,
    save_output = FALSE,
    return_products = TRUE,
    clean_tempfiles = TRUE
  )
  expect_setequal(names(result), schemes)
  lapply(result, function(x) expect_s4_class(x, "SpatRaster"))
})

test_that("hydroweight returns file path if return_products = FALSE", {
  dir <- tempdir()
  dem <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:4326")
  flow_accum <- terra::rast(
    nrows = 5,
    ncols = 5,
    vals = 1:25,
    crs = "EPSG:4326"
  )
  target_O <- terra::vect(
    matrix(c(0, 0), ncol = 2),
    type = "points",
    crs = "EPSG:4326"
  )
  uid <- "testsite"
  result <- hydroweight(
    hydroweight_dir = dir,
    dem = dem,
    flow_accum = flow_accum,
    target_O = target_O,
    target_uid = uid,
    weighting_scheme = "lumped",
    save_output = TRUE,
    return_products = FALSE,
    clean_tempfiles = TRUE
  )
  expect_type(result, "character")
  expect_true(grepl("inv_distances.zip", result))
})

test_that("hydroweight errors if target_uid is missing", {
  dir <- tempdir()
  dem <- terra::rast(nrows = 5, ncols = 5, vals = 1:25, crs = "EPSG:4326")
  flow_accum <- terra::rast(
    nrows = 5,
    ncols = 5,
    vals = 1:25,
    crs = "EPSG:4326"
  )
  target_O <- terra::vect(
    matrix(c(0, 0), ncol = 2),
    type = "points",
    crs = "EPSG:4326"
  )
  expect_error(
    hydroweight(
      hydroweight_dir = dir,
      dem = dem,
      flow_accum = flow_accum,
      target_O = target_O,
      weighting_scheme = "lumped",
      save_output = FALSE,
      return_products = TRUE,
      clean_tempfiles = TRUE
    ),
    "target_uid"
  )
})
