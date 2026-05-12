library(testthat)
skip_on_cran()

# Integration test: generate vectors and trace flowpaths works
test_that("generate_vectors and trace_flowpaths integration", {
  skip_if_not_installed("whitebox")

  dem <- ihydro::ex_data("elev_ned_30m.tif")
  streams <- ihydro::ex_data("streams.shp")
  points <- ihydro::ex_data("sites_nc.shp")

  out_gpkg <- file.path(tempdir(), "test_vectors.gpkg")
  if (file.exists(out_gpkg)) file.remove(out_gpkg)

  hydro_out <- ihydro::process_flowdir(
    dem = dem,
    threshold = 200L,
    burn_streams = streams,
    burn_depth = 5,
    depression_corr = "breach",
    return_products = TRUE,
    output_filename = out_gpkg,
    verbose = FALSE
  )

  # generate vectors
  gv <- ihydro::generate_vectors(
    input = hydro_out,
    points = points,
    site_id_col = "site_id",
    snap_distance = Inf,
    break_on_noSnap = FALSE,
    return_products = TRUE,
    temp_dir = NULL,
    verbose = FALSE
  )

  expect_s3_class(gv, "ihydro")
  layers <- ihydro::ihydro_layers(gv)
  expect_true("stream_lines" %in% layers$layer_name)
  expect_true("snapped_points" %in% layers$layer_name)

  # trace flowpaths
  tf <- ihydro::trace_flowpaths(input = gv, pwise_dist = TRUE, pwise_all_links = FALSE, return_products = TRUE)
  expect_s3_class(tf, "ihydro")
  layers2 <- ihydro::ihydro_layers(tf)
  expect_true("us_flowpaths" %in% layers2$layer_name)
  expect_true("ds_flowpaths" %in% layers2$layer_name)

  # check that us_flowpaths has expected columns
  us <- tf$us_flowpaths
  expect_true(all(c("pour_point_id", "origin_link_id") %in% names(us)))
})
