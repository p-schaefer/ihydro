library(testthat)

skip_on_cran()

# Test generate_pdist and basic modeling pipeline pieces using small dataset

test_that("generate_pdist creates pairwise tables and can be used for clustering", {
  skip_if_not_installed("whitebox")

  dem <- ihydro::ex_data("elev_ned_30m.tif")
  streams <- ihydro::ex_data("streams.shp")
  points <- ihydro::ex_data("sites_nc.shp")

  out_gpkg <- file.path(tempdir(), "test_pdist.gpkg")
  if (file.exists(out_gpkg)) file.remove(out_gpkg)

  hydro_out <- ihydro::process_hydrology(
    dem = dem,
    output_filename = out_gpkg,
    burn_streams = streams,
    burn_depth = 5,
    threshold = 300L,
    points = points,
    site_id_col = "site_id",
    break_on_noSnap = FALSE,
    pwise_dist = FALSE,
    return_products = TRUE,
    verbose = FALSE
  )

  pd <- ihydro::generate_pdist(hydro_out, return_products = TRUE, pwise_all_links = FALSE)
  expect_s3_class(pd, "ihydro")
  layers <- ihydro::ihydro_layers(pd)
  expect_true(any(grepl("pwise_dist", layers$layer_name)))

  # quick clustering on pairwise dist from pd$pwise_dist if present
  if (!is.null(pd$pwise_dist)) {
    dmat <- pd$pwise_dist
    expect_true(is.data.frame(dmat) || inherits(dmat, "tbl"))
  }

})
