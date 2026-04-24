library(testthat)
library(terra)
library(sf)
library(exactextractr)

# These tests focus on verifying that the weighted SD calculation in exactextractr
# matches manual calculations, including handling of NAs and pooling across chunks.
# Note: exactextractr currently uses 32-bit precision, so expect errors ~ 1e-7

test_that("weighted SD matches manual calculation (no NAs)", {
  vals <- c(1, 2, 3, 4, 5)
  wts <- c(1, 2, 3, 4, 5)
  # Manual weighted mean
  wm <- sum(vals * wts) / sum(wts)
  # Manual weighted population variance
  wvar_pop <- sum(wts * (vals - wm)^2) / sum(wts)
  # Manual weighted population SD
  wsd_pop <- sqrt(wvar_pop)
  # Manual weighted sample SD
  n_eff <- sum(wts)^2 / sum(wts^2)
  wsd_samp <- sqrt(wvar_pop * n_eff / (n_eff - 1))

  # Raster setup
  r <- rast(matrix(vals, nrow = 1))
  wr <- rast(matrix(wts, nrow = 1))
  v <- vect(
    matrix(c(0, 0, 0, 1, 5, 1, 5, 0), ncol = 2, byrow = TRUE),
    type = "polygons"
  )

  # exactextractr weighted SD (population)
  ee_sd_pop <- exact_extract(r, st_as_sf(v), 'weighted_stdev', weights = wr)
  expect_equal(ee_sd_pop, wsd_pop, tolerance = 1e-7)

  # exactextractr weighted SD (sample, manual conversion)
  ee_sd_samp <- ee_sd_pop * sqrt(n_eff / (n_eff - 1))
  expect_equal(ee_sd_samp, wsd_samp, tolerance = 1e-7)
})

test_that("weighted SD handles NAs correctly", {
  vals <- c(1, 2, NA, 4, 5)
  wts <- c(1, 2, 3, 4, 5)
  # Remove NAs for manual calculation
  idx <- !is.na(vals)
  vals2 <- vals[idx]
  wts2 <- wts[idx]
  wm <- sum(vals2 * wts2) / sum(wts2)
  wvar_pop <- sum(wts2 * (vals2 - wm)^2) / sum(wts2)
  wsd_pop <- sqrt(wvar_pop)
  n_eff <- sum(wts2)^2 / sum(wts2^2)
  wsd_samp <- sqrt(wvar_pop * n_eff / (n_eff - 1))

  # Raster setup
  r <- rast(matrix(vals, nrow = 1))
  wr <- rast(matrix(wts, nrow = 1))
  v <- vect(
    matrix(c(0, 0, 0, 1, 5, 1, 5, 0), ncol = 2, byrow = TRUE),
    type = "polygons"
  )

  ee_sd_pop <- exact_extract(r, st_as_sf(v), 'weighted_stdev', weights = wr)
  expect_equal(ee_sd_pop, wsd_pop, tolerance = 1e-7)

  ee_sd_samp <- ee_sd_pop * sqrt(n_eff / (n_eff - 1))
  expect_equal(ee_sd_samp, wsd_samp, tolerance = 1e-7)
})

test_that("weighted SD is NA if all values are NA", {
  vals <- c(NA, NA, NA, NA, NA)
  wts <- c(1, 2, 3, 4, 5)
  r <- rast(matrix(vals, nrow = 1))
  wr <- rast(matrix(wts, nrow = 1))
  v <- vect(
    matrix(c(0, 0, 0, 1, 5, 1, 5, 0), ncol = 2, byrow = TRUE),
    type = "polygons"
  )
  ee_sd_pop <- exact_extract(r, st_as_sf(v), 'weighted_stdev', weights = wr)
  expect_true(is.na(ee_sd_pop))
})


test_that("pooled weighted SD from chunk summaries matches direct calculation (Hmisc unbiased)", {
  set.seed(1)
  vals <- matrix(runif(20, 1, 10), nrow = 4)
  wts <- matrix(runif(20, 0.5, 2), nrow = 4)
  vals[c(11, 8)] <- NA_real_
  wts[12] <- NA_real_
  r <- rast(vals)
  wr <- rast(wts)

  v1 <- vect(
    matrix(c(0, 0, 0, 2, 2, 2, 2, 0), ncol = 2, byrow = TRUE),
    type = "polygons"
  )
  v2 <- vect(
    matrix(c(2, 0, 2, 2, 4, 2, 4, 0), ncol = 2, byrow = TRUE),
    type = "polygons"
  )
  v_all <- vect(
    matrix(
      c(
        0,
        0,
        0,
        2,
        2,
        2,
        2,
        0,
        2,
        0,
        2,
        2,
        4,
        2,
        4,
        0
      ),
      ncol = 2,
      byrow = TRUE
    ),
    type = "polygons"
  )

  # Helper for chunk stats using Hmisc denominator
  get_chunk_stats <- function(poly) {
    df <- exact_extract(
      r,
      st_as_sf(poly),
      fun = c("mean", "stdev", "count", "weighted_mean", "weighted_stdev"),
      weights = wr,
      default_weight = 0
    )
    rw <- terra::clamp(r, 1, 1)
    W <- exact_extract(
      wr,
      st_as_sf(poly),
      fun = "weighted_sum",
      weights = rw,
      default_weight = 0
    )[[1]]
    W2 <- exact_extract(
      wr^2,
      st_as_sf(poly),
      fun = "weighted_sum",
      weights = rw,
      default_weight = 0
    )[[1]]
    pop_var <- (df$weighted_stdev)^2
    var_unbiased <- pop_var * W / (W - 1)

    lumped_pop_var <- (df$stdev)^2
    lumped_var_unbiased <- lumped_pop_var * df$count / (df$count - 1)

    list(
      mean = df$weighted_mean,
      lumped_mean = df$mean,
      lumped_stdev = df$stdev,
      lumped_count = df$count,
      lumped_pop_var = lumped_pop_var,
      lumped_var_unbiased = lumped_var_unbiased,
      pop_var = pop_var,
      var_unbiased = var_unbiased,
      W = W,
      W2 = W2
    )
  }

  s1 <- get_chunk_stats(v1)
  s2 <- get_chunk_stats(v2)

  # Pooled weighted mean
  mu_pooled <- (s1$mean * s1$W + s2$mean * s2$W) / (s1$W + s2$W)

  # lumped pooled mean
  lumped_mu_pooled <- (s1$lumped_mean *
    s1$lumped_count +
    s2$lumped_mean * s2$lumped_count) /
    (s1$lumped_count + s2$lumped_count)

  # Pooled unbiased variance (Hmisc style)
  num <- (s1$W - 1) *
    s1$var_unbiased +
    (s2$W - 1) * s2$var_unbiased +
    s1$W * (s1$mean - mu_pooled)^2 +
    s2$W * (s2$mean - mu_pooled)^2
  denom <- (s1$W + s2$W) - 1
  pooled_var_unbiased <- num / denom
  pooled_sd <- sqrt(pooled_var_unbiased)

  # Pooled unbiased lumped variance
  num <- (s1$lumped_count - 1) *
    s1$lumped_var_unbiased +
    (s2$lumped_count - 1) * s2$lumped_var_unbiased +
    s1$lumped_count * (s1$lumped_mean - lumped_mu_pooled)^2 +
    s2$lumped_count * (s2$lumped_mean - lumped_mu_pooled)^2
  denom <- (s1$lumped_count + s2$lumped_count) - 1
  lumped_pooled_var_unbiased <- num / denom
  lumped_pooled_sd <- sqrt(lumped_pooled_var_unbiased)

  # Direct calculation on all values
  vals_all0 <- exact_extract(r, st_as_sf(v_all))[[1]]$value
  wts_all0 <- exact_extract(wr, st_as_sf(v_all))[[1]]$value
  ok_all <- !is.na(vals_all0) & !is.na(wts_all0) & wts_all0 > 0
  vals_all <- vals_all0[ok_all]
  wts_all <- wts_all0[ok_all]
  W_all <- sum(wts_all)
  mean_all <- sum(wts_all * vals_all) / W_all
  var_all_unbiased <- sum(wts_all * (vals_all - mean_all)^2) / (W_all - 1)
  sd_all <- sqrt(var_all_unbiased)

  # Hmisc reference
  raw_wtd_sd <- sqrt(Hmisc::wtd.var(vals_all0, wts_all0))
  raw_wtd_mean <- Hmisc::wtd.mean(vals_all0, wts_all0)
  raw_lmp_sd <- sqrt(Hmisc::wtd.var(vals_all0))
  raw_lmp_mean <- Hmisc::wtd.mean(vals_all0)

  # Check means
  expect_equal(mu_pooled, mean_all, tolerance = 1e-7)
  expect_equal(mu_pooled, raw_wtd_mean, tolerance = 1e-7)
  expect_equal(mean_all, raw_wtd_mean, tolerance = 1e-7)

  expect_equal(lumped_mu_pooled, raw_lmp_mean, tolerance = 1e-7)

  # check SDs
  expect_equal(pooled_sd, sd_all, tolerance = 1e-7)
  expect_equal(pooled_sd, raw_wtd_sd, tolerance = 1e-7)
  expect_equal(sd_all, raw_wtd_sd, tolerance = 1e-7)

  expect_equal(lumped_pooled_sd, raw_lmp_sd, tolerance = 1e-7)
})

test_that("combine_* functions match Hmisc and direct calculation using exactextractr chunk output", {
  set.seed(123)
  vals <- matrix(runif(20, 1, 10), nrow = 4)
  wts <- matrix(runif(20, 0.5, 2), nrow = 4)
  vals[c(7, 13)] <- NA_real_
  wts[5] <- NA_real_
  r <- rast(vals)
  wr <- rast(wts)

  v1 <- vect(
    matrix(c(0, 0, 0, 2, 2, 2, 2, 0), ncol = 2, byrow = TRUE),
    type = "polygons"
  )
  v2 <- vect(
    matrix(c(2, 0, 2, 2, 4, 2, 4, 0), ncol = 2, byrow = TRUE),
    type = "polygons"
  )
  v_all <- vect(
    matrix(
      c(0, 0, 0, 2, 2, 2, 2, 0, 2, 0, 2, 2, 4, 2, 4, 0),
      ncol = 2,
      byrow = TRUE
    ),
    type = "polygons"
  )

  # Get chunk stats
  get_chunk_stats <- function(poly) {
    df <- exact_extract(
      r,
      st_as_sf(poly),
      fun = c("mean", "stdev", "count", "weighted_mean", "weighted_stdev"),
      weights = wr,
      default_weight = 0
    )
    rw <- terra::clamp(r, 1, 1)
    W <- exact_extract(
      wr,
      st_as_sf(poly),
      fun = "weighted_sum",
      weights = rw,
      default_weight = 0
    )[[1]]
    lumped_var <- (df$stdev)^2
    weighted_var <- (df$weighted_stdev)^2
    list(
      mean = df$mean,
      var = lumped_var,
      count = df$count,
      weighted_mean = df$weighted_mean,
      weighted_var = weighted_var,
      weighted_sum = W
    )
  }

  s1 <- get_chunk_stats(v1)
  s2 <- get_chunk_stats(v2)

  # Pool using your functions
  lumped_mean <- combine_lumped_mean(
    c(s1$mean, s2$mean),
    c(s1$count, s2$count),
    c(s1$count, s2$count)
  )
  lumped_var <- combine_lumped_sample_var(
    c(s1$mean, s2$mean),
    c(s1$var, s2$var),
    c(s1$count, s2$count),
    c(s1$count, s2$count)
  )
  lumped_sd <- sqrt(lumped_var)

  weighted_mean <- combine_weighted_mean(
    c(s1$weighted_mean, s2$weighted_mean),
    c(s1$weighted_sum, s2$weighted_sum),
    c(s1$weighted_sum, s2$weighted_sum)
  )
  weighted_var <- combine_weighted_sample_var(
    c(s1$weighted_mean, s2$weighted_mean),
    c(s1$weighted_var, s2$weighted_var),
    c(s1$weighted_sum, s2$weighted_sum),
    c(s1$weighted_sum, s2$weighted_sum)
  )
  weighted_sd <- sqrt(weighted_var)

  # Direct calculation on all values
  vals_all0 <- exact_extract(r, st_as_sf(v_all))[[1]]$value
  wts_all0 <- exact_extract(wr, st_as_sf(v_all))[[1]]$value
  ok_all <- !is.na(vals_all0) & !is.na(wts_all0) & wts_all0 > 0
  vals_all <- vals_all0[ok_all]
  wts_all <- wts_all0[ok_all]
  W_all <- sum(wts_all)
  mean_all <- sum(wts_all * vals_all) / W_all
  var_all_unbiased <- sum(wts_all * (vals_all - mean_all)^2) / (W_all - 1)
  sd_all <- sqrt(var_all_unbiased)

  # Hmisc reference
  raw_wtd_mean <- Hmisc::wtd.mean(vals_all0, wts_all0)
  raw_wtd_sd <- sqrt(Hmisc::wtd.var(vals_all0, wts_all0))
  raw_lmp_mean <- Hmisc::wtd.mean(vals_all0)
  raw_lmp_sd <- sqrt(Hmisc::wtd.var(vals_all0))

  # Check means
  expect_equal(lumped_mean, raw_lmp_mean, tolerance = 1e-7)
  expect_equal(weighted_mean, raw_wtd_mean, tolerance = 1e-7)

  # Check SDs
  expect_equal(lumped_sd, raw_lmp_sd, tolerance = 1e-7)
  expect_equal(weighted_sd, raw_wtd_sd, tolerance = 1e-7)
  expect_equal(weighted_sd, sd_all, tolerance = 1e-7)
})
