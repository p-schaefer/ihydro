# ============================================================================
# Benchmark: fasttrib_points() vs exactextractr::exact_extract()
# ============================================================================
#
# Compares:
#   1. Speed of lumped AND weighted catchment-level summarisation
#   2. Numerical agreement of results (mean, sum, sd, weighted mean, weighted sd)
#
# Requirements:
#   - ihydro (this package)
#   - exactextractr
#   - bench (for benchmarking)
#
# Run with: source("inst/benchmarks/benchmark_fasttrib_vs_exactextract.R")
# ============================================================================

library(ihydro)
library(terra)
library(sf)
library(dplyr)
library(exactextractr)
library(bench)

future::plan("sequential") # keep it simple; adjust for parallel tests

# ── 1. Build hydrology + LOI from bundled example data ──────────────────────

ex_loc <- tempfile(pattern = "bench_ihydro_")
dir.create(ex_loc, recursive = TRUE)

dem <- ex_data("elev_ned_30m.tif")
sites <- ex_data("sites_nc.shp")
strms <- ex_data("streams.shp")

writeRaster(dem, file.path(ex_loc, "dem.tif"), overwrite = TRUE)
writeVector(sites, file.path(ex_loc, "sites.shp"), overwrite = TRUE)
writeVector(strms, file.path(ex_loc, "streams.shp"), overwrite = TRUE)

# Landcover raster (categorical)
lc <- ex_data("landuse_r.tif") |> setNames("LC")
writeRaster(lc, file.path(ex_loc, "LC.tif"), overwrite = TRUE)

# Slope raster (numeric)
whitebox::wbt_slope(
  dem = file.path(ex_loc, "dem.tif"),
  output = file.path(ex_loc, "slope.tif")
)

gpkg_hydro <- file.path(ex_loc, "hydro.gpkg")

message("Processing hydrology...")
hydro <- process_hydrology(
  dem = dem,
  output_filename = gpkg_hydro,
  threshold = 500L,
  burn_streams = file.path(ex_loc, "streams.shp"),
  burn_depth = 5,
  burn_slope_dist = 250,
  burn_slope_depth = 5,
  points = file.path(ex_loc, "sites.shp"),
  site_id_col = "site_id",
  snap_distance = NULL,
  break_on_noSnap = FALSE,
  pwise_dist = TRUE,
  pwise_all_links = FALSE,
  return_products = TRUE,
  verbose = FALSE
)

gpkg_loi <- file.path(ex_loc, "loi.gpkg")

message("Processing LOI...")
loi <- process_loi(
  dem = dem,
  num_inputs = list(slope = file.path(ex_loc, "slope.tif")),
  cat_inputs = list(landcover = file.path(ex_loc, "LC.tif")),
  output_filename = gpkg_loi,
  return_products = TRUE,
  verbose = FALSE
)

# Pick a handful of sample points
sp_ids <- as.character(hydro$snapped_points$site_id[1:5])

message("Setup complete. Sample points: ", paste(sp_ids, collapse = ", "))


# ── 2. Run fasttrib_points (lumped + iFLS weighted, mean + sd + sum) ────────

out_csv <- file.path(ex_loc, "fasttrib_result.csv")

message("Running fasttrib_points() [lumped + iFLS]...")
ft_time <- system.time({
  ft_result <- fasttrib_points(
    input = hydro,
    out_filename = out_csv,
    loi_file = loi,
    loi_cols = NULL,
    sample_points = sp_ids,
    target_o_type = "point",
    weighting_scheme = c("lumped", "iFLS"),
    loi_numeric_stats = c("mean", "sd", "sum"),
    inv_function = function(x) (x * 0.001 + 1)^-1,
    temp_dir = NULL,
    verbose = FALSE,
    backend = "tibble"
  )
})


# ── 3. Run exact_extract on the same catchments & rasters ───────────────────

# Get catchment polygons for the same points
catch_polys <- get_catchment(
  hydro,
  sample_points = sp_ids,
  verbose = FALSE
)

# Read LOI rasters from the gpkg
loi_meta <- read_ihydro(loi, "loi_meta")
loi_cols <- loi_meta$loi_var_nms
loi_rast <- rast(loi$outfile, loi_cols)

num_cols <- loi_meta$loi_var_nms[loi_meta$loi_type == "num_rast"]
cat_cols <- loi_meta$loi_var_nms[loi_meta$loi_type == "cat_rast"]

# Prepare iFLS weight raster from the iDW file created by fasttrib_points
# fasttrib_points stores weights alongside the hydro gpkg or a temp file;
# we re-derive them here with prep_weights for a clean comparison.
message("Preparing iFLS weights for exact_extract...")
idw_gpkg <- file.path(ex_loc, "weights.gpkg")
sf::write_sf(
  read_ihydro(hydro, "DEM_Extent"),
  idw_gpkg,
  layer = "DEM_Extent",
  append = TRUE,
  delete_layer = FALSE,
  delete_dsn = FALSE
)
idw_obj <- prep_weights(
  input = hydro,
  output_filename = idw_gpkg,
  sample_points = sp_ids,
  target_o_type = "point",
  weighting_scheme = "iFLS",
  inv_function = function(x) (x * 0.001 + 1)^-1,
  verbose = FALSE
)
ifls_rast <- rast(idw_obj$outfile, "iFLS")

message("Running exact_extract()...")
ee_time <- system.time({
  # ── 3a. Lumped: mean, sum, sd ──────────────────────────────────────────
  if (length(num_cols) > 0) {
    ee_mean <- exact_extract(
      subset(loi_rast, num_cols),
      catch_polys,
      fun = "mean",
      append_cols = "link_id"
    )
    colnames(ee_mean) <- c("link_id", paste0(num_cols, "_lumped_mean"))

    ee_sum <- exact_extract(
      subset(loi_rast, num_cols),
      catch_polys,
      fun = "sum",
      append_cols = "link_id"
    )
    colnames(ee_sum) <- c("link_id", paste0(num_cols, "_lumped_sum"))

    ee_sd <- exact_extract(
      subset(loi_rast, num_cols),
      catch_polys,
      fun = "stdev",
      append_cols = "link_id"
    )
    colnames(ee_sd) <- c("link_id", paste0(num_cols, "_lumped_sd"))
  } else {
    ee_mean <- ee_sum <- ee_sd <- tibble(link_id = catch_polys$link_id)
  }

  # Proportion (mean of 0/1) for categorical layers
  if (length(cat_cols) > 0) {
    ee_prop <- exact_extract(
      subset(loi_rast, cat_cols),
      catch_polys,
      fun = "mean",
      append_cols = "link_id"
    )
    colnames(ee_prop) <- c("link_id", paste0(cat_cols, "_lumped_prop"))
  } else {
    ee_prop <- tibble(link_id = catch_polys$link_id)
  }

  # ── 3b. iFLS-weighted mean ────────────────────────────────────────────
  # exact_extract weighted_mean uses coverage_fraction * weight_raster as w
  # fasttrib multiplies LOI * iFLS then does sum(LOI*w)/sum(w) — same formula
  if (length(num_cols) > 0) {
    ee_wmean <- exact_extract(
      subset(loi_rast, num_cols),
      catch_polys,
      fun = "weighted_mean",
      weights = ifls_rast,
      append_cols = "link_id"
    )
    colnames(ee_wmean) <- c("link_id", paste0(num_cols, "_iFLS_mean"))
  } else {
    ee_wmean <- tibble(link_id = catch_polys$link_id)
  }

  # iFLS-weighted proportion for categorical layers
  if (length(cat_cols) > 0) {
    ee_wprop <- exact_extract(
      subset(loi_rast, cat_cols),
      catch_polys,
      fun = "weighted_mean",
      weights = ifls_rast,
      append_cols = "link_id"
    )
    colnames(ee_wprop) <- c("link_id", paste0(cat_cols, "_iFLS_prop"))
  } else {
    ee_wprop <- tibble(link_id = catch_polys$link_id)
  }

  # ── 3c. iFLS-weighted SD ──────────────────────────────────────────────
  # fasttrib uses reliability-weighted SD:
  #   sqrt( sum(w * (x - wmean)^2) / ((N-1)/N * sum(w)) )
  # where N = count of non-zero weights and wmean = sum(w*x)/sum(w)
  # This requires a custom function in exact_extract.
  if (length(num_cols) > 0) {
    ee_wsd_list <- exact_extract(
      c(subset(loi_rast, num_cols), ifls_rast),
      catch_polys,
      fun = function(df) {
        w <- df[["iFLS"]]
        # match fasttrib: only non-NA cells matter
        result <- vapply(
          num_cols,
          function(col) {
            x <- df[[col]]
            ok <- !is.na(x) & !is.na(w)
            xo <- x[ok]
            wo <- w[ok]
            N <- sum(wo != 0)
            if (N < 2) {
              return(NA_real_)
            }
            wm <- sum(xo * wo, na.rm = TRUE) / sum(wo, na.rm = TRUE)
            sqrt(
              sum(wo * (xo - wm)^2, na.rm = TRUE) /
                ((N - 1) / N * sum(wo, na.rm = TRUE))
            )
          },
          numeric(1)
        )
        setNames(as.list(result), paste0(num_cols, "_iFLS_sd"))
      },
      append_cols = "link_id",
      summarize_df = TRUE
    )
    ee_wsd <- bind_rows(ee_wsd_list)
    if (!"link_id" %in% names(ee_wsd)) {
      ee_wsd <- bind_cols(tibble(link_id = catch_polys$link_id), ee_wsd)
    }
  } else {
    ee_wsd <- tibble(link_id = catch_polys$link_id)
  }

  # ── Combine all exact_extract results ──────────────────────────────────
  ee_result <- ee_mean |>
    left_join(ee_sum, by = "link_id") |>
    left_join(ee_sd, by = "link_id") |>
    left_join(ee_prop, by = "link_id") |>
    left_join(ee_wmean, by = "link_id") |>
    left_join(ee_wprop, by = "link_id") |>
    left_join(ee_wsd, by = "link_id") |>
    mutate(link_id = as.character(link_id))
})


# ── 4. Compare speed ───────────────────────────────────────────────────────

message("\n=== Timing ===")
message("fasttrib_points : ", round(ft_time["elapsed"], 2), "s")
message("exact_extract   : ", round(ee_time["elapsed"], 2), "s")


# ── 5. Compare numerical results ──────────────────────────────────────────

ft_df <- ft_result |>
  as_tibble() |>
  mutate(link_id = as.character(link_id)) |>
  select(link_id, any_of(colnames(ee_result)))

comparison <- inner_join(
  ft_df,
  ee_result,
  by = "link_id",
  suffix = c(".ft", ".ee")
)

shared_cols <- setdiff(
  intersect(colnames(ft_df), colnames(ee_result)),
  "link_id"
)

message("\n=== Numerical Comparison ===")
diffs <- purrr::map_dfr(shared_cols, function(col) {
  v_ft <- comparison[[paste0(col, ".ft")]]
  v_ee <- comparison[[paste0(col, ".ee")]]
  ok <- !is.na(v_ft) & !is.na(v_ee)
  if (!any(ok)) {
    return(tibble(
      column = col,
      max_diff = NA_real_,
      mean_diff = NA_real_,
      pct_diff = NA_real_,
      n = 0L
    ))
  }
  d <- abs(v_ft[ok] - v_ee[ok])
  avg_mag <- mean(abs(v_ee[ok]), na.rm = TRUE)
  tibble(
    column = col,
    max_diff = max(d),
    mean_diff = mean(d),
    pct_diff = if (avg_mag > 0) mean(d) / avg_mag * 100 else NA_real_,
    n = sum(ok)
  )
})

print(diffs, n = 50)


# ── 6. Formal micro-benchmark (optional, more precise) ────────────────────

message("\n=== Micro-benchmark (3 iterations) ===")
bm <- bench::mark(
  fasttrib = {
    fasttrib_points(
      input = hydro,
      out_filename = tempfile(fileext = ".csv"),
      loi_file = loi,
      sample_points = sp_ids,
      target_o_type = "point",
      weighting_scheme = c("lumped", "iFLS"),
      loi_numeric_stats = c("mean", "sd", "sum"),
      inv_function = function(x) (x * 0.001 + 1)^-1,
      verbose = FALSE,
      backend = "tibble"
    )
  },
  exact_extract = {
    ee_m <- exact_extract(subset(loi_rast, num_cols), catch_polys, "mean")
    ee_s <- exact_extract(subset(loi_rast, num_cols), catch_polys, "sum")
    ee_d <- exact_extract(subset(loi_rast, num_cols), catch_polys, "stdev")
    ee_p <- exact_extract(subset(loi_rast, cat_cols), catch_polys, "mean")
    ee_wm <- exact_extract(
      subset(loi_rast, num_cols),
      catch_polys,
      "weighted_mean",
      weights = ifls_rast
    )
    ee_wp <- exact_extract(
      subset(loi_rast, cat_cols),
      catch_polys,
      "weighted_mean",
      weights = ifls_rast
    )
    ee_ws <- exact_extract(
      c(subset(loi_rast, num_cols), ifls_rast),
      catch_polys,
      fun = function(df) {
        w <- df[["iFLS"]]
        vapply(
          num_cols,
          function(col) {
            x <- df[[col]]
            ok <- !is.na(x) & !is.na(w)
            xo <- x[ok]
            wo <- w[ok]
            N <- sum(wo != 0)
            if (N < 2) {
              return(NA_real_)
            }
            wm <- sum(xo * wo) / sum(wo)
            sqrt(sum(wo * (xo - wm)^2) / ((N - 1) / N * sum(wo)))
          },
          numeric(1)
        )
      }
    )
    bind_cols(ee_m, ee_s, ee_d, ee_p, ee_wm, ee_wp, ee_ws)
  },
  iterations = 3,
  check = FALSE
)

print(bm)
