#' Benchmark parallel vs. sequential exact_extract
#'
#' @param rast A terra SpatRaster.
#' @param polys An sf or SpatVector polygon object.
#' @param weights Optional weights for the extraction.
#' @param fun Summary function(s) to use.
#' @param n_workers Number of parallel workers.
#' @param max_cells_in_memory Passed to exact_extract.
#' @param n_reps Number of repetitions for timing.
#' @return Data.frame with timing and memory results.
benchmark_exact_extract <- function(
  rast,
  polys,
  weights = NULL,
  fun = "mean",
  n_workers = 2,
  chunks_per_worker = 1L,
  max_cells_in_memory = 1e7,
  n_reps = 3
) {
  td <- tempfile()
  dir.create(td)
  on.exit(suppressWarnings(unlink(td, recursive = TRUE)), add = TRUE)
  ras_file <- file.path(td, "rast.tif")
  terra::writeRaster(rast, ras_file, overwrite = TRUE)
  polys_file <- file.path(td, "polys.gpkg")
  sf::st_write(polys, polys_file, delete_dsn = TRUE, quiet = TRUE)

  weights_file <- NULL
  if (!is.null(weights)) {
    weights_file <- file.path(td, "weights.tif")
    terra::writeRaster(weights, weights_file, overwrite = TRUE)
  }

  resolve_weights <- function(weights_file) {
    if (!is.null(weights_file)) {
      return(terra::rast(weights_file))
    }
    NULL
  }

  # Sequential
  seq_time <- microbenchmark::microbenchmark(
    exactextractr::exact_extract(
      x = terra::rast(ras_file),
      y = sf::st_read(polys_file, quiet = TRUE),
      weights = resolve_weights(weights_file),
      fun = fun,
      progress = FALSE,
      max_cells_in_memory = max_cells_in_memory
    ),
    times = n_reps
  )

  # Parallel
  future::plan(future::multisession, workers = n_workers)
  n <- nrow(polys)
  idx <- split(
    seq_len(n),
    cut(seq_len(n), n_workers, labels = FALSE)
  )

  par_time <- microbenchmark::microbenchmark(
    {
      res <- future.apply::future_lapply(
        idx,
        function(i) {
          exactextractr::exact_extract(
            x = terra::rast(ras_file),
            y = sf::st_read(polys_file, quiet = TRUE)[i, ],
            weights = resolve_weights(weights_file),
            fun = fun,
            progress = FALSE,
            max_cells_in_memory = floor(
              (max_cells_in_memory / n_workers) * 0.8
            )
          )
        },
        future.globals = list(
          ras_file = ras_file,
          polys_file = polys_file,
          weights_file = weights_file,
          fun = fun,
          max_cells_in_memory = max_cells_in_memory,
          n_workers = n_workers
        ),
        future.seed = NULL,
        future.packages = c("exactextractr", "terra", "sf")
      )
      do.call(rbind, res)
    },
    times = n_reps
  )

  if (chunks_per_worker > 1) {
    idx <- split(
      seq_len(n),
      cut(seq_len(n), n_workers * chunks_per_worker, labels = FALSE)
    )

    chunk_par_time <- microbenchmark::microbenchmark(
      {
        res <- future.apply::future_lapply(
          idx,
          function(i) {
            exactextractr::exact_extract(
              x = terra::rast(ras_file),
              y = sf::st_read(polys_file, quiet = TRUE)[i, ],
              weights = resolve_weights(weights_file),
              fun = fun,
              progress = FALSE,
              max_cells_in_memory = floor(
                (max_cells_in_memory / n_workers) * 0.8
              )
            )
          },
          future.globals = list(
            ras_file = ras_file,
            polys_file = polys_file,
            weights_file = weights_file,
            fun = fun,
            max_cells_in_memory = max_cells_in_memory,
            n_workers = n_workers
          ),
          future.seed = NULL,
          future.packages = c("exactextractr", "terra", "sf"),
          future.scheduling = chunks_per_worker
        )
        do.call(rbind, res)
      },
      times = n_reps
    )

    chunk_par_time <- data.frame(
      method = c("chunked"),
      median_time = median(chunk_par_time$time / 1e9), # seconds
      min_time = min(chunk_par_time$time) / 1e9,
      max_time = max(chunk_par_time$time) / 1e9
    )
  }
  future::plan(future::sequential)

  out <- data.frame(
    method = c("sequential", "parallel"),
    median_time = c(median(seq_time$time) / 1e9, median(par_time$time) / 1e9), # seconds
    min_time = c(min(seq_time$time) / 1e9, min(par_time$time) / 1e9),
    max_time = c(max(seq_time$time) / 1e9, max(par_time$time) / 1e9)
  )

  if (chunks_per_worker > 1) {
    out <- rbind(out, chunk_par_time)
  }
  return(out)
}
