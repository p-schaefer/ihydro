#' Process and align spatial inputs (refactored to use package helpers)
#'
#' This version reuses ensure_temp_dir, set_terra_options, and restore_terra_options.
#' It adds write_for_whitebox to produce a single GeoTIFF path for downstream tools.
#' @noRd
process_input <- function(
  input = NULL,
  input_variable_names = NULL,
  align_to = NULL,
  clip_region = NULL,
  resample_type = c("bilinear", "near"),
  snap = c("near", "in", "out"),
  working_dir = NULL,
  ...
) {
  if (is.null(input)) {
    return(input)
  }

  resample_type <- match.arg(resample_type, c("bilinear", "near"))
  snap <- match.arg(snap, c("near", "in", "out"))

  # stable working dir for this call
  working_root <- ensure_temp_dir(working_dir)
  # set terra options for this call and capture old options for restore
  n_cores <- n_workers()
  old_terra_opts <- set_terra_options(
    n_cores = n_cores,
    temp_dir = working_root,
    verbose = FALSE
  )
  on.exit(
    {
      tryCatch(restore_terra_options(old_terra_opts), error = function(e) NULL)
    },
    add = TRUE
  )

  tmp_file <- function(ext = ".tif", prefix = "tmp_") {
    file.path(working_root, paste0(prefix, basename(tempfile()), ext))
  }

  # Coerce input to terra object
  if (is.character(input)) {
    if (grepl("\\.shp$", input, ignore.case = TRUE)) {
      output <- sf::read_sf(input)
      output <- terra::vect(output)
    } else if (grepl("\\.(tif|tiff)$", input, ignore.case = TRUE)) {
      output <- terra::rast(input)
    } else {
      cli::cli_abort(
        "{.arg input} must end .shp, .tif, or .tiff"
      )
    }
  } else if (inherits(input, c("RasterLayer", "PackedSpatRaster"))) {
    output <- terra::rast(input)
  } else if (inherits(input, c("PackedSpatVector", "sf", "sfc"))) {
    output <- terra::vect(input)
  } else if (inherits(input, c("SpatRaster", "SpatVector"))) {
    output <- input
  } else {
    stop("'input' is not a supported spatial type")
  }

  if (
    is.null(terra::crs(output)) ||
      is.na(terra::crs(output)) ||
      terra::crs(output) == ""
  ) {
    cli::cli_abort(
      "{.arg output} crs() is NULL or NA. Apply a CRS before continuing."
    )
  }

  orig_output <- output

  # Subset variables/layers
  if (is.null(input_variable_names)) {
    input_variable_names <- names(output)
  }
  if (any(!input_variable_names %in% names(output))) {
    cli::cli_abort(
      "Some {.arg input_variable_names} not present in {.arg input}"
    )
  }
  if (inherits(output, "SpatVector") && length(input_variable_names) > 0) {
    output <- output[, input_variable_names, drop = FALSE]
  }
  if (inherits(output, "SpatRaster") && length(input_variable_names) > 0) {
    output <- terra::subset(output, input_variable_names)
  }

  # Align to align_to if provided (reuse same working_root)
  if (!is.null(align_to)) {
    align_to <- process_input(
      input = align_to,
      working_dir = working_root
    )

    if (is.null(terra::crs(align_to)) || is.na(terra::crs(align_to))) {
      cli::cli_abort(
        "{.arg align_to} crs() is NULL or NA. Apply a CRS before continuing."
      )
    }

    # Clip / mask if requested
    if (!is.null(clip_region)) {
      align_crs <- sf::st_crs(output)
      if (!is.null(align_to)) {
        align_crs <- sf::st_crs(align_to)
      }
      cr <- process_input(
        input = clip_region,
        align_to = terra::vect(sf::st_as_sfc(
          "POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",
          crs = align_crs
        )),
        working_dir = working_root
      )
      if (inherits(output, "SpatVector")) {
        output <- terra::crop(output, cr)
      } else if (inherits(output, "SpatRaster")) {
        output <- terra::crop(
          x = output,
          y = cr,
          snap = snap,
          mask = TRUE,
          overwrite = TRUE
        )
        output <- terra::mask(output, cr, overwrite = TRUE)
      }
    }

    # Align vector -> raster or raster -> raster
    did_categorical_rasterization <- FALSE # flag to track if we did categorical rasterization for later splitting
    if (inherits(align_to, "SpatRaster")) {
      if (inherits(output, "SpatVector")) {
        output <- terra::project(output, align_to)
        if (resample_type == "near") {
          # categorical rasterization: write minimal temporary rasters
          out_layers <- list()
          rasters <- lapply(input_variable_names, function(varname) {
            var_vals <- output[[varname]][[1]]
            vals <- unique(var_vals[!is.na(var_vals)])
            for (v in vals) {
              vx <- output[var_vals == v, ]
              out_file <- tmp_file(".tif")
              rr <- terra::rasterize(
                vx,
                align_to,
                field = "",
                overwrite = TRUE,
                filename = out_file,
                wopt = list(gdal = "COMPRESS=NONE", datatype = "FLT4S"),
                ...
              )
              names(rr) <- paste0(varname, "_", v)
              out_layers <- c(out_layers, list(rr))
            }
            output <- terra::rast(out_layers)
          })
          output <- terra::rast(rasters)
          did_categorical_rasterization <- TRUE
        } else {
          rasters <- lapply(input_variable_names, function(varname) {
            out_file <- tmp_file(".tif")
            rr <- terra::rasterize(
              output,
              align_to,
              field = varname,
              overwrite = TRUE,
              filename = out_file,
              wopt = list(gdal = "COMPRESS=NONE", datatype = "FLT4S"),
              ...
            )
            names(rr) <- varname
            rr
          })
          output <- terra::rast(rasters)
        }
      }

      if (inherits(output, "SpatRaster")) {
        need_reproj <- !terra::compareGeom(
          output,
          align_to,
          lyrs = FALSE,
          crs = TRUE,
          ext = TRUE,
          rowcol = TRUE,
          res = TRUE,
          warncrs = FALSE,
          stopOnError = FALSE,
          messages = FALSE
        )
        if (need_reproj) {
          orig_lyr_names <- names(output)
          output <- terra::project(
            output,
            align_to,
            method = resample_type,
            overwrite = TRUE,
            ...
          )
          names(output) <- orig_lyr_names
        }
      }
    }

    if (inherits(align_to, "SpatVector")) {
      if (inherits(output, "SpatRaster")) {
        if (terra::is.polygons(align_to)) {
          output <- terra::as.polygons(output, ...)
        }
        if (terra::is.points(align_to)) {
          output <- terra::as.points(output, ...)
        }
        if (terra::is.lines(align_to)) {
          output <- terra::as.lines(output, ...)
        }
      }
      if (inherits(output, "SpatVector")) {
        output <- terra::project(output, align_to, ...)
      }
    }
  }

  # Split categorical rasters for resample_type = "near"
  if (
    inherits(output, "SpatRaster") &&
      resample_type == "near" &&
      !did_categorical_rasterization
  ) {
    if (terra::nlyr(output) > 1) {
      brick_list <- terra::split(output, names(output))
    } else {
      brick_list <- list(output)
      names(brick_list) <- names(output)
    }
    out_list <- lapply(names(brick_list), function(y) {
      lyr <- brick_list[[y]]
      vals <- terra::unique(as.numeric(lyr))[[1]]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 1) {
        outs <- lapply(vals, function(v) {
          out <- terra::classify(
            as.numeric(lyr),
            cbind(v, 1),
            others = NA_integer_
          )
          names(out) <- paste0(y, "_", v)
          out
        })
        terra::rast(outs)
      } else {
        lyr
      }
    })
    output <- terra::rast(out_list)
  }

  # final sanity check
  if (
    is.null(terra::crs(output)) ||
      is.na(terra::crs(output)) ||
      is.na(terra::crs(output) == "")
  ) {
    cli::cli_abort("{.arg output} crs() is NULL or NA after processing")
  }

  output
}
