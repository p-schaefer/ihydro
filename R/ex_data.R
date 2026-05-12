#' Access example data layers included in the package
#'
#' @param layer Character. Name of the example layer to load. Must be one of:
#'   "elev_ned_30m.tif", "geology.shp", "lakes.shp", "landuse_r.tif",
#'   "pointsources.shp", "sites_nc.shp", "streams.shp"
#'
#' @return A `SpatRaster` or `SpatVector` object, depending on the layer.
#'
#' @details Original source of data: https://github.com/MiKatt/openSTARS
#'
#' @examples
#' \dontrun{
#' ex_data("elev_ned_30m.tif") # example DEM raster
#' ex_data("geology.shp") # example geology vector
#' ex_data("lakes.shp") # example lakes vector
#' ex_data("landuse_r.tif") # example land use raster
#' ex_data("pointsources.shp") # example point sources vector
#' ex_data("sites_nc.shp") # example sampling sites vector
#' ex_data("streams.shp") # example stream network vector
#' }
#'
#' @export
ex_data <- function(
  layer = c(
    "elev_ned_30m.tif",
    "geology.shp",
    "lakes.shp",
    "landuse_r.tif",
    "pointsources.shp",
    "sites_nc.shp",
    "streams.shp"
  )
) {
  layer <- match.arg(layer)

  ex_files <- system.file(
    "extdata",
    "openSTARS_data.zip",
    package = "ihydro"
  ) |>
    normalizePath() |>
    utils::unzip(list = TRUE)

  ex_files <- ex_files$Name[grepl("\\.shp|\\.tif", ex_files$Name)]
  if (!layer %in% basename(ex_files)) {
    cli::cli_abort(
      c(
        "{.arg layer} must be a valid example file name.",
        "i" = "Valid options are: {.val {basename(ex_files)}}."
      )
    )
  }

  final_file <- file.path(
    "/vsizip",
    system.file("extdata", "openSTARS_data.zip", package = "ihydro"),
    layer
  )

  file_ext <- function(x) {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }

  switch(
    file_ext(layer),
    "shp" = terra::vect(final_file),
    "tif" = terra::rast(final_file)
  )
}
