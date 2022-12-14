% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_loi.R
\name{process_loi}
\alias{process_loi}
\title{Processes layers of interest (loi) for attribution of stream network}
\usage{
process_loi(
  input = NULL,
  dem = NULL,
  clip_region = NULL,
  num_inputs = list(),
  cat_inputs = list(),
  output_filename = NULL,
  variable_names = NULL,
  return_products = F,
  temp_dir = NULL,
  verbose = F,
  overwrite = T
)
}
\arguments{
\item{input}{\code{NULL} or output from \code{process_hydrology()} or \code{process_flowdir()}. If \code{NULL}, \code{dem} must be specified.}

\item{dem}{\code{NULL} or character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster. If \code{NULL}, input must be specified.}

\item{clip_region}{character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of ESRI Shapefile type or GeoTiFF type only. Region over which loi are returned. defaults to \code{dem} extent.}

\item{num_inputs}{named list containing file paths (with extensions), and/or R GIS objects to be coerced into \code{SpatRaster}. Layers of interest to be summarized numerically (e.g., with mean, SD, etc.)}

\item{cat_inputs}{named list containing file paths (with extensions), and/or R GIS objects to be coerced into \code{SpatRaster}. Layers of interest to be summarized categorical (e.g., with with proportions or weighted proportions in the catchment)}

\item{output_filename}{\code{NULL} or character (full file path with extensions). If \code{input} is provided, loi are saved in the same .zip output folder regadless of \code{output_filename}. If \code{input} not provided, file path to write resulting .zip file.}

\item{variable_names}{named list containing variables from \code{num_inputs} and \code{cat_inputs} to include in output. If not specied, all variables are used.}

\item{return_products}{logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.}

\item{temp_dir}{character. File path for intermediate products; these are deleted once the function runs successfully.}

\item{verbose}{logical.}
}
\value{
If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
}
\description{
Processes layers of interest (loi) for attribution of stream network
}
