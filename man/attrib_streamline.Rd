% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attrib_streamline.R
\name{attrib_streamline}
\alias{attrib_streamline}
\title{Generate and attribute stream line and points along the lines}
\usage{
attrib_streamline(
  input,
  extra_attr = c("link_slope", "cont_slope", "USChnLn_To", "Elevation", "StOrd_Hack",
    "StOrd_Str", "StOrd_Hort", "StOrd_Shr"),
  points = NULL,
  site_id_col = NULL,
  snap_distance = 10L,
  break_on_noSnap = T,
  return_products = F,
  temp_dir = NULL,
  compress = F,
  verbose = F
)
}
\arguments{
\item{input}{resulting object from \code{generate_subbasins()}}

\item{extra_attr}{character. One or more of c("link_slope", "cont_slope", "USChnLn_To", "Elevation", "StOrd_Hack", "StOrd_Str", "StOrd_Hort", "StOrd_Shr"). Optional attributes to add to stream vector outputs.}

\item{points}{character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/points.shp"), or any GIS data object that will be converted to spatial points. Points representing sampling locations.}

\item{site_id_col}{character. Variable name in \code{points} that corresponds to unique site identifiers. This column will be included in all vector geospatial analysis products. Note, if multiple points have the same \code{site_id_col}, their centroid will be used and returned; if multiple points overlap after snapping, only the first is used.}

\item{snap_distance}{integer. Maximum distance which points will be snapped to stream lines in map units}

\item{break_on_noSnap}{logical. Should the function stop if any points are not snapped to a stream segment (i.e., are beyon snap_distance)}

\item{return_products}{logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.}

\item{temp_dir}{character. File path for temporary file storage, If \code{NULL}, \code{tempfile()} will be used}

\item{compress}{logical. Should output rasters be compressed, slower but more space efficient.}

\item{verbose}{logical.}
}
\value{
If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
}
\description{
Generate and attribute stream line and points along the lines
}
