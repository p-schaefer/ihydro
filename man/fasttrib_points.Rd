% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fasttrib_points.R
\name{fasttrib_points}
\alias{fasttrib_points}
\title{Quickly attribute stream segments/sampling points with layers of interest (loi)}
\usage{
fasttrib_points(
  input,
  out_filename,
  loi_file = NULL,
  loi_cols = NULL,
  iDW_file = NULL,
  store_iDW = F,
  sample_points = NULL,
  link_id = NULL,
  target_o_type = c("point", "segment_point", "segment_whole"),
  weighting_scheme = c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO"),
  loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
  inv_function = function(x) {
     (x * 0.001 + 1)^-1
 },
  temp_dir = NULL,
  verbose = F,
  backend = c("tibble", "SQLite", "data.table"),
  SQLdb_path = ":memory:"
)
}
\arguments{
\item{input}{output from \code{process_hydrology()} (if \code{process_loi()} was not run on \code{process_hydrology()}, \code{loi_file} must be specified)}

\item{out_filename}{Output file name, ending in '.csv'.}

\item{loi_file}{filepath of \code{process_loi()} output (optional, will search for loi layers in 'input' if not specified).}

\item{loi_cols}{character or NULL. Names of loi layers to include in summary. If NULL, all layers used.}

\item{iDW_file}{filepath of \code{prep_weights()} output (optional, will search for weight layers in 'input' if not specified).}

\item{store_iDW}{logical, if 'iDW_file' is NULL, weights will be added to 'input'}

\item{sample_points}{character or NULL. IDs of unique station identifiers provided in 'site_id_col' of \code{generate_vectors()} to calculate attributes for.}

\item{link_id}{character or NULL. 'link_id's of reaches to calculate attributes for.}

\item{target_o_type}{character. One of: c("point","segment_point","segment_whole"). Target for"iFLO", and "HAiFLO" weighting schemes. "Point" represents the sampling point on the stream, "segment_point" represents the upstream segment of the sampling points, and "segment_whole" will target the entire reach, regardless of where sampling occurred.}

\item{weighting_scheme}{character. One or more weighting schemes: c("lumped","iFLO", "iFLS", "HAiFLO", "HAiFLS")}

\item{loi_numeric_stats}{character. One or more of c("mean", "sd", "median", "min", "max", "sum"). Only distance-weighted versions of mean and SD are returned for all weighting schemes except lumped.}

\item{inv_function}{function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.}

\item{temp_dir}{character. File path for intermediate products; these are deleted once the function runs successfully.}

\item{verbose}{logical.}

\item{backend}{character. One of "data.table" or "SQLite"}

\item{SQLdb_path}{character. If backend = "SQLite" directory to store temporary SQL tables, or ":memory:" for an in-memory SQL table}
}
\value{
A data.frame of weighted attributes for the requested areas
}
\description{
Quickly attribute stream segments/sampling points with layers of interest (loi)
}
