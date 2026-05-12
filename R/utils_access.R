#' List layers stored in an ihydro GeoPackage
#'
#' Queries the GeoPackage database to enumerate every raster, vector, and
#' plain-table layer it contains, along with a logical grouping label
#' (`data_group`). This is useful for inspecting what has been computed so
#' far and for verifying that all required layers are present before running
#' downstream functions.
#'
#' @details
#' Layers are classified into the following groups:
#' \describe{
#'   \item{`hydro`}{DEM-derived rasters and stream network vectors
#'     (e.g., `dem_final`, `stream_links`, `Subbasins_poly`).}
#'   \item{`flow_path`}{Upstream / downstream lookup tables (`us_flowpaths`,
#'     `ds_flowpaths`).}
#'   \item{`pwise_dist`}{Pairwise distance matrices (`fcon_pwise_dist`,
#'     `funcon_pwise_dist`).}
#'   \item{`sample_points`}{Original and snapped sampling-point layers.}
#'   \item{`loi`}{Layers of interest (landscape rasters processed by
#'     [process_loi()]).}
#'   \item{`iDW`}{Inverse distance weight rasters (`iFLS`, `HAiFLO`, etc.).}
#'   \item{`meta`}{Metadata tables (`site_id_col`, `loi_meta`, etc.).}
#' }
#'
#' @param input An object of class `ihydro`, as returned by [as.ihydro()],
#'   [process_flowdir()], or any other ihydro function.
#'
#' @return A [tibble][tibble::tibble] with columns:
#' \describe{
#'   \item{`layer_name`}{Name of the layer in the GeoPackage.}
#'   \item{`data_type`}{One of `"Raster"`, `"Vector"`, or `"Table"`.}
#'   \item{`data_group`}{Logical grouping label (see **Details**).}
#' }
#'
#' @seealso [as.ihydro()] to create an `ihydro` object.
#'
#' @examples
#' \dontrun{
#' hydro <- as.ihydro("my_watershed.gpkg")
#'
#' # List everything
#' ihydro_layers(hydro)
#'
#' # Filter to just the hydrology layers
#' ihydro_layers(hydro) |>
#'   dplyr::filter(data_group == "hydro")
#' }
#'
#' @export
ihydro_layers <- function(input) {
  .check_ihydro(input)

  con <- DBI::dbConnect(RSQLite::SQLite(), input$outfile)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Raster layers
  rast_lyrs <- dplyr::tbl(con, "gpkg_contents") |>
    dplyr::collect() |>
    dplyr::filter(grepl("gridded", data_type)) |>
    dplyr::transmute(layer_name = table_name, data_type = "Raster")

  # Vector layers (identified by rtree spatial index)
  all_tables <- DBI::dbListTables(con)

  vect_lyrs <- tibble::tibble(
    layer_name = gsub(
      "^rtree_|_geom$",
      "",
      all_tables[
        grepl("^rtree_", all_tables) &
          grepl("_geom$", all_tables)
      ]
    ),
    data_type = "Vector"
  )

  # Plain tables (everything else)
  exclude <- grepl("^rtree_|^gpkg_|^sqlite_|_geom$", all_tables) |
    all_tables %in% rast_lyrs$layer_name
  tbl_lyrs <- tibble::tibble(
    layer_name = all_tables[!exclude],
    data_type = "Table"
  ) |>
    dplyr::filter(
      !layer_name %in% vect_lyrs$layer_name,
      !layer_name %in% rast_lyrs$layer_name,
    )

  out <- dplyr::bind_rows(rast_lyrs, vect_lyrs, tbl_lyrs) |>
    dplyr::mutate(
      data_group = dplyr::case_when(
        grepl("^dem_", layer_name) ~ "hydro",
        grepl("iFLS|HAiFLS|iFLO|HAiFLO", layer_name) ~ "iDW",
        layer_name %in%
          c(
            "site_id_col",
            "DEM_Extent",
            "loi_meta",
            "unnest_catchment"
          ) ~ "meta",
        layer_name %in%
          c("snapped_points", "original_points") ~ "sample_points",
        layer_name %in%
          c(
            "Elevation",
            "stream_links",
            "stream_lines",
            "stream_points",
            "Subbasins_poly",
            "Catchment_poly",
            "stream_points_attr",
            "stream_links_attr"
          ) ~ "hydro",
        layer_name %in% c("ds_flowpaths", "us_flowpaths") ~ "flow_path",
        layer_name %in%
          c("funcon_pwise_dist", "fcon_pwise_dist") ~ "pwise_dist",
        TRUE ~ "loi"
      )
    ) |>
    dplyr::arrange(data_type, data_group)

  out <- as.data.frame(out)

  return(out)
}


#' Extract layers stored in an ihydro GeoPackage
#'
#' Reads a specified layer from the ihydro GeoPackage and returns it as an R object (see return value for details).
#'
#' @param ihydro An object of class `ihydro`, as returned by [as.ihydro()],
#'   [process_flowdir()], or any other ihydro function.
#' @param layer Character. Name of the layer to read from the GeoPackage. Use [ihydro_layers()] to see what layers are available.
#' @param collect Logical. Whether to use [dplyr::collect()] to retieve data into a local tibble. Only relevant when `layer` is a table.
#'
#' @return The requested layer, read from the GeoPackage. The type of the returned object depends on the layer:
#' - Raster layers are returned as `SpatRaster` objects.
#' - Vector layers are returned as `SpatVector` objects.
#' - Plain tables are returned as `tibble` data frames.
#'
#'
#' @seealso [ihydro_layers()] to create an `ihydro` object.
#'
#' @examples
#' \dontrun{
#' hydro <- as.ihydro("my_watershed.gpkg")
#'
#' # List everything
#' ihydro_layers(hydro)
#'
#' # Return the stream lines
#' read_ihydro(hydro, "stream_lines")
#' }
#'
#' @export
read_ihydro <- function(ihydro, layer, collect = TRUE) {
  stopifnot(is.character(layer), length(layer) == 1L)
  .check_ihydro(ihydro)
  .check_ihydro_layer(ihydro, layer)

  layer_list <- ihydro_layers(ihydro)

  layer_sel <- dplyr::filter(layer_list, layer_name == layer)

  if (nrow(layer_sel) == 0) {
    if (layer == "site_id_col") {
      return(tibble::tibble(site_id_col = NA_character_)[F, ])
    }
  }
  out <- switch(
    layer_sel$data_type,
    "Raster" = .ihydro_extract_rast(ihydro, layer = layer),
    "Vector" = sf::read_sf(ihydro$outfile, layer = layer),
    "Table" = .ihydro_extract_tbl(ihydro, layer, collect),
    cli::cli_abort("Unknown data type for layer {.val {layer}}")
  )

  if (inherits(out, c("data.table", "sf"))) {
    out <- dplyr::mutate(
      out,
      dplyr::across(
        tidyselect::any_of("link_id"),
        ~ as.character(.x)
      )
    )
  }

  return(out)
}


#' Create an ihydro object from a GeoPackage file
#'
#' Loads an existing ihydro GeoPackage (`.gpkg`) file and wraps it in an
#' `ihydro` S3 object. The resulting object stores the file path and can be
#' passed to any other function in the package (e.g., [generate_vectors()],
#' [trace_flowpaths()], [process_loi()]).
#'
#' `as_ihydro()` is provided as a snake_case alias for convenience.
#'
#' @param file Character. File path to a GeoPackage file (`.gpkg`). The file
#'   must already exist on disk.
#'
#' @return An object of class `ihydro`, which is a list with element `outfile`
#'   containing the normalised path to the GeoPackage.
#'
#' @seealso [ihydro_layers()] to inspect the contents of an `ihydro` object,
#'   [process_flowdir()] or [process_hydrology()] to create a new GeoPackage
#'   from a DEM.
#'
#' @examples
#' \dontrun{
#' # Open an existing ihydro GeoPackage
#' my_hydro <- as.ihydro("path/to/hydrology.gpkg")
#'
#' # Inspect the layers stored inside
#' ihydro_layers(my_hydro)
#'
#' # The snake_case alias works identically
#' my_hydro <- as_ihydro("path/to/hydrology.gpkg")
#' }
#'
#' @export
as.ihydro <- function(file) {
  file <- normalizePath(file, mustWork = FALSE)
  if (!file.exists(file)) {
    cli::cli_abort("{.arg file} does not exist: {.path {file}}")
  }
  if (!grepl("\\.gpkg$", file)) {
    cli::cli_abort("{.arg file} must end in {.val .gpkg}.")
  }
  structure(list(outfile = file), class = "ihydro")
}

#' @export
#' @rdname as.ihydro
as_ihydro <- as.ihydro
