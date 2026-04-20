ihydro: Integrated hydrology tools for environmental science
================

- [ihydro
  <img src="man/figures/logo.png" align="right" width="280"/>](#ihydro-)
  - [1 Introduction](#1-introduction)
  - [2 System setup and installation](#2-system-setup-and-installation)
  - [3 Prepare DEM and Sampling Points for
    analysis](#3-prepare-dem-and-sampling-points-for-analysis)
    - [3.1 Generate toy terrain dataset and sampling
      points](#31-generate-toy-terrain-dataset-and-sampling-points)
    - [3.2 Generate flow direction/accumulation geospatial analysis
      products with
      `process_flowdir()`](#32-generate-flow-directionaccumulation-geospatial-analysis-products-with-process_flowdir)
    - [3.3 Generate vector geospatial analysis products with
      `generate_vectors()`](#33-generate-vector-geospatial-analysis-products-with-generate_vectors)
    - [3.4 Create lookup tables of flow-directions with
      `trace_flowpaths()`](#34-create-lookup-tables-of-flow-directions-with-trace_flowpaths)
    - [3.5 Generate complete upstream catchment areas with
      `get_catchment()`](#35-generate-complete-upstream-catchment-areas-with-get_catchment)
    - [3.6 Generate entire hydrology workflow in one
      step.](#36-generate-entire-hydrology-workflow-in-one-step)
  - [4 Examining Pairwise Instream
    Distances](#4-examining-pairwise-instream-distances)
  - [5 Process layers of interest](#5-process-layers-of-interest)
    - [5.1 Separate from hydrology
      products](#51-separate-from-hydrology-products)
    - [5.2 Together with hydrology
      products](#52-together-with-hydrology-products)
  - [6 Calculate (weighted) spatial
    summaries:](#6-calculate-weighted-spatial-summaries)
    - [6.1 At specific sampling points](#61-at-specific-sampling-points)
    - [6.2 At all sampled points](#62-at-all-sampled-points)
    - [6.3 Across entire study area
      `fasttrib_points()`](#63-across-entire-study-area-fasttrib_points)
  - [7 Example Modelling](#7-example-modelling)
    - [7.1 Train a Random Forest Model](#71-train-a-random-forest-model)
  - [8 Future Plans](#8-future-plans)
  - [9 References](#9-references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# ihydro <img src="man/figures/logo.png" align="right" width="280"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

## 1 Introduction

Aquatic environmental scientists are often interested in relating
landscape features to observed responses (e.g., fish size, water
temperature, invertebrate community composition, etc.) in streams,
rivers and lakes. The computational workflow for conducting these
investigations is complex. Simply describing how water flows and
accumulates across the landscape can be a challenge itself, but for
aquatic scientists it is only the first step. The stream network must
then be extracted from the landscape, and reaches (a.k.a. segments;
i.e., stretches of river between two confluences) identified and given
unique identifiers. These reaches must then be attributed to be
informative (e.g., slope, stream order, upstream channel length, etc.);
and upstream-downstream connectivity between reaches established.

Typically, observed responses are measured at discrete sampling
locations along the stream network. If the location information is to be
preserved (i.e., if samples are upstream and downstream of a particular
effluent outflow), they must be incorporated into the network. This is
done by splitting stream lines and catchments at these points.

Once that is done, landscape features of interest (e.g., land-cover,
geology, climate, etc.) must be related to the reach (or sampling
points). This alone can be complex as the spatial configuration of these
factors relative to flow direction and accumulation can be important
(Peterson ***et al.***, 2011). The ***ihydro*** package uses
***hydroweight*** (Kielstra ***et al.*** 2021) to calculate these
attributes.

The complexity of this workflow can be a rate limiting step in the
scope, content, quality, and applicability of investigations by aquatic
environmental scientists. The ***ihydro*** package offers tools and
workflows to simplify these complex steps. It is capable of handling all
above computation steps, leaving researchers the task of identifying
layers of interest and modeling with (potentially) large numbers of
predictors.

The ***ihydro*** package also implements descriptions of spatial
autocorrelation among sites. Due to the linear nature of flow in streams
and rivers, autocorrelation tends to be asymmetric, with sites generally
more similar when they are flow connected, than not flow connected.
***ihydro*** produces tables that can be transformed into an asymmetric
matrices that describes the relationships between sampling points based
on instream distances and/or the proportions shared catchments.
Proportion of shared upstream catchment (rather than in-stream distance)
are a more relevant measure of spatial autocorrelation in streams
because it accounts for differences in catchment areas between points,
as well as spatial proximity. For example, if water chemistry samples
are taken from a large 6th order stream, and a upstream small 1st order
tributary we would expect the small tributary to have only a small
impact on the larger stream (despite its close physical proximity).
Hence, autocorrelation should be low because the tributary does not
contribute much flow to the larger stream. Using in-stream distances
alone may misrepresent this pattern.

***ihydro*** stores its geospatial products in a geopackage file for
ease of retrieval and extraction in external software (i.e. QGIS). Many
***ihydro*** functions can be run in parallel for increased speed (if
enough memory is available). The functions are also quick at removing
internal intermediate files to keep hard drives from filling up too
fast.

[Back to top](#1-introduction)

## 2 System setup and installation

*WhiteboxTools* and *whitebox* are required for ***ihydro***. See
[whiteboxR](https://github.com/giswqs/whiteboxR) or below for
installation.

``` r
## Follow instructions for whitebox installation accordingly
## devtools::install_github("giswqs/whiteboxR") # For development version
## whitebox is now available on CRAN
# install.packages("whitebox")

library(whitebox)

if (F) {
  install_whitebox()
  # Possible warning message:
  # ------------------------------------------------------------------------
  # Could not find WhiteboxTools!
  # ------------------------------------------------------------------------
  #
  # Your next step is to download and install the WhiteboxTools binary:
  #     > whitebox::install_whitebox()
  #
  # If you have WhiteboxTools installed already run `wbt_init(exe_path=...)`':
  #    > wbt_init(exe_path='/home/user/path/to/whitebox_tools')
  #
  # For whitebox package documentation, ask for help:
  #    > ??whitebox
  #
  # For more information visit https://giswqs.github.io/whiteboxR/
  #
  # ------------------------------------------------------------------------
}
```

[Back to top](#1-introduction)

## 3 Prepare DEM and Sampling Points for analysis

### 3.1 Generate toy terrain dataset and sampling points

Begin by bringing in the digital elevation model.

``` r
## Load libraries
# remotes::install_github("p-schaefer/ihydro)
library(ihydro)
library(tmap)
library(furrr)
library(whitebox)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)

tmap_options(component.autoscale = FALSE)

# Many function in 'ihydro' can be run in parallel internally.

# Parallelization is done through the future package, so all parallel backends
# should be supported i.e.,:
plan(multisession(workers = 4))

## Generate ex_loc as a temporary directory
ex_loc <- tempdir()

## Use example data provided by ihydro
toy_dem <- ex_data("elev_ned_30m.tif")
sites <- ex_data("sites_nc.shp")

## copy to ex_loc if needed by your workflow
terra::writeRaster(toy_dem, file.path(ex_loc, "toy_dem.tif"), overwrite = TRUE)
terra::writeVector(sites, file.path(ex_loc, "sites.shp"), overwrite = TRUE)

plot(toy_dem, main = "Elevation")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

[Back to top](#1-introduction)

### 3.2 Generate flow direction/accumulation geospatial analysis products with `process_flowdir()`

The DEM must be processed in a way to remove depressions. Whitebox
offers methods for breaching and filling DEMs to remove depressions.
These may be applied before had, or through the `depression_corr`
argument in `process_flowdir()`.

Another factor to consider at this step is whether to burn a stream
vector into the DEM; ***ihydro*** implements a simplified stream burning
method that lowers the evation along the DEM by `burn_depth` meters. See
[here](https://proceedings.esri.com/library/userconf/proc99/proceed/papers/pap802/p802.htm#Trois)
for more detailed approaches.

``` r
# Outputs of all functions are always saved to output_filename by default,
# but can be included in function return with return_products=T (note this can
# be very slow for large regions)
output_filename_hydro <- file.path(ex_loc, "Processed_Hydrology.gpkg")

# Generates d8 flow direction and accumulation, extracts streams at a specified
# flow accumulation threshold
hydro_out <- process_flowdir(
  dem = toy_dem,
  burn_streams = ex_data("streams.shp"),
  burn_depth = 5,
  burn_slope_dist = 500,
  burn_slope_depth = 5,
  min_length = 3,
  depression_corr = "breach",
  threshold = 100L,
  return_products = T,
  output_filename = output_filename_hydro,
  temp_dir = NULL,
  verbose = T
)

# List files present in gpkg
ihydro_layers(hydro_out)
#>         layer_name data_type data_group
#> 1        dem_final    Raster      hydro
#> 2           dem_d8    Raster      hydro
#> 3     dem_accum_d8    Raster      hydro
#> 4 dem_accum_d8_sca    Raster      hydro
#> 5   dem_streams_d8    Raster      hydro
#> 6        Elevation    Raster        loi
#> 7       DEM_Extent    Vector       meta

# hydro_out$outfile # -> This is the full file path of the resulting geopackage (gpkg) file

# remaining outputs only present when `return_products` == T
# hydro_out$dem_final.tif # -> final dem after stream burning and depression correction
# hydro_out$dem_d8.tif # -> d8 flow direction
# hydro_out$dem_accum_d8.tif # -> d8 flow accumulation (cells)
# hydro_out$dem_accum_d8_sca.tif # -> d8 flow accumulation (specific catchment areas)
# hydro_out$dem_streams_d8.tif # -> extracted streams at specified `threshold`

# if `return_products` == F, all produces are only available in the gpkg file.
# terra and sf allow access to files directly in the gpkg file, whitebox
# requires them to be extracted to a folder

flow_accum <- read_ihydro(hydro_out, "dem_accum_d8")

plot(log10(flow_accum))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="" width="100%" />

To create an ihydro object from an existing geopackage (i.e., created
from a previous R run):

``` r
hydro_out_new <- as_ihydro(output_filename_hydro)
ihydro_layers(hydro_out_new)
#>         layer_name data_type data_group
#> 1        dem_final    Raster      hydro
#> 2           dem_d8    Raster      hydro
#> 3     dem_accum_d8    Raster      hydro
#> 4 dem_accum_d8_sca    Raster      hydro
#> 5   dem_streams_d8    Raster      hydro
#> 6        Elevation    Raster        loi
#> 7       DEM_Extent    Vector       meta
```

[Back to top](#1-introduction)

### 3.3 Generate vector geospatial analysis products with `generate_vectors()`

This function combines `generate_subbasins()` and `attrib_streamline()`
to produce subcatchment polygons, stream lines, stream links (downstream
most points along each subcatchment), stream points (stream line
represented as individual points).

Typically, observed responses are measured at discrete sampling
locations along the stream network. If the location information is to be
preserved (i.e., if samples are upstream and downstream of a particular
effluent outflow), the sampling points must be incorporated into the
network. This is done by splitting stream lines and subcatchments at
these points.

``` r
hydro_out <- generate_vectors(
  input = hydro_out,
  points = file.path(ex_loc, "sites.shp"), # These are optional sampling locations
  site_id_col = "site_id", # Column name in points layer that corresponds to
  #                      # unique IDs that will be available in data products
  snap_distance = 250L, # points that are more than 100m from closest stream are excluded
  break_on_noSnap = F, # default is to stop when any points don't snap, this will ignore that
  return_products = T,
  temp_dir = NULL,
  verbose = T
)

ihydro_layers(hydro_out)
#>            layer_name data_type    data_group
#> 1           dem_final    Raster         hydro
#> 2              dem_d8    Raster         hydro
#> 3        dem_accum_d8    Raster         hydro
#> 4    dem_accum_d8_sca    Raster         hydro
#> 5      dem_streams_d8    Raster         hydro
#> 6           Elevation    Raster           loi
#> 7   stream_links_attr     Table         hydro
#> 8  stream_points_attr     Table         hydro
#> 9         site_id_col     Table          meta
#> 10     Subbasins_poly    Vector         hydro
#> 11       stream_lines    Vector         hydro
#> 12       stream_links    Vector         hydro
#> 13      stream_points    Vector         hydro
#> 14         DEM_Extent    Vector          meta
#> 15    original_points    Vector sample_points
#> 16     snapped_points    Vector sample_points

# Several important columns are used throughout the vector layers:
# `link_id` - identifies reaches/segments (stream length between two confluences)
# `trib_id` - identifies streams/tributaries, with the shortest channel getting
#           # a new trib_id at a confluence, and longest channel retaining the original ID
# `uslink_id` and `dslink_id`  - columns identify upstream and downstream links
# `ustrib_id` and `dstrib_id`  - columns identify upstream and downstream tributaries (longest
#  continuous flow-paths)

# If 'points' are provided, this function modifies the vector outputs by inserting
# links, and splitting lines/subbasins at `points`. Inserted points are given "link_id'
# values that increment with decimal places from downstream to upstream directions
# (i.e., 10.0 remains the pour point for the segment, and 10.1, 10.2,... identify
# sample points in an upstream direction).

# New layers added by `generate_vectors()`
# hydro_out$subbasins      # -> polygon subbasins attributed with `link_id` and reach
#                          #  contributing area (in m^2)
# hydro_out$stream_lines   # -> line vectors attributed with `link_id`
# hydro_out$points         # -> point vectors along lines identifying 'nodes' (confluence
#                          #  points), vs 'links' segments joining 'nodes', and also
#                          #  attributed with `link_id`, `trib_id`, upstream
#                          #  and downstream link and trib IDS and a number of extra attributes
# hydro_out$links          # -> point vector representing pour-points for subbasins,
#                          #  attributed with `link_id` and extra attributes
# hydro_out$snapped_points # -> provided points snapped to nearest segment. Any points
#                          #  beyond snap distance are removed, with a warning if
#                          #  break_on_noSnap == F.

tm_shape(hydro_out$subbasins) +
  tm_polygons(
    fill = "link_id",
    fill.scale = tm_scale(values = "viridis"),
    fill_alpha = 0.2,
    fill.legend = tm_legend_hide()
  ) +
  tm_shape(hydro_out$stream_lines) +
  tm_lines(
    col = "blue",
    col_alpha = 0.5,
    col.legend = tm_legend_hide(),
    lwd = 3
  ) +
  tm_shape(hydro_out$links) +
  tm_dots(
    fill = "yellow",
    fill.legend = tm_legend_hide(),
    size = 0.15,
    col = "black",
    col_alpha = 1
  )
```

<img src="man/figures/README-unnamed-chunk-6-1.png" alt="" width="100%" />

``` r
# Optional step, inserts sampling points into stream vectors, splitting subbasins
# and lines at sampling points, additional links inserted at sampling points as well

tm_shape(rast(hydro_out$dem_final.tif)) +
  tm_raster(
    col.scale = tm_scale(values = "cividis"),
    col.legend = tm_legend_hide()
  ) +
  tm_shape(hydro_out$stream_lines) +
  tm_lines(
    col = "red",
    col_alpha = 0.5,
    col.legend = tm_legend_hide(),
    lwd = 2
  ) +
  tm_shape(hydro_out$snapped_points %>% mutate(site_id = as.numeric(site_id))) +
  tm_dots(
    fill = "darkgray",
    fill.legend = tm_legend_hide(),
    size = 0.35,
    col = "black",
    col_alpha = 1,
    border.lwd = 1
  ) +
  tm_shape(read_sf(file.path(ex_loc, "sites.shp"))) +
  tm_dots(
    fill = "black",
    fill.legend = tm_legend_hide(),
    size = 0.35,
    col = "black",
    col_alpha = 1,
    border.lwd = 1
  )
```

<img src="man/figures/README-unnamed-chunk-7-1.png" alt="" width="100%" />

[Back to top](#1-introduction)

### 3.4 Create lookup tables of flow-directions with `trace_flowpaths()`

To more efficiently generate catchments, look-ups are created that
identify all upstream and downstream links originating from each link.

``` r
hydro_out <- trace_flowpaths(
  input = hydro_out,
  pwise_dist = T, # This will calculate all downstream flow-connected pairwise distances
  pwise_all_links = F, # This will calculate all flow un-connected pairwise distances (can be very time consuming)
  temp_dir = NULL,
  return_products = T,
  verbose = T
)

ihydro_layers(hydro_out)
#>            layer_name data_type    data_group
#> 1           dem_final    Raster         hydro
#> 2              dem_d8    Raster         hydro
#> 3        dem_accum_d8    Raster         hydro
#> 4    dem_accum_d8_sca    Raster         hydro
#> 5      dem_streams_d8    Raster         hydro
#> 6           Elevation    Raster           loi
#> 7        ds_flowpaths     Table     flow_path
#> 8        us_flowpaths     Table     flow_path
#> 9   stream_links_attr     Table         hydro
#> 10 stream_points_attr     Table         hydro
#> 11        site_id_col     Table          meta
#> 12    fcon_pwise_dist     Table    pwise_dist
#> 13  funcon_pwise_dist     Table    pwise_dist
#> 14     Subbasins_poly    Vector         hydro
#> 15       stream_lines    Vector         hydro
#> 16       stream_links    Vector         hydro
#> 17      stream_points    Vector         hydro
#> 18         DEM_Extent    Vector          meta
#> 19    original_points    Vector sample_points
#> 20     snapped_points    Vector sample_points

head(hydro_out$us_flowpaths)
#> # A tibble: 6 × 2
#>   pour_point_id origin_link_id
#>   <chr>         <chr>         
#> 1 1             1             
#> 2 10            10            
#> 3 100           100           
#> 4 1000          207           
#> 5 1000          1000          
#> 6 1000          233

head(hydro_out$ds_flowpaths %>% arrange(destination_link_id))
#> # A tibble: 6 × 2
#>   destination_link_id origin_link_id
#>   <chr>               <chr>         
#> 1 1                   1             
#> 2 10                  10            
#> 3 100                 100           
#> 4 1000                207           
#> 5 1000                1000          
#> 6 1000                233

# Below we demonstrate how to extract the upstream and downstream flowpaths for specific points

us_800 <- hydro_out$us_flowpaths %>%
  filter(pour_point_id == "800") # get all upstream link_ids from link_id 800
ds_101.1 <- hydro_out$ds_flowpaths %>%
  filter(origin_link_id == "101.1") # get all downstream link_ids from
#                                   # link_id 101.1 (this corresponds with site_id = 62)

lines_out <- hydro_out$stream_lines %>%
  filter(
    link_id %in%
      us_800$origin_link_id |
      link_id %in% ds_101.1$destination_link_id
  )
sub_out <- hydro_out$subbasins %>%
  filter(
    link_id %in%
      us_800$origin_link_id |
      link_id %in% ds_101.1$destination_link_id
  )

tm_shape(sub_out) +
  tm_polygons(
    fill = "white",
    fill_alpha = 0.2,
    fill.legend = tm_legend_hide()
  ) +
  tm_shape(lines_out) +
  tm_lines(
    col = "blue",
    col_alpha = 0.5,
    col.legend = tm_legend_hide(),
    lwd = 3
  ) +
  tm_shape(hydro_out$links %>% filter(link_id %in% c("800", "101.1"))) +
  tm_dots(
    fill.legend = tm_legend_hide(),
    size = 0.45,
    col = "black",
    col_alpha = 1,
    border.lwd = 1
  ) +
  tm_layout(frame = FALSE)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" alt="" width="100%" />

[Back to top](#1-introduction)

### 3.5 Generate complete upstream catchment areas with `get_catchment()`

Once lookups are established, catchments can easily be retrieved

``` r
subbasin_catchments <- get_catchment(
  # retrieve catchment for an arbitrary reach
  input = hydro_out,
  sample_points = NULL,
  link_id = "838"
)

point_catchments <- get_catchment(
  # retrieve catchments at specified sampling points
  input = hydro_out,
  sample_points = c("1", "25"),
  link_id = NULL
)

tm_shape(bind_rows(subbasin_catchments, point_catchments)) +
  tm_polygons(
    col = "white",
    fill_alpha = 0.2,
    fill.legend = tm_legend_hide(),
    lwd = 4
  ) +
  tm_shape(hydro_out$stream_lines) +
  tm_lines(
    col = "blue",
    col_alpha = 0.5,
    col.legend = tm_legend_hide(),
    lwd = 2
  )
```

<img src="man/figures/README-unnamed-chunk-9-1.png" alt="" width="100%" />

[Back to top](#1-introduction)

### 3.6 Generate entire hydrology workflow in one step.

The entire workflow above can be accomplished with a single function:

``` r
output_filename_hydro_sparse <- file.path(
  ex_loc,
  "Processed_Hydrology_sparse.gpkg"
)

# In this case we will use a higher stream initiation threshold to speed up calculation.

hydro_out_sparse <- process_hydrology(
  dem = toy_dem,
  output_filename = output_filename_hydro_sparse,
  # burn_streams=file.path(ex_loc,"streams.shp"),
  burn_depth = 5,
  depression_corr = "breach",
  threshold = 500L,
  points = hydro_out$snapped_points,
  site_id_col = "site_id",
  break_on_noSnap = F,
  pwise_dist = T,
  pwise_all_links = F,
  return_products = F,
  temp_dir = NULL,
  verbose = F
)

ihydro_layers(hydro_out_sparse)
#>            layer_name data_type    data_group
#> 1           dem_final    Raster         hydro
#> 2              dem_d8    Raster         hydro
#> 3        dem_accum_d8    Raster         hydro
#> 4    dem_accum_d8_sca    Raster         hydro
#> 5      dem_streams_d8    Raster         hydro
#> 6           Elevation    Raster           loi
#> 7        ds_flowpaths     Table     flow_path
#> 8        us_flowpaths     Table     flow_path
#> 9   stream_links_attr     Table         hydro
#> 10 stream_points_attr     Table         hydro
#> 11        site_id_col     Table          meta
#> 12    fcon_pwise_dist     Table    pwise_dist
#> 13  funcon_pwise_dist     Table    pwise_dist
#> 14     Subbasins_poly    Vector         hydro
#> 15       stream_lines    Vector         hydro
#> 16       stream_links    Vector         hydro
#> 17      stream_points    Vector         hydro
#> 18         DEM_Extent    Vector          meta
#> 19    original_points    Vector sample_points
#> 20     snapped_points    Vector sample_points

# Since we didn't return the products, we'll verify the outputs exist in the gpkg file

tm_shape(read_sf(hydro_out_sparse$outfile, "Subbasins_poly")) +
  tm_polygons(col = "white", fill_alpha = 0.2, fill.legend = tm_legend_hide()) +
  tm_shape(read_sf(hydro_out_sparse$outfile, "stream_lines")) +
  tm_lines(
    col = "blue",
    col_alpha = 0.3,
    col.legend = tm_legend_hide(),
    lwd = 2
  ) +
  tm_shape(read_sf(hydro_out_sparse$outfile, "stream_links")) +
  tm_dots(
    col.legend = tm_legend_hide(),
    size = 0.2,
    col = "black",
    col_alpha = 1,
    border.lwd = 1
  )
```

<img src="man/figures/README-unnamed-chunk-10-1.png" alt="" width="100%" />

[Back to top](#1-introduction)

## 4 Examining Pairwise Instream Distances

For more complete and thorough treatment on spatial autocorrelation in
stream systems, see [Zimmerman and Hoef
(2007)](https://www.fs.usda.gov/rm/boise/AWAE/projects/NationalStreamInternet/downloads/17ZimmermanVerHoef_TheTorgeg%20amForFluvialVariography.pdf).

In case we didn’t request pairwise distances in `trace_flowpaths()`, we
can calculate them separately with `generate_pdist()`. Below we
calculate a table that summaries relevant pairwise distance measures
between reaches.

``` r
# If pairwise distances were not requested with trace_flowpaths(), they can be added
# at any point after using generate_pdist()

hydro_out <- generate_pdist(
  input = hydro_out,
  return_products = T,
  pwise_all_links = T # This will calculate pairwise distances between all stream links
  #                 #  which can be very slow for dense stream networks
)

# New data available in `fcon_pwise_dist` and 'funcon_pwise_dist'
# *_pwise_dist # -> tables of downstream path lengths between each pair of link_ids,
#                   #    with flow-connected in-stream distances (directed_path_length),
#                   #    flow-unconnected in-stream distances (undirected_path_length),
#                   #    and proportions of shared catchments (prop_shared_catchment)
#                   #    and log-transformed catchment proportions (prop_shared_logcatchment)

ihydro_layers(hydro_out)
#>            layer_name data_type    data_group
#> 1           dem_final    Raster         hydro
#> 2              dem_d8    Raster         hydro
#> 3        dem_accum_d8    Raster         hydro
#> 4    dem_accum_d8_sca    Raster         hydro
#> 5      dem_streams_d8    Raster         hydro
#> 6           Elevation    Raster           loi
#> 7        ds_flowpaths     Table     flow_path
#> 8        us_flowpaths     Table     flow_path
#> 9   stream_links_attr     Table         hydro
#> 10 stream_points_attr     Table         hydro
#> 11        site_id_col     Table          meta
#> 12    fcon_pwise_dist     Table    pwise_dist
#> 13  funcon_pwise_dist     Table    pwise_dist
#> 14     Catchment_poly    Vector         hydro
#> 15     Subbasins_poly    Vector         hydro
#> 16       stream_lines    Vector         hydro
#> 17       stream_links    Vector         hydro
#> 18      stream_points    Vector         hydro
#> 19         DEM_Extent    Vector          meta
#> 20    original_points    Vector sample_points
#> 21     snapped_points    Vector sample_points

# funcon_pwise_dist is only available if pwise_all_links==T

hydro_out$pwise_dist <- bind_rows(
  read_ihydro(hydro_out, "fcon_pwise_dist") %>%
    mutate(dist_type = "Flow Connected"),
  read_ihydro(hydro_out, "funcon_pwise_dist") %>%
    mutate(dist_type = "Flow Unconnected")
)

head(hydro_out$pwise_dist)
#> # A tibble: 6 × 7
#>   origin destination directed_path_length prop_shared_catchment
#>   <chr>  <chr>                      <dbl>                 <dbl>
#> 1 1      1                           84.9                 1    
#> 2 10     10                         447.                  1    
#> 3 100    100                        392.                  1    
#> 4 100    1115                       639.                  0.546
#> 5 100    1105                      1009.                  0.324
#> 6 100    1129                      1154.                  0.177
#> # ℹ 3 more variables: prop_shared_logcatchment <dbl>,
#> #   undirected_path_length <dbl>, dist_type <chr>

p1 <- hydro_out$stream_lines %>%
  mutate(link_id = as.character(link_id)) %>%
  left_join(
    hydro_out$pwise_dist %>%
      filter(dist_type == "Flow Connected") %>%
      filter(origin == "200"),
    by = c("link_id" = "destination")
  ) %>%
  tm_shape() +
  tm_lines(
    col = "directed_path_length",
    col.scale = tm_scale_continuous(values = "viridis"),
    col_alpha = 1,
    lwd = 2
  ) +
  tm_shape(hydro_out$links %>% filter(link_id == "200")) +
  tm_dots(
    fill.legend = tm_legend_hide(),
    size = 0.45,
    col = "black",
    col_alpha = 1,
    border.lwd = 1
  ) +
  tm_title("Flow connected path log-length from 200") +
  tm_layout(legend.outside = TRUE)

p2 <- hydro_out$stream_lines %>%
  mutate(link_id = as.character(link_id)) %>%
  left_join(
    hydro_out$pwise_dist %>%
      filter(origin == "200") %>%
      filter(dist_type == "Flow Unconnected"),
    by = c("link_id" = "destination"),
    multiple = "all"
  ) %>%
  tm_shape() +
  tm_lines(
    col = "undirected_path_length",
    col_alpha = 1,
    lwd = 2,
    col.scale = tm_scale_continuous(values = "viridis")
  ) +
  tm_shape(hydro_out$links %>% filter(link_id == "200")) +
  tm_dots(
    fill.legend = tm_legend_hide(),
    size = 0.45,
    col = "black",
    col_alpha = 1,
    border.lwd = 1
  ) +
  tm_title("Flow Unconnected path log-length from 99") +
  tm_layout(legend.outside = TRUE)

tmap_arrange(p1, p2, ncol = 1)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" alt="" width="100%" />

Below. we take these distances and visualize them as heatmaps to
illustrate how pairwise distance matrices can be used to represent
spatial relationships among sites.

``` r
# Examine relationships among sampled sites.

# Get link_id for sampled points
sel_link_id <- hydro_out$snapped_points

# filter long table to selected sites, and convert to wide format

# This table describes the flow-connected in-stream distances at destination points (rows)
# from upstream origin points (columns)
dmat <- hydro_out$pwise_dist %>%
  filter(
    origin %in% sel_link_id$link_id & destination %in% sel_link_id$link_id
  ) %>%
  filter(origin != destination) %>%
  select(
    -prop_shared_catchment,
    -undirected_path_length,
    -prop_shared_logcatchment,
    -dist_type
  ) %>%
  rename(link_id = origin) %>%
  mutate(
    directed_path_length = ifelse(
      directed_path_length == 1 | is.na(directed_path_length),
      0,
      directed_path_length
    )
  ) %>%
  distinct() %>%
  filter(directed_path_length != 0) %>%
  pivot_wider(
    names_from = destination,
    values_from = directed_path_length,
    values_fill = 0
  ) %>%
  data.frame(check.names = F) %>%
  tibble::column_to_rownames("link_id") %>%
  log1p()

head(dmat)
#>         1077.1   1034.1    999.1    945.1    1048.1    997.1   1045.2   1045.1
#> 101.1  8.15971 8.684756 9.033572 9.366995  9.719346 0.000000 0.000000 0.000000
#> 1022.1 0.00000 8.267308 8.758113 9.177639  9.590057 7.776448 0.000000 0.000000
#> 1032.1 0.00000 0.000000 0.000000 0.000000  8.827211 0.000000 7.915997 8.240093
#> 1033.1 0.00000 0.000000 8.133272 8.813156  9.364137 0.000000 0.000000 0.000000
#> 1034.1 0.00000 0.000000 7.988078 8.742223  9.323869 0.000000 0.000000 0.000000
#> 1042.1 0.00000 9.575440 9.733542 9.912962 10.131851 9.464666 0.000000 0.000000
#>          1047.1   1069.1    996.1 1068.1    971.1 1033.1 1022.1 1051.1 1072.1
#> 101.1  0.000000 0.000000 0.000000 0.0000 0.000000      0      0      0      0
#> 1022.1 0.000000 0.000000 0.000000 0.0000 0.000000      0      0      0      0
#> 1032.1 0.000000 0.000000 0.000000 0.0000 0.000000      0      0      0      0
#> 1033.1 0.000000 0.000000 0.000000 0.0000 0.000000      0      0      0      0
#> 1034.1 0.000000 0.000000 0.000000 0.0000 0.000000      0      0      0      0
#> 1042.1 7.495116 7.954495 8.687005 8.8755 9.270075      0      0      0      0
#>        1065.1 1090.1 1088.1 924.1 785.1 700.1 868.1 904.1 936.1 889.1
#> 101.1       0      0      0     0     0     0     0     0     0     0
#> 1022.1      0      0      0     0     0     0     0     0     0     0
#> 1032.1      0      0      0     0     0     0     0     0     0     0
#> 1033.1      0      0      0     0     0     0     0     0     0     0
#> 1034.1      0      0      0     0     0     0     0     0     0     0
#> 1042.1      0      0      0     0     0     0     0     0     0     0

# This table describes the proportions of shared catchments at destination points (rows)
# from upstream origin points (columns)
dmat2 <- hydro_out$pwise_dist %>%
  filter(
    origin %in% sel_link_id$link_id & destination %in% sel_link_id$link_id
  ) %>%
  filter(origin != destination) %>%
  select(
    -directed_path_length,
    -undirected_path_length,
    -prop_shared_catchment,
    -dist_type
  ) %>%
  rename(link_id = origin) %>%
  mutate(
    prop_shared_logcatchment = ifelse(
      prop_shared_logcatchment == 1 | is.na(prop_shared_logcatchment),
      1,
      prop_shared_logcatchment
    )
  ) %>%
  distinct() %>%
  filter(prop_shared_logcatchment != 0) %>%
  pivot_wider(
    names_from = destination,
    values_from = prop_shared_logcatchment,
    values_fill = 0
  ) %>%
  data.frame(check.names = F) %>%
  tibble::column_to_rownames("link_id")

head(dmat2)
#>           1077.1    1034.1     999.1     945.1    1048.1     997.1    1045.2
#> 101.1  0.8858744 0.7676318 0.7585807 0.7486138 0.7316749 0.0000000 0.0000000
#> 1022.1 0.0000000 0.8714326 0.8611575 0.8498429 0.8306134 0.8827067 0.0000000
#> 1032.1 0.0000000 0.0000000 0.0000000 0.0000000 0.7380410 0.0000000 0.9014871
#> 1033.1 0.0000000 0.0000000 0.8148296 0.8041237 0.7859287 0.0000000 0.0000000
#> 1034.1 0.0000000 0.0000000 0.9882090 0.9752250 0.9531586 0.0000000 0.0000000
#> 1042.1 0.0000000 0.7807891 0.7715828 0.7614451 0.7442158 0.7908905 0.0000000
#>           1045.1    1047.1    1069.1     996.1    1068.1     971.1 1033.1
#> 101.1  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000      0
#> 1022.1 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000      0
#> 1032.1 0.8933912 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000      0
#> 1033.1 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000      0
#> 1034.1 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000      0
#> 1042.1 0.0000000 0.9531981 0.8823503 0.8387889 0.8220535 0.8050773      0
#>        1022.1 1051.1 1072.1 1065.1 1090.1 1088.1 924.1 785.1 700.1 868.1 904.1
#> 101.1       0      0      0      0      0      0     0     0     0     0     0
#> 1022.1      0      0      0      0      0      0     0     0     0     0     0
#> 1032.1      0      0      0      0      0      0     0     0     0     0     0
#> 1033.1      0      0      0      0      0      0     0     0     0     0     0
#> 1034.1      0      0      0      0      0      0     0     0     0     0     0
#> 1042.1      0      0      0      0      0      0     0     0     0     0     0
#>        936.1 889.1
#> 101.1      0     0
#> 1022.1     0     0
#> 1032.1     0     0
#> 1033.1     0     0
#> 1034.1     0     0
#> 1042.1     0     0

# Here we multiply the matrices (using the proportions of shared catchments as a rough weighting scheme)
# calculate manhattan distances and generate a heatmap
(dmat * dmat2) %>%
  dist("man") %>%
  as.matrix() %>%
  heatmap()
```

<img src="man/figures/README-unnamed-chunk-12-1.png" alt="" width="100%" />

Using these relationships, we can perform a clustering analysis to
identify groups of sites with potentially high spatial autocorrelation.
These groups could for instance be used for cross-validation purposes.

``` r
# Using the above approach, we create 5 groups of spatially proximate points
km <- (dmat * dmat2) %>%
  dist("man") %>%
  hclust() %>%
  cutree(k = 5)

gps <- tibble::enframe(km, "link_id", "group") %>%
  mutate(link_id = as.numeric(link_id))

point_groups <- hydro_out$snapped_points %>%
  left_join(gps) %>%
  filter(!is.na(group))

# These can be used for cross-validation purposes to see how well the models extrapolate outside of
# sampled areas

tm_shape(hydro_out$subbasins) +
  tm_polygons(
    fill = "white",
    fill_alpha = 0.2,
    fill.legend = tm_legend_hide()
  ) +
  tm_shape(hydro_out$stream_lines) +
  tm_lines(
    col = "blue",
    col_alpha = 0.3,
    col.legend = tm_legend_hide(),
    lwd = 2
  ) +
  tm_shape(point_groups) +
  tm_dots(
    fill = "group",
    fill.scale = tm_scale(values = "hcl.dark2"),
    size = 0.45
  ) +
  tm_layout(legend.outside = TRUE)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" alt="" width="100%" />

Finally, we’ll examine relationships among our response variable in
terms of in-stream distances, contrasting flow connected sites vs, flow
unconnected sites. Here we see greater differences between responses
with increasing in-stream distance. We expect greater similarity among
flow-connected than flow-unconnected sites, but don’t see it here. This
may be a product of this being an artificial data set.

``` r
# get response variables
response_table <- file.path(ex_loc, "sites.shp") %>%
  read_sf() %>%
  as_tibble() %>%
  select(site_id, value) %>%
  left_join(
    hydro_out$snapped_points %>% # join by snapped sites to get link_id values
      as_tibble() %>%
      select(site_id, link_id)
  ) %>%
  mutate(link_id = as.character(link_id))

# Combine pairwise data with values and examine spatial relationships.
dmat <- hydro_out$pwise_dist %>%
  filter(
    origin %in% sel_link_id$link_id & destination %in% sel_link_id$link_id
  ) %>%
  # Add response values for origin and destination points into the table
  left_join(
    response_table %>% rename(origin_value = value),
    by = c("origin" = "link_id")
  ) %>%
  left_join(
    response_table %>% rename(destination_value = value),
    by = c("destination" = "link_id")
  ) %>%
  mutate(value_diff = sqrt((origin_value - destination_value)^2)) %>% # Calculate the squared difference
  filter(origin != destination) %>%
  select(
    -dist_type,
    -site_id.x,
    -site_id.y,
    -origin_value,
    -destination_value
  ) %>%
  pivot_longer(
    c(
      directed_path_length,
      undirected_path_length,
      prop_shared_catchment,
      prop_shared_logcatchment
    ),
    names_to = "dist_type",
    values_to = "dist"
  ) %>%
  filter(!is.na(dist)) %>%
  mutate(
    `Distance Type` = case_when(
      dist_type == "directed_path_length" ~ "Flow Connected",
      dist_type == "undirected_path_length" ~ "Flow Unconnected",
      dist_type == "prop_shared_catchment" ~ "Shared Catchment",
      dist_type == "prop_shared_logcatchment" ~ "Shared log-Catchment",
    )
  ) %>%
  distinct()

dmat
#> # A tibble: 7,968 × 6
#>    origin destination value_diff dist_type                  dist `Distance Type`
#>    <chr>  <chr>            <dbl> <chr>                     <dbl> <chr>          
#>  1 101.1  1077.1               1 directed_path_length    3.50e+3 Flow Connected 
#>  2 101.1  1077.1               1 undirected_path_length  3.50e+3 Flow Unconnect…
#>  3 101.1  1077.1               1 prop_shared_catchment   1.84e-1 Shared Catchme…
#>  4 101.1  1077.1               1 prop_shared_logcatchme… 8.86e-1 Shared log-Cat…
#>  5 101.1  1034.1               2 directed_path_length    5.91e+3 Flow Connected 
#>  6 101.1  1034.1               2 undirected_path_length  5.91e+3 Flow Unconnect…
#>  7 101.1  1034.1               2 prop_shared_catchment   1.86e-2 Shared Catchme…
#>  8 101.1  1034.1               2 prop_shared_logcatchme… 7.68e-1 Shared log-Cat…
#>  9 101.1  999.1                2 directed_path_length    8.38e+3 Flow Connected 
#> 10 101.1  999.1                2 undirected_path_length  8.38e+3 Flow Unconnect…
#> # ℹ 7,958 more rows

require(ggplot2)
ggplot(
  dmat %>%
    filter(
      dist_type %in% c("directed_path_length", "undirected_path_length")
    ) %>%
    filter(dist > 0),
  aes(x = dist, y = value_diff, col = `Distance Type`)
) +
  geom_point(bins = 50) +
  geom_smooth(method = "gam", se = F) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylab("Pairwise Value Difference") +
  xlab("Pairwise Distance (m)") +
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~`Distance Type`)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" alt="" width="100%" />

Finally, we’ll do a similar comparison but using percent of shared
catchments. Here we expect pairwise differences to decrease as the
percent of shared catchments increases.

``` r
ggplot(
  dmat %>%
    filter(
      dist_type %in% c("prop_shared_catchment", "prop_shared_logcatchment")
    ) %>%
    filter(dist > 0),
  aes(x = dist, y = value_diff, colour = `Distance Type`)
) +
  geom_point() +
  geom_smooth(method = "gam", se = F) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylab("Pairwise Value Difference") +
  xlab("Percent of shared catchments (Flow Connected only)") +
  scale_x_continuous(labels = scales::percent) +
  facet_wrap(~`Distance Type`, scales = "free_x")
```

<img src="man/figures/README-unnamed-chunk-15-1.png" alt="" width="100%" />

[Back to top](#1-introduction)

## 5 Process layers of interest

### 5.1 Separate from hydrology products

Layers of interest (loi) represent the features on the landscape we are
interested in relating to the observed responses in the stream. These
can include: land-cover, climate, geology, soils, NDVI, slope, etc.

``` r
## Predictors from openSTARS
ex_data("landuse_r.tif") %>%
  setNames("LC") %>%
  writeRaster(file.path(ex_loc, "LC.tif"), overwrite = T)

writeVector(
  ex_data("geology.shp"),
  file.path(ex_loc, "geology.shp"),
  overwrite = T
)
writeVector(
  ex_data("pointsources.shp"),
  file.path(ex_loc, "pointsources.shp"),
  overwrite = T
)

landuse_r_path <- file.path(ex_loc, "LC.tif")
geology_path <- file.path(ex_loc, "geology.shp")
pointsources_path <- file.path(ex_loc, "pointsources.shp")

read_sf(pointsources_path) %>%
  mutate(pointsource = "pontsrc") %>%
  st_buffer(60) %>%
  write_sf(file.path(ex_loc, "pointsources.shp"), overwrite = T)

pointsources_path <- file.path(ex_loc, "pointsources.shp")


# Numeric Raster

wbt_slope(
  dem = file.path(ex_loc, "toy_dem.tif"),
  output = file.path(ex_loc, "slope.tif")
)

# Combine loi layers
output_filename_loi <- file.path(ex_loc, "Processed_loi.gpkg")

# This function standardizes numeric and categorical loi layers.

loi_combined <- process_loi(
  dem = toy_dem,
  num_inputs = list(
    # Can be given as a mixture of input types (file paths, or any sf or terra format)
    slope = file.path(ex_loc, "slope.tif")
  ),
  cat_inputs = list(
    # Can be given as a mixture of input types (file paths, or any sf or terra format)
    landcover = landuse_r_path,
    geology = geology_path,
    pointsources = pointsources_path
  ),
  variable_names = list(
    # any unlisted inputs will be used in their entirety
    geology = "GEO_NAME", # names listed here will subset those attributes or layers from the inputs
    pointsources = "pontsrc"
  ),
  output_filename = output_filename_loi,
  return_products = T,
  temp_dir = NULL,
  verbose = T
)

ihydro_layers(loi_combined)
#>         layer_name data_type data_group
#> 1    GEO_NAME.CZam    Raster        loi
#> 2    GEO_NAME.CZbg    Raster        loi
#> 3    GEO_NAME.CZfg    Raster        loi
#> 4     GEO_NAME.CZg    Raster        loi
#> 5    GEO_NAME.CZig    Raster        loi
#> 6    GEO_NAME.CZlg    Raster        loi
#> 7    GEO_NAME.CZve    Raster        loi
#> 8      GEO_NAME.Km    Raster        loi
#> 9             LC_1    Raster        loi
#> 10            LC_2    Raster        loi
#> 11            LC_3    Raster        loi
#> 12            LC_4    Raster        loi
#> 13            LC_5    Raster        loi
#> 14            LC_6    Raster        loi
#> 15            LC_7    Raster        loi
#> 16 pontsrc_pontsrc    Raster        loi
#> 17           slope    Raster        loi
#> 18        loi_meta     Table       meta
#> 19      DEM_Extent    Vector       meta
```

All categorical layers have been transformed to rasters with 1
indicating presence, and NA for absence. The layers have been rescaled
and projected to match the DEM

``` r
plot(loi_combined$cat_inputs, type = "classes", col = "darkgreen")
```

<img src="man/figures/README-unnamed-chunk-17-1.png" alt="" width="100%" />

All numeric layers have been rescaled and projected to match the DEM.

``` r
plot(loi_combined$num_inputs, type = "continuous")
```

<img src="man/figures/README-unnamed-chunk-18-1.png" alt="" width="100%" />

[Back to top](#1-introduction)

### 5.2 Together with hydrology products

Instead of processing the loi separately (as was done above), they can
instead be added to the completed workflow, and added to the existing
gpkg file for convenient file storage/organization.

``` r
# In this case, we will use our previously calculated loi results, but if `process_loi`
# is run with an input parameter specified, the loi rasters will be added to the
# output. This can make for convenient data storage.

hydro_out_sparse <- process_loi(
  input = hydro_out_sparse,
  num_inputs = list(
    # Can be given as a mixture of input types (file paths, or any sf or terra format)
    slope = file.path(ex_loc, "slope.tif")
  ),
  cat_inputs = list(
    # Can be given as a mixture of input types (file paths, or any sf or terra format)
    landcover = landuse_r_path,
    geology = geology_path,
    pointsources = pointsources_path
  ),
  variable_names = list(
    # any unlisted inputs will be used in their entirety
    geology = "GEO_NAME", # names listed here will subset those attributes or layers from the inputs
    pointsources = "pontsrc"
  ),
  return_products = F, # these layers can get large, and it is generally not advised to return them into R
  temp_dir = NULL,
  verbose = F
)

ihydro_layers(hydro_out_sparse)
#>            layer_name data_type    data_group
#> 1           dem_final    Raster         hydro
#> 2              dem_d8    Raster         hydro
#> 3        dem_accum_d8    Raster         hydro
#> 4    dem_accum_d8_sca    Raster         hydro
#> 5      dem_streams_d8    Raster         hydro
#> 6           Elevation    Raster           loi
#> 7       GEO_NAME.CZam    Raster           loi
#> 8       GEO_NAME.CZbg    Raster           loi
#> 9       GEO_NAME.CZfg    Raster           loi
#> 10       GEO_NAME.CZg    Raster           loi
#> 11      GEO_NAME.CZig    Raster           loi
#> 12      GEO_NAME.CZlg    Raster           loi
#> 13      GEO_NAME.CZve    Raster           loi
#> 14        GEO_NAME.Km    Raster           loi
#> 15               LC_1    Raster           loi
#> 16               LC_2    Raster           loi
#> 17               LC_3    Raster           loi
#> 18               LC_4    Raster           loi
#> 19               LC_5    Raster           loi
#> 20               LC_6    Raster           loi
#> 21               LC_7    Raster           loi
#> 22    pontsrc_pontsrc    Raster           loi
#> 23              slope    Raster           loi
#> 24       ds_flowpaths     Table     flow_path
#> 25       us_flowpaths     Table     flow_path
#> 26  stream_links_attr     Table         hydro
#> 27 stream_points_attr     Table         hydro
#> 28           loi_meta     Table          meta
#> 29        site_id_col     Table          meta
#> 30    fcon_pwise_dist     Table    pwise_dist
#> 31  funcon_pwise_dist     Table    pwise_dist
#> 32     Subbasins_poly    Vector         hydro
#> 33       stream_lines    Vector         hydro
#> 34       stream_links    Vector         hydro
#> 35      stream_points    Vector         hydro
#> 36         DEM_Extent    Vector          meta
#> 37    original_points    Vector sample_points
#> 38     snapped_points    Vector sample_points
```

[Back to top](#1-introduction)

## 6 Calculate (weighted) spatial summaries:

New research is showing that the presence of particular features on the
landscape is not always sufficient to predict the in-stream response to
those features, and the location of those features relative to the
locations of streams and areas of higher flow accumulation is very
important (Peterson ***et al.*** 2011).

***ihydro*** provides two functions for calculating weighted spatial
summaries: `attrib_points()` and `fasttrib_points()`. `attrib_points()`
uses the [hydroweight](https://github.com/GLFC-WET/hydroweight) package
to calculate weighted spatial summaries of supplied loi layers. It can
be used to examine the resulting distance-weighted rasters and
distance-weighted loi layers. However, `attrib_points()` is slow, so
`fasttrib_points()` is available for improved performance for larger
datasets.

### 6.1 At specific sampling points

The `attrib_points()` function uses the
[hydroweight](https://github.com/GLFC-WET/hydroweight) package and can
calculate and return weighted attributes for any points identified
either through supplied sampling points and/or arbitrary link_ids. The
returned products can be helpful in visualizing how the final weighted
attributes are derived.

``` r
attrib_points_time_small <- system.time(
  final_attributes_sub_slow <- attrib_points(
    input = hydro_out,
    out_filename = file.path(tempdir(), "attrib_points1.csv"),
    loi_file = output_filename_loi,
    loi_cols = NULL,
    sample_points = c("1", "25", "80"),
    link_id = NULL,
    clip_region = NULL,
    target_o_type = c("point"),
    weighting_scheme = c("lumped", "iFLO", "iFLS", "HAiFLO", "HAiFLS"),
    loi_numeric_stats = c(
      "distwtd_mean",
      "distwtd_sd",
      "mean",
      "sd",
      "median",
      "min",
      "max",
      "sum"
    ),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir = NULL,
    return_products = T,
    verbose = T
  )
)


final_attributes_sub_slow
#> # A tibble: 3 × 98
#>   link_id site_id status   products          slope_lumped_mean slope_lumped_sd
#>   <chr>   <chr>   <chr>    <list>                        <dbl>           <dbl>
#> 1 1099.1  1       Complete <named list [15]>              3.06           1.45 
#> 2 1048.1  25      Complete <named list [15]>              3.52           2.24 
#> 3 555.1   80      Complete <named list [15]>              2.38           0.962
#> # ℹ 92 more variables: slope_lumped_median <dbl>, slope_lumped_min <dbl>,
#> #   slope_lumped_max <dbl>, slope_lumped_sum <dbl>, slope_iFLO_mean <dbl>,
#> #   slope_iFLO_sd <dbl>, slope_iFLS_mean <dbl>, slope_iFLS_sd <dbl>,
#> #   slope_HAiFLO_mean <dbl>, slope_HAiFLO_sd <dbl>, slope_HAiFLS_mean <dbl>,
#> #   slope_HAiFLS_sd <dbl>, LC_1_lumped_prop <dbl>, LC_2_lumped_prop <dbl>,
#> #   LC_3_lumped_prop <dbl>, LC_4_lumped_prop <dbl>, LC_5_lumped_prop <dbl>,
#> #   LC_6_lumped_prop <dbl>, LC_7_lumped_prop <dbl>, …
```

We can access the weighting layers and weighted attribute layers (if
return_products==T) of the `attrib_points()` output.

``` r
plot(
  rast(
    list(
      rast(final_attributes_sub_slow$products[[2]]$iFLS) %>%
        setNames("iFLS Weighting"),
      log10(rast(final_attributes_sub_slow$products[[2]]$HAiFLO)) %>%
        setNames("log10-HAiFLO Weighting"),
      rast(final_attributes_sub_slow$products[[2]]$iFLS.loi_dist_rast_num) %>%
        setNames("iFLS Weighted Slope"),
      log10(rast(
        final_attributes_sub_slow$products[[2]]$HAiFLO.loi_dist_rast_cat
      )[[1]]) %>%
        setNames("log10-HAiFLO Weighted Landcover Class 1")
    )
  ),
  col = viridis::viridis(101),
  axes = F
)
```

<img src="man/figures/README-unnamed-chunk-21-1.png" alt="" width="100%" />

The `fasttrib_points()` function is faster for larger data sets, more
loi, and more sampling points. For very small datasets `attrib_points()`
may be faster.

``` r
fasttrib_points_time_small <- system.time(
  final_attributes_sub <- fasttrib_points(
    input = hydro_out,
    out_filename = file.path(tempdir(), "sample_points_wgtattr.csv"),
    loi_file = output_filename_loi, # specify loi file, if NULL, function will look for loi in 'input'
    loi_cols = NULL, # Specify loi columns to use, if NULL, all will be used
    iDW_file = NULL, # Leaving this as NULL will look for iDW in 'input' and calculate any not available
    store_iDW = T, # This will save the distance weights to the input or iDW_file if it is specified
    sample_points = c("1", "25", "80"),
    link_id = NULL,
    target_o_type = c("point"),
    weighting_scheme = c("lumped", "iFLS", "HAiFLS", "iFLO", "HAiFLO"), #
    loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    # temp_dir=NULL,
    verbose = T
  )
)

ihydro_layers(hydro_out)
#>                      layer_name data_type    data_group
#> 1                     dem_final    Raster         hydro
#> 2                        dem_d8    Raster         hydro
#> 3                  dem_accum_d8    Raster         hydro
#> 4              dem_accum_d8_sca    Raster         hydro
#> 5                dem_streams_d8    Raster         hydro
#> 6                          iFLS    Raster           iDW
#> 7                        HAiFLS    Raster           iDW
#> 8  HAiFLO_unn_group6JYyR87N7Gtc    Raster           iDW
#> 9  HAiFLO_unn_groupkkcqlRUD3oKk    Raster           iDW
#> 10   iFLO_unn_group6JYyR87N7Gtc    Raster           iDW
#> 11   iFLO_unn_groupkkcqlRUD3oKk    Raster           iDW
#> 12                    Elevation    Raster           loi
#> 13                 ds_flowpaths     Table     flow_path
#> 14                 us_flowpaths     Table     flow_path
#> 15            stream_links_attr     Table         hydro
#> 16           stream_points_attr     Table         hydro
#> 17                  site_id_col     Table          meta
#> 18                target_o_meta     Table          meta
#> 19              fcon_pwise_dist     Table    pwise_dist
#> 20            funcon_pwise_dist     Table    pwise_dist
#> 21               Catchment_poly    Vector         hydro
#> 22               Subbasins_poly    Vector         hydro
#> 23                 stream_lines    Vector         hydro
#> 24                 stream_links    Vector         hydro
#> 25                stream_points    Vector         hydro
#> 26                   DEM_Extent    Vector          meta
#> 27              original_points    Vector sample_points
#> 28               snapped_points    Vector sample_points

attr1 <- final_attributes_sub_slow %>%
  dplyr::select(any_of(colnames(final_attributes_sub)), -link_id, -site_id)

attr2 <- final_attributes_sub %>%
  dplyr::select(
    any_of(colnames(final_attributes_sub_slow)),
    -link_id,
    -site_id
  ) %>%
  dplyr::select(any_of(colnames(final_attributes_sub)))

stopifnot(all.equal(attr1, attr2))

final_attributes_sub
#> # A tibble: 3 × 97
#>   link_id site_id status   slope_lumped_mean slope_lumped_sd slope_lumped_min
#>   <chr>   <chr>   <chr>                <dbl>           <dbl>            <dbl>
#> 1 1099.1  1       Complete              3.06           1.45           0.213  
#> 2 1048.1  25      Complete              3.52           2.24           0.00216
#> 3 555.1   80      Complete              2.38           0.962          0.400  
#> # ℹ 91 more variables: slope_lumped_max <dbl>, slope_lumped_median <dbl>,
#> #   slope_lumped_sum <dbl>, slope_iFLS_mean <dbl>, slope_HAiFLS_mean <dbl>,
#> #   slope_iFLO_mean <dbl>, slope_HAiFLO_mean <dbl>, slope_HAiFLO_sd <dbl>,
#> #   slope_HAiFLS_sd <dbl>, slope_iFLO_sd <dbl>, slope_iFLS_sd <dbl>,
#> #   LC_1_lumped_prop <dbl>, LC_1_iFLS_prop <dbl>, LC_1_HAiFLS_prop <dbl>,
#> #   LC_1_iFLO_prop <dbl>, LC_1_HAiFLO_prop <dbl>, LC_2_lumped_prop <dbl>,
#> #   LC_2_iFLS_prop <dbl>, LC_2_HAiFLS_prop <dbl>, LC_2_iFLO_prop <dbl>, …
```

[Back to top](#1-introduction)

### 6.2 At all sampled points

``` r
final_attributes <- fasttrib_points(
  input = hydro_out,
  out_filename = file.path(tempdir(), "sample_points_wgtattr.csv"),
  loi_file = output_filename_loi,
  loi_cols = NULL,
  iDW_file = NULL,
  store_iDW = T,
  sample_points = hydro_out$snapped_points$site_id, # here we specify all sampling points
  link_id = NULL,
  target_o_type = c("point"),
  weighting_scheme = c("lumped", "iFLS", "iFLO", "HAiFLO", "HAiFLS"),
  loi_numeric_stats = c("mean", "sd", "min", "max"),
  inv_function = function(x) {
    (x * 0.001 + 1)^-1
  },
  temp_dir = NULL,
  verbose = F
)


final_attributes
#> # A tibble: 60 × 95
#>    link_id site_id status   slope_lumped_mean slope_lumped_sd slope_lumped_min
#>    <chr>   <chr>   <chr>                <dbl>           <dbl>            <dbl>
#>  1 1114.1  41      Complete              2.43            1.58          0.00253
#>  2 1099.1  1       Complete              3.06            1.45          0.213  
#>  3 101.1   62      Complete              3.09            1.81          0.0364 
#>  4 143.1   26      Complete              2.78            1.37          0.152  
#>  5 1083.1  32      Complete              2.60            1.49          0.0404 
#>  6 1065.1  4       Complete              3.40            1.74          0.105  
#>  7 1090.1  43      Complete              2.95            1.66          0.00253
#>  8 1077.1  63      Complete              3.65            2.13          0.0364 
#>  9 1047.1  7       Complete              2.77            1.76          0.0323 
#> 10 1042.1  5       Complete              2.56            1.61          0.0323 
#> # ℹ 50 more rows
#> # ℹ 89 more variables: slope_lumped_max <dbl>, slope_iFLS_mean <dbl>,
#> #   slope_HAiFLS_mean <dbl>, slope_iFLO_mean <dbl>, slope_HAiFLO_mean <dbl>,
#> #   slope_HAiFLO_sd <dbl>, slope_HAiFLS_sd <dbl>, slope_iFLO_sd <dbl>,
#> #   slope_iFLS_sd <dbl>, LC_1_lumped_prop <dbl>, LC_1_iFLS_prop <dbl>,
#> #   LC_1_HAiFLS_prop <dbl>, LC_1_iFLO_prop <dbl>, LC_1_HAiFLO_prop <dbl>,
#> #   LC_2_lumped_prop <dbl>, LC_2_iFLS_prop <dbl>, LC_2_HAiFLS_prop <dbl>, …
```

[Back to top](#1-introduction)

### 6.3 Across entire study area `fasttrib_points()`

In order to make predictions across the landscape, we will need to
calculate our attributes across the landscape as well. We leave
`sample_points` and `link_id` as NULL to predict across all reaches. At
this point, we may also consider switching our target_o parameter from
the sampling point (as was done above) to the entire reach by setting
`target_o_type`=“segment_whole”. This will calculate all target_o
weighting schemes to be the entire reach. This may be more conservative
for predicting beyond sampling points as it integrates landscape factors
across the length of the whole reach.

``` r
# This function is optional. It will save the hydroweight rasters in either a separate gpkp,
# or the input pkg.

# These can be large and take up a lot of room, so consider whether they are needed before storing.
# Having the raster available does decrease the computation time of fasttrib_points().

dw_time <- system.time(
  Processed_HW <- prep_weights(
    input = hydro_out,
    output_filename = file.path(tempfile(), "Processed_HW_out.gpkg"),
    weighting_scheme = c("iFLS", "iFLO", "HAiFLO", "HAiFLS"),
    target_o_type = c("point"),
    write_strategy = "batched"
  )
)


fasttrib_points_time_big <- system.time(
  final_attributes_all <- fasttrib_points(
    input = hydro_out,
    out_filename = "sample_points_wgtattr.csv",
    loi_file = output_filename_loi,
    loi_cols = NULL,
    iDW_file = file.path(tempfile(), "iDW_temp.gpkg"), # to use previously calculated weights file.path(tempfile(),"Processed_HW_out.gpkg")
    store_iDW = T,
    sample_points = NULL, # by specifying neither sample_points nor link_id
    link_id = NULL, #   we get all reaches
    target_o_type = c("segment_point"),
    weighting_scheme = c("lumped", "iFLS", "iFLO", "HAiFLO", "HAiFLS"),
    loi_numeric_stats = c("mean", "sd", "min", "max"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir = NULL,
    verbose = F
  )
)

final_attributes_all
#> # A tibble: 1,251 × 95
#>    link_id site_id status   slope_lumped_mean slope_lumped_sd slope_lumped_min
#>    <chr>   <chr>   <chr>                <dbl>           <dbl>            <dbl>
#>  1 7       <NA>    Complete              3.97            1.88           0.178 
#>  2 5       <NA>    Complete              4.81            2.39           0.385 
#>  3 1196    <NA>    Complete              3.94            1.92           0.0656
#>  4 1200    <NA>    Complete              4.11            2.39           0.119 
#>  5 1203    <NA>    Complete              3.59            2.21           0.0577
#>  6 1       <NA>    Complete              4.34            2.32           0.330 
#>  7 1195    <NA>    Complete              3.60            2.18           0.117 
#>  8 2       <NA>    Complete              3.55            1.70           0.317 
#>  9 1197    <NA>    Complete              4.29            2.59           0.162 
#> 10 32      <NA>    Complete              4.71            2.79           0.277 
#> # ℹ 1,241 more rows
#> # ℹ 89 more variables: slope_lumped_max <dbl>, slope_iFLS_mean <dbl>,
#> #   slope_HAiFLS_mean <dbl>, slope_iFLO_mean <dbl>, slope_HAiFLO_mean <dbl>,
#> #   slope_HAiFLO_sd <dbl>, slope_HAiFLS_sd <dbl>, slope_iFLO_sd <dbl>,
#> #   slope_iFLS_sd <dbl>, LC_1_lumped_prop <dbl>, LC_1_iFLS_prop <dbl>,
#> #   LC_1_HAiFLS_prop <dbl>, LC_1_iFLO_prop <dbl>, LC_1_HAiFLO_prop <dbl>,
#> #   LC_2_lumped_prop <dbl>, LC_2_iFLS_prop <dbl>, LC_2_HAiFLS_prop <dbl>, …

attrib_points_time_big <- system.time(
  final_attributes_slow <- attrib_points(
    input = hydro_out,
    out_filename = file.path(tempdir(), "attrib_points1.csv"),
    loi_file = output_filename_loi,
    loi_cols = NULL,
    sample_points = NULL,
    link_id = NULL,
    clip_region = NULL,
    target_o_type = c("segment_point"),
    weighting_scheme = c("lumped", "iFLO", "iFLS", "HAiFLO", "HAiFLS"),
    loi_numeric_stats = c(
      "distwtd_mean",
      "distwtd_sd",
      "mean",
      "sd",
      "min",
      "max"
    ),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir = NULL,
    return_products = F,
    verbose = F
  )
)

final_attributes_slow
#> # A tibble: 1,251 × 96
#>    link_id site_id status   products slope_lumped_mean slope_lumped_sd
#>    <chr>   <chr>   <chr>    <list>               <dbl>           <dbl>
#>  1 7       <NA>    Complete <NULL>                3.97            1.88
#>  2 5       <NA>    Complete <NULL>                4.81            2.39
#>  3 1196    <NA>    Complete <NULL>                3.94            1.92
#>  4 1200    <NA>    Complete <NULL>                4.11            2.39
#>  5 1203    <NA>    Complete <NULL>                3.59            2.21
#>  6 1       <NA>    Complete <NULL>                4.34            2.32
#>  7 1195    <NA>    Complete <NULL>                3.60            2.18
#>  8 2       <NA>    Complete <NULL>                3.55            1.70
#>  9 1197    <NA>    Complete <NULL>                4.29            2.59
#> 10 32      <NA>    Complete <NULL>                4.71            2.79
#> # ℹ 1,241 more rows
#> # ℹ 90 more variables: slope_lumped_min <dbl>, slope_lumped_max <dbl>,
#> #   slope_iFLO_mean <dbl>, slope_iFLO_sd <dbl>, slope_iFLS_mean <dbl>,
#> #   slope_iFLS_sd <dbl>, slope_HAiFLO_mean <dbl>, slope_HAiFLO_sd <dbl>,
#> #   slope_HAiFLS_mean <dbl>, slope_HAiFLS_sd <dbl>, LC_1_lumped_prop <dbl>,
#> #   LC_2_lumped_prop <dbl>, LC_3_lumped_prop <dbl>, LC_4_lumped_prop <dbl>,
#> #   LC_5_lumped_prop <dbl>, LC_6_lumped_prop <dbl>, LC_7_lumped_prop <dbl>, …

pmap(
  list(
    list(
      "attrib_points()",
      "fasttrib_points()",
      "attrib_points()",
      "fasttrib_points()"
    ),
    list(
      attrib_points_time_small,
      fasttrib_points_time_small,
      attrib_points_time_big,
      fasttrib_points_time_big
    ),
    list(
      final_attributes_sub_slow %>%
        select(
          -any_of("products"),
          -any_of("link_id"),
          -any_of("site_id"),
          -any_of("status")
        ),
      final_attributes_sub %>%
        select(
          -any_of("products"),
          -any_of("link_id"),
          -any_of("site_id"),
          -any_of("status")
        ),
      final_attributes_slow %>%
        select(
          -any_of("products"),
          -any_of("link_id"),
          -any_of("site_id"),
          -any_of("status")
        ),
      final_attributes_all %>%
        select(
          -any_of("products"),
          -any_of("link_id"),
          -any_of("site_id"),
          -any_of("status")
        )
    )
  ),
  function(.z, .x, .y) {
    paste0(
      .z,
      " took ",
      round(.x[[3]] / 60, 2),
      " min to calculate for ",
      nrow(.y),
      " reaches with ",
      ncol(.y),
      " attributes using ",
      nbrOfWorkers(),
      " cores."
    )
  }
)
#> [[1]]
#> [1] "attrib_points() took 0.39 min to calculate for 3 reaches with 94 attributes using 4 cores."
#> 
#> [[2]]
#> [1] "fasttrib_points() took 0.35 min to calculate for 3 reaches with 94 attributes using 4 cores."
#> 
#> [[3]]
#> [1] "attrib_points() took 37.02 min to calculate for 1251 reaches with 92 attributes using 4 cores."
#> 
#> [[4]]
#> [1] "fasttrib_points() took 9.98 min to calculate for 1251 reaches with 92 attributes using 4 cores."

paste0(
  round(dw_time[[3]] / 60, 2),
  " min to calculate distance weights for ",
  nrow(final_attributes_all),
  " reaches using ",
  nbrOfWorkers(),
  " cores."
)
#> [1] "4.9 min to calculate distance weights for 1251 reaches using 4 cores."
```

[Back to top](#1-introduction)

## 7 Example Modelling

### 7.1 Train a Random Forest Model

We’ll finish off with a brief demonstration of how to use the geospatial
products calculated above to build, and validate a predictive model.
First we’ll create a combined dataset:

``` r
# As we will ultimately be predicting to our sparse landscape, we will only keep
# autocorrelation variables that are relevant to the sparse landscape. This would
# not necessarily be required if predicting to the more dense landscape.

# get response variables
response_table <- file.path(ex_loc, "sites.shp") %>%
  read_sf() %>%
  as_tibble() %>%
  select(site_id, value) %>%
  left_join(
    hydro_out$links %>%
      as_tibble() %>%
      select(site_id, link_id) %>%
      mutate(site_id = as.character(site_id))
  ) %>%
  mutate(link_id = as.character(link_id)) %>%
  filter(!is.na(link_id))

head(response_table)
#> # A tibble: 6 × 3
#>   site_id value link_id
#>   <chr>   <dbl> <chr>  
#> 1 1           1 1099.1 
#> 2 4           1 1065.1 
#> 3 5           1 1042.1 
#> 4 7           1 1047.1 
#> 5 8           1 1069.1 
#> 6 11          2 996.1

# Columns for spatial cross-validation
# Here we will use a matrix of directed path lengths to perform the
# spatial cross-validation:
clust_data <- hydro_out$pwise_dist %>%
  filter(origin %in% response_table$link_id) %>%
  filter(origin != destination) %>%
  select(
    -prop_shared_catchment,
    -undirected_path_length,
    -prop_shared_logcatchment,
    -dist_type
  ) %>%
  rename(link_id = origin) %>%
  mutate(
    directed_path_length = ifelse(
      directed_path_length == 1 | is.na(directed_path_length),
      0,
      directed_path_length
    )
  ) %>%
  distinct() %>%
  filter(directed_path_length != 0) %>%
  pivot_wider(
    names_from = destination,
    values_from = directed_path_length,
    values_fill = 0
  ) %>%
  data.frame(check.names = F) %>%
  tibble::column_to_rownames("link_id") %>%
  log1p() %>%
  tibble::rownames_to_column("link_id") %>%
  as_tibble() %>%
  rename_with(.cols = c(everything(), -link_id), .fn = ~ paste0("CLUST_", .))

head(clust_data)
#> # A tibble: 6 × 238
#>   link_id CLUST_101 CLUST_1128 CLUST_1087 CLUST_1063 CLUST_1077.1 CLUST_1077
#>   <chr>       <dbl>      <dbl>      <dbl>      <dbl>        <dbl>      <dbl>
#> 1 101.1        7.71       7.76       7.98       8.02         8.16       8.28
#> 2 1022.1       0          0          0          0            0          0   
#> 3 1032.1       0          0          0          0            0          0   
#> 4 1033.1       0          0          0          0            0          0   
#> 5 1034.1       0          0          0          0            0          0   
#> 6 1042.1       0          0          0          0            0          0   
#> # ℹ 231 more variables: CLUST_1041 <dbl>, CLUST_1094 <dbl>, CLUST_1059 <dbl>,
#> #   CLUST_1119 <dbl>, CLUST_1034.1 <dbl>, CLUST_1034 <dbl>, CLUST_989 <dbl>,
#> #   CLUST_1073 <dbl>, CLUST_1009 <dbl>, CLUST_969 <dbl>, CLUST_1036 <dbl>,
#> #   CLUST_977 <dbl>, CLUST_999.1 <dbl>, CLUST_999 <dbl>, CLUST_1014 <dbl>,
#> #   CLUST_963 <dbl>, CLUST_944 <dbl>, CLUST_1025 <dbl>, CLUST_988 <dbl>,
#> #   CLUST_941 <dbl>, CLUST_968 <dbl>, CLUST_945.1 <dbl>, CLUST_945 <dbl>,
#> #   CLUST_974 <dbl>, CLUST_1140 <dbl>, CLUST_949 <dbl>, CLUST_992 <dbl>, …

# Combine the data into a single dataset.
comb_data <- response_table %>%
  left_join(
    final_attributes %>% mutate(site_id = as.character(site_id)) # replaces missing proportions with 0's
  ) %>%
  left_join(
    hydro_out$pwise_dist %>%
      filter(origin %in% response_table$link_id) %>%
      select(origin, destination, prop_shared_catchment) %>%
      filter(origin != destination) %>%
      rename(link_id = destination) %>%
      mutate(origin = paste0("Prop_catch_", origin)) %>%
      mutate(
        prop_shared_catchment = ifelse(
          prop_shared_catchment == 1 | is.na(prop_shared_catchment),
          0,
          prop_shared_catchment
        )
      ) %>%
      distinct() %>%
      filter(prop_shared_catchment != 0) %>%
      pivot_wider(
        names_from = origin,
        values_from = prop_shared_catchment,
        values_fill = 0
      )
  ) %>%
  left_join(clust_data) %>%
  filter(!is.na(link_id)) %>%
  mutate(across(starts_with("Prop_catch_"), ~ ifelse(is.na(.), 0, .)))

head(comb_data)
#> # A tibble: 6 × 393
#>   site_id value link_id status   slope_lumped_mean slope_lumped_sd
#>   <chr>   <dbl> <chr>   <chr>                <dbl>           <dbl>
#> 1 1           1 1099.1  Complete              3.06            1.45
#> 2 4           1 1065.1  Complete              3.40            1.74
#> 3 5           1 1042.1  Complete              2.56            1.61
#> 4 7           1 1047.1  Complete              2.77            1.76
#> 5 8           1 1069.1  Complete              3.25            1.90
#> 6 11          2 996.1   Complete              3.32            2.24
#> # ℹ 387 more variables: slope_lumped_min <dbl>, slope_lumped_max <dbl>,
#> #   slope_iFLS_mean <dbl>, slope_HAiFLS_mean <dbl>, slope_iFLO_mean <dbl>,
#> #   slope_HAiFLO_mean <dbl>, slope_HAiFLO_sd <dbl>, slope_HAiFLS_sd <dbl>,
#> #   slope_iFLO_sd <dbl>, slope_iFLS_sd <dbl>, LC_1_lumped_prop <dbl>,
#> #   LC_1_iFLS_prop <dbl>, LC_1_HAiFLS_prop <dbl>, LC_1_iFLO_prop <dbl>,
#> #   LC_1_HAiFLO_prop <dbl>, LC_2_lumped_prop <dbl>, LC_2_iFLS_prop <dbl>,
#> #   LC_2_HAiFLS_prop <dbl>, LC_2_iFLO_prop <dbl>, LC_2_HAiFLO_prop <dbl>, …
```

Then, we’ll follow the tidymodels workflow from here:
<https://www.tidymodels.org/>.

``` r
# install.packages("tidymodels")
require(tidymodels)
# install.packages("recipes")
require(recipes)
# install.packages("ranger")
require(ranger)

# Define Model - tune 3 main parameters
rf_mod <- rand_forest(
  trees = tune(),
  mtry = tune(),
  min_n = tune()
) %>%
  set_engine(
    "ranger",
    keep.inbag = TRUE,
    quantreg = TRUE,
    splitrule = "extratrees",
    num.random.splits = 25,
  ) %>%
  set_mode("regression")

# Setup recipes, define column roles, and preprocessing steps
recip <- recipe(
  x = comb_data %>%
    select(
      -starts_with("CLUST_"),
      -contains("site_id"),
      -contains("link_id"),
      -contains("pour_point_id"),
      -contains("status")
    )
) %>%
  update_role(c(everything()), new_role = "predictor") %>%
  update_role(value, new_role = "outcome") %>%
  step_zv(starts_with("Prop_catch_")) %>% # remove variables with zero variance
  step_nzv(
    contains("lumped"),
    contains("iFLS"),
    contains("iFLO")
  ) # remove variables with near zero variance

# Setup Cross-Validation Strategies
set.seed(1234)
cv_strats <- list(
  standard = vfold_cv(comb_data, v = 5), # standard random leave group_out cross-validation
  spatial = clustering_cv(
    comb_data, # using spatial information to leave out groups of nearby stations
    cluster_function = "hclust", # hclust did a decent job clustering sites previously
    vars = colnames(comb_data)[grepl("CLUST_", colnames(comb_data))],
    v = 5
  )
)

# Map hold-out data for each cross-validation fold
spatial_cv <- map2(
  cv_strats$spatial$splits,
  1:length(cv_strats$spatial$splits),
  function(x, y) {
    tm_shape(hydro_out$stream_lines) +
      tm_lines(
        col = "blue",
        col_alpha = 0.3,
        col.legend = tm_legend_hide(),
        lwd = 2
      ) +
      tm_shape(hydro_out$links[
        hydro_out$links$link_id %in% assessment(x)$link_id,
      ]) +
      tm_dots(
        col.legend = tm_legend_hide(),
        size = 0.3,
        col = "black",
        col_alpha = 1,
        border.lwd = 1
      ) +
      tm_title(paste0("Spatial CV", y)) +
      tm_layout()
  }
)

standard_cv <- map2(
  cv_strats$standard$splits,
  1:length(cv_strats$standard$splits),
  function(x, y) {
    tm_shape(hydro_out$stream_lines) +
      tm_lines(
        col = "blue",
        col_alpha = 0.3,
        col.legend = tm_legend_hide(),
        lwd = 2
      ) +
      tm_shape(hydro_out$links[
        hydro_out$links$link_id %in% assessment(x)$link_id,
      ]) +
      tm_dots(
        col.legend = tm_legend_hide(),
        size = 0.3,
        col = "black",
        col_alpha = 1,
        border.lwd = 1
      ) +
      tm_title(paste0("Standard CV", y)) +
      tm_layout()
  }
)

tmap_arrange(c(spatial_cv, standard_cv), ncol = 5)
```

<img src="man/figures/README-unnamed-chunk-26-1.png" alt="" width="100%" />

We see the spatial cross-validation hold-out sites cluster nicely
together, whereas the standard cross-validation are randomly distributed
across the region. Depending on the purpose of the model, we may favour
one type of cross-validation over another.

``` r
# Setup final workflow
wf <- workflow() %>%
  add_model(rf_mod) %>%
  add_recipe(recip)

# Run cross-validation strategies

# tune_grid() can only work in parallel with the doParallel package
plan(sequential)
#
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
set.seed(1234)

par_info <- hardhat::extract_parameter_set_dials(rf_mod) %>%
  finalize(bake(prep(recip, comb_data), comb_data)) %>%
  update(trees = dials::trees(range = c(100L, 1000L))) # Must increase the default number of trees to get quantile predictions

final_out <- map(
  cv_strats,
  ~ tune_grid(
    wf,
    resamples = .,
    metrics = metric_set(mae, rmse, rsq),
    param_info = par_info,
    grid = 200 # Choose 200 random hyper-parameter configurations
  )
)

stopCluster(cl)

# We expect the standard cross-validation to have higher accuracy than spatial
# cross-validation because the spatial autocorrelation variables should allow more accurate
# predictions from nearby sites. Conversely, the spatial cross-validation models
# show the best possible accuracy for predicting beyond the spatial extents of
# Where the data was collected.

map_dfr(
  final_out,
  show_best,
  n = 5,
  metric = "mae",
  .id = "Cross-validation strategy"
)
#> # A tibble: 10 × 10
#>    `Cross-validation strategy`  mtry trees min_n .metric .estimator  mean     n
#>    <chr>                       <int> <int> <int> <chr>   <chr>      <dbl> <int>
#>  1 standard                       18   145     5 mae     standard    2.23     5
#>  2 standard                       66   231     5 mae     standard    2.25     5
#>  3 standard                       91   271     8 mae     standard    2.25     5
#>  4 standard                       59   574     3 mae     standard    2.25     5
#>  5 standard                       83   339     3 mae     standard    2.25     5
#>  6 spatial                       107   136    22 mae     standard    2.51     5
#>  7 spatial                       112   656    11 mae     standard    2.52     5
#>  8 spatial                       104   724     5 mae     standard    2.53     5
#>  9 spatial                       106   932    16 mae     standard    2.53     5
#> 10 spatial                       113   760    17 mae     standard    2.53     5
#> # ℹ 2 more variables: std_err <dbl>, .config <chr>

map_dfr(
  final_out,
  show_best,
  n = 5,
  metric = "rsq",
  .id = "Cross-validation strategy"
)
#> # A tibble: 10 × 10
#>    `Cross-validation strategy`  mtry trees min_n .metric .estimator  mean     n
#>    <chr>                       <int> <int> <int> <chr>   <chr>      <dbl> <int>
#>  1 standard                      113   760    17 rsq     standard   0.327     5
#>  2 standard                      119   484    12 rsq     standard   0.326     5
#>  3 standard                       91   271     8 rsq     standard   0.326     5
#>  4 standard                      120   588    17 rsq     standard   0.324     5
#>  5 standard                       99   204    15 rsq     standard   0.322     5
#>  6 spatial                         1   208     7 rsq     standard   0.160     5
#>  7 spatial                       120   588    17 rsq     standard   0.142     5
#>  8 spatial                        90   167     3 rsq     standard   0.138     5
#>  9 spatial                       119   484    12 rsq     standard   0.134     5
#> 10 spatial                       121   692    23 rsq     standard   0.133     5
#> # ℹ 2 more variables: std_err <dbl>, .config <chr>
```

Next we build the final model. We will use the hyper-parameters that
yielded the highest mean absolute error (mae) prediction accuracy using
standard cross-validation. The spatial cross-validation may provide
better predictions as the distance from sampled areas increases, but
with so few samples, it is difficult to evaluate that aspect of the
model.

``` r
# We will use mae to select best metrics overall as it is better suited
# for model selection in this context than r^2

best_tunes <- map(final_out, select_best, metric = "mae")

best_tunes
#> $standard
#> # A tibble: 1 × 4
#>    mtry trees min_n .config          
#>   <int> <int> <int> <chr>            
#> 1    18   145     5 pre0_mod035_post0
#> 
#> $spatial
#> # A tibble: 1 × 4
#>    mtry trees min_n .config          
#>   <int> <int> <int> <chr>            
#> 1   107   136    22 pre0_mod183_post0

# Final ranger results
final_model <- finalize_workflow(wf, best_tunes$standard) %>%
  fit(comb_data) %>%
  extract_fit_engine()

# Overall observed vs expected leaves room for improvement, but is acceptable
# given the small data size and limited predictor consideration

tibble(
  Observed = comb_data$value,
  Predicted = predict(
    final_model,
    data = bake(prep(recip), new_data = comb_data)
  )$predictions
) %>%
  ggplot(aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(0, 10)) +
  theme_bw()
```

<img src="man/figures/README-unnamed-chunk-28-1.png" alt="" width="100%" />

Finally, we will use the model to predict across the whole landscape.
First we’ll assemble a table to predict from:

``` r
prediction_data <- final_attributes_all %>%
  left_join(
    # Add spatial autocorrelation variables
    hydro_out$pwise_dist %>%
      filter(origin %in% response_table$link_id) %>%
      select(origin, destination, prop_shared_catchment) %>%
      filter(origin != destination) %>%
      rename(link_id = destination) %>%
      mutate(origin = paste0("Prop_catch_", origin)) %>%
      mutate(
        prop_shared_catchment = ifelse(
          prop_shared_catchment == 1 | is.na(prop_shared_catchment),
          0,
          prop_shared_catchment
        )
      ) %>%
      distinct() %>%
      filter(prop_shared_catchment != 0) %>%
      pivot_wider(
        names_from = origin,
        values_from = prop_shared_catchment,
        values_fill = 0
      ),
    by = c("link_id")
  ) %>%
  mutate(across(starts_with("Prop_catch_"), ~ ifelse(is.na(.), 0, .))) %>%
  tidyr::drop_na(c(everything(), -site_id))

prediction_data
#> # A tibble: 1,251 × 155
#>    link_id site_id status   slope_lumped_mean slope_lumped_sd slope_lumped_min
#>    <chr>   <chr>   <chr>                <dbl>           <dbl>            <dbl>
#>  1 7       <NA>    Complete              3.97            1.88           0.178 
#>  2 5       <NA>    Complete              4.81            2.39           0.385 
#>  3 1196    <NA>    Complete              3.94            1.92           0.0656
#>  4 1200    <NA>    Complete              4.11            2.39           0.119 
#>  5 1203    <NA>    Complete              3.59            2.21           0.0577
#>  6 1       <NA>    Complete              4.34            2.32           0.330 
#>  7 1195    <NA>    Complete              3.60            2.18           0.117 
#>  8 2       <NA>    Complete              3.55            1.70           0.317 
#>  9 1197    <NA>    Complete              4.29            2.59           0.162 
#> 10 32      <NA>    Complete              4.71            2.79           0.277 
#> # ℹ 1,241 more rows
#> # ℹ 149 more variables: slope_lumped_max <dbl>, slope_iFLS_mean <dbl>,
#> #   slope_HAiFLS_mean <dbl>, slope_iFLO_mean <dbl>, slope_HAiFLO_mean <dbl>,
#> #   slope_HAiFLO_sd <dbl>, slope_HAiFLS_sd <dbl>, slope_iFLO_sd <dbl>,
#> #   slope_iFLS_sd <dbl>, LC_1_lumped_prop <dbl>, LC_1_iFLS_prop <dbl>,
#> #   LC_1_HAiFLS_prop <dbl>, LC_1_iFLO_prop <dbl>, LC_1_HAiFLO_prop <dbl>,
#> #   LC_2_lumped_prop <dbl>, LC_2_iFLS_prop <dbl>, LC_2_HAiFLS_prop <dbl>, …
```

And finally, we can predict across the landscape. We’ll use the range of
quantile predictions to estimate uncertainty in the model predictions.
The final map we produce here shows sampled points at their measured
value, stream segments are coloured according to their predicted values,
and the width of the stream segment corresponds to the uncertainty in
the prediction.

``` r
prediction_tbl <- tibble(
  link_id = prediction_data$link_id
) %>%
  bind_cols(
    predict(
      final_model,
      data = bake(prep(recip), new_data = prediction_data),
      type = "quantiles",
      quantiles = c(0.25, 0.5, 0.75)
    )$predictions %>%
      as_tibble() %>%
      setNames(c("p25", "p50", "p75"))
  ) %>%
  mutate(
    `Uncertainty` = p75 - p25,
    `Predicted` = p50
  )

prediction_tbl
#> # A tibble: 1,251 × 6
#>    link_id   p25   p50   p75 Uncertainty Predicted
#>    <chr>   <dbl> <dbl> <dbl>       <dbl>     <dbl>
#>  1 7           2     3     5           3         3
#>  2 5           3     4     5           2         4
#>  3 1196        1     3     4           3         3
#>  4 1200        1     3     4           3         3
#>  5 1203        2     4     7           5         4
#>  6 1           3     6     8           5         6
#>  7 1195        2     4     8           6         4
#>  8 2           1     4     8           7         4
#>  9 1197        4     7     8           4         7
#> 10 32          3     4     8           5         4
#> # ℹ 1,241 more rows

# Since we only have predictions for entire stream segments
# This will merge and breaks in stream lines
Streams <- read_sf(hydro_out$outfile, "stream_lines") %>%
  mutate(link_id = floor(link_id)) %>%
  group_by(link_id) %>%
  summarize(geom = sf::st_union(geom)) %>%
  ungroup() %>%
  mutate(link_id = as.character(link_id))

Streams <- Streams %>%
  left_join(prediction_tbl)

Points <- read_sf(hydro_out$outfile, "snapped_points") %>%
  mutate(link_id = as.character(link_id)) %>%
  left_join(response_table %>% select(-link_id)) %>%
  filter(!is.na(link_id)) %>%
  mutate(Observed = value)

tm_shape(Streams) +
  tm_lines(
    col = "Predicted",
    col.scale = tm_scale_intervals(
      values = "viridis",
      breaks = c(0, 2, 4, 6, 8, 10)
    ),
    lwd = "Uncertainty",
    lwd.scale = tm_scale_continuous(values.scale = 3) # numeric multiplier moved into values.scale
  ) +
  tm_shape(Points) +
  tm_dots(
    fill = "Observed",
    fill.scale = tm_scale_intervals(
      values = "viridis",
      breaks = c(0, 2, 4, 6, 8, 10)
    ),
    size = 0.75
  ) +
  tm_layout(legend.outside = TRUE)
```

<img src="man/figures/README-unnamed-chunk-30-1.png" alt="" width="100%" />

We can also summarize the percent of stream lengths that fall within
predicted value classes:

``` r
Streams %>%
  mutate(strm_len = as.numeric(st_length(geom))) %>%
  as_tibble() %>%
  drop_na() %>%
  mutate(
    `Prediction Groups` = cut(
      Predicted,
      breaks = c(0, 2, 4, 6, 8, 10),
      include.lowest = TRUE
    )
  ) %>%
  mutate(total_strm_len = sum(strm_len)) %>%
  group_by(`Prediction Groups`) %>%
  summarize(
    `% of Total Stream Length` = (sum(strm_len) / total_strm_len[[1]] * 100)
  ) %>%
  distinct() %>%
  ggplot(aes(
    x = "",
    y = `% of Total Stream Length`,
    fill = `Prediction Groups`
  )) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_viridis_d() +
  theme_bw()
```

<img src="man/figures/README-unnamed-chunk-31-1.png" alt="" width="100%" />

Next we will explore the final model in more detail. The
[DALEX](https://github.com/ModelOriented/DALEX) package is a great tool
for evaluating machine learning models. A supplemental ebook is
available with great detailed examples and explanations of all the
features available in this package ([Biecek and Burzykowski,
2020](https://ema.drwhy.ai/)).

``` r
# install.packages("DALEX)
library(DALEX)

# Prepare some datasets using the tidymodels recipe we specified for the model.
# We will prepare a full dataset to allow inferences about the full study region,
# as well as some spatial subsets to explore things more locally.

# Full dataset
final_data <- bake(prep(recip, comb_data), comb_data)

# One of the spatial validation folds in the upper stream
subset1_data <- assessment(cv_strats$spatial$splits[[1]])
subset1_data <- bake(prep(recip, comb_data), subset1_data)

# One of the spatial validation folds in the lower stream
subset2_data <- assessment(cv_strats$spatial$splits[[3]])
subset2_data <- bake(prep(recip, comb_data), subset2_data)

# Create some explainer objects
explainer_full <- explain(
  final_model,
  data = final_data %>% select(-value),
  y = final_data %>% select(value),
  label = "Full-Dataset"
)
#> Preparation of a new explainer is initiated
#>   -> model label       :  Full-Dataset 
#>   -> data              :  60  rows  123  cols 
#>   -> data              :  tibble converted into a data.frame 
#>   -> target variable   :  Argument 'y' was a data frame. Converted to a vector. (  WARNING  )
#>   -> target variable   :  60  values 
#>   -> predict function  :  yhat.ranger  will be used (  default  )
#>   -> predicted values  :  No value for predict function target column. (  default  )
#>   -> model_info        :  package ranger , ver. 0.18.0 , task regression (  default  ) 
#>   -> predicted values  :  numerical, min =  1.285862 , mean =  4.560791 , max =  8.224023  
#>   -> residual function :  difference between y and yhat (  default  )
#>   -> residuals         :  numerical, min =  -2.901839 , mean =  0.005875479 , max =  2.586552  
#>   A new explainer has been created!

explainer_sub1 <- explain(
  final_model,
  data = subset1_data %>% select(-value),
  y = subset1_data %>% select(value),
  label = "Subset 1"
)
#> Preparation of a new explainer is initiated
#>   -> model label       :  Subset 1 
#>   -> data              :  9  rows  123  cols 
#>   -> data              :  tibble converted into a data.frame 
#>   -> target variable   :  Argument 'y' was a data frame. Converted to a vector. (  WARNING  )
#>   -> target variable   :  9  values 
#>   -> predict function  :  yhat.ranger  will be used (  default  )
#>   -> predicted values  :  No value for predict function target column. (  default  )
#>   -> model_info        :  package ranger , ver. 0.18.0 , task regression (  default  ) 
#>   -> predicted values  :  numerical, min =  1.285862 , mean =  2.90092 , max =  6.272874  
#>   -> residual function :  difference between y and yhat (  default  )
#>   -> residuals         :  numerical, min =  -1.276092 , mean =  -0.01203065 , max =  2.194828  
#>   A new explainer has been created!

explainer_sub2 <- explain(
  final_model,
  data = subset2_data %>% select(-value),
  y = subset2_data %>% select(value),
  label = "Subset 2"
)
#> Preparation of a new explainer is initiated
#>   -> model label       :  Subset 2 
#>   -> data              :  22  rows  123  cols 
#>   -> data              :  tibble converted into a data.frame 
#>   -> target variable   :  Argument 'y' was a data frame. Converted to a vector. (  WARNING  )
#>   -> target variable   :  22  values 
#>   -> predict function  :  yhat.ranger  will be used (  default  )
#>   -> predicted values  :  No value for predict function target column. (  default  )
#>   -> model_info        :  package ranger , ver. 0.18.0 , task regression (  default  ) 
#>   -> predicted values  :  numerical, min =  1.819195 , mean =  5.633678 , max =  8.224023  
#>   -> residual function :  difference between y and yhat (  default  )
#>   -> residuals         :  numerical, min =  -2.901839 , mean =  0.2754127 , max =  2.167701  
#>   A new explainer has been created!

explainer_list <- list(
  `Full-Dataset` = explainer_full,
  `Subset 1` = explainer_sub1,
  `Subset 2` = explainer_sub2
)

imp_full <- lapply(explainer_list, feature_importance, B = 100, type = "ratio")
# Setting type = "ratio" means variable importance

plot(
  imp_full,
  max_vars = 7, # just show the top 7 predictors for each model
  bar_width = 4
)
```

<img src="man/figures/README-unnamed-chunk-32-1.png" alt="" width="100%" />

The above plots show the relative importance of individual predictor
variables expressed as a percent loss in model accuracy (measured as
RMSE) when that variable is permuted. A greater loss in predictive
accuracy implies that variable is more important in model accuracy.

As we expect, the relative importance of predictors varies depending on
the region we are investigating. However, with this many predictors, it
is difficult to make confident statements about the overall importance
of broader groups of predictors (i.e., the relative importance of
Land-cover overall vs geology; or the relative importance of local
affects measured through HAIFLO vs upstream effects measured through
HAIFLS). We can make these inferences by defining predictor groups.

``` r
# Define predictor variable groups
predictor_groups <- list(
  Autocorrelation = colnames(final_data)[grepl(
    "Prop_catch_",
    colnames(final_data)
  )],
  Landcover = colnames(final_data)[grepl("^LC_", colnames(final_data))],
  Geology = colnames(final_data)[grepl("^GEO_NAME_", colnames(final_data))],
  Slope = colnames(final_data)[grepl("^slope_", colnames(final_data))]
)

weighting_groups <- list(
  `lumped` = colnames(final_data)[grepl("lumped", colnames(final_data))],
  `iFLS` = colnames(final_data)[grepl("_iFLS", colnames(final_data))],
  `iFLO` = colnames(final_data)[grepl("_iFLO", colnames(final_data))],
  `HAiFLS` = colnames(final_data)[grepl("HAiFLS", colnames(final_data))],
  `HAiFLO` = colnames(final_data)[grepl("HAiFLO", colnames(final_data))]
)

imp_pred_gp <- map(
  explainer_list,
  ~ feature_importance(
    .,
    label = paste0(.$label, " Predictor Groups"),
    B = 100,
    variable_groups = predictor_groups,
    type = "ratio"
  )
)

imp_weighting_gp <- map(
  explainer_list,
  ~ feature_importance(
    .,
    label = paste0(.$label, " Weighting Groups"),
    B = 100,
    variable_groups = weighting_groups,
    type = "ratio"
  )
)

# install.packages("cowplot)

cowplot::plot_grid(
  plot(imp_pred_gp, bar_width = 7, title = NULL, subtitle = ""),
  plot(imp_weighting_gp, bar_width = 7, title = NULL, subtitle = ""),
  ncol = 2
)
```

<img src="man/figures/README-unnamed-chunk-33-1.png" alt="" width="100%" />

These plots show the relative importance of different predictors and
weighting schemes for the three subsets of data.

There are many more ways to explore model results with the DALEX
package, but we’ll leave that for you to explore.

## 8 Future Plans

This package is considered in early alpha, and bugs are expected. There
is lots of room for improvements in many of the functions. The core
functions may still change as the package develops and matures.

- One big area of improvement will be the handling of the output gpkg
  file structure.
- output gpkg may be replaced by a spatialite/rasterlite database file.
- Many functions can be cleaned and optimized further.

## 9 References

Kielstra, B., Mackereth, R., Melles, S., & Emilson E. (2021).
hydroweight: Inverse distance-weighted rasters and landscape attributes
(v1.0.0). Zenodo. <https://doi.org/10.5281/zenodo.4728559>

Lindsay, J.B. (2016). Whitebox GAT: A case study in geomorphometric
analysis. Computers & Geosciences, 95: 75-84.
<https://doi.org/10.1016/j.cageo.2016.07.003>

Peterson, E. E., Sheldon, F., Darnell, R., Bunn, S. E., & Harch, B. D.
(2011). A comparison of spatially explicit landscape representation
methods and their relationship to stream condition. Freshwater Biology,
56(3), 590–610. <https://doi.org/10.1111/j.1365-2427.2010.02507.x>

Przemyslaw Biecek and Tomasz Burzykowski. (2020). Explanatory Model
Analysis. <https://pbiecek.github.io/ema/>

R Core Team (2022). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>

Wickham, H., Bryan, J. (2021). R Packages. 2nd edition.
<https://r-pkgs.org/>.

Wu, Q. (2020). whitebox: ‘WhiteboxTools’ R Frontend. R package version
1.4.0. <https://github.com/giswqs/whiteboxR>

Zimmerman, Dale L., and Jay M. Ver Hoef. (2017).”The Torgegram for
fluvial variography: characterizing spatial dependence on stream
networks. Journal of Computational and Graphical Statistics 26.2
253-264. <https://doi.org/10.1080/10618600.2016.1247006>

[Back to top](#1-introduction)
