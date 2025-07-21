# README – bioclimatic data

## Overview

The folder contains spatially explicit bioclimatic predictor variables used for modeling virtual species distributions. The data are derived from WorldClim version 2.1 at a spatial resolution of 30 arc-seconds and have been specifically clipped to selected Australian regions.

### Subdirectory Structure

`climate/wc2.1_country/`:
Contains processed bioclimatic raster layers, downloaded from the WorldClim dataset and masked to Australia.

## File Description

`AUS_wc2.1_30s_bio.tif`:
A GeoTIFF raster stack containing all 19 bioclimatic variables. The resolution is 30 arc-seconds (approximately 1 km²). The data covers all of Australia, and is avialable in **EPSG:4326** (WGS84).

## Data Generation

The file `AUS_wc2.1_30s_bio.tif` was generated using the R script `R/01_virtualspecies.R`.

**References:**

Fick, S.E., and Hijmans, R.J. (2017). WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, 37(12), 4302–4315. https://doi.org/10.1002/joc.5086

Climate data source: WorldClim 2.1 – [https://www.worldclim.org/data/worldclim21.html](https://www.worldclim.org/data/worldclim21.html)

**Software used:**

R version 4.4.2
geodata version 0.6-2 
terra version 1.8-5
