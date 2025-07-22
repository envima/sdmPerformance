# README – Data Directory

## Overview

The `data/` directory contains all data used for simulation, modeling, and evaluation in this study. This dataset consists of environmental variables, virtual species distributions, presence-absence points.

## Directory Structure

The directory is organized as follows:

```
data/
├── climate/          # Raw bioclimatic variables
├── gadm/             # Administrative boundaries
├── virtualSpecies/   # Simulated virtual species
├── paRaster/         # Binarized presence-absence rasters
├── PA/               # Sampled presence-absence data
├── run2/             # Output directory (e.g., models, maps, evaluation results)
├── variables.tif     # Environmental predictor stack used for modeling
├── bg.gpkg           # Background points
```

See the README files of the subdirectories for more information.

## File Descriptions

### `variables.tif`

The code used to create the environmenal variables is available in script `R/001_virtualspecies.R`. A multi-layer raster stack containing four bioclimatic predictors:

* `bio_1` – Annual Mean Temperature
* `bio_3` – Isothermality (BIO2/BIO7 × 100)
* `bio_7` – Temperature Annual Range (BIO5–BIO6)
* `bio_12` – Annual Precipitation

**Source**: [WorldClim v2.1](https://www.worldclim.org/) via `geodata::worldclim_country()`
**Projection**: GDA94 / Australian Albers (EPSG:3577)
**Extent**: Cropped and masked to New South Wales, Victoria, and the Australian Capital Territory

### `bg.gpkg`

A GeoPackage containing 10,000 background points randomly sampled across the extent of `variables.tif`. Used for model fitting in presence-background based modeling approaches.

These points were randomly sampled across the extent of the environmental predictor raster (`variables.tif`) and annotated with environmental variables values, simulated absence labels, and folds for cv.

Note: While the column names reflect different spatial cross-validation strategies (random, KNNDM, clusters, block1, block2), these assignments were generated using uniform random sampling (folds 1–6). This is due to the fact that background points do not exhibit sampling bias, unlike presence-absence samples.

**Format**: `sf` object with point geometries
**Coordinate Reference System**: GDA94 / Australian Albers (EPSG:3577)

**Attributes**:

| Column     | Description                                                  |
| ---------- | ------------------------------------------------------------ |
| `Real`     | Always 0 for background points |
| `Observed` | Always 0 for background points |
| `bio_1`    | Annual Mean Temperature (°C × 10)                            |
| `bio_3`    | Isothermality                                                |
| `bio_7`    | Temperature Annual Range                                     |
| `bio_12`   | Annual Precipitation (mm)                                    |
| `random`   | Random fold assignment (1–6)                                 |
| `KNNDM`    | Random fold assignment (1–6)      |
| `clusters` | Random fold assignment (1–6)                                |
| `block1`   | Random fold assignment (1–6)            |
| `block2`   | Random fold assignment (1–6)          |
| `geom`     | Simple feature geometry column with point coordinates        |


**Source**: Sampled using `predicts::backgroundSample()`. The code used to create the background points is available in script `R/002_dataPartition.R`.

**Software used:**

R version 4.4.2
sf version 1.0-20 
predicts version 0.1-17
terra version 1.8-5 



