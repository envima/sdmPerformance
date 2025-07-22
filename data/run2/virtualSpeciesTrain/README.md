# README – virtualSpeciesTrain

This directory contains 600 training datasets for species distribution modeling, derived from presence–absence samples (`data/PA`) and annotated with environmental covariates and multiple (spatial) fold separation assignments. The datasets were generated with the script `R/002_dataPartition.R`. These are the input files for model training.

Each file is stored in GeoPackage (`.gpkg`) format and contains a single dataset representing one combination of species, sample size, and replicate as well as five different fold separations of the points. File names follow the structure:

```
VSXX_YYY_ZZ.gpkg
```

Where:

* `VSXX`: Virtual species identifier (e.g., `VS01`)
* `YYY`: Number of presence–absence samples (e.g., `40`)
* `ZZ`: Replicate index (1–10)


### Format

Each `.gpkg` file contains an `sf` object with point geometries and 12 attributes:

| Column     | Description                                                                   |
| ---------- | ----------------------------------------------------------------------------- |
| `Real`     | Ground-truth presence (1) or absence (0) label                                |
| `Observed` | Observed label (same as Real in this case)                                  |
| `KNNDM`    | Fold assignment from K-nearest neighbor distance matching                     |
| `random`   | Random fold assignment (1–6)                                                  |
| `block1`   | Spatial block fold (hexagonal, \~300 km)                                      |
| `block2`   | Spatial block fold (square, \~100 km)                                         |
| `clusters` | Environmental clustering-based fold assignment                                |
| `bio_1`    | Annual Mean Temperature (°C × 10), from bioclimatic variables                 |
| `bio_3`    | Isothermality (BIO3)                                                          |
| `bio_7`    | Temperature Annual Range (BIO7)                                               |
| `bio_12`   | Annual Precipitation (mm) (BIO12)                                             |
| `geom`     | Geometry column with point coordinates (EPSG:3577; GDA94 / Australian Albers) |


### Dependencies

Each dataset was generated using:

* Presence–absence data: `data/PA/`
* Spatial folds: `data/run2/folds/`
* Environmental predictors: `data/variables.tif` (see `data/`)
* Script: `002_dataPartition.R`

**Software used:**

R version 4.4.2
terra version 1.8-5
sf version 1.0-20    
