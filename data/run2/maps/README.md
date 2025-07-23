# README – `maps`

## Overview

This directory contains raster prediction outputs from species distribution models (SDMs) trained on virtual species data. Each file is a GeoTIFF (`.tif`) representing model-predicted environmental suitability across Australia for a single modeling configuration (species × sampling strategy × algorithm × data split × replicate).

All raster files were generated during the execution of the script `002_modeling.R`, which implements the modeling workflow adapted from Valavi et al. (2023). These files correspond one-to-one with the model objects saved in the `models` directory and are used in subsequent performance evaluations.

## File Structure and Naming Convention

Each file follows this standardized naming pattern:

```
<species>_<partitioning>_<model>_testData<fold>_points<sample_size>_replicates<replicate>.tif
```

### Components:

* `species`: Identifier for the virtual species (e.g., `VS01`–`VS10`)
* `partitioning`: Spatial data partitioning method (`block1`, `block2`, `clusters`, `KNNDM`, or `random`)
* `model`: Modeling algorithm (`BRT`, `RF`, `GAM`, `Lasso`, or `Maxent`)
* `testData`: Integer from 1 to 6, indicating the test data fold
* `points`: Number of presence points used for model training (e.g., `40`, `80`, etc.)
* `replicates`: Replicate number (1 to 5)

### Example:

```
VS01_block1_BRT_testData1_points80_replicates1.tif
```

This file is the prediction output from a BRT model trained on species `VS01`, using the `block1` partitioning scheme, trained with 80 presence points in replicate 1, with fold 1 used as test data.

## Format and Spatial Reference

* **File format**: GeoTIFF (`.tif`)
* **Resolution**: Matches input predictor layers (see `variables.tif`)
* **Coordinate reference system (CRS)**: Same as the original bioclimatic raster layers (typically WGS 84, EPSG:4326)

Each raster encodes predicted environmental suitability (continuous values from 0 to 1) at each pixel based on the trained SDM.

## Dependencies

These prediction rasters are generated from:

* Predictor raster stack (`variables.tif`)
* Virtual species training datasets (`virtualSpeciesTrain/*.gpkg`)
* Background points (`bg.gpkg`)
* Modeling functions (`trainSpeciesDistributionModel.R`)
* Model configurations from the `002_modeling.R` script

## Use

These files serve as the spatial prediction layer for evaluating model performance and generalizability. Evaluation results based on these predictions are stored in the `results` folder and summarized in `results.RDS`.