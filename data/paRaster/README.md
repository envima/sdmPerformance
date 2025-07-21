# README – Presence-Absence Raster

## Overview

The folder contains binary presence–absence (PA) raster layers and corresponding metadata objects for ten virtual species (VS01–VS10). These data were derived from  suitability maps in folder `virtualSpecies`.

## File Description

`VS01.tif` to `VS10.tif`:
Raster layers representing binary presence (1) and absence (0) of each virtual species.

  * CRS: **EPSG:3577** (GDA94 / Australian Albers)
  * Data type: integer (0 = absence, 1 = presence)

`VS01.RDS` to `VS10.RDS`:
Serialized R objects containing detailed information about the presence–absence conversion process. Each object is a list with the following components:

| Element                     | Description                                                                                                       |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `suitab.raster`             | Packed `SpatRaster` of the original suitability raster used for PA conversion.                                    |
| `probability.of.occurrence` | Packed `SpatRaster` representing the estimated probability of species occurrence derived from suitability values. |
| `PA.conversion`             | Named character vector providing parameters of the PA conversion:                                                 |
|                             | - `conversion.method`: Method used ("probability")                                                                |
|                             | - `probabilistic.method`: Logistic function                                                                       |
|                             | - `alpha` and `beta`: Parameters of the logistic function applied                                                 |
|                             | - `species.prevalence`: Target prevalence value enforced during conversion                                        |
| `pa.raster`                 | Packed `SpatRaster` containing the final presence–absence map (binary values 0/1).                                |


## Data Generation

The files were generated using the R script `R/01_virtualspecies.R` by applying the function `convertToPA()` from the `virtualspecies` package to each environmental suitability raster. Species-specific prevalence values were predefined as follows (Grimmet et al. 2020):

| Species | Prevalence |
| ------- | ---------- |
| VS01    | 0.35       |
| VS02    | 0.34       |
| VS03    | 0.33       |
| VS04    | 0.29       |
| VS05    | 0.26       |
| VS06    | 0.21       |
| VS07    | 0.15       |
| VS08    | 0.12       |
| VS09    | 0.11       |
| VS10    | 0.05       |

## References

Grimmet, L., Whitsed, R., Horta, A. (2020). Presence-only species distribution models are sensitive to sample prevalence: Evaluating models using spatial prediction stability and accuracy metrics. *Ecological Modelling*. [https://doi.org/10.1016/j.ecolmodel.2020.109194](https://doi.org/10.1016/j.ecolmodel.2020.109194)

Leroy, B., Meynard, C.N., Bellard, C., Courchamp, F. (2015). virtualspecies, an R package to generate virtual species distributions. Ecography. [https://doi.org/10.1111/ecog.01388](https://doi.org/10.1111/ecog.01388)

**Software used:**

R version 4.4.2
virtualspecies version 1.6
terra version 1.8-5