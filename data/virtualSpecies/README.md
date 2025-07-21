# README – Virtual Species Data

## Overview

The folder contains environmental suitability maps and full metadata objects for ten virtual species (VS01–VS10). These species were generated using the `virtualspecies` R package based on ecological response functions applied to bioclimatic variables. The virtual species were originally created by Grimmet et al. 2020.

## File Description

`VS01.tif` to `VS10.tif`:
Raster files representing continuous habitat suitability values for each virtual species. All rasters are in **EPSG:3577** (GDA94 / Australian Albers). Pixel values range from 0 (unsuitable) to 1 (highly suitable).

`VS01.RDS` to `VS10.RDS`:
Serialized R objects containing full metadata for each virtual species, including the suitability raster, the response functions, and the parameters used in their construction.
These R objects of class `virtualspecies` (a list with `S3` and `S4` components) contain:

`approach`:
A character string indicating the generation method (`"response"`, i.e., based on ecological response curves).

`details`:
A list documenting how suitability was computed:
    * `variables`: The four bioclimatic variables used (`bio_1`, `bio_3`, `bio_7`, `bio_12`)
    * `formula`: The multiplicative combination of response functions used to compute suitability (`bio_1 * bio_3 * bio_7 * bio_12`)
    * `rescale.each.response`: Whether each response curve was individually rescaled to \[0,1]
    * `rescale`: Whether the final suitability raster was rescaled to \[0,1]
    * `parameters`: A named list of parameter settings for each variable.

`suitab.raster`:
A `terra::SpatRaster` object containing the final suitability, with the same spatial characteristics as the corresponding `.tif` file.

## Data Generation

All files were generated using the R script `R/01_virtualspecies.R`. The virtual species were constructed with `virtualspecies::generateSpFromFun()` using a Gaussian response functions (`dnorm`) for four selected bioclimatic predictors: `bio_1` (annual mean temperature), `bio_3` (isothermality), `bio_7` (temperature annual range), and `bio_12` (annual precipitation).

## References

Grimmet, L., Whitsed, R., Horta, A. (2020). Presence-only species distribution models are sensitive to sample prevalence: Evaluating models using spatial prediction stability and accuracy metrics. *Ecological Modelling*. [https://doi.org/10.1016/j.ecolmodel.2020.109194](https://doi.org/10.1016/j.ecolmodel.2020.109194)

Leroy, B., Meynard, C.N., Bellard, C., Courchamp, F. (2015). virtualspecies, an R package to generate virtual species distributions. Ecography. [https://doi.org/10.1111/ecog.01388](https://doi.org/10.1111/ecog.01388)

**Software used:**

R version 4.4.2
virtualspecies version 1.6
terra version 1.8-5