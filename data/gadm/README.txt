# README – Administrative Boundary Data

## Overview

The folder contains spatial administrative boundary data used for masking and subsetting bioclimatic predictors to specific regions within Australia. The data are sourced from the GADM database (version 4.1).

## File Description

`gadm41_AUS_1_pk.rds`:
An R serialized object (RDS file) containing level 1 administrative boundaries for Australia (i.e., states and territories). The data are provided in the native GADM coordinate reference system (WGS84, **EPSG:4326**).

## Data Generation

The file `gadm41_AUS_1_pk.rds` was generated using the R script `R/01_virtualspecies.R`.

## References

GADM (2022). Database of Global Administrative Areas, version 4.1. [https://gadm.org](https://gadm.org)

**Software used:**

R version 4.4.2
geodata version 0.6-2

