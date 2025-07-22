# README – Creation of Spatial Folds

This directory contains spatial cross-validation fold separations for each sampled virtual species dataset (`data/PA`). The folds were generated with the script `002_dataPartition.R` and are stored in `.RDS` format. A total of 2,400 files are present.

Each file encodes a single fold separation object for a given points dataset and partitioning strategy. File names follow the format:

```
VSXX_YYY_ZZ_FOLDS.RDS
```

Where:

* `VSXX` is the species identifier (e.g., `VS10`)
* `YYY` is the number of sampled points (e.g., `400`)
* `ZZ` indicates the replicate number (e.g., `1`)
* `FOLDS` is the fold creation strategy (`KNNDM`, `random`, `clusters`, `block1`, `block2`)

These folds are used to assign presence–absence points to one of six cross-validation folds (`k = 6`), supporting spatially explicit model evaluation.

### Partitioning strategies

* **KNNDM**: K-nearest neighbor distance matching, as implemented in the `CAST::knndm()` function. This method clusters spatial points based on multivariate similarity using environmental covariates.
* **random**: Simple random assignment of fold numbers from 1 to 6.
* **block1**: Spatial blocking via `blockCV::cv_spatial()` with hexagonal blocks of \~300 km diameter.
* **block2**: Spatial blocking via `blockCV::cv_spatial()` using square tiles of \~100 km extent.
* **clusters**: Environmental clustering using `blockCV::cv_cluster()` which groups spatial units based on environmental similarity.

### Format

Each `.RDS` file is an object specific to the corresponding method:

* `block1`, `block2`, `clusters`: List objects with slot `$folds_ids` mapping each data point to a fold.
* `KNNDM`: An object of class `"knndm"` with a `$clusters` slot.
* `random` folds are not stored independently but as part of the presence-absence `.gpkg` in `virtualSpeciesTrain/`.

### **Dependencies**

Each fold was generated using:

* Presence-absence data from: `data/PA/`
* Environmental predictors from: `data/variables.tif` (see `data/` folder)
* Code used: `002_dataPartition.R`

### References

Linnenbrink, J., Mila, C., Ludwig, M., Meyer, H. (2024). **kNNDM CV: k-fold nearest-neighbour distance matching cross-validation for map accuracy estimation.** *Geoscientific Model Development*. 17(15), 5897–5912. [https://doi.org/10.5194/gmd-17-5897-2024](https://doi.org/10.5194/gmd-17-5897-2024)

Valavi, R., Elith, J., Lahoz-Monfort, J. J., & Guillera-Arroita, G. (2019). **blockCV: An R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models**. *Methods in Ecology and Evolution*, 10(2), 225–232. [https://doi.org/10.1111/2041-210X.13107](https://doi.org/10.1111/2041-210X.13107)


**Software used:**

R version 4.4.2
sf version 1.0-20     
terra version 1.8-5   
CAST version 1.0.3    
blockCV version 3.1-5

