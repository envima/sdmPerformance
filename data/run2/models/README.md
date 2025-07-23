# README – Models

This directory contains 45,000 species distribution model (SDM) objects trained across multiple virtual species, partitioning strategies, algorithms, and spatial folds. Each model was generated using the script `R/003_modeling.R`.

All models are saved in `.RDS` format and represent fully trained models.

### File Naming Convention

```
VSXX_PARTITION_ALGO_testDataFOLD_pointsNN_replicatesRR.RDS
```

Where:

* `VSXX` — Virtual species identifier (e.g., `VS01`)
* `PARTITION` — Partitioning strategy used for spatial folds (e.g., `block1`, `KNNDM`)
* `ALGO` — SDM algorithm used (e.g., `BRT`, `RF`, `GAM`, `Lasso`, `Maxent`)
* `FOLD` — Index of the held-out test fold (`1`–`6`)
* `NN` — Number of sampled presence–absence points (e.g., `40`)
* `RR` — Replicate number (`1`–`5`)

Each `.RDS` file stores a fitted model object.

### Model Algorithms

Five modeling algorithms were implemented:

* **BRT** — Boosted Regression Trees
* **RF** — Random Forest
* **GAM** — Generalized Additive Models
* **Lasso** — L1-regularized Generalized Linear Models
* **Maxent** — Maximum Entropy modeling

Model selection, training, and evaluation were standardized across algorithms using the function `trainSpeciesDistributionModel()` from `R/functions/`. For each model on of the test folds is withheld for later model evalation.

---

### Dependencies

Each model was generated using:

* Species occurrence and background data: `data/run2/virtualSpeciesTrain/`
* Spatial folds: column values within those GeoPackages (see `README – virtualSpeciesTrain`)
* Environmental predictors: `data/variables.tif`
* Background data: `data/bg.gpkg`
* Code: `002_modeling.R`, and helper functions:

  * `trainSpeciesDistributionModel.R`

---

### Output Content

Each `.RDS` file is a serialized R object, typically of class:

* `"gbm"` for BRT models
* `"randomForest"` for RF
* `"gam"` for GAM
* `"glmnet"` for Lasso
* `"MaxEnt"` or custom S4 class for Maxent (via `dismo` or `maxnet`)


**Software used:**

R version 4.4.1

Platform: x86_64-pc-linux-gnu

Running under: Ubuntu 24.04.2 LTS

randomForest version 4.7-1.2 

ENMeval version 2.0.5.2

glmnet version 4.1-8

CAST version 1.0.3

rJava version 1.0-11

caret version 7.0-1

gbm version 2.2.2

sf version 1.0-19

terra version 1.8-54

gam version 1.22-5       

Custom scripts from Valavi et al. (2023) benchmark [https://osf.io/puk8v](https://osf.io/puk8v)

---

### References

* Valavi, R., Elith, J., & Guillera-Arroita, G. (2023). A comprehensive benchmark for SDM evaluation. *Global Ecology and Biogeography*, 32(9), 1413–1429. [https://doi.org/10.1111/geb.13639](https://doi.org/10.1111/geb.13639)
* Valavi, R., Elith, J., Lahoz-Monfort, J. J., & Guillera-Arroita, G. (2019). blockCV: An R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models. *Methods in Ecology and Evolution*, 10(2), 225–232. [https://doi.org/10.1111/2041-210X.13107](https://doi.org/10.1111/2041-210X.13107)


maxent