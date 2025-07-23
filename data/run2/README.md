# README – run2 Directory

This directory stores all processed data and outputs of the model training.

## Directory Structure

```
run2/
├── folds/                  # Spatial partitioning files for virtual species (.RDS)
├── maps/                   # Raster predictions (.tif) for trained SDMs
├── models/                 # Fitted model objects (.RDS)
├── results/                # Evaluation metrics on test data (.RDS)
├── virtualSpeciesTrain/    # GeoPackage files of training data
└── results.RDS             # Aggregated evaluation results (data frame)
```

## File Descriptions

### `results.RDS`

This file aggregates all individual evaluation files from the `results/` folder into a single `data.frame` with 124,081 observations and 24 variables. The file has been created with the script `R/004_evaluation.R`

**Structure:**

| **Column Name**    | **Type**    | **Description**                                                                                                               |
| ------------------ | ----------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `metric`           | `numeric`   | Composite evaluation index; reflects overall model performance across multiple metrics.                 |
| `AUC`              | `numeric`   | Area Under the Receiver Operating Characteristic Curve. |
| `COR`              | `numeric`   | Pearson correlation coefficient.                    |
| `Spec`             | `numeric`   | Specificity: proportion of true absences correctly predicted (true negative rate).                                            |
| `Sens`             | `numeric`   | Sensitivity: proportion of true presences correctly predicted (true positive rate).                                           |
| `Kappa`            | `numeric`   | Cohen’s Kappa statistic: agreement between predicted and observed occurrences corrected for chance.                           |
| `PCC`              | `numeric`   | Percent Correctly Classified.                                                              |
| `TSS`              | `numeric`   | True Skill Statistic: `Sensitivity + Specificity − 1`, ranging from -1 to +1.                                                 |
| `stability`        | `numeric`   | Prediction stability index across the six replicates of the spatial prediction.         |
| `PRG`              | `numeric`   | AUC-PRG: Area Under the Precision-Recall-Gain curve.                                |
| `MAE`              | `numeric`   | Mean Absolute Error.                            |
| `BIAS`             | `numeric`   | Mean prediction bias (signed error); measures tendency to over- or underpredict presence.                                     |
| `noPresencePoints` | `numeric`   | Number of presence records in the test fold.                                  |
| `trueCor`          | `numeric`   | Spatial Pearson correlation between predicted raster and the known (true) virtual species suitability raster.                       |
| `model`            | `character` | SDM algorithm used. One of: `Lasso`, `RF`, `Maxent`, `BRT`, `GAM`.                 |
| `size`             | `character` | Fold partitioning strategy: one of `random`, `block1`, `block2`, `KNNDM`, or `clusters`.                                      |
| `testData`         | `integer`   | Identifier for the held-out fold in 6-fold cross-validation (values: 1–6).                                                    |
| `method`           | `character` | Evaluation method label (e.g., `PA`, `PAA`, `PBG`, etc.) specifying the test dataset used.                          |
| `corTrainTest`     | `numeric`   | Spatial correlation between training and test points based on bio_1. For each training point, the mean of the five nearest bio_1 values in the test set is computed, and Pearson’s correlation is calculated with the training bio_1 values. Serves as a proxy for environmental similarity across spatial folds.      |
| `moranTrain`       | `numeric`   | Moran’s I spatial autocorrelation value for presence points in the training data.                                             |
| `moranTest`        | `numeric`   | Moran’s I spatial autocorrelation value for presence points in the test data.                                                 |
| `replicate`        | `character` | Replicate number (1–5).                  |
| `points`           | `character` | Number of presence–absence points sampled. One of: `40`, `80`, `120`, `160`, `200`, `400`, etc.                       |
| `species`          | `character` | Virtual species identifier in the format `VSXX`, where `XX` ranges from `01` to `10`.                                         |


This dataset serves as the primary output for downstream statistical analyses.