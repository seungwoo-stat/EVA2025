# EVA 2025

This repo contains team **SHSmultiscale**'s source code for the 14th International Conference on Extreme Value Analysis (EVA 2025) [Data Challenge](https://eva2025.unc.edu/data-challenge/).

## Overview

- `SUMMER.rds` and `WINTER.rds` are datasets we have used in this data challenge. The datasets provided by the data challenge organizers were converted to 121440 x 25 (`SUMMER.rds`) and 119460 x 25 (`WINTER.rds`) numeric matrices, each consisting of the values for summer months (May--Oct.) and winter months (Nov.--Apr.). Each row represents a specific date, and each column represents one of the 25 cells.

- `Functions.R` contains all the functions for conducting the extremal PCA. We modified the code from the Supplemental Content B of Rohrbeck and Cooley (2023).

- Run `m-bandwidth-select.R` to find the optimal number of PCs (m) and the bandwidth of the spherical kernel. Note that the computation is executed in parallel across a number of cores, whose number is set to 6. Modify the variable `cl` if needed.

- Run `EPCA.R` to get point estimates of the target quantities and point estimates after bootstrapping the given dataset. Note that the computation is executed in parallel across a number of cores, whose number is set to `min(detectCores()-2, 15)`. See the `Functions.R` file.

- Run `Confidence-interval.R` to compute the 95% bootstrap confidence interval for each of the target quantities.


## References
- Rohrbeck, C. and Cooley, D. (2023). [Simulating flood event sets using extremal principal components](https://doi.org/10.1214/22-AOAS1672). *Ann. Appl. Stat.* **17**(2) 1333--1352.
