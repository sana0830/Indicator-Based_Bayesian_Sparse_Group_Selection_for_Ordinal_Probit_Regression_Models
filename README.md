# Indicator-Based Bayesian Sparse Group Selection for Ordinal Probit Regression Models

This repository contains R scripts and R Markdown files for implementing indicator-based Bayesian selection methods for ordinal probit regression models. The main method is a Bayesian sparse group selection approach, which selects important groups and important variables within selected groups.

---

## Folder Structure

- `R/`  
  Contains the main functions used for model fitting and diagnostics.

  - `fit_sparse_group.R`  
    Proposed Bayesian sparse group selection method.  
    Main function: `fit_ordinal_sparse_group()`.

  - `fit_sparse_group_fast.R`  
    Accelerated version of the proposed sparse group selection method.  
    Main function: `fit_ordinal_sparse_group_fast()`.

  - `fit_baseline_selection.R`  
    Baseline methods for comparison.  
    Main function: `fit_baseline_selection()`.

    If `group_list = list()`, ordinary variable selection is performed.  
    If `group_list` is provided, group selection is performed.

  - `utils.R`  
    Utility functions shared by different fitting methods.

  - `mcse_diagnostics.R`  
    MCSE-based convergence diagnostic functions.

- `simulation_data/`  
  Contains simulated datasets and fitting results.

- `simulation_data_fitting.Rmd`  
  R Markdown file for simulation studies.

- `real_data_fitting_template.Rmd`  
  R Markdown template for empirical data analysis.  
  The empirical dataset is not included due to confidentiality restrictions.

- `performance_metrics.R`  
  Functions for evaluating variable selection, group selection, and prediction performance.

---

## Data Format

All fitting functions assume that the input data are stored as an `n x (p + 1)` data frame or matrix.

- The first `p` columns are predictors.
- The last column is the ordinal response variable.
- The response variable must be coded as consecutive integers:

```r
1, 2, ..., K
```

---

## Group Structure

For sparse group selection and group selection, the group structure is specified by `group_list`.

```r
group_list <- list(
  c(1, 2, 3),
  c(4, 5),
  c(6, 7, 8)
)
```

Each element contains the column indices of predictors in one group.

If `group_list = list()`, ordinary variable selection is performed.

---

## How to Run

Install the required packages:

```r
install.packages(c("MASS", "msm"))
```

Load the functions:

```r
source("R/utils.R")
source("R/mcse_diagnostics.R")
source("R/fit_sparse_group.R")
source("R/fit_sparse_group_fast.R")
source("R/fit_baseline_selection.R")
source("R/performance_metrics.R")
```

Run the proposed sparse group selection method:

```r
fit <- fit_ordinal_sparse_group(
  fitting_data = train_data,
  iter = 30000,
  iter_inn = 2000,
  a = 1,
  b = 10,
  theta = 0.5,
  rho = 0.5,
  seed_num = 1,
  group_list = group_list,
  use_mcse = TRUE
)
```

Run the accelerated sparse group selection method:

```r
fit_fast <- fit_ordinal_sparse_group_fast(
  fitting_data = train_data,
  iter = 30000,
  a = 1,
  b = 10,
  theta = 0.5,
  rho = 0.5,
  seed_num = 1,
  group_list = group_list,
  use_mcse = TRUE
)
```

Run the baseline variable selection method:

```r
fit_var <- fit_baseline_selection(
  fitting_data = train_data,
  iter = 30000,
  a = 1,
  b = 10,
  theta = 0.5,
  seed_num = 1,
  group_list = list(),
  use_mcse = TRUE
)
```

Run the baseline group selection method:

```r
fit_group <- fit_baseline_selection(
  fitting_data = train_data,
  iter = 30000,
  a = 1,
  b = 10,
  theta = 0.5,
  seed_num = 1,
  group_list = group_list,
  use_mcse = TRUE
)
```

---

## Output

The fitting functions return a list containing posterior samples after burn-in.

Common output elements include:

- `gamma_record`: posterior samples of variable indicators.
- `beta_record`: posterior samples of regression coefficients.
- `tau_record`: posterior samples of cutpoints.
- `eta_record`: posterior samples of group indicators, if group structure is used.
- `burn_in`: final burn-in iteration.
- `last_iter`: last MCMC iteration used.
- `mcse_iteration`: iteration at which the MCSE stopping criterion is reached.

---

## Methods Included

This repository includes three types of methods:

1. **Sparse group selection**  
   Selects important groups and important variables within selected groups.

2. **Accelerated sparse group selection**  
   Uses a faster group-level update and incremental within-group updates.

3. **Baseline selection methods**  
   Uses `fit_baseline_selection()` for both ordinary variable selection and group selection.  
   The setting `group_list = list()` gives ordinary variable selection, while providing `group_list` gives group selection.

---

## Notes

- The ordinal response variable must be placed in the last column of `fitting_data`.
- The response variable must be coded as `1, ..., K`.
- The main proposed method is implemented in `R/fit_sparse_group.R`.
- The accelerated version is implemented in `R/fit_sparse_group_fast.R`.
- The baseline comparison methods are implemented in `R/fit_baseline_selection.R`.
