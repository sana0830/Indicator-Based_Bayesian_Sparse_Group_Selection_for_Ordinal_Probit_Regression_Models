# Simulation Data

This folder contains simulated datasets and simulation results used in the simulation study.

The results are organized by method, tuning setting, simulation background, and parameter version.

## Folder structure

- `sparse_group/`: Results from the proposed Bayesian sparse group selection method.
- `group_selection/`: Results from the Bayesian group selection method used for comparison.
- `sparse_group_fast/`: Results from the accelerated version of the proposed sparse group selection method.

## Sparse group

The `sparse_group/` folder contains both tuning results and fixed-parameter results.

- `tuning/vary_a_fixed_b10/`: Results obtained by varying `a` while fixing `b = 10`.
- `tuning/vary_b_fixed_a1/`: Results obtained by varying `b` while fixing `a = 1`.
- `fixed_parameter/`: Results obtained after selecting a fixed hyperparameter setting.

For `vary_a_fixed_b10/`, only `background_1/` is included, and it is divided into three parameter versions.

For `vary_b_fixed_a1/`, `background_1/` contains `param_version_1/` and `param_version_2/`, while `background_2/` is not further divided by parameter version.

For fixed-parameter results, both `background_1/` and `background_2/` are divided into three parameter versions.

## Group selection

The `group_selection/` folder contains results under the fixed setting `a = 1` and `b = 10`.

Both `background_1/` and `background_2/` are divided into three parameter versions.

## Sparse group fast

The `sparse_group_fast/` folder contains results for the accelerated version.

Only two cases are included:

- `background_1/`: Simulation Background 1 with parameter version 1.
- `background_2/`: Simulation Background 2 with parameter version 2.
