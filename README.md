# Dynamic Policy Evaluation and Learning

This repository contains the simulation and empirical application code used to reproduce the results for the paper *Dynamic Policy Evaluation and Learning with Spatio-temporal Data*. 

## Quick Start

Clone the repository, open a terminal in the repository root, and install the required R packages listed below. The fastest reproducibility check uses the prepared simulation plotting data:

```bash
Rscript R/simulation/plots/policy_learning_plot.R
Rscript R/simulation/plots/policy_evaluation_plot.R
```

This writes the simulation figures under `outputs/plots/simulation/`. To inspect the computation for a single Monte Carlo replicate, run either simulation entry point with the arguments documented below. The empirical application is a separate, substantially heavier workflow that uses fixed public inputs from Harvard Dataverse plus two restricted local inputs and is described in [Application Analysis](#application-analysis).

## Reproducibility Scope

| Goal | Entry point | Inputs | Output |
|---|---|---|---|
| Recreate paper simulation figures | `R/simulation/plots/*.R` | Prepared objects in `data/simulation/` | PDFs in `outputs/plots/simulation/` |
| Run one policy-learning replicate | `R/simulation/policy_learning_simulation.R` | Command-line seed/model/sample-size arguments | RDS file in `outputs/results/simulation/policy_learning/` |
| Run one policy-evaluation replicate | `R/simulation/policy_evaluation_simulation.R` | Command-line seed/model/sample-size arguments | RDS file in `outputs/results/simulation/policy_evaluation/` |
| Rebuild the empirical application | `R/application/run_application_analysis.R` | Two restricted local CSV files (not distributed), fixed Dataverse files, and `geocausal` package data | Intermediate data, result objects, figures, and console tables |

Prepared plotting objects are included so that readers can recreate the reported simulation figures without rerunning and aggregating every Monte Carlo replicate. The repository also contains saved application results and figures for comparison with a fresh run.

## Folder Structure

```text
Dynamic-Policy-Evaluation-and-Learning/
  R/
    application/
      run_application_analysis.R
      functions/
        application_utilities.R
    simulation/
      policy_learning_simulation.R
      policy_evaluation_simulation.R
      functions/
        data_generating_process.R
        estimate_policy_value.R
        estimate_propensity_score.R
        estimator_transforms.R
        optimize_policy.R
        perform_test_only.R
        policy_learning_test.R
        true_optimal_policy.R
        true_policy_value.R
      plots/
        policy_learning_plot.R
        policy_evaluation_plot.R
  data/
    application/
      EXTERNAL_DATA_INVENTORY.md
    simulation/
      policy_learning_plot_data.rds
      policy_evaluation_plot_data.rds
  outputs/
    results/
    plots/
      application/
      simulation/
```

The simulation files in `data/simulation/` are prepared plotting data. They let the plot scripts reproduce the paper figures directly without requiring readers to rerun all simulation replicates or use the internal combine scripts.

The application workflow rebuilds its district-week analysis data from two restricted local application inputs plus fixed Harvard Dataverse file IDs. The two local inputs are not currently distributed in this repository while permission to release them is being confirmed.

## Requirements

The scripts were checked with R 4.2.0. Required R packages for simulations:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "mvtnorm", "scales"))
```

The application script additionally uses spatial, plotting, and table packages:

```r
install.packages(c(
  "cowplot", "data.table", "dplyr", "forcats", "geepack", "ggplot2",
  "httr", "kableExtra", "lubridate", "nngeo", "purrr", "readr",
  "rlang", "sf", "sfheaders", "spatstat", "spatstat.geom", "tibble",
  "tidytext", "tidyr", "viridis", "zoo"
))
```

The application also uses package datasets and functions from `geocausal`; the replication environment used `geocausal` version `0.3.4`. Base/recommended R packages such as `splines` and `stats` are used without separate installation.

Run all commands from the repository root. The code was syntax-checked with R 4.2.0; newer package releases can change numerical or graphical output slightly.


## Reproduce Simulation Plots

Policy learning plots:

```bash
Rscript R/simulation/plots/policy_learning_plot.R
```

Policy evaluation plots:

```bash
Rscript R/simulation/plots/policy_evaluation_plot.R
```

The scripts write PDFs to:

```text
outputs/plots/simulation/policy_learning/
outputs/plots/simulation/policy_evaluation/
```

## Run One Simulation Replicate

Policy learning:

```bash
Rscript R/simulation/policy_learning_simulation.R <sim> <model> <time_points> <n>
```

Example:

```bash
Rscript R/simulation/policy_learning_simulation.R 1 1 50 10
```

Policy evaluation:

```bash
Rscript R/simulation/policy_evaluation_simulation.R <sim> <model> <time_points> <n>
```

Example:

```bash
Rscript R/simulation/policy_evaluation_simulation.R 1 1 50 10
```

Single-replicate simulation outputs are saved under:

```text
outputs/results/simulation/
```

The simulation scripts are intended to show the reproducible computation for individual replicates and to provide code that readers can adapt. The full paper figures were generated from many replicates. Those already-combined plotting datasets are provided in `data/simulation/`.


## Application Analysis

The cleaned application workflow has one entry point:

```text
R/application/run_application_analysis.R
```

The application script expects these two local inputs before running:

```text
data/application/Maaws_PublicData.csv
data/application/CSM_filter.csv
```

Both input files are included in this repository. `CSM_filter.csv` is a minimal, US-only subset of the public troop-density source; see `data/application/EXTERNAL_DATA_INVENTORY.md` for provenance and the source DOI.


