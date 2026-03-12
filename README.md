# algorithmDiagnostics

Diagnostic plots for evaluating clinical algorithms for conditional biases across subgroups.

## Installation

```r
# install.packages("devtools")
devtools::install_github("AartikSarma/algorithm_diagnostics")
```

## Usage

```r
library(algorithmDiagnostics)

result <- conditional_bias_plot(
  data = my_data,
  dependent_vars = c("outcome1", "outcome2"),
  independent_vars = c("predictor1"),
  grouping_vars = c("race", "sex"),
  n_tiles = 10,
  conf_level = 0.95
)
print(result)
```

## Overview

`conditional_bias_plot()` generates multi-panel figures that show:

- **Percentile-level scatter**: raw outcome averages at each percentile of the predictor, stratified by group
- **Quantile-aggregated estimates**: point estimates with confidence intervals (Wilson for binary outcomes, t-based for continuous)
- **LOESS smoothers**: trend lines for each subgroup
- **Beeswarm strips**: distribution of observations across predictor percentiles by group

Panels are tiled across dependent variables (rows), independent variables (columns), and grouping variables (vertical sections) using patchwork.
