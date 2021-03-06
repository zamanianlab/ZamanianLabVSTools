
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ZamanianLabVSTools

<!-- badges: start -->
<!-- badges: end -->

## Overview

`ZamanianLabVSTools` contains a set of utilities used during virtual
screening/docking analyses, particularly during evaluation of docking
results. Modeling of the data is performed in the
[`tidymodels`](https://github.com/tidymodels) framework.

## Installation

You can install the development version of ZamanianLabVSTools from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wheelern/ZamanianLabVSTools")
```

## Usage

The primary contribution of this package is during fit evaluation, where
some non-standard metrics are used. These have been developed
specifically for virtual screening and seek to solve the “early
recognition” problem. These metrics can be used to evaluate model
performance on hold-out data. Here’s an example using the output from a
docking run with [GNINA](https://github.com/gnina/gnina), along with a
few custom features; these features were then use for fitting with a
random forest:

``` r
library(ZamanianLabVSTools)
library(dplyr)

data("vs_rf_predictions")

set <- yardstick::metric_set(ef, bedroc, rie)

vs_rf_predictions %>% 
  set(truth, .pred_TRUE)
#> # A tibble: 3 × 3
#>   .metric .estimator .estimate
#>   <chr>   <chr>          <dbl>
#> 1 ef      binary         4.04 
#> 2 bedroc  binary         0.115
#> 3 rie     binary         1.97
```

## References

BEDROC: Truchon J-F, Bayly CI. Evaluating virtual screening methods:
good and bad metrics for the “early recognition” problem. J Chem Inf
Model. 2007 Mar;47(2):488–508.

RIE: Sheridan RP, Singh SB, Fluder EM, Kearsley SK. Protocols for
bridging the peptide to nonpeptide gap in topological similarity
searches. J Chem Inf Comput Sci. 2001 Sep;41(5):1395–406.

## Acknowledgments

Much of the code for `bedroc()` and `rie()` was adapted from
[`enrichvs`](https://rdrr.io/cran/enrichvs/)
