# Formal Privacy Feasibility Appendix

This repository contains code and data for a feasibility study of formally private methods for histograms, means, quantiles, and linear regression by Andres Felipe Barrientos, Aaron R. Williams, Joshua Snoke, and Claire McKay Bowen. 

This repository contains extended results for the tax examples and CPS examples. However, this repository only contains replication materials for the CPS code for privacy and data access reasons. 

## Table of Contents

## Repository Contents

The `.R` scripts in `analyses/` iterate the tests with different methods, applications, and data. `03` is summary statistics for the CPS, and `04` is regression methods for the CPS. These scripts were mostly run on AWS instances of varying sizes, though the scripts contain non-iterated examples that are commented out. 

* `analysis/03_cps_histograms.R`
* `analyses/03_cps_mean-income.R`
* `analyses/03_cps_mean-earned-income.R`
* `analyses/03_cps_quantiles.R`
* `analyses/04_cps_female-regression.R`
* `analyses/04_cps_male-regression.R`

The iteration scripts write results and summaries of results to the `results/` folder. The following five documents contain documentation about the applications and summarize the results from `results/`:

* `analyses/03_cps-summary-stats.Rmd`
* `analyses/04_cps_female-regression.Rmd`
* `analyses/04_cps_male-regression.Rmd`

The results of those documents can be viewed here (for the version on the `master` branch):

* [03 CPS Summary Statistics](https://urbaninstitute.github.io/formal-privacy-comp-appendix/analyses/03_cps-summary-stats)
* [04 CPS Regression Analysis (Female)](https://urbaninstitute.github.io/formal-privacy-comp-appendix/analyses/04_cps_female-regression)
* [04 CPS Regression Analysis (Male)](https://urbaninstitute.github.io/formal-privacy-comp-appendix/analyses/04_cps_male-regression)

Most of the methods have wrapper functions like `lm_*` and `quantile_*` in `R/`. The wrapper functions provide a consistent interface around methods with varying styles in R and Python. These methods are in `dp-code/`. 

`R/prep_*.R` contain code for setting up the example analyses. 

## Replication

### Environments

* The histogram iteration programs do not use parallel computing. 
* The quantile iteration programs do not use parallel computing because they reference python code. The "smooth" method is very slow. 
* The mean programs use 14 workers and need a computer like an AWS c5.4xlarge.
* The regression programs use 34 workers need and need a computer like an AWS c5.9xlarge.

The programs rely on R packages from CRAN and GitHub, and Python packages from PIP, Conda, and GitHub. `config.sh` sets up the computing environment including the packages and Conda environments. 

## License

## Contact

Please contact [Aaron R. Williams](awilliams@urban.org) with questions. 
