# Formal Privacy Feasibility Study

This repository contains code and data for a feasibility study of formally private methods for histograms, means, quantiles, and linear regression by Andres Felipe Barrientos, Aaron R. Williams, Joshua Snoke, and Claire McKay Bowen. 

## Table of Contents

* [Background](#background)
* [Repository Contents](#repository-contents)
* [Data](#data)
* [Environments](#environments)
* [Contact](#contact)

## Background

This repository contains code and data for A Feasibility Study of Differentially Private Summary Statistics and Regression Analyses with Evaluations on Administrative and Survey Data ([JASA](https://arxiv.org/abs/2110.12055](https://www.tandfonline.com/doi/full/10.1080/01621459.2023.2270795))). The paper evaluates the results of several differentially private methods for a range of statistics on two case studies. 

One case study uses confidential tax microdata that is not available for reproduction. The second case study uses public microdata from the 1994 to 1996 Current Population Survey Annual Social and Economic Supplements (CPS ASEC). 

## Repository Contents

The `.R` scripts in `analyses/` iterate the tests with different methods, applications, and data. `01` is summary statistics for the SOI, `02` is regression methods for the SOI, `03` is summary statistics for the CPS, and `04` is regression methods for the CPS. These scripts were mostly run on AWS instances of varying sizes, though the scripts contain non-iterated examples that are commented out. 

* `analyses/01_soi_histograms.R`
* `analyses/01_soi_means.R`
* `analyses/01_soi_quantiles.R`
* `analyses/02_soi_regression.R`
* `analysis/03_cps_histograms.R`
* `analyses/03_cps_mean-income.R`
* `analyses/03_cps_mean-earned-income.R`
* `analyses/03_cps_quantiles.R`
* `analyses/04_cps_female-regression.R`
* `analyses/04_cps_male-regression.R`

The iteration scripts write results and summaries of results to the `results/` folder. The following five documents contain documentation about the applications and summarize the results from `results/`:

* `analyses/01_soi-puf-summary-stats.Rmd`
* `analyses/02_soi-puf-regression.Rmd`
* `analyses/03_cps-summary-stats.Rmd`
* `analyses/04_cps_female-regression.Rmd`
* `analyses/04_cps_male-regression.Rmd`

The results of those documents can be viewed here (for the version on the `master` branch):

* [01 PUF Summary Statistics](https://ui-research.github.io/formal-privacy-comp/analyses/01_soi-puf-summary-stats.html)
* [02 PUF Regression Analysis](https://ui-research.github.io/formal-privacy-comp/analyses/02_soi-puf-regression.html)
* [03 CPS Summary Statistics](https://ui-research.github.io/formal-privacy-comp/analyses/03_cps-summary-stats.html)
* [04 CPS Regression Analysis (Female)](https://ui-research.github.io/formal-privacy-comp/analyses/04_cps_female-regression.html)
* [04 CPS Regression Analysis (Male)](https://ui-research.github.io/formal-privacy-comp/analyses/04_cps_male-regression.html)

Most of the methods have wrapper functions like `lm_*` and `quantile_*` in `R/`. The wrapper functions provide a consistent interface around methods with varying styles in R and Python. These methods are in `dp-code/`. 

`tax-code/` and `R/prep_*.R` contain code for setting up the example analyses. 

## Data

The CPS example uses data pulled from [IPUMS CPS](https://cps.ipums.org/cps/). 

> Sarah Flood, Miriam King, Renae Rodgers, Steven Ruggles, J. Robert Warren and Michael Westberry. Integrated Public Use Microdata Series, Current Population Survey: Version 10.0 [dataset]. Minneapolis, MN: IPUMS, 2022. https://doi.org/10.18128/D030.V10.0

The data are loaded using `library(ipumsr)` and the functions `read_ipums_ddi()` and `read_ipums_micro()`. These functions use the `*.xml` and `*.dat.gz` files in the `data/` directory to load the CPS microdata into R. 

`cps_00008.xml` is used for quantiles and means. `cps_00009.xml` is used for histograms. `cps_00012.xml` and `cps_00013.xml` are used for linear regressions. 

## Environments

* The histogram iteration programs do not use parallel computing. 
* The quantile iteration programs do not use parallel computing because they reference python code. The "smooth" method is very slow. 
* The mean programs use 14 workers and need a computer like an AWS c5.4xlarge.
* The regression programs use 34 workers and need a computer like an AWS c5.9xlarge.

The programs rely on R packages from CRAN and GitHub, and Python packages from PIP, Conda, and GitHub. `config.sh` sets up the computing environment including the packages and Conda environments. `library(synthpuf)` is cloned by `config.sh` but needs to be built by opening the package in RStudio and building. 

## Contact

Please contact [Aaron R. Williams](awilliams@urban.org) with questions. 
