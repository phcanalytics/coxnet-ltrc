# Penalized regression for left-truncated and right-censored survival data
This repository contains code for a [paper](https://www.medrxiv.org/content/10.1101/2021.02.09.21251373v1) on estimation and evaluation of penalized survival models with high dimensional left-truncated and right-censored (LTRC) survival data. All code for the manuscript can be executed by running the file [`main.R`](main.R). All functions required are available within the `R` directory. An example analysis is also provided as described [below](#Example).

## Data availability 
The analysis was performed with the Flatiron Health and Foundation Medicine Clinico-Genomic Database (CGDB). A single analysis dataset named `data.rds` is required for the real-world data application, but we are unfortuntely not permitted to share it. 

## Running the simulations
Although we are not permitted to share the data, the simulations can still be run by using the simulation settings stored in `sim_settings.rds` and runing the following `R` code:

```{r}
source("01-setup.R")
source("02-simulation.R")
```

Note, however, that these simulations are slightly different than what is reported in the paper because we report results from simulations in which the the features of the input matrix based on the CGDB data were randomly drawn from replacement from the raw data; in the executable version of the simulations, the (standardized) CGDB based features are drawn from a multivariate normal distribution. Furthermore, parallel computing is used to speed up the simulations, so the parallelization settings specified in `01-setup.R` may need to be modified for a users computing environment.

## Example
Since the simulations can be difficult to follow, we also have an example `R` Markdown file (HTML output [here](https://phcanalytics.github.io/coxnet-ltrc/example.html)) showing how to use penalized survival models for prediction with left-truncated and right-censored data (LTRC). Note that the analysis is based on simulated dataset `simdata.rds` created using the `simdata.R` script.  

## Dependencies
R package dependencies are managed through the [`renv`](https://rstudio.github.io/renv/articles/renv.html) package. You can view all packages and their versions in the lockfile [`renv.lock`](renv.lock). All required packages and the appropriate versions can be installed with `renv::restore()`. 
