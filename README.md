# Penalized regression for left-truncated and right-censored survival data
This repository contains code for a paper on estimation and evaluation of penalized survival models with high dimensional left-truncated and right-censored (LTRC) survival data. All code for the manuscript is contained in the file [`main.R`](main.R). All functions used in `main.R` are available within the `R` directory.

The analysis was performed with the Flatiron Health and Foundation Medicine Clinico-Genomic Database (CGDB). A single analysis dataset named `data.rds` is read into `main.R`, but we are unfortuntely not permitted to share it. 

R package dependencies are managed through the [`renv`](https://rstudio.github.io/renv/articles/renv.html) package. You can view all packages and their versions in the lockfile [`renv.lock`](renv.lock). All required packages and the appropriate versions can be installed with `renv::restore()`. 

