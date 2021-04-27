# R Packages -------------------------------------------------------------------
library("data.table")
library("doParallel")
library("doRNG") # Reproducible parallel foreach loops
library("dplyr") 
library("foreach")
library("glmnet")
library("ggplot2")
library("ggpubr")
library("impute") # From bioconductor
library("knitr")
library("parallel")
library("purrr")
library("rngtools") # Required by dorng
library("rsample")
library("survival")
library("tidyr") # (>= 1.0.0) for pivot_longer() and pivot_wider()
# library("xfun") # We will use xfun::cache_rds() below, but not attach the package here
library("xtable")

# Custom R functions -----------------------------------------------------------
source("R/plot_patient_followup.R")
source("R/impute_genes.R")
source("R/make_xy.R")
source("R/calibrate_sim.R")
source("R/run_sim.R")
source("R/adjust_surv.R")
source("R/tidycoef.R")
source("R/concordance.R")
source("R/calibrate.R")

# Settings ---------------------------------------------------------------------
set.seed(77)
theme_set(theme_bw())
center_title <- function() theme(plot.title = element_text(hjust = 0.5))
N_SIMS <- 200

# Caching
RERUN_CACHE <- TRUE # Set to TRUE to rerun all cached results
if (!dir.exists("cache")){
  dir.create("cache")
}

# Parallel
PARALLEL <- TRUE
if (PARALLEL) {
  cl <- parallel::makeCluster(20, setup_strategy = "sequential", outfile = "simout")
  registerDoParallel(cl)
}

# Model formula ----------------------------------------------------------------
# This formula is used for both the simulations and the real-world data application
f_small <- formula(
  ~ PracticeType + index_date_year + Race + age_at_dx + RE_TP53 +
    CN_TP53 + SV_TP53 + SV_KRAS + SV_EGFR
)
vars_small <- c("Practice type: community", "Index year", "African American",
                "Other race", "White", "Age", "TP53 RE",
                "TP53 CN", "TP53 SV", "KRAS SV", "EGFR SV")

# Followup plot ----------------------------------------------------------------
p_followup <- plot_patient_followup()
ggsave("figs/followup.pdf", p_followup, height = 5, width = 7)