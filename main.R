# Setup ------------------------------------------------------------------------
rm(list = ls())

# Packages
library("data.table")
library("doParallel")
library("dplyr") 
library("foreach")
library("ggplot2")
library("ggpubr")
library("impute") # From bioconductor
library("knitr")
library("parallel")
library("purrr")
library("rsample")
library("survival")
library("tidyr") # (>= 1.0.0) for pivot_longer() and pivot_wider()
library("xfun") # (>= 0.13.0) for cache_rds()
library("xtable")

# Sourced R files
source("R/plot_patient_followup.R")
source("R/impute_genes.R")
source("R/make_xy.R")
source("R/calibrate_sim.R")
source("R/run_sim.R")

# Settings
set.seed(77)
theme_set(theme_bw())
center_title <- function() theme(plot.title = element_text(hjust = 0.5))
N_SIMS <- 40

# Caching
RERUN_CACHE <- TRUE # Set to TRUE to rerun all cached results
if (!dir.exists("cache")){
  dir.create("cache")
}

# Parallel
PARALLEL <- TRUE
if (PARALLEL) {
  cl <- parallel::makeCluster(20, setup_strategy = "sequential")
  registerDoParallel(cl)
}

# Store numbers to use in text of model documentation
txt <- list() # List for in-line text statistics

# Followup plot ----------------------------------------------------------------
p_followup <- plot_patient_followup()
ggsave("figs/followup.pdf", p_followup, height = 5, width = 7)

# Load data --------------------------------------------------------------------
data <- readRDS("data.rds")

# Split data -------------------------------------------------------------------
# Using rsample package
data_split <- initial_split(data, prop = .75)
train_data <- training(data_split)
test_data <- testing(data_split)

# Clean up
rm(data_split)

# Impute genomic data ----------------------------------------------------------
train_test_imputed <- impute_genes(train_data, test_data)
train_data <- train_test_imputed$train
test_data <- train_test_imputed$test
selected_genes <- train_test_imputed$selected_genes
  
# Clean up
rm(train_test_imputed)

# Model formulas ---------------------------------------------------------------
# Small p formula
f_small <- formula(
  ~ PracticeType + index_date_year + Race + age_at_dx + RE_TP53 +
    CN_TP53 + SV_TP53 + SV_KRAS + SV_EGFR
)
vars_small <- c("Practice type: community", "Index year", "African American",
                "Other race", "White", "Age", "TP53 RE",
                "TP53 CN", "TP53 SV", "KRAS SV", "EGFR SV")

# Big p formula 
f_big <- as.formula(paste0(
  "~ PracticeType + index_date_year + Race + age_at_dx + `",
  paste(selected_genes ,collapse = "`+`"),
  "`")
)
vars_big <- c(vars_small[1:6], 
              gsub("_", ":", selected_genes))

# XY test/train data -----------------------------------------------------------
# Small p
train_small <- make_xy(train_data, f_small)
test_small <- make_xy(test_data, f_small)

# Big p
train_big <- make_xy(train_data, f_big)
test_big <- make_xy(test_data, f_big)

# Clean up
rm(train_data)
rm(test_data)

# Calibrate simulation ---------------------------------------------------------
sim_settings <- calibrate_sim(f_small, data = data)

# Cumulative hazard plot comparing different parametric models
# with Kaplan-Meier estimator
ggsave("figs/sim_calibration_cumhaz.pdf", 
       sim_settings$os_comparisons$cumhaz_plot,
       height = 5, width = 7)

# Example simulated data -------------------------------------------------------
x_sim <- sim_x(n_pats = 5000, sim_settings, p_bin = 10)
params <- set_params(x_sim, sim_settings, dist = "weibullPH")
simdata <- sim_survdata(x_sim, params)
simdata_summary <- summarize_simdata(simdata, save = TRUE, name = "p10")

# Distribution of entry times --------------------------------------------------
p_left_trunc <- ggplot(data, aes(x = entry_days_dx)) +
  geom_histogram(binwidth = 60, colour = "white") + 
  xlab("Days between diagnosis and FMI test") + 
  ylab("Number of patients") +
  coord_cartesian(xlim = c(0, 1500))
ggsave("figs/left_trunc_hist.pdf", p_left_trunc, height = 5, width = 7)

# Number of deaths -------------------------------------------------------------
n_deaths <- sum(train_big$y[, "status"])
max_p <- n_deaths/15 # Rule of thumb for max number of predictors in Cox model
