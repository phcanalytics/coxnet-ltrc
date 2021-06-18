library("data.table")
source("sim.R") # Contains functions used for the simulation

# Simulate both predictors and survival times based on the CGDB
sim_settings <- readRDS("sim_settings.rds")
x_sim <- sim_x(n_pats = 1000, sim_settings, p_bin = 0)
params <- set_params(x_sim, sim_settings, dist = "weibullPH")
simdata <- sim_survdata(X = x_sim, params = params, n_pats = 10000)

# Save
simdata[, f_hidden := NULL]
saveRDS(simdata, "simdata.rds")