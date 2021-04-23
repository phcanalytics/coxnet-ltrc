# Calibrate simulation ---------------------------------------------------------
if (file.exists("sim_settings.rds")) { # Interested readers can run this
  sim_settings <- readRDS("sim_settings.rds")
} else { # This setting is run for the paper
  data <- readRDS("data.rds")
  sim_settings <- calibrate_sim(f_small, data = data, store_x = TRUE)
}

# Save cumulative hazard plot
ggsave("figs/sim_calibration_cumhaz.pdf", 
       sim_settings$cumhaz_plot,
       height = 5, width = 10)

# Example simulated data -------------------------------------------------------
# No predictors (intercept only model)
params <- set_params(sim_settings = sim_settings, dist = "weibullPH")
simdata <- sim_survdata(params = params, n_pats = 10000)
simdata_summary_int <- summarize_simdata(simdata, km_rwd = sim_settings$os_km$rc,  
                                         save = TRUE, name = "int")

# Now include predictors (and hence create informative censoring when
# using the Kaplan-Meier estimator). The design matrix X and parameters will 
# be used in the simulations that follow
x_sim <- sim_x(n_pats = 5000, sim_settings, p_bin = 0)
params <- set_params(x_sim, sim_settings, dist = "weibullPH")
simdata <- sim_survdata(x_sim, params)
simdata_summary_p11 <- summarize_simdata(simdata, km_rwd = sim_settings$os_km$rc,
                                         save = TRUE, name = "p11")

# Run simulation for unpenalized Cox model -------------------------------------
sim_coxph_p21 <- xfun::cache_rds({
  run_sim(n_sims = N_SIMS, x = x_sim, params = params, method = "coxph")
}, file = "sim_coxph_p21.rds", rerun = RERUN_CACHE)
sim_coxph_21_summary <- summarize_sim(sim_coxph_p21, save = TRUE,
                                      model_name = "coxph_p21")
rm(sim_coxph_p21)

# Run simulation for lasso model with lambda = 0 and small p -------------------
# run_sim1(simdata, method = "coxnet", lambda = 0) # For debugging
sim_coxlasso_lambda0_p21 <- xfun::cache_rds({
  run_sim(n_sims = N_SIMS, x = x_sim, params = params, 
          method = "coxnet", lambda = c(1, 0))
}, file = "sim_coxlasso_lambda0_p21.rds", rerun = RERUN_CACHE)
coxlasso_lambda0_p21_summary <- summarize_sim(
  sim_coxlasso_lambda0_p21, save = TRUE, 
  model_name = "coxlasso_lambda0_p21"
)
rm(sim_coxlasso_lambda0_p21)

# Run simulation for lasso model with small p ----------------------------------
sim_coxlasso_p21 <- xfun::cache_rds({
  run_sim(n_sims = N_SIMS, x = x_sim, params = params, 
          method = "coxnet")
}, file = "sim_coxlasso_p21.rds", rerun = RERUN_CACHE)
sim_coxlasso_p21_summary <- summarize_sim(sim_coxlasso_p21, save = TRUE, 
                                          model_name = "coxlasso_p21")
rm(sim_coxlasso_p21)

# Run simulation for ridge model with small p ----------------------------------
sim_coxridge_p21 <- xfun::cache_rds({
  run_sim(n_sims = N_SIMS, x = x_sim, params = params, 
          method = "coxnet", alpha = 0)
}, file = "sim_coxridge_p21.rds", rerun = RERUN_CACHE)
sim_coxridge_p21_summary <- summarize_sim(sim_coxridge_p21, save = TRUE, 
                                          model_name = "coxridge_p21")
rm(sim_coxridge_p21)

# Run simulation for lasso model with large p ----------------------------------
x_sim <- sim_x(n_pats = 5000, sim_settings, p_bin = 1000)
params <- set_params(x_sim, sim_settings, dist = "weibullPH")

sim_coxlasso_p1011 <- xfun::cache_rds({
  run_sim(n_sims = N_SIMS, x = x_sim, params = params, 
          method = "coxnet")
}, file = "sim_coxlasso_p1011.rds", rerun = RERUN_CACHE)
sim_coxlasso_p1011_summary <- summarize_sim(sim_coxlasso_p1011, save = TRUE, 
                                            model_name = "coxlasso_p1011")
rm(sim_coxlasso_p1011)

# Calibration plot comparing small and large p simulation ----------------------
sim_coxlasso_calplot <- ggarrange(
  sim_coxlasso_p21_summary$calibration_plots$complete +
    ggtitle("(A) Small p") + center_title(),
  sim_coxlasso_p1011_summary$calibration_plots$complete +
    ggtitle("(B) Large p") + center_title(),
  nrow = 2, common.legend = TRUE, legend = "bottom"
)
ggsave("figs/sim_coxlasso_calibration_complete.pdf", sim_coxlasso_calplot, 
       width = 7, height = 9)