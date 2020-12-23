# Simulate design matrix -------------------------------------------------------
sim_vars <- function(n_pats, K){
  
  # Build correlation matrix
  rmat <- matrix(stats::rnorm(K * K), K, K)
  cov_mat <- rmat %*% t(rmat)
  corr_mat <- cov_mat/sqrt(diag(cov_mat) %*% t(diag(cov_mat))) # e ~ N(0, 1)
  
  # Simulate correlated binomial data using latent probit model
  # ystar = xB + e
  probs <- stats::runif(K, .2, .8) 
  sim <- mvtnorm::rmvnorm(n_pats, stats::qnorm(probs), sigma = corr_mat)
  sim <- ifelse(sim >= 0, 1, 0)
  
  # Return
  colnames(sim) <- paste0("binvar_", 1:K)
  attr(sim, "type") <- rep("Binary", K)
  return(sim)
}

#' @export
sim_x <- function(n_pats = 1000, sim_settings, p_bin) {
  # RWD
  rows_to_sample <- sample.int(n = nrow(sim_settings$xy$x), 
                               size = n_pats, replace = TRUE)
  X <- sim_settings$xy$x[rows_to_sample, ]
  n_rwdvars <- ncol(X)
  rwd_type <- rep("RWD", ncol(X))
  attr(X, "type") <- rwd_type
  
  # Simulated binary variables
  if (p_bin > 0){
    X_sim <- sim_vars(n_pats, K = p_bin)
    X <- cbind(X, X_sim)
    attr(X, "type") <- c(rwd_type, attr(X_sim, "type"))
  } 
  return(X)
}

# Get parameters of model given design matrix ----------------------------------
set_beta <- function(X, sim_settings, dist) {
  
  basefit <- sim_settings$os_base[[dist]]
  coxfit = sim_settings$os_cox
  
  type <- attr(X, "type")
  beta <- rep(NA, length(type))
  names(beta) <- colnames(X)
  
  # Coefficients for RWD variables
  rwd_index <- which(type == "RWD")
  name <- basefit$dlist$location
  pos <- which(basefit$dlist$pars == name)
  beta[rwd_index] <- stats::coef(coxfit)
  
  ## Binary variables
  which_bin <- which(type == "Binary")
  p_bin <- length(which_bin)
  if (p_bin > 0){
    is_zero_beta_bin <- purrr::rbernoulli(p_bin, .5)
    beta_bin <- ifelse(is_zero_beta_bin,
                       0,
                       stats::runif(p_bin, -.25, .25))
    beta[which_bin] <- beta_bin
  }
  return(beta)
}

#' @export
set_params <- function(X = NULL, sim_settings,
                       dist = c("weibullPH", "survspline", "exp")) {
  dist <- match.arg(dist)
  if (!is.null(X)) beta <- set_beta(X, sim_settings, dist)
  
  list(
    beta = beta,
    os_base = sim_settings$os_base[[dist]],
    rc_base = sim_settings$rc_base[[dist]] 
  )
}  

# Simulate survival data from parameters and design matrix ---------------------
predict_location <- function(intercept, beta = NULL, X = NULL, inv_transforms){
  if (!is.null(X)) {
    xb <- X %*% as.matrix(beta)
  } else{
    xb <- 0
  }
  lp <- intercept + xb - mean(xb)
  res <- inv_transforms(lp)
  return(res)
}

surv_rng <- function(n, beta, basefit, X){
  
  # Location parameter
  loc_pos <- which(basefit$dlist$pars == basefit$dlist$location)
  loc_name <- basefit$dlist$pars[loc_pos]
  loc <- predict_location(intercept = c(basefit$res.t[loc_name, "est"]), 
                          beta, X = X, 
                          inv_transforms = basefit$dlist$inv.transforms[[loc_pos]])
  loc <- list(loc)
  names(loc) <- loc_name
  
  
  # Ancillary parameters
  if (length(basefit$dlist$pars) > 1){
    anc_pos <- which(basefit$dlist$pars != basefit$dlist$location)
    anc_names <- basefit$dlist$pars[anc_pos]
    n_anc <- length(anc_names)
    anc <- vector(mode = "list", length = n_anc)
    names(anc) <- anc_names
    for (j in 1:n_anc){
      inv_transforms <- basefit$dlist$inv.transforms[[anc_pos[j]]]
      anc[[j]] <- inv_transforms(basefit$res.t[anc_names[j], "est"])
    }
    pars <- c(loc, anc)
  } else{
    pars <- loc
  }
  
  # Simulate
  if (basefit$dlist$name == "survspline"){
    anc_mat <- matrix(unlist(anc), ncol = length(anc), 
                      nrow = n, byrow = TRUE)
    gamma <- cbind(loc[[1]], anc_mat)
    res <- rsurvspline(n = n, gamma = gamma, knots = basefit[["knots"]])
    
  } else{
    random_fun <- get(paste0("r", basefit$dlist$name) , asNamespace("flexsurv"))
    res <- do.call(random_fun, c(list(n = n), pars))
  }
  return(res)
}

entry_time_rng <- function(n, p_zero = .2, pos_mean = 1.6, pos_median = 1){
  if (pos_mean < pos_median){
    stop("'pos_mean' must be greater than or equal to 'pos_median'.")
  }
  # Two part model
  ## Zeros
  zero <- purrr::rbernoulli(n, p_zero) 
  
  ## Positive
  pos_mu <- log(pos_median)
  pos_sigma <- sqrt(2 * (log(pos_mean) - pos_mu)) 
  pos_time <- stats::rlnorm(n, meanlog = pos_mu, sdlog = pos_sigma)
  
  ## Combine
  time <- ifelse(zero == 1, 0, pos_time)
  return(time)
}

#' @export
sim_survdata <- function(X = NULL, params, n_pats = 5000){
  if (!is.null(X)) {
    n_pats <- nrow(X)
  }
  
  # Latent survival
  death_time <- surv_rng(n = n_pats, beta = params$beta,
                         basefit = params$os_base,
                         X = X)/365.25
  death_time <- pmin(death_time, 100) # Maximum survival
  
  # Right censoring time
  censored_time <- surv_rng(n = n_pats,  beta = params$beta,
                            basefit = params$rc_base, 
                            X = X)/365.25
  
  # Observed survival time
  event_time <- pmin(censored_time, death_time)
  dead <- ifelse(event_time < censored_time, 1, 0)
  
  # Entry time
  entry_time <- entry_time_rng(n = n_pats)
  hidden <- 1 * (entry_time > event_time)
  
  # Create a dataset of survival times
  data <- data.table(patient_id = 1:n_pats,
                     death_time = death_time, 
                     entry_time = entry_time,
                     event_time = event_time,
                     dead = dead,
                     hidden = hidden, 
                     X)
  data[, event_time := ifelse(hidden == 0 & event_time - entry_time < .01,
                              entry_time + .01,
                              event_time)] # Fix potential computational errors
  data[, f_hidden := factor(hidden, 
                            labels = c("No", "Yes"),
                            levels = c(0, 1))]
  setattr(data, "x_vars", colnames(X))
  setattr(data, "start", "entry_time")
  setattr(data, "stop", "event_time")
  setattr(data, "status", "dead")
  return(data[, ])
}

# Summarize simulated survival data --------------------------------------------
plot_latent_surv <- function(data){
  ggplot(data, aes(x = death_time)) +
    geom_histogram(binwidth = 1/12, color = "black", fill = "white") +
    xlab("True survival time (in days)") +
    ylab("Count")
}

plot_latent_entry <- function(data, stratify = FALSE){
  if (!stratify){
    p <- ggplot(data = data, aes(x = entry_time)) +
      geom_histogram(binwidth = 1/12, color = "white", fill = "black", 
                     alpha = .8) +
      xlab("Entry time") + ylab("Count") +
      coord_cartesian(xlim = c(0, 10))
  } else{
    p <- ggplot() +
      geom_histogram(data = data[f_hidden == "No"],
                     binwidth = 1/12, color = "white", alpha = .8,
                     aes(x = entry_time, fill = "No")) +
      geom_histogram(data = data[f_hidden == "Yes"],
                     binwidth = 1/12, color = "white", alpha = .8,
                     aes(x = entry_time, y = -..count.., fill = "Yes")) +
      xlab("Entry times") +
      ylab("Count") +
      scale_fill_discrete("Left truncated") +
      scale_y_continuous(labels = function(x) abs(x)) +
      theme(legend.position = "bottom") +
      coord_cartesian(xlim = c(0, 10))
  }
  return(p)
}

fit_km <- function(data){
  kmfits <- list(
    # True survival
    list(fit = survival::survfit(survival::Surv(death_time) ~ 1, data = data),
         label = "True survival",
         est_sample = "Complete"),
    
    # Right censored adjusted survival
    ## Complete sample (pretty close to truth)
    list(fit = survival::survfit(survival::Surv(event_time, dead) ~ 1, 
                                 data = data),
         label = "Right censored adjusted survival",
         est_sample = "Complete"),
    
    ## Observed sample (biased upwards)
    list(fit = survival::survfit(survival::Surv(event_time, dead) ~ 1, 
                                 data = data[hidden == 0]),
         label = "Right censored adjusted survival",
         est_sample = "Observed"),
    
    # Left-truncated and right-censored adjusted survival (pretty close to truth)
    list(fit = survival::survfit(survival::Surv(entry_time, event_time, dead) ~ 1, 
                                 data = data[hidden == 0]),
         label = "Left truncated and right censored adjusted survival",
         est_sample = "Observed")
  )
  return(kmfits)
}

plot_km <- function(fits){
  n_fits <- length(fits)
  surv_list <- vector(mode = "list", length = n_fits)
  for (i in 1:n_fits){
    surv_list[[i]] <- data.table(label = fits[[i]]$label,
                                 est_sample = fits[[i]]$est_sample,
                                 time = fits[[i]]$fit$time,
                                 prob = fits[[i]]$fit$surv)
  }
  surv_df <- rbindlist(surv_list)
  ggplot(surv_df,
         aes(x = time, y = prob, col = label, linetype = est_sample)) + 
    geom_line() + 
    scale_linetype_discrete(name = "Estimation sample") +
    scale_color_manual(name = "Estimate",
                       values = c("True survival" = "black",
                                  "Right censored adjusted survival" = "red",
                                  "Left truncated and right censored adjusted survival" = "blue")) +
    xlab("Time") + ylab("Survival probability") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(ncol = 1),
           linetype = guide_legend(ncol = 1))
}

summarize_km <- function(fits, probs = c(.1, .25, .5, .75, .9)){
  n_fits <- length(fits)
  q_mat <- matrix(NA, nrow = n_fits, ncol = length(probs))
  label <- est_sample <- rep(NA, n_fits)
  colnames(q_mat) = paste0("q", probs)
  for (i in 1:n_fits){
    label[i] <- fits[[i]]$label
    est_sample[i] <- fits[[i]]$est_sample
    q_mat[i, ] <- stats::quantile(fits[[i]]$fit, probs)$quantile
  }
  res <- data.frame(Estimate = label,
                    "Estimation sample" = est_sample,
                    q_mat,
                    check.names = FALSE)
  return(res)
}

#' @export
summarize_simdata <- function(data, save = FALSE, ...) {
  # Descriptive statistics for simulated survival data
  latent_surv_plot <- plot_latent_surv(data) + coord_cartesian(xlim = c(0, 10))
  prop_rc <- 1 - mean(data$dead) 
  prop_truncated <-  mean(data$hidden)
  latent_entry_plot <- plot_latent_entry(data, stratify = FALSE)
  stratified_entry_plot <- plot_latent_entry(data, stratify = TRUE)
  
  # Kaplan-Meier estimates
  km_fits <- fit_km(data)
  km_plot <- plot_km(km_fits) + coord_cartesian(xlim = c(0, 10))
  km_quantiles <- summarize_km(km_fits)
  
  # Return and optionally save
  res <- list(
    latent_surv_plot = latent_surv_plot,
    prop_rc = prop_rc,
    prop_truncated = prop_truncated,
    latent_entry_plot = latent_entry_plot,
    stratified_entry_plot = stratified_entry_plot,
    km_plot = km_plot,
    km_quantiles = km_quantiles
  )
  
  if (save) save_simdata(res, ...)
  return(res)
}

save_simdata <- function(object, name = NULL) {
  if (!is.null(name)) {
    name <- paste0(name, "_")
  }
  fig_path <- paste0("figs/sim_", name)
  tbl_path <- paste0("tables/sim_", name)
  output_path <- paste0("output/sim_", name)
  
  # Save
  write.csv(object$prop_rc, paste0(tbl_path, "prop_rc.csv"))
  write.csv(object$prop_truncated, paste0(tbl_path, "prop_truncated.csv"))
  ggsave(paste0(fig_path, "km_plot.pdf"), 
         object$km_plot, 
         height = 5, width = 7)
  write.csv(object$km_quantiles, paste0(tbl_path, "km_quantiles.csv"))
  saveRDS(object$data, paste0(output_path, "data.rds"))
}

# Run Monte Carlo simulation to test methods -----------------------------------
make_train_test <- function(data){
  data_split <- rsample::initial_split(data, prop = .75)
  train_data <- rsample::training(data_split)[hidden == 0]
  test_data_complete <- rsample::testing(data_split) 
  test_data_observed <- test_data_complete[hidden == 0]  
  
  make_xy <- function(data, x_vars, start, stop, status){
    x <- as.matrix(data[, x_vars, with = FALSE])
    y <- data.frame(start = data[[start]],
                    stop = data[[stop]],
                    status = data[[status]])
    return(list(x = x, y = as.matrix(y),
                patient_id = data$patient_id))
  }
  x_vars <- attr(data, "x_vars")
  start <- attr(data, "start")
  stop <- attr(data, "stop")
  status <- attr(data, "status")
  train <- make_xy(train_data, x_vars, start, stop, status)
  test_complete <- make_xy(test_data_complete, x_vars, start, stop, status)
  test_observed <- make_xy(test_data_observed, x_vars, start, stop, status)
  attr(test_complete, "sample") <- "complete"
  attr(test_observed, "sample") <- "observed"
  attr(train, "split") <- "train"
  attr(test_complete, "split") <- "test"
  attr(test_observed, "split") <- "test"
  return(list(train = train, test_complete = test_complete,
              test_observed = test_observed))
}

fit_model <- function(train, left_truncation,
                      method = c("coxph", "coxnet"),
                      parallel,alpha, lambda){
  if (is.null(alpha)) alpha <- 1
  method <- match.arg(method)
  
  # Fit
  if (method == "coxph"){ # coxph()
    if (left_truncation == "Yes"){
      train_df <- data.frame(cbind(train$x, train$y))
      fit <- survival::coxph(survival::Surv(start, stop, status) ~ ., 
                             data = train_df)
    } else{
      train_df <- data.frame(cbind(train$x, train$y[, -1]))
      fit <- survival::coxph(survival::Surv(stop, status) ~ ., data = train_df)
    }
  } else{ # cv.glmnet()
    if (is.null(lambda)) {
      n_folds <- 10
    } else{
      n_folds <- 3
    }
    if (left_truncation == "Yes"){
      train_y <- survival::Surv(time = train$y[, 1], time2 = train$y[, 2], 
                                event = train$y[, 3])
      fit <- glmnet::cv.glmnet(x = train$x, y = train_y, 
                               standardize = FALSE, parallel = parallel,
                               alpha = alpha, lambda = lambda, 
                               nfolds = n_folds, family = "cox")
    } else{
      train_y <- adjust_Surv(train$y, left_trunc = FALSE)
      fit <- glmnet::cv.glmnet(x = train$x, y = train_y, 
                               standardize = FALSE, parallel = parallel,
                               alpha = alpha, lambda = lambda,
                               nfolds = n_folds, family = "cox")
    }
    fit$x <- train$x; fit$y <- train_y
    if (!is.null(lambda)){
      fit$lambda.min <- lambda[length(lambda)] 
    }
  }
  return(fit)
}

fit_models <- function(train, x_vars,
                       method = c("coxnet", "coxph"),
                       parallel, alpha, lambda){
  fit_std <- fit_model(train = train, left_truncation = "No",
                       method = method,
                       parallel, alpha = alpha, lambda = lambda)
  fit_lt <- fit_model(train = train, left_truncation = "Yes",
                      method = method,
                      parallel, alpha = alpha, lambda = lambda)
  return(dplyr::tibble(
    left_truncation = c("No", "Yes"),
    fit = list(fit_std, fit_lt)
  ))  
}

make_coef_tbl <- function(beta, fits) {
  if (inherits(fits$fit[[1]], "coxph")) {
    betahat <- sapply(fits$fit, stats::coef)
  } else{
    betahat <- sapply(fits$fit, function(x) stats::coef(x, s = "lambda.min")[, 1])
  }
  tbl <- data.table(
    Variable = names(beta),
    Truth = beta, 
    `RC-adjusted` = betahat[, 1],
    `LTRC-adjusted` = betahat[, 2]
  )
  return(tbl)
}

predict_xb <- function(fit, new_xy, x_vars = NULL){
  if (inherits(fit, "coxph")){
    newdata <- data.frame(cbind(new_xy$x, new_xy$y))
    return(stats::predict(fit, newdata))
    
  } else{
    return(stats::predict(fit, new_xy$x, s = "lambda.min"))
  }
}

predict_surv <- function(fit, new_xy){
  if (inherits(fit, "coxph")){ 
    survival::survfit(fit, newdata = data.frame(new_xy$x), se.fit = FALSE)
  } else {
    survival::survfit(fit, newx = new_xy$x, s = "lambda.min", 
                      x = fit$x, y = fit$y)
  }
}

compute_concordance <- function(xb, new_xy, left_truncation){
  newdata <- data.frame(cbind(new_xy$x, new_xy$y))
  test_sample <- attr(new_xy, "sample")
  if (left_truncation == "Yes" & test_sample == "observed") {
    f <- survival::Surv(start, stop, status) ~ xb
  } else{
    f <- survival::Surv(stop, status) ~ xb
  }  
  return(survival::concordance(f, newdata, reverse = TRUE))
}

eval_model_scenarios <- function(fits, train_test){
  scenarios <- dplyr::tibble(
    left_truncation = rep(c("No", "Yes"), 2),
    test_sample = rep(c("Observed", "Complete"), each = 2),
    concordance = vector(mode = "list", length = 4),
    survprobs_train = vector(mode = "list", length = 4),
    survprobs_test = vector(mode = "list", length = 4)
  )
  for (i in 1:nrow(scenarios)){
    # Arguments
    ## Fit
    fit_i <- fits %>% 
      dplyr::filter(left_truncation == scenarios$left_truncation[i]) %>%
      dplyr::pull(fit) %>% .[[1]]
    
    ## Test set
    if (scenarios[i, "test_sample", drop = TRUE] == "Observed"){
      test_i <- train_test$test_observed
    } else{
      test_i <- train_test$test_complete
    }
    
    ## Left truncation
    lt_i <- scenarios$left_truncation[i]
    
    # Evaluate
    xb <- predict_xb(fit_i, new_xy = test_i)
    scenarios$survprobs_train[[i]] <- predict_surv(fit_i, new_xy = train_test$train)
    scenarios$survprobs_test[[i]] <- predict_surv(fit_i, new_xy = test_i)
    scenarios$concordance[[i]] <- compute_concordance(xb, test_i, lt_i)
  }
  
  # Return
  return(scenarios)
}

make_cindex_tbl <- function(scenarios){
  x <- cbind(scenarios[, c("left_truncation", "test_sample")],
             concordance = sapply(scenarios$concordance, 
                                  function (y) y$concordance))
  colnames(x) <- c("Left truncation adjustment", "Test sample",
                   "C-index")
  return(x)
}

make_quantile_tbl <- function(scenarios){
  
  quantile_fun <- function(x, probs = c(.1, .25, .5, .75, .9)){
    q <- stats::quantile(x, probs, conf.int = FALSE)
    q_mean <- apply(q, 2, mean, na.rm = TRUE)
    names(q_mean) <- paste0("q", probs)
    return(q_mean)
  } 
  x <- cbind(scenarios[, c("left_truncation", "test_sample")],
             t(sapply(scenarios$survprobs_test, quantile_fun)))
  return(x)
}

make_surv_tbl <- function(scenarios, train_test, 
                          data_split = c("Test", "Train")){
  
  surv_tbl_fun <- function(x, train_test, data_split){
    if (data_split == "Train"){
      survprobs <- x$survprobs_train
      newdata <- train_test$train
    } else{
      survprobs <- x$survprobs_test
      if (x$test_sample == "Complete"){
        newdata <- train_test$test_complete
      } else{
        newdata <- train_test$test_observed
      }
    }
    patient_id <- newdata$patient_id
    n_patients <- ncol(survprobs$surv)
    t <- seq(0, 10, .01)
    survfit_summary <- summary(survprobs, times = t)
    n_times <- length(survfit_summary$time)
    res <- data.table(test_sample = x$test_sample,
                      left_truncation = x$left_truncation,
                      patient_id = rep(patient_id, each = n_times),
                      time = rep(survfit_summary$time, n_patients),
                      probs = c(survfit_summary$surv)
    )
    res <- cbind(res, newdata$y[rep(1:nrow(newdata$y), each = n_times), ])
    return(res)
  }
  data_split <- match.arg(data_split)
  x <- rbindlist(apply(scenarios, 1, surv_tbl_fun, 
                       train_test = train_test, data_split = data_split))
  return(x)
} 

calibrate_simsurv <- function(scenarios, train_test){
  u <- c(.5, 1, 1.5, 2, 2.5, 3.0)
  n_groups <- 10
  
  # Complete sample
  y_complete <- survival::Surv(time = train_test$test_complete$y[, "stop"],
                               event = train_test$test_complete$y[, "status"])
  scenarios_complete <- scenarios %>% dplyr::filter(test_sample == "Complete")
  cal_complete <- calibrate(scenarios_complete %>% dplyr::pull(survprobs_test),
                            y = y_complete, times = u, n_groups = n_groups,
                            group_df = scenarios_complete %>% 
                              dplyr::select(left_truncation))
  
  # Observed data
  ## No truncation adjustment
  y_obs <- survival::Surv(time = train_test$test_observed$y[, "stop"],
                          event = train_test$test_observed$y[, "status"])
  scenarios_observed <- scenarios %>% dplyr::filter(test_sample == "Observed")
  cal_observed <-calibrate(scenarios_observed %>% dplyr::pull(survprobs_test),
                           y = y_obs, times = u, n_groups = n_groups,
                           group_df = scenarios_observed %>% 
                             dplyr::select(left_truncation))
  
  ## Left truncation adjustment
  y_obs_lt <- survival::Surv(time = train_test$test_observed$y[, "start"],
                             time2 = train_test$test_observed$y[, "stop"],
                             event = train_test$test_observed$y[, "status"])
  cal_observed_lt <- calibrate(scenarios_observed %>% dplyr::pull(survprobs_test),
                               y = y_obs_lt, times = u, n_groups = n_groups,
                               group_df = scenarios_observed %>% 
                                 dplyr::select(left_truncation))
  
  
  # Return
  res <- rbind(
    data.table(cal_complete)[, obs_method := "complete"],
    data.table(cal_observed)[, obs_method := "observed"],
    data.table(cal_observed_lt)[, obs_method := "observed_lt"]
  )
  return(res)
}

plot_calibration_sim <- function() {
  calplot_complete <- autoplot(cal_complete, col = "left_truncation") +
    scale_colour_discrete("Left truncation adjustment") +
    theme(legend.position = "bottom") 
  
  calplot_observed <- autoplot(cal_observed, col = "left_truncation") +
    scale_colour_discrete("Left truncation adjustment") +
    theme(legend.position = "bottom") 
  
  calplot_observed_lt <- autoplot(cal_observed_lt, col = "left_truncation") +
    scale_colour_discrete("Left truncation adjustment") +
    theme(legend.position = "bottom") 
}

run_sim1 <- function(data, method = "coxph", parallel = FALSE,
                     alpha = 1, lambda = NULL) {
  # Modeling
  train_test <- make_train_test(data)
  fits <- fit_models(train = train_test$train, method = method,
                     parallel = parallel, alpha = alpha, 
                     lambda = lambda)
  scenarios <- eval_model_scenarios(fits, train_test)
  cindex <- make_cindex_tbl(scenarios)
  quantiles <- make_quantile_tbl(scenarios)
  surv <- make_surv_tbl(scenarios, train_test)
  
  ## Survival calibration
  calibration <- calibrate_simsurv(scenarios, train_test)
  
  # Return
  res <- list(data = data[, .(event_time, entry_time, dead, hidden)],
              coef = make_coef_tbl(params$beta, fits),
              cindex = cindex,
              quantiles = quantiles,
              calibration = calibration)
  return(res)
}

#' @import foreach
#' @importFrom doRNG %dorng%
#' @export
run_sim <- function(n_sims = 1, x, params,
                    method = c("coxph", "coxnet"),
                    parallel = FALSE,
                    alpha = 1, lambda = NULL){
  ptm <- proc.time()
  method <- match.arg(method)
  if (file.exists("simout")) {
    msg <- paste0("Running simulation with 'method' = ", method, 
                  " , 'p' = ", ncol(x), ", and ",
                  "'lambda' = ", paste(lambda, collapse=", "))
    cat(msg, file= "simout", append = TRUE, sep = "\n")
  }
  
  # Main loop
  combine_parallel <- function(x, ...) {
    inner_fun <- function(i) {
      c(x[[i]], lapply(list(...), function(y) y[[i]]))
    }
    lapply(seq_along(x), inner_fun)
  }
  out <- foreach(i = 1:n_sims, 
                 .export = c("sim_survdata", "surv_rng", "predict_location",
                             "entry_time_rng", "predict_xb",
                             "make_train_test", "fit_models", 
                             "eval_model_scenarios", "make_cindex_tbl", "make_quantile_tbl", 
                             "make_surv_tbl", "make_coef_tbl", "fit_model", "adjust_Surv", "predict_surv", 
                             "compute_concordance", "calibrate_simsurv", "calibrate", 
                             "calibrate.list", "calibrate.survfit", "run_sim1"),
                 .packages = c("data.table", "dplyr"), # Used in functions without namespaces (pipes for dplyr)
                 .combine = "combine_parallel",
                 .errorhandling = "pass",
                 .init = rep(list(list()), 5)
  ) %dorng% {
    if (i %% 5 == 0) print(i)
    data <- sim_survdata(X = x, params = params)
    run_sim1(data = data, method = method, parallel = parallel,
             alpha = alpha, lambda = lambda)
  }
  run_time <- proc.time() - ptm

  # Combine results
  names(out) <- c("surv_data", "coef", "cindex", "quantiles", "calibration")
  res <- list(
    surv_data = rbindlist(out$surv_data, idcol = "sim"),
    coef = rbindlist(out$coef, idcol = "sim"),
    quantiles = rbindlist(out$quantiles, idcol = "sim"),
    cindex = rbindlist(out$cindex, idcol = "sim"),
    calibration = rbindlist(out$calibration, idcol = "sim"),
    run_time = run_time,
    n_patients = nrow(x),
    p = ncol(x),
    n_sims = n_sims
  )
  return(res)
}

#' @export
summarize_sim <- function(object, save = FALSE, ...) {
  res <- list()
  
  # Helper function
  summarize_dist <- function(data, var, by = NULL) {
    data[, .(mean = mean(get(var)),
             min = min(get(var)),
             q01 = quantile(get(var), .01),
             q10 = quantile(get(var), .1),
             q25 = quantile(get(var), .25),
             median = median(get(var)),
             q75 = quantile(get(var), .75),
             q90 = quantile(get(var), .9),
             q99 = quantile(get(var), .99),
             max = max(get(var))),
         by = by][, ]
  }

  # Survival data
  ## Proportion censored
  object$surv_data[, rc := 1 - dead]
  prop_rc_dist <- object$surv_data[, .(prop = mean(rc)), by = "sim"]
  res$prop_rc <- summarize_dist(prop_rc_dist, "prop")
  
  # Proportion truncated
  prop_lt_dist <- object$surv_data[, .(prop = mean(hidden)), by = "sim"]
  res$prop_lt <- summarize_dist(prop_lt_dist, "prop")
  
  # Distribution of event times
  res$event_time <- summarize_dist(object$surv_data, "event_time")
  
  # Coefficients
  res$coef <- object$coef[, lapply(.SD, mean, na.rm = TRUE),
                          .SDcols = c("Truth", "RC-adjusted", 
                                      "LTRC-adjusted"),
                          by = Variable]
  
  # Quantiles
  q_cols <- colnames(object$quantiles)[grep("q0", colnames(object$quantiles))]
  res$quantiles <- object$quantiles[, lapply(.SD, mean),
                                    .SDcols = q_cols,
                                    by = c("left_truncation", "test_sample")]
  
  # C-index
  res$cindex <- summarize_dist(object$cindex, "C-index",
                               by = c("Left truncation adjustment",
                                      "Test sample"))
  
  # Calibration
  calibration_data <- object$calibration[, 
                                         .(pred = mean(pred), obs = mean(obs)),
                                         by = c("time", "interval", "left_truncation", "obs_method")]
  calibration_data <- split(calibration_data, by = "obs_method")
  calibration_plots <- vector(mode = "list", length = length(calibration_data))
  names(calibration_plots) <- unlist(lapply(calibration_data, 
                                            function (x) unique(x$obs_method)))
  for (i in 1:length(calibration_plots)) {
    setattr(calibration_data[[i]], "class", c("calibrate", "data.frame"))
    calibration_plots[[i]] <- autoplot(calibration_data[[i]], 
                                       col = "left_truncation") +
      scale_colour_discrete("Left truncation adjustment") +
      theme(legend.position = "bottom") 
  }
  res$calibration_plots <- calibration_plots
  
  # Run time
  res$run_time <- data.frame(t(c(object$run_time)),
                             n_patients = object$n_patients,
                             n_sims = object$n_sims)
  
  # Return
  if (save) save_sim(res, ...)
  return(res)
}

#' @export
save_sim <- function(object, model_name){
  fig_path <- paste0("figs/sim_", model_name, "_")
  tbl_path <- paste0("tables/sim_", model_name, "_")
  output_path <- paste0("output/sim_", model_name, "_")
  
  # Tables
  ## CSV
  write.csv(object$prop_lt, paste0(tbl_path, "prop_lt.csv"))
  write.csv(object$prop_rc, paste0(tbl_path, "prop_rc.csv"))
  write.csv(object$event_time, paste0(tbl_path, "event_time.csv"))
  write.csv(object$coef, paste0(tbl_path, "coef.csv"))
  write.csv(object$quantiles, paste0(tbl_path, "quantiles.csv"))
  write.csv(object$cindex, paste0(tbl_path, "cindex.csv"))
  write.csv(object$run_time, paste0(tbl_path, "run_time.csv"),
            row.names = FALSE)
  
  ## LaTeX
  print(xtable(
    object$cindex[, c("Left truncation adjustment", "Test sample", 
                      "mean"), with = FALSE]),
    include.rownames = FALSE, include.colnames = FALSE,
    only.contents = TRUE, sanitize.text.function = identity,
    file = paste0(tbl_path, "cindex.txt"
    ))
  
  # Figures
  ggsave(paste0(fig_path, "calibration_complete.pdf"), 
         object$calibration_plots$complete, 
         height = 5, width = 7)
  ggsave(paste0(fig_path, "calibration_observed_lt.pdf"), 
         object$calibration_plots$observed_lt, 
         height = 5, width = 7)
  ggsave(paste0(fig_path, "calibration_observed.pdf"), 
         object$calibration_plots$observed, 
         height = 5, width = 7)
}