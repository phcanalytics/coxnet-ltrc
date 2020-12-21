fit_survregs <- function(f, data){
  dist <- c("exp", "weibullPH")
  pretty_dist <- c("Exponential", "Weibull")
  fits <- vector(mode = "list", length = length(dist)  + 1)
  names(fits) <- c(dist, "survspline")
  for (i in 1:length(dist)){
    fits[[i]] <- try(flexsurv::flexsurvreg(f, data = data, dist = dist[i])) 
    fits[[i]]$pretty_dist <- pretty_dist[i]
  }
  fits$survspline <- flexsurv::flexsurvspline(f, data = data, k = 1)
  fits$survspline$pretty_dist <- "Spline"
  return(fits)
}

plot_cumhaz <- function(data){
  # Mnual line and colors
  model_names <- unique(data$model)
  n_models <-length(model_names)
  parametric_model_names <- model_names[model_names != "Kaplan-Meier"]
  colors <-  c("black", viridis::viridis_pal()(8))
  linetypes <- c(1:6, 2:4)
  names(linetypes) <- names(colors) <- c("Kaplan-Meier", parametric_model_names)
  
  # Plot
  ggplot(data, aes(x = time/365.25, y = est, col = model,
                   linetype = model)) +
    geom_line() + xlab("Time") + ylab("Cumulative hazard") +
    scale_color_manual(name = "", values = colors) +
    scale_linetype_manual(name = "", values = linetypes) 
}

compare_survregs <- function(fits, km){
  n_fits <- length(fits)
  ic <- data.frame(NA, nrow = n_fits, ncol = 2)
  colnames(ic) <- c("Model", "AIC", "BIC")
  cumhaz <- vector(mode = "list", length = n_fits)
  for (i in 1:n_fits){
    # Information criteria
    ic[i, "Model"] <- fits[[i]]$pretty_dist
    ic[i, "AIC"] <- fits[[i]]$AIC
    ic[i, "BIC"] <- stats::BIC(fits[[i]])
    
    # Cumulative hazard
    cumhaz[[i]] <- summary(fits[[i]], type = "cumhaz", 
                           ci = FALSE, tidy = TRUE)
    cumhaz[[i]]$model <- fits[[i]]$pretty_dist
  }
  cumhaz <- do.call("rbind", cumhaz)
  
  # Add Kaplan-Meier
  cumhaz <-  rbind(
    cumhaz,
    data.frame(time = km$time,
               est = km$cumhaz,
               model = "Kaplan-Meier")
  )
  
  # Plot cumulative hazard
  cumhaz_plot <- plot_cumhaz(cumhaz)
  
  # Return
  return(list(ic = ic, cumhaz = cumhaz,
              cumhaz_plot = cumhaz_plot))
}

add_coef <- function(fit, coef){
  pars <- fit$dlist$pars
  n_pars <- length(pars)
  pars_loc <- which(pars == fit$dlist$location)
  pars_coef <- vector(mode = "list", length = n_pars)
  names(pars_coef) <- pars
  for (j in 1:n_pars){
    if (j == pars_loc){
      intercept <- fit$res.t[pars[j], "est"] - sum(coef)
      pars_coef[[j]] <- c("intercept" = intercept, coef)
    } else{
      pars_coef[[j]] <- c("intercept" = fit$res.t[pars[j], "est"])
    }
  }
  return(pars_coef)
}

make_survregs <- function(survregfits, coxfit){
  n_fits <- length(survregfits)
  fits <- vector(mode = "list", length = n_fits)
  for (i in 1:n_fits){
    fits[[i]]$coefs <- add_coef(survregfits[[i]], stats::coef(coxfit))
    fits[[i]]$knots <- survregfits[[i]]$knots  
    fits[[i]]$dlist <- survregfits[[i]]$dlist
    browser()
  }
  names(fits) <- names(survregfits)
  return(fits)
}

#' @export
calibrate_sim <- function(f, data){
  xy <- make_xy(data, f)
  
  make_f <- function(f_lhs, f_rhs){
    return(stats::update.formula(f_rhs, f_lhs))
  }
  
  # Overall survival
  f_os <- survival::Surv(entry_days_dx, os_days_dx, died) ~ .
  os_km <- survival::survfit(make_f(f_os, ~1), data) 
  os_base <- fit_survregs(make_f(f_os, ~1), data)
  os_cox <- survival::coxph(make_f(f_os, f), data)
  os_comparisons <- compare_survregs(os_base, os_km)
  
  # Right censoring
  f_rc <- survival::Surv(entry_days_dx, os_days_dx, 1 - died) ~ .
  rc_km <- survival::survfit(make_f(f_rc, ~1), data) 
  rc_base <- fit_survregs(make_f(f_rc, ~1), data)
  rc_cox <- survival::coxph(make_f(f_rc, f), data) 
  rc_comparisons <- compare_survregs(rc_base, rc_km)
  
  # Return
  res <- list(os_km = os_km,
              os_base = os_base,
              os_cox = os_cox,
              os_comparisons = os_comparisons,
              rc_base = rc_base,
              rc_cox = rc_cox,
              rc_comparisons = rc_comparisons,
              xy = xy,
              formula = f)
  return(res)
}