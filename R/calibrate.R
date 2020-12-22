#' Calibrate survival model
#' 
#' Calibrate survival models by comparing predicted survival
#' probabilities with pseudo-observed probabilities.
#' 
#' @param object An object of the appropriate class.
#' @param ... Additional arguments to pass to other methods.
#' @return A [data.frame] with class `calibrate` comparing predicted survival
#' probabilities to pseudo-observed probabilitie at specified time points.
#' @export
calibrate <- function (object, ...) {
  UseMethod("calibrate", object)
}

#' @name calibrate
#' @param y A [survival::Surv] object used to validate the model.
#' @param times The time points at which predicted survival probabilities
#' are compared to pseudo-observed probabilities. 
#' @param n_groups Number of groups to break predicted survival probabilities
#' into. 
#' @export
calibrate.survfit <- function(object, times, y, n_groups = 10){
  if (!inherits(y, "Surv")){
    stop("'y' must be an object of class 'Surv'.")
  }
  survfit_summary <- summary(object, times = times) 
  ymat <- as.matrix(y)
  u <- survfit_summary[["time"]]
  survhat <- survfit_summary[["surv"]]
  if (ncol(survhat) != nrow(y)){
    stop(paste0("The number of rows in 'y' must equal the number of columns ",
                "in 'x$surv'; that is, the number of observations in 'x' must ",
                "equal the number of observations in 'y'."))
  }
  surv <- data.table(u = rep(u, ncol(survhat)),
                     pred = c(survhat))
  surv <- cbind(surv,
                ymat[rep(1:nrow(ymat), each = length(u)), ])
  
  
  # (1) Cutpoints of predicted S(u|x)
  ntile <- function(x, n){
    cut(x, 
        breaks = quantile(x, 
                          probs = seq(0, 1, length = n + 1), 
                          na.rm = TRUE,
                          type = 2),
        include.lowest = TRUE,
        labels = FALSE)
  }
  if (n_groups > 1){
    #NOTE: ntile() and dplyr::ntile() should be equivalent but ntile() seems to cause more errors
    # when n + 1 quantiles cannot be recreated.
    surv[, interval := dplyr::ntile(pred, n = n_groups), by = "u"] 
  } else{
    surv[, interval := 1]
  }
  
  # (2) Average predicted S(u|x) in each interval
  surv_mean <- surv[, .(pred = mean(pred),
                        n_patients = .N),
                    by = c("u", "interval")]
  
  # (3) Compute KM for patients in each interval
  if(ncol(y) == 3){
    kmfit <- survival::survfit(survival::Surv(start, stop, status) ~ survival::strata(interval), 
                               data = surv)
  } else{
    kmfit <- survival::survfit(survival::Surv(time, status) ~ survival::strata(interval), 
                               data = surv)
  }
  kmfit_summary <- summary(kmfit, times = u, extend = TRUE)
  if (n_groups > 1){
    strata <- as.integer(gsub("survival::strata(interval)=interval=", "",
                              kmfit_summary$strata, fixed = TRUE))
  } else{
    strata <- 1
  }
  kmfit_df <-  data.frame(
    interval = strata,
    time = kmfit_summary$time,
    obs = kmfit_summary$surv
  )
  setnames(surv_mean, "u", "time")
  surv_mean <- merge(surv_mean, kmfit_df, 
                     by = c("time", "interval"), 
                     all.x = TRUE)
  
  # Return
  res <- data.frame(surv_mean)
  class(res) <- c("calibrate", class(res))
  return(res)
}

#' @name calibrate
#' @param group_df A data frame of grouping variables. Must have the same
#' number of rows as `object`,
#' @export
calibrate.list <- function(object, times, y, n_groups = 10, group_df){
  n_objects <- length(object)
  res <- vector(mode = "list", n_objects)
  for (i in 1:n_objects){
    res[[i]] <- calibrate(object[[i]], times = times, y = y, n_groups = n_groups)
    res[[i]] <- cbind(res[[i]], group_df[i, ])
  }
  res <- rbindlist(res)
  res <- data.frame(res)
  class(res) <- c("calibrate", class(res))
  return(res)
}

#' Plot survival calibration curves
#' 
#' A [ggplot2::autoplot()] method for creating survival calibration curves.
#' 
#' @param object An object of class [calibrate].
#' @return A [ggplot] object.
#' @export
autoplot.calibrate <- function(object, colour = NULL){
  object$f_time <- factor(object$time,
                          levels = object$time,
                          labels = paste0("Time = ", object$time))
  p <- ggplot(object)
  if (is.null(colour)){
    aes <- aes_string(x = "pred", y = "obs", label = "interval")
  } else{
    aes <- aes_string(x = "pred", y = "obs", col = colour, label = "interval")
  }
  p <- p +
    aes +
    geom_point(size = 2) +
    geom_abline(slope = 1) +
    facet_wrap(~f_time) +
    scale_shape_discrete(name = "Model") +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    xlab("Predicted survival probability") + 
    ylab("Observed survival probability") 
  return(p)
}