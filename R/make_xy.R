#' Standardize regression variables
#' 
#' Standardize variables by subtracting the mean and dividing by the 
#' standard deviation.
#' @param x A matrix. 
#' @export
standardize <- function(x){
  mx <- colMeans(x)
  sx <- sqrt(apply(x, 2, stats::var))
  x <- scale(x, mx, sx)
  x[, sx == 0] <- 0 
  return(x)
}

#' @export
make_xy <- function(data, formula){
  # x
  x <- data %>% 
    stats::model.matrix(formula, data = .)  %>%
    .[, colnames(.) != "(Intercept)"] # Intercept is absorbed by baseline hazard
  colnames(x) <- gsub(" ", "_", colnames(x)) 
  x <- standardize(x) 
  
  # y
  indices <- as.integer(rownames(x))
  surv_data <- data %>% 
    dplyr::slice(indices) %>%
    dplyr::select(entry_days_dx, os_days_dx, died)
  y <- survival::Surv(time = surv_data$entry_days_dx/365.25,
                      time2 = surv_data$os_days_dx/365.25,
                      event = surv_data$died)
  # Return
  return(list(x = x, y = y))
}