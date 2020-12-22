#' Tidy coefficient output
#' 
#' Return coefficients from a model in a tidy format.
#' 
#' @param object The fitted model.
#' @param varnames Variable names for coefficients. Must be a character vector
#' with length equal to the number of coefficients in the model.
#' @param ... Additional arguments to pass to other methods. Currently unused.
#' @return A [tibble] with class `tidycoef` containing the model coefficients.
#' @export
tidycoef <- function (object, ...) {
  UseMethod("tidycoef", object)
}

#' @name tidycoef
#' @export
tidycoef.coxph <- function(object, varnames = NULL, ...){
  res <- broom::tidy(object)
  if (is.null(varnames)){
    varnames <- res$term
  }
  res <- res %>% 
    mutate(variable = varnames,
           hr = exp(estimate))
  class(res) <- c("tidycoef", class(res))
  return(res)
}

#' @name tidycoef
#' @export
tidycoef.cv.glmnet <- function(object, varnames, ...){
  est <- predict(object, type = "coef", s = "lambda.min", exact = FALSE)
  res <- tibble(variable = varnames,
                estimate = est[, 1]) %>%
    mutate(hr = exp(estimate))
  class(res) <- c("tidycoef", class(res))
  return(res)
}

#' @name tidycoef
#' @export
tidycoef.list <- function(object, varnames, ...){
  res <- purrr::map_df(object, tidycoef, varnames = varnames,
                       .id = "model") %>%
    group_by(model)
  class(res) <- c("tidycoef", class(res))
  return(res)
}