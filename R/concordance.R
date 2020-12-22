#' Compute concordance
#' 
#' Compute the concordance statistic from a [cv.glmnet] object.
#' 
#' @param object An object of class [cv.glmnet].
#' @param newx Matrix of new values for `x` at which predictions are to be made.
#' @param newy A new value of `y`. See the description of `y` in [glmnet] for 
#' Cox models.
#' @param s Values of the penalty parameter `lambda` at which predictions are
#' required.
#' @param ... Further arguments to pass to [survival::concordance.formula()].
#' @return An object of class [survival::concordance].
#' @importFrom survival concordance
#' @method concordance cv.glmnet
#' @export
concordance.cv.glmnet <- function(object, newx, newy, 
                                    s = c("lambda.min", "lambda.1se"), ...){
  if (!is.Surv(newy)) stop("'newy' must be a 'Surv' object.")
  xb <- stats::predict(object, newx = newx, s = s)
  concordance(object = newy ~ xb, ...)
}