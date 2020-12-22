#' Adjust survival object
#' 
#' Adjust a survival object so that it is in the right format given the assumptions
#' made about left truncation. If there is no left truncation adjustment,
#' then the object `y` will be reduced to 2 columns. If there is left truncation,
#' no adjustment is required.
#' @param y An object of class [Surv][survival::Surv] with 3 columns.
#' @param left_trunc `TRUE` if the model adjusts for left truncation and `FALSE`
#' otherwise.
#' @return A [Surv][survival::Surv] object.
#' @export
adjust_Surv <- function(y, left_trunc){
  if (!left_trunc){
    y <- as.matrix(y)
    y <- survival::Surv(time = y[, "stop"],
                        event = y[, "status"])
  } else{
    return(y)
  }
}