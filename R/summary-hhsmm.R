#' @export
#'
summary.hhsmm <- function(object,...) {
  cat("\n Starting distribution = \n")
  print(object$model$init,2)	
  cat("\n Transition matrix = \n")
  print(object$model$transition,2)
  if(!all(!object$model$semi)) cat("\nSojourn distribution parameters = \n")
  if(!all(!object$model$semi)) print(object$model$sojourn)
  cat("\n Emission distribution parameters = \n")
  print(object$model$parms.emission)
}