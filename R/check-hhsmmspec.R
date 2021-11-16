.check.hhsmmspec <- function(object) {
  if(is.null(object$dens.emission)) stop("No emission density function provided!")
  if(is.null(object$init)) stop("No initial distribution specified!")
  if(is.null(object$transition)) stop("No initial distribution specified!")
  if(is.null(object$parms.emission)) stop("No emission parameters specified!")
  if(!all(!object$semi)) if(is.null(object$sojourn)) stop("No sojourn distribution specified!")
  if(!is.null(object$sojourn$d)) if(NCOL(object$sojourn$d)!=nrow(object$transition)) stop("Inconsistent sojourn d")
  if(length(object$init)!=NROW(object$transition))    stop('length(init)!=NROW(transition)')
  if(NROW(object$transition)!=NCOL(object$transition)) stop('NROW(transition)!=NCOL(transition)')
  for(j in 1:object$J) if(object$semi[j] & object$transition[j,j]!=0) stop('Semi-markov states must have diagonal zero transition elements')
}
