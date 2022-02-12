#' Random data generation from the mixture of Gaussian linear (Markov-switching) autoregressive models for hhsmm model
#'
#' Generates vectors of observations
#' from mixture of Gaussian linear (Markov-switching) autoregressive
#' model in a specified state and using the 
#' parameters of a specified model
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param j a specified state
#' @param model a \code{\link{hhsmmspec}} model
#' @param x the previous x vector as the covariate of the autoregressive model
#'
#' @return a random matrix of observations from mixture of Gaussian linear 
#' (Markov-switching) autoregressive model
#'
#' 
#' @export
#'
rmixar <- function(j, model, x){
	dx = length(x)
	if (!is.null(dim(x))) if (nrow(x) > 1)	stop("x must be a vector not a matrix!")
	if (length(model$parms.emission$mix.p[[j]]) > 1){
		k = length(model$parms.emission$mix.p[[j]])
		u = runif(1)
		pc = cumsum(c(0, model$parms.emission$mix.p[[j]]))
		p = length(model$parms.emission$intercept[[j]][[1]])
		q = length(model$parms.emission$coefficient[[j]][[1]])
		r = length(model$parms.emission$csigma[[j]][[1]])
		if (p != dx | q != p | r != q) stop("The length of x and model parameters are not consistent!")
		for (i in 1:k) {
			if ((u >= pc[i]) & (u < pc[i + 1])) {
				if (p > 1)
 					y = rmvnorm(1, mean = model$parms.emission$intercept[[j]][[i]] +
						x %*% t(model$parms.emission$coefficient[[j]][[i]]),
                             sigma = model$parms.emission$csigma[[j]][[i]])
				else
 					y = rnorm(1, model$parms.emission$intercept[[j]][[i]] +
						x %*% t(model$parms.emission$coefficient[[j]][[i]]),
                             sqrt(model$parms.emission$csigma[[j]][[i]]))
			}
		}#for i 
	} else {
		 p = length(model$parms.emission$intercept[[j]])
		 if(p > 1)
			y = rmvnorm(1, mean = model$parms.emission$intercept[[j]] +
				x %*% t(model$parms.emission$coefficient[[j]]),
                    sigma = model$parms.emission$csigma[[j]])
		 else
 			y = rnorm(1, model$parms.emission$intercept[[j]] +
				x %*% t(model$parms.emission$coefficient[[j]]),
                    sqrt(model$parms.emission$csigma[[j]]))	
	}#if else 
	y
}
