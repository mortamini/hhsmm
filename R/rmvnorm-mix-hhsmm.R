#' Random data generation from the mixture of multivariate normals for hhsmm model
#'
#' Generates a vector of observations from mixture multivariate 
#' normal distribution in a specified state and using the 
#' parameters of a specified model
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param j a specified state
#' @param model a \code{\link{hhsmmspec}} model
#'
#' @return a random vector of observations from mixture of multivariate 
#' normal distributions
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- c(FALSE, TRUE, FALSE)
#' P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), nrow = J, 
#' byrow = TRUE)
#' par <- list(mu = list(list(7, 8), list(10, 9, 11), list(12, 14)),
#' sigma = list(list(3.8, 4.9), list(4.3, 4.2, 5.4), list(4.5, 6.1)),
#' mix.p = list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.5, 0.5)))
#' sojourn <- list(shape = c(0, 3, 0), scale = c(0, 10, 0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' x = rmixmvnorm(1, model)
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @export
#'
rmixmvnorm <- function(j, model)
{
	if (length(model$parms.emission$mix.p[[j]]) > 1) {
		k = length(model$parms.emission$mix.p[[j]])
		u = runif(1)
		pc = cumsum(c(0, model$parms.emission$mix.p[[j]]))
		p = length(model$parms.emission$mu[[j]][[1]])
		for (i in 1:k) {
			if ((u >= pc[i]) & (u < pc[i + 1])) {
				if (p > 1)
 					x = rmvnorm(1, mean = model$parms.emission$mu[[j]][[i]],
                             sigma = model$parms.emission$sigma[[j]][[i]])
				else
 					x = rnorm(1, model$parms.emission$mu[[j]][[i]],
                             sqrt(model$parms.emission$sigma[[j]][[i]]))
			}
		}#for i 
	} else {
		 p = length(model$parms.emission$mu[[j]])
		 if(p > 1)
		 	x = rmvnorm(1, mean = model$parms.emission$mu[[j]],
                     sigma = model$parms.emission$sigma[[j]])
		 else
		 	x = rnorm(1, model$parms.emission$mu[[j]],
                      sqrt(model$parms.emission$sigma[[j]]))
	}#if else 
	x
}
