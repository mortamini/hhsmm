#' pdf of the mixture of multivariate normals for hhsmm
#'
#' The probability density function of a mixture multivariate normal 
#' for a specified observation vector, a specified state and a specified 
#' model's parameters
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat, \email{aftbayat@@gmail.com}
#'
#' @param x an observation vector or matrix
#' @param j a specified state between 1 to nstate
#' @param model a hhsmmspec model
#'
#' @return the probability density function value
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- c(FALSE, TRUE, FALSE)
#' P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), 
#' nrow = J, byrow = TRUE)
#' par <- list(mu = list(list(7, 8),list(10, 9, 11), list(12, 14)),
#' sigma = list(list(3.8, 4.9), list(4.3, 4.2, 5.4), list(4.5, 6.1)),
#' mix.p = list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.5, 0.5)))
#' sojourn <- list(shape = c(0, 3, 0), scale = c(0, 10, 0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, 
#' remission = rmixmvnorm)
#' p = dmixmvnorm(train$x, 1, model)
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @export
#'
dmixmvnorm <- function(x, j, model){
	x = as.matrix(x)
	d = ncol(x)
	if (length(model$parms.emission$mix.p[[j]]) > 1) {
		dens = rep(0, nrow(x))
		k = length(model$parms.emission$mix.p[[j]])
		for (i in 1:k) {
			if (!is.na(model$parms.emission$mix.p[[j]][i]))
			if (d > 1) {
 				dens = dens + model$parms.emission$mix.p[[j]][i] * 
					dmvnorm(x, mean = model$parms.emission$mu[[j]][[i]],
                               sigma = model$parms.emission$sigma[[j]][[i]])
			} else {
				dens = dens + model$parms.emission$mix.p[[j]][i] * 
					dnorm(x, model$parms.emission$mu[[j]][[i]],
                               sqrt(model$parms.emission$sigma[[j]][[i]]))
			}
		}#for i 
	} else {
		if (d > 1) {
		 	dens = dmvnorm(x, mean = model$parms.emission$mu[[j]],
                     sigma = model$parms.emission$sigma[[j]])
		} else {
		 	dens = dnorm(x, model$parms.emission$mu[[j]],
                     sqrt(model$parms.emission$sigma[[j]]))
		}
	}#if else 
	dens
}
