#' Random data generation from the mixture of Gaussian linear (Markov-switching) models for hhsmm model
#'
#' Generates vectors of covariate and response observations
#' from mixture of Gaussian linear (Markov-switching) models in a specified state and using the 
#' parameters of a specified model
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param j a specified state
#' @param model a \code{\link{hhsmmspec}} model
#' @param covar either a function which generates the covariate vector or a list containing the following items:
#' \itemize{
#' \item \code{mean}{ the mean vector of covariates (to be generated from multivariate normal distribution)}
#' \item \code{cov}{ the variance-covariance matrix of covariates (to be generated from multivariate normal distribution)}
#'}
#' @param ... additional arguments of the \code{covar} function 
#' 
#' @return a random matrix of observations from mixture of Gaussian linear (Markov-switching) models,
#' in which the first columns are associated with the responses 
#' and the last columns are associated with the covariates
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- rep(FALSE, 3)
#' P <- matrix(c(0.5, 0.2, 0.3, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), nrow = J, 
#' byrow = TRUE)
#' par <- list(intercept = list(3, list(-10, -1), 14),
#' coefficient = list(-1, list(1, 5), -7),
#' csigma = list(1.2, list(2.3, 3.4), 1.1),
#' mix.p = list(1, c(0.4, 0.6), 1))
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixlm, semi = semi)
#'
#' #use the covar as the list of mean and 
#' #variance of the normal distribution 
#'
#' train1 <- simulate(model, nsim = c(20, 30, 42, 50), seed = 1234, 
#' remission = rmixlm, covar = list(mean = 0, cov = 1))
#' plot(train1$x[,1] ~ train1$x[,2], col = train1$s, pch = 16, 
#' xlab = "x", ylab = "y")
#'
#' #use the covar as the runif function 
#' #to generate one covariate from standard uniform distribution 
#'
#' train2 <- simulate(model, nsim = c(20, 30, 42, 50), seed = 1234, 
#' remission = rmixlm, covar = runif, 1)
#' plot(train2$x[,1] ~ train2$x[,2], col = train2$s, pch = 16, 
#' xlab = "x", ylab = "y")

#' 
#' @references
#' Kim, C. J., Piger, J. and Startz, R. (2008). Estimation of Markov 
#' regime-switching regression models with endogenous switching. 
#' Journal of Econometrics, 143(2), 263-273.
#'
#' @export
#'
rmixlm <- function(j, model, covar, ...){
	if (is.list(covar)) {
		if (length(covar) != 2) stop("covar must be a list including mean and cov!")
		if (is.null(names(covar)))	 names(covar) <- c("mean","cov")
		dx = length(covar$mean)
		if (dx > 1){
			if (is.null(dim(covar$cov))){
				stop("covar$cov must be a var-covar matrix!")
			} else {
				if (ncol(covar$cov) != nrow(covar$cov)) stop("covar$cov must be a square matrix!")
				if (!isSymmetric(covar$cov)) stop("covar$cov must be a symmetric matrix!")
				if (ncol(covar$cov) != dx) stop("The dim of var-covar matrix and the mean vector are not the same!")
			}
	 		x = rmvnorm(1, mean = covar$mean,
    		       sigma = covar$cov)
		} else {
			if (length(covar$cov) > dx) stop("The dim of var-covar matrix and the mean vector are not the same!")
 			x = rnorm(1, covar$mean, sqrt(covar$cov))	
		}
	} else if (is.function(covar)) {
		x = as.vector(covar(...))
	} else stop("covar must be either a list of mean and cov or a function!")
	if (length(model$parms.emission$mix.p[[j]]) > 1) {
		k = length(model$parms.emission$mix.p[[j]])
		u = runif(1)
		pc = cumsum(c(0, model$parms.emission$mix.p[[j]]))
		p = length(model$parms.emission$intercept[[j]][[1]])
		for (i in 1:k) {
			if ((u >= pc[i]) & (u < pc[i + 1])){
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
		 if (p > 1)
			y = rmvnorm(1, mean = model$parms.emission$intercept[[j]] +
				x %*% t(model$parms.emission$coefficient[[j]]),
                    sigma = model$parms.emission$csigma[[j]])
		 else
 			y = rnorm(1, model$parms.emission$intercept[[j]] +
				x %*% t(model$parms.emission$coefficient[[j]]),
                    sqrt(model$parms.emission$csigma[[j]]))	
	}#if else 
	out = cbind(y, x)
	if(is.null(names(x))) colnames(out)[-1] <- paste0("X",1:length(x))
	out
}
