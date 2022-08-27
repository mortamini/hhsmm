#' Random data generation from the Gaussian additive (Markov-switching) model for hhsmm model
#'
#' Generates vectors of covariate and response observations
#' from the Gaussian additive (Markov-switching) model, using B-Splines in a specified state and using the 
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
#' }
#' @param ... additional arguments of the \code{covar} function 
#' 
#' @return a random matrix of observations from Gaussian additive (Markov-switching) model,
#' in which the first columns are associated with the responses 
#' and the last columns are associated with the covariates
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- rep(FALSE, 3)
#' P <- matrix(c(0.5, 0.2, 0.3, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), nrow = J, 
#' byrow = TRUE)
#' par <- list(intercept = list(-21, -83, 33),
#' coef = list(array(c(1, 8, 52, 27, 38), dim = c(5, 1, 1)),  
#' array(c(99, 87, 94, 77, 50), dim = c(5, 1, 1)), 
#' array(c(-1, -8, -40, -22, -28), dim = c(5, 1, 1))),
#' sigma = list(0.2, 0.4, 0.1))
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dnorm_additive_reg, semi = semi)
#' train <- simulate(model, nsim = 70, seed = 1234, 
#' remission = raddreg, covar = list(mean = 0, cov = 1))
#' plot(train$x[, 1] ~ train$x[, 2], col = train$s, pch = 16, 
#' xlab = "x", ylab = "y")
#' 
#' @references
#' Langrock, R., Adam, T., Leos-Barajas, V., 
#' Mews, S., Miller, D. L., and Papastamatiou, Y. P. (2018). 
#' Spline-based nonparametric inference in general state-switching models. 
#' Statistica Neerlandica, 72(3), 179-200.
#'
#' @export
#'
raddreg <- function(j, model, covar, ...){
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
			sx = rmvnorm(100, mean = covar$mean,
    		       sigma = covar$cov)
			grid = apply(sx, 2, function(t) seq(min(t), max(t), length = 100))
	 		x = rmvnorm(1, mean = covar$mean,
    		       sigma = covar$cov)
		} else {
			if (length(covar$cov) > dx) stop("The dim of var-covar matrix and the mean vector are not the same!")
 			sx = rnorm(100, covar$mean, sqrt(covar$cov))
			grid = as.matrix(seq(min(sx), max(sx), length = 100)	)
 			x = rnorm(1, covar$mean, sqrt(covar$cov))	
		}
	} else if (is.function(covar)) {
		x = as.vector(covar(...))
		sx = sapply(1:100, function(t) as.vector(covar(...)))
		if (length(x)>1) {
			grid = apply(sx, 1, function(t) seq(min(t), max(t), length = 100))
		} else {
			grid = as.matrix(seq(min(sx), max(sx), length = 100)	)
		}
	} else stop("covar must be either a list of mean and cov or a function!")
	K = dim(model$parms.emission$coef[[1]])[1]
	dy = length(model$parms.emission$intercept[[j]])
	dx = length(x)
  	basis = lapply(1:dx, function(i) bSpline(c(x[i],grid[, i]),
                 df = K)[1,])
	if (dy > 1){
    		cmean = sapply(1:dy, function(p) 
      		model$parms.emission$intercept[[j]][p] + 
        			sum(sapply(1:dx, function(i)
                    basis[[i]] %*%
                       as.matrix(model$parms.emission$coef[[j]][, i, p]))))
    		ccov = model$parms.emission$sigma[[j]]
		y = rmvnorm(1, mean = cmean, sigma = ccov)
  	} else {
    		cmean = model$parms.emission$intercept[[j]] +
      		sum(sapply(1:dx, function(i)
                    basis[[i]] %*%
                       as.matrix(model$parms.emission$coef[[j]][, i,])))
    		csd = as.vector(sqrt(model$parms.emission$sigma[[j]]))
		y = rnorm(1, cmean, csd)	
	}
	out = cbind(y, x)
	if(is.null(names(x))) colnames(out)[-1] <- paste0("X",1:length(x))
	out
}
