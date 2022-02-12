#' the M step function of the EM algorithm
#'
#' The M step function of the EM algorithm for the mixture 
#' of splines nonparametric density estimator 
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, 
#' Reza Salehian,  \email{reza.salehian@@ut.ac.ir}
#'
#' @param x the observation matrix
#' @param wt the state probabilities matrix (number of observations 
#' times number of states)
#' @param K the degrees of freedom for the B-spline, default is \code{K=5}
#' @param lambda0 the initial value of the smoothing parameter, default is \code{lambda0=0.5}
#' 
#' @return list of emission (nonparametric mixture of splines) parameters:
#' (\code{coef})
#'
#' @examples
#' K <- 5
#' x <- rmvnorm(100, rep(0, 2), matrix(c(4, 2, 2, 3), 2, 2))
#' wt <- matrix(rep(1, 100), 100, 1)
#' emission = nonpar_mstep(x, wt, K = K)
#' coef <- emission$coef[[1]]
#' x_axis <- seq(min(x[, 1]), max(x[, 1]), length.out = 100)
#' y_axis <- seq(min(x[, 2]), max(x[, 2]), length.out = 100)
#' f1 <- function(x, y) { 
#'   data = matrix(c(x, y), ncol = 2)
#'   tmpmodel = list(parms.emission = emission)
#' 	 dnonpar(data, 1, tmpmodel)
#' }
#' z1 <- outer(x_axis, y_axis, f1)
#' f2 <- function(x, y) { 
#'   data = matrix(c(x, y), ncol = 2)
#'   dmvnorm(data, rep(0, 2), matrix(c(4, 2, 2, 3), 2, 2))
#' }
#' z2 <- outer(x_axis, y_axis, f2)
#' par(mfrow = c(1, 2))
#' persp(x_axis, y_axis, z1, theta = -60, phi = 45, col = rainbow(50))
#' persp(x_axis, y_axis, z2, theta = -60, phi = 45, col = rainbow(50))
#' 
#' @references
#' Langrock, R., Kneib, T., Sohn, A., & DeRuiter, S. L. (2015).
#' Nonparametric inference in hidden Markov 
#' models using P-splines. \emph{Biometrics}, 71(2), 520-528.
#'
#' @importFrom psych tr
#' @importFrom cpr btensor
#' @importFrom utils object.size
#' 
#' @export
#'
nonpar_mstep = function(x, wt, K = 5, lambda0 = 0.5)
{
  nstate <- ncol(wt)
  emission <- list(coef = list(), lambda = numeric(nstate))
  lambda <- numeric(nstate)
  d <- ncol(x)
  n <- nrow(x)
  tryCatch(
         {
	a<-matrix(0, nrow = n, ncol = K^d)
	if(object.size(a) > 1.8e+9) 
		warning("The dimension of the data or the degree of the spline is large! 
		This will result in a very slow progress!")
	rm(a)
  },
        error = function(cond) {
	stop("The dimension of the data or the degree of the spline is too large!
		There is no enough memory for fitting! Try another emission distribution.")

  })
  basis = btensor(lapply(1:d, function(i) x[, i]),
                  df = K, bknots = lapply(1:d, 
                     function(i) c(min(x[, i]) - 0.01,
                                 max(x[, i]) + 0.01)))
  for (j in 1:nstate) {
    	lambda[j] <- lambda0
    	mloglike_lambda0 <- function(beta) {
		dbeta <- diff(beta,differences = 2)
      	omega <- exp(beta) / sum(exp(beta))
      	loglike <- t(wt[, j]) %*% log(basis %*% omega) -
          		lambda0 / 2 * sum(dbeta ^ 2)
      	return(- loglike) 
    }	
    start <- runif(K^d)
    suppressWarnings(fit <- nlm(mloglike_lambda0, start, hessian = T))
    H_lambda0 <- - fit$hessian
    difference <- 1; eps <- 1e-6
    cntr <- 1
    beta_hat <- list(rep(1, K))
    while (difference > eps) {
      	mloglike <- function(beta){
			dbeta <- diff(beta,differences = 2)
        		omega <- exp(beta) / sum(exp(beta))
        		inf_index <- which(is.infinite(log(basis %*% omega)))
        		loglike <- t(wt[, j]) %*% log(basis %*% omega) -
          			lambda[j] / 2 * sum(dbeta ^ 2)
        		return(- loglike) 
     	}	
      	start <- runif(K ^ d)
      	suppressWarnings(fit <- nlm(mloglike, start, hessian = T))
      	H <- - fit$hessian
      	beta_hat[[cntr + 1]] <- fit$estimate
      	df_lambda <- tr(ginv(H) %*% H_lambda0)
		dbeta <- diff(beta_hat[[cntr+1]],differences = 2)
      	lambda[j] <- (df_lambda - d) / (sum(dbeta ^ 2))
      	difference <- sum(beta_hat[[cntr+1]] - beta_hat[[cntr]])
      	cntr <- cntr + 1
    	}
    	emission$coef[[j]] <- exp(beta_hat[[cntr]]) / sum(exp(beta_hat[[cntr]]))
    emission$lambda[j] <- lambda[j]
  }# for j
  emission
}