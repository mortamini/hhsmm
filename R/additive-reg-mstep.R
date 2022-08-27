#' the M step function of the EM algorithm
#'
#' The M step function of the EM algorithm for the Gaussian linear (Markov-switching) regression 
#' as the emission distribution using the 
#' responses and covariates matrices and the estimated weight vectors
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, 
#' Reza Salehian,  \email{reza.salehian@@ut.ac.ir}
#'
#' @param x the observation matrix
#' @param wt the state probabilities matrix (number of observations 
#' times number of states)
#' @param control the parameters to control the M-step function. 
#' The simillar name is chosen with that of \code{\link{dnorm_additive_reg}}, 
#' to be used in \code{...} argument of the \code{\link{hhsmmfit}} function.
#' Here, it contains the following items:
#' \itemize{
#' \item \code{K}{ the degrees of freedom for the B-spline, default is \code{K=5}}
#' \item \code{lambda0}{ the initial value of the smoothing parameter, default is \code{lambda0=0.01}}
#' \item \code{resp.ind}{ a vector of the column numbers of \code{x} which contain response variables. 
#' The default is 1, which means that the first column of \code{x} is the univariate 
#' response variable}
#'}
#' 
#' @return list of emission (nonparametric mixture of splines) parameters:
# \code{coef}, \code{intercept} and \code{sigma}
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
#' train <- simulate(model, nsim = c(20, 30, 42, 50), seed = 1234, 
#' remission = rmixlm, covar = list(mean = 0, cov = 1))
#' clus = initial_cluster(train = train, nstate = 3, nmix = NULL,
#' ltr = FALSE, final.absorb = FALSE, verbose = TRUE, regress = TRUE)
#' initmodel = initialize_model(clus = clus ,mstep = additive_reg_mstep,
#' dens.emission = dnorm_additive_reg, sojourn = NULL, semi = rep(FALSE, 3),
#' M = max(train$N),verbose = TRUE)
#' fit1 = hhsmmfit(x = train, model = initmodel, mstep = additive_reg_mstep,
#' M = max(train$N))
#' plot(train$x[, 1] ~ train$x[, 2], col = train$s, pch = fit1$yhat, 
#' xlab = "x", ylab = "y")
#' text(0,30, "colors are real states",col="red")
#' text(0,28, "characters are predicted states")
#' pred <- addreg_hhsmm_predict(fit1, train$x[, 2], 5)
#' yhat1 <- pred[[1]]
#' yhat2 <- pred[[2]]
#' yhat3 <- pred[[3]]
#' 
#' lines(yhat1[order(train$x[, 2])]~sort(train$x[, 2]),col = 2)
#' lines(yhat2[order(train$x[, 2])]~sort(train$x[, 2]),col = 1)
#' lines(yhat3[order(train$x[, 2])]~sort(train$x[, 2]),col = 3)
#'
#' @references
#' Langrock, R., Adam, T., Leos-Barajas, V., 
#' Mews, S., Miller, D. L., and Papastamatiou, Y. P. (2018). 
#' Spline-based nonparametric inference in general state-switching models. 
#' Statistica Neerlandica, 72(3), 179-200.
#' 
#' @importFrom magic adiag
#' @importFrom splines2 bSpline
#' 
#' @export
additive_reg_mstep <- function(x, wt, control = list(K = 5, lambda0 = 0.01, resp.ind = 1)) 
{
	defcon <- list(K = 5, lambda0 = 0.01, resp.ind = 1)
  	control <- modifyList(defcon, control)
  	K <- control$K
  	lambda0 <- control$lambda0
	resp.ind <- control$resp.ind
  	nstate = ncol(wt)
  	wt <- wt / rowSums(wt)
  	n = nrow(x)
  	y = as.matrix(x[, resp.ind])
  	x = as.matrix(x[1:nrow(y), - resp.ind])
  	dx = ncol(x)
  	dy = ncol(y)
  	basis = lapply(1:dx, function(i) bSpline(x[,i], df = K,
                         Boundary.knots = c(min(x[, i]) - 0.01,
                                            max(x[, i]) + 0.01)))
  	z <- Reduce(cbind, basis)
  	z <- cbind(1, z)
	emission = list(coef = list(), intercept = list(),
                  sigma =  list())
    D2 = diff(diag(K),differences = 2)
    D2 = t(D2) %*% D2
    lambda = array(lambda0, dim = c(nstate, dx, dy))
  	for (j in 1:nstate) {
		emission$sigma[[j]] = cov(y)
    		S = list()
    		for (i in 1:dy) {
      		S[[i]] = adiag(0)
      		for (k in 1:dx) {
        			S[[i]] = adiag(S[[i]], lambda[j, k, i] * D2)
     	 	}
    		}
    		W = diag(wt[, j])
    		coef = sapply(1:dy, function(i)
                    ginv(t(z) %*% W %*% z + S[[i]]) %*%
                    	t(z) %*% W %*% y[,i])
    		emission$intercept[[j]] = coef[1,]
    		emission$coef[[j]] = array(coef[-1,], dim = c(K, dx, dy))
		residuals = y - sapply(1:dy, function(p)
      		emission$intercept[[j]][p] + rowSums(sapply(1:dx, function(i)
        		basis[[i]] %*% emission$coef[[j]][, i, p])))
		swt = wt[, j]/sum(wt[,j])
    		residuals = sapply(1:n, function(i) sqrt(swt[i]) * residuals[i,])
    		residuals = matrix(residuals, n, dy, byrow=TRUE)
    		emission$sigma[[j]] = t(residuals) %*% (residuals)/(1-sum(swt^2))
    		for(i in 1:dx){
      		F = n * ginv(t(basis[[i]]) %*% basis[[i]])
      		for(p in 1:dy){
        			lambda[j, i, p] = 2 /
          		((t(emission$coef[[j]][, i, p]) %*%
              		D2 %*% emission$coef[[j]][, i, p])/
					emission$sigma[[j]][p, p] +
             			sum(diag(F %*% D2) / n))
      		}
    		}
    		S = list()
    		for(i in 1:dy) {
      		S[[i]] = adiag(0)
      		for (k in 1:dx) {
        			S[[i]] = adiag(S[[i]], lambda[j, k, i] * D2)
      		}
    		}
    		coef = sapply(1:dy, function(i)
                    ginv(t(z) %*% W %*% z + S[[i]]) %*%
                    	t(z) %*% W %*% y[,i])
    		emission$intercept[[j]] = coef[1,]
    		emission$coef[[j]] = array(coef[-1,], dim = c(K, dx, dy))
		residuals = y - sapply(1:dy, function(p)
      		emission$intercept[[j]][p] + rowSums(sapply(1:dx, function(i)
        		basis[[i]] %*% emission$coef[[j]][, i, p])))
    		residuals = sapply(1:n, function(i) sqrt(swt[i]) * residuals[i,])
    		residuals = matrix(residuals, n, dy, byrow=TRUE)
    		emission$sigma[[j]] = t(residuals) %*% (residuals)/(1-sum(swt^2)) 
  	}# for j
  	emission
}# end of function
