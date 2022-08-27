#' the M step function of the EM algorithm
#'
#' The M step function of the EM algorithm for the mixture 
#' of Gaussian linear (Markov-switching) regressions as the emission distribution using the 
#' responses and covariates matrices and the estimated weight vectors
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param x the observation matrix including responses and covariates
#' @param wt1 the state probabilities matrix (number of observations 
#' times number of states)
#' @param wt2 the mixture components probabilities list (of length 
#' nstate) of matrices (number of observations times number of 
#' mixture components)
#' @param resp.ind a vector of the column numbers of \code{x} which contain response variables. 
#' The default is 1, which means that the first column of \code{x} is the univariate 
#' response variable 
#'
#' @return list of emission (mixture of Gaussian linear regression models) parameters:
#' (\code{intercept}, \code{coefficients}, \code{csigma} (conditional covariance) and \code{mix.p})
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
#' clus = initial_cluster(train = train, nstate = 3, nmix = c(1, 2, 1),
#' ltr = FALSE, final.absorb = FALSE, verbose = TRUE, regress = TRUE)
#' initmodel = initialize_model(clus = clus ,mstep = mixlm_mstep,
#' dens.emission = dmixlm, sojourn = NULL, semi = rep(FALSE, 3),
#' M = max(train$N),verbose = TRUE)
#' fit1 = hhsmmfit(x = train, model = initmodel, mstep = mixlm_mstep,
#' M = max(train$N))
#' plot(train$x[, 1] ~ train$x[, 2], col = train$s, pch = 16, 
#' xlab = "x", ylab = "y")
#' abline(fit1$model$parms.emission$intercept[[1]],
#' fit1$model$parms.emission$coefficient[[1]], col = 1)
#' abline(fit1$model$parms.emission$intercept[[2]][[1]],
#' fit1$model$parms.emission$coefficient[[2]][[1]], col = 2)
#' abline(fit1$model$parms.emission$intercept[[2]][[2]],
#' fit1$model$parms.emission$coefficient[[2]][[2]], col = 2)
#' abline(fit1$model$parms.emission$intercept[[3]],
#' fit1$model$parms.emission$coefficient[[3]], col = 3)
#' 
#' @references
#' Kim, C. J., Piger, J. and Startz, R. (2008). Estimation of Markov 
#' regime-switching regression models with endogenous switching. 
#' Journal of Econometrics, 143(2), 263-273.
#' 
#' @export
#'
mixlm_mstep <- function(x, wt1, wt2, resp.ind = 1) 
{
  	emission <- list(mix.p = list(), intercept = list(), 
		coefficients = list(), csigma = list())
	nstate <- ncol(wt1)
	nmix <- c()
	y <- as.matrix(x[, resp.ind])
	x <- as.matrix(x[, - resp.ind])
	dx <- ncol(x)
	dy <- ncol(y)
	x <- x[1:nrow(y), ]
	x <- cbind(1, x)
  	for (j in 1:nstate) {
		nmix[j] <- dim(wt2[[j]])[2]
		if (nmix[j] > 1) {
			emission$intercept[[j]] <- list()
			emission$coefficients[[j]] <- list()
			emission$csigma[[j]] <- list()
			emission$mix.p[[j]] <- rep(0, nmix[j])
			for (i in 1:nmix[j]) {
    				tmp <- cov.mix.wt(y, wt1[, j],wt2[[j]][, i])
				upcov <- tmp$cov
				wt <- wt1[, j] * wt2[[j]][, i]
				W <- diag(wt)
				beta <- ginv(t(x) %*% W %*% x) %*% t(x) %*% W %*% y 
				res <- y - x %*% beta
				wt <- wt / sum(wt)
				res <- sqrt(wt) * res
    				emission$intercept[[j]][[i]] <- as.matrix(beta[1,])
    				emission$coefficients[[j]][[i]] <- as.matrix(beta[-1,])
    				rescov <- t(res) %*% res /(1 - sum(wt^2))
				if(det(as.matrix(rescov)) > det(as.matrix(upcov))) rescov <- upcov
				emission$csigma[[j]][[i]] <- rescov
				emission$mix.p[[j]][i] <- mean(wt2[[j]][, i])
			}
			emission$mix.p[[j]] <- emission$mix.p[[j]] / sum(emission$mix.p[[j]])
		} else {
    			tmp <- cov.mix.wt(y, wt1[, j], wt2[[j]][, 1])
			upcov <- tmp$cov
			wt <- wt1[, j] 
			W <- diag(wt)
			beta <- ginv(t(x) %*% W %*% x) %*% t(x) %*% W %*% y 
			wt <- wt / sum(wt)
			res <- y - x %*% beta
			res <- sqrt(wt) * res
    			emission$intercept[[j]] <- as.matrix(beta[1,])
    			emission$coefficients[[j]] <- as.matrix(beta[-1,])
   			rescov <- t(res) %*% res /(1 - sum(wt^2))
			if(det(as.matrix(rescov)) > det(as.matrix(upcov))) rescov <- upcov
			emission$csigma[[j]] <- rescov
			emission$mix.p[[j]] <- 1
		}#if else 
	}# for j
  	emission
}
