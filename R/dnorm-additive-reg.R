#' pdf of the Gaussian additive (Markov-switching) model for hhsmm
#'
#' The probability density function of a Gaussian additive (Markov-switching) model
#' for a specified observation vector, a specified state and a specified 
#' model's parameters
#' 
#' 
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, 
#' Reza Salehian,  \email{reza.salehian@@ut.ac.ir}
#'
#' @param x the observation matrix including responses and covariates
#' @param j a specified state between 1 to nstate
#' @param model a hhsmmspec model
#' @param control the parameters to control the density function. 
#' The simillar name is chosen with that of \code{\link{additive_reg_mstep}}, 
#' to be used in \code{...} argument of the \code{\link{hhsmmfit}} function.
#' Here, it contains the following items:
#' \itemize{
#' \item \code{K}{ the degrees of freedom for the B-spline, default is \code{K=5}}
#' \item \code{resp.ind}{ a vector of the column numbers of \code{x} which contain response variables. 
#' The default is 1, which means that the first column of \code{x} is the univariate 
#' response variable}
#'}
#'
#' @return the probability density function value
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
#' @export
#' 
dnorm_additive_reg <- function(x, j, model, control = list(K = 5, resp.ind = 1))
{
	defcon <- list(K = 5, resp.ind = 1)
  	control <- modifyList(defcon, control)
  	K <- control$K
	resp.ind <- control$resp.ind
	y = as.matrix(x[, resp.ind])
	x = as.matrix(x[, -resp.ind])
	dx = ncol(x)
	dy = ncol(y)
	n = nrow(x)
	dens = rep(0, n)
  	K = dim(model$parms.emission$coef[[j]])[1]
  	basis = lapply(1:dx, function(i) bSpline(x[, i],
                      df = K, Boundary.knots = c(min(x[, i]) - 0.01,
                                           max(x[, i]) + 0.01)))
	if (dy > 1){
    		cmean = sapply(1:dy, function(p) 
      		model$parms.emission$intercept[[j]][p] + 
        			rowSums(sapply(1:dx, function(i)
                    basis[[i]] %*%
                       as.matrix(model$parms.emission$coef[[j]][, i, p]))))
    		ccov = model$parms.emission$sigma[[j]]
    		dens = dmvnorm(y, mean = cmean, sigma = ccov)
  	} else {
    		cmean = model$parms.emission$intercept[[j]] +
      		rowSums(sapply(1:dx, function(i)
                    basis[[i]] %*%
                       as.matrix(model$parms.emission$coef[[j]][, i,])))
    		csd = as.vector(sqrt(model$parms.emission$sigma[[j]]))
    		dens = dnorm(y, cmean, csd)
	}
  	dens
}
