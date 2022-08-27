#' predicting the response values for the regime switching model 
#'
#' This function computes the predictions of the response variable for the Gaussian linear (Markov-switching) regression 
#' model for different states for any observation matrix of the covariates 
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param object a fitted model of class \code{"hhsmm"} estimated by \code{hhsmmfit}
#' @param x the observation matrix of the covariates
#' @param K the degrees of freedom for the B-spline
#' 
#' @return list of predictions of the response variable
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
#' @importFrom splines2 bSpline
#' 
#' @export
addreg_hhsmm_predict <- function(object, x, K){
	p <- dim(object$model$parms.emission$coef[[1]])[2]
	J <- object$model$J
	if(is.null(dim(x))){
		x <- as.matrix(x)
		q <- ncol(x)
		if(q != p | p > 1) stop("dimension mismatch for x, 
			perhaps you have entered just one sample! You need to rbind
			the new sample with the train data set for calculating 
			the correct B-spline basis function values.") 
	}
	basis = lapply(1:p, function(i) 
		bSpline(x[, i], df = K,
		Boundary.knots = c(min(x[, i]) - 0.01,
		max(x[, i]) + 0.01)))
	pred = lapply(1:J, function(j)
		object$model$parms.emission$intercept[[j]][1] +
		rowSums(as.matrix(sapply(1:p, function(i)
		basis[[i]] %*% object$model$parms.emission$coef[[j]][,i,1]))))
	return(pred)
}
