#' pdf of the mixture of Gaussian linear (Markov-sweaching) models for hhsmm
#'
#' The probability density function of a mixture Gaussian linear (Markov-sweaching) models
#' for a specified observation vector, a specified state and a specified 
#' model's parameters
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat, \email{aftbayat@@gmail.com}
#'
#' @param x the observation matrix including responses and covariates
#' @param j a specified state between 1 to nstate
#' @param model a hhsmmspec model
#' @param resp.ind a vector of the column numbers of \code{x} which contain response variables. 
#' The default is 1, which means that the first column of \code{x} is the univariate 
#' response variable 
#'
#' @return the probability density function value
#'
#' @examples
#'J <- 3
#'initial <- c(1,0,0)
#'semi <- rep(FALSE,3)
#'P <- matrix(c(0.5, 0.2, 0.3, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), nrow = J, byrow=TRUE)
#'par <- list(intercept = list(3,list(-10,-1),14),
#'coefficient = list(-1,list(1,5),-7),
#'csigma = list(1.2,list(2.3,3.4),1.1),
#'mix.p = list(1,c(0.4,0.6),1))
#'model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#'dens.emis = dmixlm, semi = semi)
#'train <- simulate(model, nsim = c(20,30,42,50), seed = 1234, remission = rmixlm, 
#'covar.mean=0, covar.cov=1)
#'plot(train$x[,1]~train$x[,2],col=train$s,pch=16)
#'clus = initial_cluster(train=train,nstate=3,nmix=c(1,2,1),ltr=FALSE,
#'final.absorb=FALSE,verbose=TRUE,regress=TRUE)
#'initmodel = initialize_model(clus=clus,mstep = mixlm_mstep,
#'dens.emission = dmixlm, sojourn=NULL, semi=rep(FALSE,3),M=max(train$N),
#'verbose=TRUE)
#'fit1 = hhsmmfit(x = train, model = initmodel, mstep = mixlm_mstep,
#'M = max(train$N), maxit = 100, lock.transition = FALSE, 
#'lock.d = FALSE, lock.init=TRUE, graphical = FALSE,verbose = TRUE)
#'abline(fit1$model$parms.emission$intercept[[1]],
#'fit1$model$parms.emission$coefficient[[1]],col=1)
#'abline(fit1$model$parms.emission$intercept[[2]][[1]],
#'fit1$model$parms.emission$coefficient[[2]][[1]],col=2)
#'abline(fit1$model$parms.emission$intercept[[2]][[2]],
#'fit1$model$parms.emission$coefficient[[2]][[2]],col=2)
#'abline(fit1$model$parms.emission$intercept[[3]],
#'fit1$model$parms.emission$coefficient[[3]],col=3)
#'
#' @references
#' Kim, C. J., Piger, J. and Startz, R. (2008). Estimation of Markov 
#' regime-switching regression models with endogenous switching. 
#' Journal of Econometrics, 143(2), 263-273.
#' 
#' @export
#'
dmixlm <- function(x, j, model, resp.ind = 1){
	y = as.matrix(x[,resp.ind])
	x = as.matrix(x[,-resp.ind])
	dx = ncol(x)
	dy = ncol(y)
	if(length(model$parms.emission$mix.p[[j]])>1){
		dens=rep(0,nrow(x))
		k = length(model$parms.emission$mix.p[[j]])
		for(i in 1:k){
			if(!is.na(model$parms.emission$mix.p[[j]][i]))
			if(dy>1){
				cmean = as.matrix(as.vector(model$parms.emission$intercept[[j]][[i]])*rep(1,nrow(x)))
					+x%*%t(model$parms.emission$coefficients[[j]][[i]])
				ccov = model$parms.emission$csigma[[j]][[i]]
 				dens = dens + model$parms.emission$mix.p[[j]][i]*dmvnorm(y,mean = cmean,
                               sigma = ccov)
			} else {
				cmean = as.vector(as.matrix(as.vector(model$parms.emission$intercept[[j]][[i]])*rep(1,nrow(x)))
					+x%*%t(model$parms.emission$coefficients[[j]][[i]]))
				csd = as.vector(sqrt(model$parms.emission$csigma[[j]][[i]]))
				dens = dens + model$parms.emission$mix.p[[j]][i]*dnorm(y,cmean,csd)
			}
		}#for i 
	} else {
		if(dy>1){
			cmean = as.matrix(as.vector(model$parms.emission$intercept[[j]])*rep(1,nrow(x)))
					+x%*%t(model$parms.emission$coefficients[[j]])
			ccov = model$parms.emission$csigma[[j]]
		 	dens = dmvnorm(x,mean = cmean,
                     sigma = ccov)
		} else {
			cmean = as.vector(as.matrix(as.vector(model$parms.emission$intercept[[j]])*rep(1,nrow(x)))
					+x%*%t(model$parms.emission$coefficients[[j]]))
			csd = as.vector(sqrt(model$parms.emission$csigma[[j]]))
		 	dens = dnorm(x,cmean,csd)
		}
	}#if else 
	dens
}
