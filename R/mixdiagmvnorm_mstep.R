#' the M step function of the EM algorithm
#'
#' The M step function of the EM algorithm for the mixture 
#' of multivariate normals with diagonal covariance matrix
#' as the emission distribution using the 
#' observation matrix and the estimated weight vectors
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param x the observation matrix
#' @param wt1 the state probabilities matrix (number of observations 
#' times number of states)
#' @param wt2 the mixture components probabilities list (of length 
#' nstate) of matrices (number of observations times number of 
#' mixture components)
#'
#' @return list of emission (mixture multivariate normal) parameters:
#' (\code{mu}, \code{sigma} and \code{mix.p}), where \code{sigma} is a diagonal matrix
#'
#' @examples
#' data(CMAPSS)
#' n = nrow(CMAPSS$train$x)
#' wt1 <- matrix(runif(3 * n), nrow = n, ncol = 3)
#' wt2 <- list()
#' for(j in 1:3) wt2[[j]] <- matrix(runif(5 * n), nrow = n, ncol = 5)
#' emission <- mixdiagmvnorm_mstep(CMAPSS$train$x, wt1, wt2)
#'
#' 
#' @export
#'
mixdiagmvnorm_mstep <- function(x, wt1, wt2) 
{
  	emission <- list(mix.p = list(), mu = list(), sigma = list())
	nstate <- ncol(wt1)
	nmix <- c()
  	for (j in 1:nstate) {
		nmix[j] <- dim(wt2[[j]])[2]
		if (nmix[j] > 1) {
			emission$mu[[j]] <- list()
			emission$sigma[[j]] <- list()
			emission$mix.p[[j]] <- rep(0, nmix[j])
			for (i in 1:nmix[j]) {	
    				tmp <- cov.mix.wt(x, wt1[, j], wt2[[j]][, i])
    				emission$mu[[j]][[i]] <- tmp$center
    				emission$sigma[[j]][[i]] <- diag(diag(tmp$cov))
				emission$mix.p[[j]][i] <- tmp$pmix
			}
			emission$mix.p[[j]] <- emission$mix.p[[j]] / sum(emission$mix.p[[j]])
		} else {
    			tmp <- cov.mix.wt(x, wt1[, j], wt2[[j]][, 1])
    			emission$mu[[j]] <- tmp$center
    			emission$sigma[[j]] <- diag(diag(tmp$cov))
			emission$mix.p[[j]] <- tmp$pmix
		}#if else 
	}# for j
  	emission
}
