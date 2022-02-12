#' the score of new observations  
#'
#' computes the score (log-likelihood) of new observations using a trained model 
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param xnew a new single observation, observation matrix 
#' or a list of the class \code{\link{hhsmmdata}} containing $x and $N elements
#' @param fit a fitted model using the \code{\link{hhsmmfit}} function
#' 
#' @return the vector of scores (log-likelihood) of \code{xnew}
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- c(FALSE, TRUE, FALSE)
#' P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), nrow = J, 
#' byrow = TRUE)
#' par <- list(mu = list(list(7, 8), list(10, 9, 11), list(12, 14)),
#' sigma = list(list(3.8, 4.9), list(4.3, 4.2, 5.4), list(4.5, 6.1)),
#' mix.p = list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.5, 0.5)))
#' sojourn <- list(shape = c(0, 3, 0), scale = c(0, 10, 0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, 
#' remission = rmixmvnorm)
#' test <- simulate(model, nsim = c(5, 4, 6, 7), seed = 1234, 
#' remission = rmixmvnorm)
#' clus = initial_cluster(train, nstate = 3, nmix = c(2, 2, 2), ltr = FALSE,
#' final.absorb = FALSE, verbose = TRUE)
#' semi <- c(FALSE, TRUE, FALSE)
#' initmodel1 = initialize_model(clus = clus, sojourn = "gamma",
#' M = max(train$N), semi = semi)
#' fit1 = hhsmmfit(x = train, model = initmodel1, M = max(train$N),
#' maxit = 100, lock.transition = FALSE, lock.d = FALSE, lock.init = FALSE,
#' graphical = FALSE)
#' score(test, fit1)
#' 
#' @export
#'
score <- function(xnew, fit) 
{
  	if (mode(xnew) == "numeric" | mode(xnew) == "integer") {
		if (is.null(dim(xnew))) {
			N = nrow(xnew <- t(as.matrix(xnew)))
		} else {
	    		N = nrow(xnew <- as.matrix(xnew))    
		}
  	} else {
    		N = xnew$N
    		xnew = as.matrix(xnew$x)
  	}
	Nc = c(0, cumsum(N))
	score = c()
	for (i in 1:length(N)) {
		for (j in (Nc[i] + 1):(Nc[i + 1])) {
			xx = matrix(xnew[j, ], 1, ncol(xnew))
			suppressWarnings(score <- c(score, hhsmmfit(xx, fit$model, fit$mstep, maxit = 1, 
				M = fit$M, verbose = FALSE)$loglik))
		}
	}
	score
}