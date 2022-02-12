#' prediction of state sequence for hhsmm
#'
#' Predicts the state sequence of a fitted hidden hybrid Markov/semi-Markov model estimated by
#' \code{\link{hhsmmfit}} for a new (test) data of class \code{"hhsmmdata"} with an optional prediction of the 
#' residual useful lifetime (RUL) for a left to right model
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @seealso \code{\link{predict.hhsmmspec}}
#'
#' @param object a fitted model of class \code{"hhsmm"} estimated by \code{hhsmmfit}
#' @param newdata a new (test) data of class \code{"hhsmmdata"},
#' which also can contain missing values (NA or NaN)
#' @param future number of future states to be predicted 
#' @param method the prediction method with two options:
#' \itemize{
#' \item \code{"viterbi"}{ (default) uses the Viterbi algorithm for prediction}
#' \item \code{"smoothing"}{ uses the smoothing algorithm for prediction}
#' }
#' @param RUL.estimate logical. if TRUE the residual useful lifetime (RUL) of a left to right model, as well as 
#' the prediction interval will also be predicted (default is FALSE)
#' @param confidence the method for obtaining the prediction interval of the RUL, with two cases:
#' \itemize{
#' \item \code{"max"}{ (default) the maximum probability as the point predict and the high probability critical
#' values as the lower and upper bounds}
#' \item \code{"mean"}{ the mean value as the point predict and the normal confidence lower and upper bounds as the 
#' prediction interval}
#' }
#' @param conf.level the confidence level of the prediction interval (default 0.95) 
#' @param ... additional parameters for the dens.emission and mstep functions
#'
#' @return a list containing the following items:
#' \itemize{
#' \item \code{x}{ the observation sequence}
#' \item \code{s}{ the predicted state sequence}
#' \item \code{N}{ the vector of sequence lengths}
#' \item \code{p}{ the state probabilities }
#' \item \code{RUL}{ the point predicts of the RUL}
#' \item \code{RUL.low}{ the lower bounds for the prediction intervals of the RUL}
#' \item \code{RUL.up}{ the upper bounds for the prediction intervals of the RUL}
#' }
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
#' train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, remission = rmixmvnorm)
#' test <-  simulate(model, nsim = c(7, 3, 3, 8), seed = 1234, remission = rmixmvnorm)
#' clus = initial_cluster(train, nstate = 3, nmix = c(2, 2, 2),ltr = FALSE,
#' final.absorb = FALSE, verbose = TRUE)
#' semi <- c(FALSE, TRUE, FALSE)
#' initmodel1 = initialize_model(clus = clus,sojourn = "gamma",
#' M = max(train$N), semi = semi)
#' fit1 = hhsmmfit(x = train, model = initmodel1, M = max(train$N),
#' maxit = 100, lock.transition = FALSE, lock.d = FALSE, lock.init = FALSE,
#' graphical = FALSE)
#' yhat1 <- predict(fit1, test)
#'
#' @references
#' Guedon, Y. (2005). Hidden hybrid Markov/semi-Markov chains. 
#' \emph{Computational statistics and Data analysis}, 49(3), 663-688.
#'
#' OConnell, J., & Hojsgaard, S. (2011). Hidden semi Markov 
#' models for multiple observation sequences: The mhsmm package 
#' for R. \emph{Journal of Statistical Software}, 39(4), 1-22.
#'
#' @export
#'
predict.hhsmm <- function(object, newdata, future = 0, method = "viterbi", 
	RUL.estimate = FALSE, confidence = "max", conf.level = 0.95, ...) 
{
  if (missing(newdata)){ stop("newdata missing!") } else { x <- newdata }
  J <- object$J
  semi <- object$model$semi
  m <- -1e300
  if (mode(x) == "numeric" | mode(x) == "integer") {
    warning('x is a primitive vector.  Assuming single sequence.')
    x0 <- x
    N <- NROW(x)
    NN <- c(0, N)
    if (N < 1) stop("N less than one")
  } else {
    N <- x$N
    NN <- cumsum(c(0, x$N))
    x0 <- as.matrix(x$x)
  }
  statehat <- integer(NROW(x))
  statehat <- NA
  if (method == "viterbi") {
    M <- nrow(object$model$d)
    a <- object$model$transition  
    loga <- as.double(log(object$model$transition))
    loga[loga == -Inf] <- m
    loga[loga == Inf] <- 1e300
    loga[is.na(loga) | is.nan(loga)] <- m
    logstart = as.double(log(object$model$init))
    logstart[logstart == -Inf | is.nan(logstart)] <- m
    logstart[logstart == Inf | is.nan(logstart)] <- 1e300
    logstart[is.na(logstart) | is.nan(logstart)] <- m
    d <- apply(object$model$d, 2, function(x) x / sum(x))
    D <- apply(d, 2, function(x) rev(cumsum(rev(x))))
    d <- log(d)
    d[d == -Inf] <- m
    d[d == Inf] <- 1e300
	d[d == 0] <- m
    d[is.na(d) | is.nan(d)] <- m
    D <- log(D)
    D[D == -Inf] <- m
    D[D == Inf] <- 1e300
    D[is.na(D) | is.nan(D)] <- m
    loglik <- 0
	RUL = RUL.up = RUL.low = c()
    for (i in 1:length(N)) {
    		if (NCOL(x0) == 1){
		   	xi <- as.matrix(x0[(NN[i] + 1):NN[i + 1]])
      	} else {
			xi <- as.matrix(x0[(NN[i] + 1):NN[i + 1], ])
	  	}
		if (anyNA(xi) | any(is.nan(xi))) {
			b <- log(.densComputeMiss(xi, object$model, ...))
		} else {
		   	b <- log(unlist(sapply(1:J,function(state) object$f(xi, state, object$model, ...))))
		}# if else missing 
     	b[b == -Inf] <- m
      	b[b == Inf] <- 1e300
      	b[is.na(b) | is.na(b)] <- m
      	tmp = .C("viterbi", a = loga, pi = logstart,
               p = as.double(b), d = as.double(d),
               D = as.double(D), timelength = as.integer(N[i]),
               J = as.integer(J), M = as.integer(rep(M,J)),
               alpha = double(N[i] * J), statehat = integer(N[i]),
               psi_state0 = integer(N[i] * J), 
			  psi_time0 = integer(N[i] * J),
               semi = as.double(semi),PACKAGE = "hhsmm")
		pmat = matrix(tmp$alpha, ncol = object$J)
		pmat <- exp(pmat) / rowSums(exp(pmat))
      	loglik <- loglik+max(tmp$alpha[N[i] * (1:J)])
		statehat[(NN[i] + 1):NN[i + 1]] <- tmp$statehat + 1
		if (RUL.estimate) {
      		alpha =  exp(tmp$alpha[N[i] * (1:J)] - max(tmp$alpha[N[i] * (1:J)]))
	  		delta_bar = alpha / sum(alpha)
      		statehat_seq = tmp$statehat + 1
	  		if (J %in% statehat_seq) {
				RUL[i] = 0
				RUL.low[i] <- RUL.up[i] <- 0
	  		} else {
				state_hat_next <- statehat_seq[N[i]]
				duration_hat = matrix(nrow = J , ncol = N[i])
				duration_hat[, 1] <- rep(1, J)
				for (t in 2:N[i]) {
					alpha_t  <- exp(tmp$alpha[t - 1 + (0:(J - 1)) * N[i]] - max(tmp$alpha[t - 1 + (0:(J-1)) * N[i]]))
					delta_t <- alpha_t / sum(alpha_t) 
					for (j in 1:J) duration_hat[j, t] = delta_t[j] * duration_hat[j, t - 1]
				}
				duration_hat = duration_hat[, N[i]]
				RUL[i] <- RUL.low[i] <- RUL.up[i] <- 0
				if (confidence == "mean") {
					duration_center = apply(object$model$d, 2, function(t) weighted.mean(1:M, t))
					sd_duration = apply(object$model$d, 2, function(t) sqrt(cov.wt(data.frame(1:M), t)$cov))
					lower_duration = duration_center - qnorm(0.5 + conf.level / 2) * sd_duration
					upper_duration = duration_center + qnorm(0.5 + conf.level / 2) * sd_duration
				} else if (confidence == "max") {
					duration_center = apply(object$model$d, 2, which.max)
					lower_duration = apply(object$model$d, 2, function(t) min(which(cumsum(t) >= (1 - conf.level) / 2)))
					upper_duration = apply(object$model$d, 2, function(t) max(which(cumsum(t) <= 1 - conf.level / 2)))
					lower_duration[!is.finite(lower_duration)] <- 0
				} else stop("confidence type must be mean or max")
				while (state_hat_next != J) {
					d_tilde_max <- sum((duration_center - duration_hat) * delta_bar)
					d_tilde_low <- sum((lower_duration - duration_hat) * delta_bar)
					d_tilde_up <- sum((upper_duration - duration_hat) * delta_bar)
					RUL[i] <- RUL[i] + d_tilde_max
					RUL.low[i] <- RUL.low[i] + d_tilde_low
					RUL.up[i] <- RUL.up[i] + d_tilde_up
      		   		delta_bar_next <- as.vector(t(a) %*% delta_bar)
					state_hat_next <- which.max(delta_bar_next)
					delta_bar <- delta_bar_next
					dhat <- 0
				}# while
		  	}# if else 
		}else{
			RUL[i] <- RUL.low[i] <- RUL.up[i] <- NA
		}# if else RUL.estimate
		if (future > 0) {
      		alpha <-  exp(tmp$alpha[N[i] * (1:J)] - max(tmp$alpha[N[i] * (1:J)]))
	  		delta_bar <- alpha / sum(alpha)
      		statehat_seq <- tmp$statehat + 1
			state_hat_next <- statehat_seq[N[i]]
			pmatf = rep(0, J)
			for (fu in 1:future) {
     		   	delta_bar_next <- as.vector(t(a) %*% delta_bar)
				state_hat_next <- which.max(delta_bar_next)
				delta_bar <- delta_bar_next
				pmatf <- rbind(pmatf, delta_bar)
				statehat <- c(statehat, state_hat_next)
			}# for fu 
			pmatf <- pmatf[-1, ]
			rownames(pmatf) <- NULL
			pmat <- rbind(pmat, pmatf)
		}# if furure > 0
	} # for i 
    ans <- list(x = x, s = statehat, N = N, p = pmat,
		loglik = loglik, RUL = RUL, RUL.up = RUL.up, RUL.low = RUL.low)
  } else if (method == "smoothing") {
	M <- nrow(object$model$d)    
    	m <- object$model
    	m$dens.emission <- object$f
    tmp <- hhsmmfit(x, m, object$mstep, M = M, maxit = 1, 
	lock.transition = TRUE, lock.d = TRUE, 
	lock.init = TRUE, graphical = FALSE,
	verbose = FALSE, ...)
	RUL <- RUL.up <- RUL.low <- c()
	a <- object$model$transition
	pmat <- matrix(tmp$estep_variables$gamma, ncol = object$J)
	for (i in 1:length(N)) {
		gamma <- tmp$estep_variables$gamma[N[i] * (1:J)]
      	statehat_seq <- tmp$yhat[(NN[i] + 1):NN[i + 1]]
		if (RUL.estimate) { 
	  		if (J %in% statehat_seq) {
				RUL[i] <- 0
				RUL.low[i] <- RUL.up[i] <- 0
	  		} else {
				state_hat_next <- statehat_seq[N[i]]
				duration_hat <- matrix(nrow = J, ncol = N[i])
				duration_hat[,1] <- rep(1, J)
				for (t in 2:N[i]) {
					for (j in 1:J) duration_hat[j,t] <- tmp$estep_variables$gamma[t - 1 + (j - 1) * N[i]] * duration_hat[j, t - 1]
				}
				duration_hat <- duration_hat[,N[i]]
				delta_bar <- gamma / sum(gamma)
				RUL[i] <- RUL.low[i] <- RUL.up[i] <- 0
				if (confidence == "mean") {
					duration_center <- apply(object$model$d, 2, function(t) weighted.mean(1:M, t))
					sd_duration <- apply(object$model$d, 2, function(t) sqrt(cov.wt(data.frame(1:M), t)$cov))
					lower_duration <- duration_center - qnorm(0.5 + conf.level / 2) * sd_duration
					upper_duration <- duration_center + qnorm(0.5 + conf.level / 2) * sd_duration
				} else if (confidence == "max") {
					duration_center <- apply(object$model$d, 2, which.max)
					lower_duration <- apply(object$model$d, 2, function(t) min(which(cumsum(t) >= (1 - conf.level) / 2)))
					upper_duration <- apply(object$model$d, 2, function(t) max(which(cumsum(t) <= 1 - conf.level / 2)))
					lower_duration[!is.finite(lower_duration)] <- 0
				} else stop("confidence type must be mean or max")
				while (state_hat_next != J) {
					d_tilde_max <- sum((duration_center - duration_hat) * delta_bar)
					d_tilde_low <- sum((lower_duration - duration_hat) * delta_bar)
					d_tilde_up <- sum((upper_duration - duration_hat) * delta_bar)
					RUL[i] <- RUL[i] + d_tilde_max
					RUL.low[i] <- RUL.low[i] + d_tilde_low
					RUL.up[i] <- RUL.up[i] + d_tilde_up
    		     		delta_bar_next <- as.vector(t(a) %*% delta_bar)
					state_hat_next <- which.max(delta_bar_next)
					delta_bar <- delta_bar_next
					dhat <- 0
				}# while 
		  	}# if elese 
		}else{
			RUL[i] <- RUL.low[i] <- RUL.up[i] <- NA
		}
		if (future > 0) {
			gamma <- tmp$estep_variables$gamma[N[i] * (1:J)]
	  		delta_bar <- gamma / sum(gamma)
      		statehat_seq <- tmp$yhat[(NN[i] + 1):NN[i + 1]]
			state_hat_next <- statehat_seq[N[i]]
			pmatf <- rep(0, J)
			for (fu in 1:future) {
     		   	delta_bar_next <- as.vector(t(a) %*% delta_bar)
				state_hat_next <- which.max(delta_bar_next)
				delta_bar <- delta_bar_next
				pmatf <- rbind(pmatf, delta_bar)
				tmp$yhat <- c(tmp$yhat, state_hat_next)
			}# for fu 
			pmatf <- pmatf[-1, ]
			rownames(pmatf) <- NULL
			pmat <- rbind(pmat, pmatf)
		}# if furure > 0
	}# for i 
    	ans <- list(x = x$x, s = tmp$yhat, N = x$N, p = pmat,
		RUL = RUL, RUL.low = RUL.low, RUL.up = RUL.up)
  }else stop(paste("Unavailable prediction method",method))
  ans
}
