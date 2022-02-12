#' Create hhsmm data of lagged time series
#'
#' Creates a data of class \code{hhsmmdata} containing lagged time series
#' which can be used for fitting auto-regressive hidden hybrid Makrov/semi-Markov 
#' model (AR-HHSMM)
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param data a data of class \code{hhsmmdata} containing a multivariate and multi-state time series
#' @param lags a positive integer which is the number of lags to be calculated
#'
#' @return a data of class \code{hhsmmdata} containing lagged time series
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- rep(FALSE, 3)
#' P <- matrix(c(0.5, 0.2, 0.3, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), nrow = J, 
#' byrow = TRUE)
#' par <- list(intercept = list(0.1, -0.1, 0.2),
#' coefficient = list(-0.6, 0.7, -0.5),
#' csigma = list(5.5, 4, 3.5), mix.p = list(1, 1, 1))
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixlm, semi = semi)
#' train <- simulate(model, nsim = c(50, 60, 84, 100), seed = 1234, 
#' autoregress = TRUE)
#' laggedtrain = lagdata(train)
#' 
#' @export
#'
lagdata<-function(data, lags = 1) 
{
	if (!("hhsmmdata" %in% class(data))) stop("data must be of class hhsmmdata")
	newdata <- data
	Nc <- c(0, cumsum(data$N))
	Nc2 <- c(0, cumsum(data$N - lags))
	newx <- matrix(nrow = nrow(data$x) - lags * length(data$N),
		ncol = ncol(data$x) * (1 + lags))
	for (i in 1:length(data$N)) {
		tmpx <- as.matrix(data$x[(Nc[i] + 1):Nc[i + 1],])
		tmpnewx <- cbind(tmpx[-(1:lags), ], tmpx[-c(0:(lags - 1), nrow(tmpx)), ])
		if (lags > 1) {
			for (l in 2:lags) {
				tmpnewx <- cbind(tmpnewx, tmpx[-c(0:(lags - l), (nrow(tmpx):(nrow(tmpx) - l + 1))), ])
			}
		}
		if (!is.null(data$s)) newdata$s[(Nc2[i] + 1):Nc2[i + 1]] <- data$s[(Nc[i] + 1):(Nc[i + 1] - lags)]
		newx[(Nc2[i] + 1):Nc2[i + 1], ] <- tmpnewx
	}
	if (!is.null(data$s)) newdata$s <- newdata$s[1:(sum(data$N) - lags * length(data$N))]
	newdata$x <- newx[1:(sum(data$N) - lags * length(data$N)), ]
	newdata$N <- data$N - lags
	newdata
}