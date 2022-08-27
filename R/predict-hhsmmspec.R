#' prediction of state sequence for hhsmm
#'
#' Predicts the state sequence of a hidden hybrid Markov/semi-Markov model 
#' for a new (test) data of class \code{"hhsmmdata"} with an optional prediction of the 
#' residual useful lifetime (RUL) for a left to right model
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @seealso \code{\link{predict.hhsmm}}
#'
#' @param object a hidden hybrid Markov/semi-Markov model 
#' @param newdata a new (test) data of class \code{"hhsmmdata"} 
#' @param method the prediction method with two options:
#' \itemize{
#' \item \code{"viterbi"}{ (default) uses the Viterbi algorithm for prediction}
#' \item \code{"smoothing"}{ uses the smoothing algorithm for prediction}
#' }
#' @param M maximum duration in states
#' @param ... additional parameters of the function \code{\link{predict.hhsmm}}
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
#' test <-  simulate(model, nsim = c(5, 3, 3, 8), seed = 1234, remission = rmixmvnorm)
#' clus = initial_cluster(train, nstate = 3, nmix = c(2, 2, 2), ltr = FALSE,
#' final.absorb = FALSE, verbose = TRUE)
#' semi <- c(FALSE, TRUE, FALSE)
#' initmodel1 = initialize_model(clus = clus, sojourn = "gamma", M = max(train$N), semi = semi)
#' yhat1 <- predict(initmodel1, test)
#'
#' @references
#' Guedon, Y. (2005). Hidden hybrid Markov/semi-Markov chains. 
#' \emph{Computational statistics and Data analysis}, 49(3), 663-688.
#'
#' OConnell, J., & Hojsgaard, S. (2011). Hidden semi Markov 
#' models for multiple observation sequences: The mhsmm package 
#' for R. \emph{Journal of Statistical Software}, 39(4), 1-22.
#'
#' @useDynLib hhsmm viterbi
#'
#' @export
#'
predict.hhsmmspec <- function(object, newdata, method = "viterbi", M = NA, ...)
 {
  if (inherits(newdata, "hhsmmdata")) NN = newdata$N
  else NN = length(newdata)
  if(is.na(M)) M = max(NN)
  .check.hhsmmspec(object)
  model <- .build_d(object, M)
  object2 <- list(model = model, J = model$J, f = model$dens.emission)
  class(object2) <- "hhsmm"
  predict.hhsmm(object2, newdata, method = method, ...)    
}
