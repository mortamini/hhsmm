#' hhsmm specification
#'
#' Specify a model of class \code{"hhsmmspec"} using the model parameters 
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param init vector of initial probabilities
#' @param transition the transition matrix
#' @param parms.emission the parameters of the emission distribution
#' @param sojourn the sojourn distribution, which is one of the following cases:
#' \itemize{
#' \item a list containing \code{d}, which is a \code{nobs} 
#' (number of observations) times 
#' \code{nstates} (number of states) matrix of probabilities,
#' and \code{type = "nonparametric"} for non-parametric sojourn distribution
#' \item a list containing the parameters 
#' \code{mu}, \code{shift} and \code{size} of a shifted 
#' negative binomial distribution, for each semi-Markovian state, and 
#' \code{type = "nbinom"} for negative binomial sojourn distribution
#' \item a list containing the parameters
#' \code{shape} and \code{shift} of a shifted logarithmic 
#' distribution, for each semi-Markovian state, and 
#' \code{type = "logarithmic"} for logarithmic sojourn distribution
#' \item a list containing the parameters 
#' \code{lambda} and \code{shift} of the 
#' shifted poisson distribution, for each semi-Markovian state, and 
#' \code{type = "poisson"} for Poisson sojourn distribution
#' \item a list containing the parameters
#' \code{shape} and \code{scale} of the gamma distribution, 
#' for each semi-Markovian state, and 
#' \code{type = "gamma"} for gamma sojourn distribution
#' \item a list containing the parameters 
#' \code{shape} and \code{scale} of the Weibull distribution, 
#' for each semi-Markovian state, and 
#' \code{type = "weibull"} for Weibull sojourn distribution
#' \item a list containing the parameters
#' \code{meanlog} and \code{sdlog} of the log-normal distribution, 
#' for each semi-Markovian state, 
#' and \code{type = "lnorm"} for log-normal sojourn distribution
#' }
#' @param dens.emission the probability density function of the emission
#' @param remission the random sample generation from the emission distribution 
#' @param mstep the M step function for the EM algorithm 
#' @param semi a logical vector of length nstate: the TRUE associated states are considered as semi-markov
#'
#' @return a model of class \code{"hhsmmspec"}
#'
#' @examples
#' init = c(1, 0)
#' transition = matrix(c(0, 1, 1, 0), 2, 2)
#' parms.emission = list(mix.p = list(c(0.5, 0.5), 1),
#' 				mu = list(list(c(1, 2), c(5, 1)), c(2, 7)),
#'               sigma = list(list(diag(2), 2 * diag(2)), 0.5 * diag(2)))
#' sojourn = list(lambda = 1, shift = 5, type = "poisson")
#' dens.emission = dmixmvnorm
#' remission = rmixmvnorm
#' mstep = mixmvnorm_mstep
#' semi = rep(TRUE,2)
#' model = hhsmmspec(init, transition, parms.emission, sojourn, 
#' dens.emission, remission, mstep, semi)
#'
#' @export
#'
hhsmmspec <- function(init, transition, parms.emission, sojourn = NULL,
  dens.emission, remission = NULL, mstep = NULL, semi = NULL) 
{
  if (is.null(dens.emission)) stop("dens.emission not specified")
  if (length(init)!=NROW(transition)) stop('length(init)!=NROW(transition)')
  if (NROW(transition)!=NCOL(transition)) stop('NROW(transition)!=NCOL(transition)')
  if (is.null(semi)) semi <- which(diag(transition) == 0)
  if (!all(!semi)) if(is.null(sojourn$type)) stop("Sojourn distribution type not specified.")
  if (!all(!semi)) if(all(sojourn$type!=c("nonparametric", "gamma", "poisson", "lnorm", "weibull",
		"nbinom","logarithmic"))) stop(paste("Invalid sojourn type specified (",sojourn$type,")"))
  ans = list(J = length(init), init = init, transition = transition, parms.emission = parms.emission,
		sojourn = sojourn, remission = remission, dens.emission = dens.emission,
		mstep = mstep, semi = semi)
  class(ans) <- c('hhsmmspec')
  .check.hhsmmspec(ans)
  ans  
}

