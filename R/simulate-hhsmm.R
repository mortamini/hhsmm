#' Simulation of data from hhsmm model 
#'
#' Simulates a data set of class \code{"hhsmmdata"} using a \code{\link{hhsmmspec}} model
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat, \email{aftbayat@@gmail.com}
#'
#' @param object a \code{\link{hhsmmspec}} model 
#' @param nsim a vector of sequence lengths (might be of length 1)
#' @param seed a random seed to be set
#' @param remission a random emission generation function (default = \code{\link{rmixmvnorm}})
#' @param ... additional parameters of the \code{remission} function
#'
#' @return a list of class \code{"hsmm.data"} containing the following items:
#' \itemize{
#' \item \code{s}{ the vector of states}
#' \item \code{x}{ observation matrix}
#' \item\code{N}{ vector of sequence lengths}
#' }
#'
#' @examples
#' J <- 3
#' initial <- c(1,0,0)
#' semi <- c(FALSE,TRUE,FALSE)
#' P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), nrow = J, byrow=TRUE)
#' par <- list(mu = list(list(7,8),list(10,9,11),list(12,14)),
#' sigma = list(list(3.8,4.9),list(4.3,4.2,5.4),list(4.5,6.1)),
#' mix.p = list(c(0.3,0.7),c(0.2,0.3,0.5),c(0.5,0.5)))
#' sojourn <- list(shape = c(0,3,0), scale = c(0,8,0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(8,5,5,10), seed = 1234, remission = rmixmvnorm)
#'
#' @importFrom stats simulate
#'
#' @export
#'
simulate.hhsmmspec <- function(object, nsim, seed=NULL, remission=rmixmvnorm,...){
	right.truncate=left.truncate=0
	if(!is.null(seed)) set.seed(seed)
	if(is.null(remission) & is.null(object$remission)) stop("remission not specified")
	if(!is.null(remission)) object$remission=remission
	for(j in 1:object$J) if(object$semi[j] & object$transition[j,j]!=0) stop('Semi-markov states must have diagonal zero transition elements')
	if(length( nsim)==1) M =  nsim else M = max( nsim)
	s0 = .simulate_markov(object$init, object$transition,  nsim)
	if (object$sojourn$type == "poisson") {
      	rsojourn <- function(ii){ 
			if(object$semi[ii]) .rpois.hhsmm(1,object$sojourn$lambda[ii],object$sojourn$shift[ii])
			else{ 
				if(object$transition[ii,ii]<1) rgeom(1, 1-object$transition[ii,ii])
				else 1
				}
		}
    	}else if (object$sojourn$type == "gamma") {
    		rsojourn <- function(ii){ 
				if(object$semi[ii])  trunc(rgamma(1,shape=object$sojourn$shape[ii],scale=object$sojourn$scale[ii]))+1
				else{ 
					if(object$transition[ii,ii]<1) rgeom(1, 1-object$transition[ii,ii])
					else 1
				}
		} 
   	} else if (object$sojourn$type == "lnorm") {
    		rsojourn <- function(ii){
			if(object$semi[ii])  trunc(rlnorm(1,object$sojourn$meanlog[ii],object$sojourn$sdlog[ii]))+1
			else{ 
				if(object$transition[ii,ii]<1) rgeom(1, 1-object$transition[ii,ii])
				else 1
			}
		}
	} else if (object$sojourn$type == "logarithmic") {
      	rsojourn <- function(ii){
			if(object$semi[ii]) .rlog(1,object$sojourn$shape[ii])
			else{ 
				if(object$transition[ii,ii]<1) rgeom(1, 1-object$transition[ii,ii])
				else 1
			}
		}        
	} else if (object$sojourn$type == "nbinom") {
      	rsojourn <- function(ii){
				if(object$semi[ii]) .rnbinom.hhsmm(1,object$sojourn$shape[ii])
				else{ 
					if(object$transition[ii,ii]<1) rgeom(1, 1-object$transition[ii,ii])
					else 1
				}
		}   
    } else if (object$sojourn$type == "nonparametric") {
    		rsojourn <- function(ii){
			if(object$semi[ii]) sample(1:nrow(object$sojourn$d), 1, prob = object$sojourn$d[,ii])
			else{ 
				if(object$transition[ii,ii]<1) rgeom(1, 1-object$transition[ii,ii])
				else 1
			}
		} 
    	} else stop("Unknown sojourn distribution")
	u = sapply(s0,rsojourn)
    	s1 = rep(s0, u)[(left.truncate + 1):(sum(u) - right.truncate)]
    	x = as.matrix(sapply(s1,function(i) remission(i,object,...)))
	if(length( nsim)> 1){
  		N =  nsim
  		tmp = cumsum(c(1,N))
  		for(i in 1:(length(tmp)-1)) N[i] = sum(u[tmp[i]:(tmp[i+1]-1)])
  		if(NCOL(x)>1) ret = list(s=s1,x=t(x),N=N)
  		else ret = list(s=s1,x=x,N=N)
	} else {
    		ret = list(s = s1, x = x, N = NROW(x))
	}
    	class(ret) <- "hhsmmdata"
    	ret
}
