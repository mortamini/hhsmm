#' the M step function of the EM algorithm
#'
#' The M step function of the EM algorithm for the mixture 
#' of multivariate normals as the emission distribution with 
#' missing values using the observation matrix and the estimated 
#' weight vectors
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param x the observation matrix which can contain missing values (NA or NaN)
#' @param wt1 the state probabilities matrix (number of observations 
#' times number of states)
#' @param wt2 the mixture components probabilities list (of length 
#' nstate) of matrices (number of observations times number of 
#' mixture components)
#' @param par the parameters of the model in the previous step of 
#' the EM algorithm. For initialization of the model when the data 
#' is initially imputed, \code{par} can be NULL
#'
#' @return list of emission (mixture multivariate normal) parameters:
#' (\code{mu}, \code{sigma} and \code{mix.p})
#'
#' @examples
#' data(CMAPSS)
#' n = nrow(CMAPSS$train$x)
#' wt1 = matrix(runif(3*n),nrow=n,ncol=3)
#' wt2 = list()
#' for(j in 1:3) wt2[[j]] = matrix(runif(5*n),nrow=n,ncol=5)
#' emission = miss_mixmvnorm_mstep(CMAPSS$train$x, wt1, wt2, par=NULL)
#'
#' @import CMAPSS
#'
#' @export
#'
miss_mixmvnorm_mstep <- function(x, wt1, wt2, par) {
	if(anyNA(x) | any(is.nan(x))){
  		emission <- list(mix.p=list() ,mu = list(), sigma = list())
		nstate = ncol(wt1)
		nmix = c()
		missed = apply(x,1,function(t) which(is.na(t)|is.nan(t)))
		means = secm = list()
		d = ncol(x)
  		for(j in 1:nstate) {
			nmix[j] = dim(wt2[[j]])[2]
			if(nmix[j]>1){
				emission$mu[[j]]=list()
				emission$sigma[[j]]=list()
				emission$mix.p[[j]]=rep(0,nmix[j])
				means[[j]] = secm[[j]] = list()
				for(i in 1:nmix[j]){	
					means[[j]][[i]]=sapply(1:length(missed), function(ii){
						l = missed[[ii]]
						if(length(l)==0){NA}else{
							if(length(l) == d){
								par$mu[[j]][[i]][l]
							}else{
								par$sigma[[j]][[i]][l,-l]%*%ginv(par$sigma[[j]][[i]][-l,-l])%*%(x[ii,-l]-par$mu[[j]][[i]][-l])+par$mu[[j]][[i]][l]	
							}
						}
					})
					secm[[j]][[i]]=sapply(1:length(missed), function(ii){
						l = missed[[ii]]
						if(length(l)==0){NA}else{
							if(length(l) == d){
								par$sigma[[j]][[i]][l,l] + par$mu[[j]][[i]][l]%*%t(par$mu[[j]][[i]][l])
							}else{
								par$sigma[[j]][[i]][l,l] - par$sigma[[j]][[i]][l,-l]%*%ginv(par$sigma[[j]][[i]][-l,-l])%*%par$sigma[[j]][[i]][-l,l]+
									par$mu[[j]][[i]][l]%*%t(par$mu[[j]][[i]][l])+par$sigma[[j]][[i]][l,-l]%*%ginv(par$sigma[[j]][[i]][-l,-l])%*%(x[ii,-l]-par$mu[[j]][[i]][-l])%*%
									t(x[ii,-l]-par$mu[[j]][[i]][-l])%*%ginv(par$sigma[[j]][[i]][-l,-l])%*%par$sigma[[j]][[i]][-l,l]+
									2*par$sigma[[j]][[i]][l,-l]%*%ginv(par$sigma[[j]][[i]][-l,-l])%*%(x[ii,-l]-par$mu[[j]][[i]][-l])%*%t(par$mu[[j]][[i]][l])
							}
						}
					})
					tmp.model1 = list(parms.emission = par)
					tmp.model2 = tmp.model1
					tmp.model2$parms.emission$mix.p[[j]][-i]=0
					f = dmixmvnorm
					xr = x
					for(ii in 1:nrow(xr)) xr[ii,is.na(xr[ii,])|is.nan(xr[ii,])] = means[[j]][[i]][[ii]]
					w = f(xr,j,tmp.model2)/f(xr,j,tmp.model1)
					w[is.nan(w)] = 1e-300
					wt2[[j]][,i][is.na(wt2[[j]][,i])|is.nan(wt2[[j]][,i])] = w[is.na(wt2[[j]][,i])|is.nan(wt2[[j]][,i])]
    					tmp <- cov.miss.mix.wt(x, means[[j]][[i]], secm[[j]][[i]], wt1[, j],wt2[[j]][, i])
    					emission$mu[[j]][[i]] <- tmp$center
    					emission$sigma[[j]][[i]] <- tmp$cov
					emission$mix.p[[j]][i] <-tmp$pmix
				}
				emission$mix.p[[j]]=emission$mix.p[[j]]/sum(emission$mix.p[[j]])
			}else{
				means[[j]]=sapply(1:length(missed), function(ii){
					l = missed[[ii]]
					if(length(l)==0){NA}else{
						if(length(l) == d){
							par$mu[[j]][l]	
						}else{
							par$sigma[[j]][l,-l]%*%ginv(par$sigma[[j]][-l,-l])%*%(x[ii,-l]-par$mu[[j]][-l])+par$mu[[j]][l]	
						}
					}
				})
				secm[[j]]=sapply(1:length(missed), function(ii){
					l = missed[[ii]]
					if(length(l)==0){NA}else{
						if(length(l) == d){
							par$sigma[[j]][l,l] + par$mu[[j]][l]%*%t(par$mu[[j]][l])
						}else{
							par$sigma[[j]][l,l] - par$sigma[[j]][l,-l]%*%ginv(par$sigma[[j]][-l,-l])%*%par$sigma[[j]][-l,l]+
								par$mu[[j]][l]%*%t(par$mu[[j]][l])+par$sigma[[j]][l,-l]%*%ginv(par$sigma[[j]][-l,-l])%*%(x[ii,-l]-par$mu[[j]][-l])%*%
								t(x[ii,-l]-par$mu[[j]][-l])%*%ginv(par$sigma[[j]][-l,-l])%*%par$sigma[[j]][-l,l]+
								2*par$sigma[[j]][l,-l]%*%ginv(par$sigma[[j]][-l,-l])%*%(x[ii,-l]-par$mu[[j]][-l])%*%t(par$mu[[j]][l])
						}
					}
				})
    				tmp <- cov.miss.mix.wt(x, means[[j]], secm[[j]], wt1[, j],wt2[[j]][, 1])
    				emission$mu[[j]] <- tmp$center
    				emission$sigma[[j]] <- tmp$cov
				emission$mix.p[[j]] <- tmp$pmix
			}#if else 
		}# for j
	}else{
		emission = mixmvnorm_mstep(x,wt1,wt2)
	}
  	emission
}
