#' initial estimation of the model parameters for a specified emission distribution 
#'
#' Provides the initial estimates of the model parameters of a specified emission 
#' distribution characterized by the \code{mstep} function, for an initial clustering 
#' obtained by \code{\link{initial_cluster}}
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param clus an initial clustering obtained by \code{initial_cluster}
#' @param mstep the mstep function of the EM algorithm with an style simillar to that of \code{\link{mixmvnorm_mstep}}
#' @param verbose logical. if TRUE the outputs will be printed 
#' @param ... additional parameters of the \code{mstep} function
#' 
#' @return a list containing the following items:
#' \itemize{
#' \item \code{emission}{ list the estimated parameterers of the emission distribution}
#' \item \code{leng}{ list of waiting times in each state for each sequence}
#' \item \code{clusters}{ the exact clusters of each observation (available if \code{ltr}=FALSE)}
#' \item \code{nmix}{ the number of mixture components (a vector of positive (non-zero) integers of length \code{nstate})}
#' \item \code{ltr}{ logical. if TRUE a left to right hidden hybrid Markovian/semi-Markovianmodel is assumed}
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
#' sojourn <- list(shape = c(0,3,0), scale = c(0,10,0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(10,8,8,18), seed = 1234, remission = rmixmvnorm)
#' clus = initial_cluster(train,nstate=3,nmix=c(2,2,2),ltr=FALSE,
#' final.absorb=FALSE,verbose=TRUE)
#' par = initial_estimate(clus,verbose=TRUE)
#'
#' @export
#'
initial_estimate<-function(clus,mstep=mixmvnorm_mstep,verbose=FALSE,...){
		mix.clus = clus$mix.clus
		state.clus = clus$state.clus
		ltr = clus$ltr
		final.absorb = clus$final.absorb
		nmix = clus$nmix
		clust.X = clus$clust.X
		num.units = length(clust.X)
		p = ncol(clust.X[[1]][[1]])
		nstate = length(nmix)
		if(verbose) cat("Intitial estimation .... \n")
		artx = matrix(1,10,p)
		artwt1 = matrix(1,10,nstate)
		artwt2 = list()
		for(j in 1:nstate) artwt2[[j]] = matrix(1,10,nmix[j])
		artem = tryCatch({mstep(artx,artwt1,artwt2,...)},
			error=function(e){stop("mstep function is not suitable!")})
		leng = list()
		if(any(nmix>1) & !("mix.p" %in% names(artem))) stop("mstep function is not suitable for mixture emissions!")
		Tx = list()
		for(j in 1:nstate){
			for(mp in 1:num.units){
				if(!all(is.na(clust.X[[mp]][[j]]))){
					Tx[[j]]= clust.X[[mp]][[j]]
					break
				}
			}
			if(num.units>1) for(m in setdiff(1:num.units,mp)) Tx[[j]]= rbind(Tx[[j]],clust.X[[m]][[j]])
			Tx[[j]] = as.matrix(Tx[[j]][apply(Tx[[j]],1,function(xx) !any(is.na(xx))),])
			if(verbose) cat("State ",j," estimation \n")
			leng[[j]]= 0
			for(m in 1:num.units){
				if(j < nstate | !final.absorb)
					if(ltr)
						if(!all(is.na(clust.X[[m]][[j]]))){	
							leng[[j]][m]=nrow(clust.X[[m]][[j]])
						}else{
							leng[[j]][m]=0
						}
				if(verbose) .progress(x=m,max=num.units)
			}# for m
			leng[[j]][is.na(leng[[j]])]=0
			mixind = which(names(artem)=="mix.p")
 			if(nmix[j]>1){
				for(k in 1:nmix[j]){
					Dmat = as.matrix(Tx[[j]][mix.clus[[j]] == k,])
					if(ncol(Dmat) == 1 & p > 1) Dmat = t(Dmat) 
					wt1 = as.matrix(rep(1/nrow(Dmat),nrow(Dmat)))
					wt2 = list(as.matrix(rep(1,nrow(Dmat))))
					em = mstep(Dmat,wt1,wt2,...)
					for(l in setdiff(1:length(artem),mixind)){
						artem[[l]][[j]][[k]] = em[[l]][[1]]
						if(is.null(dim(artem[[l]][[j]][[k]]))){
							if(length(artem[[l]][[j]][[k]])==p) names(artem[[l]][[j]][[k]])<-colnames(clust.X[[1]][[1]])
						}else{
							if(ncol(artem[[l]][[j]][[k]])==p) colnames(artem[[l]][[j]][[k]])<-colnames(clust.X[[1]][[1]])
						}
					}#for l
					if(verbose) cat("Mixture component ",k," estimation\n")
					artem[[mixind]][[j]][k] = sum(mix.clus[[j]] == k)/length(mix.clus[[j]])
				}# for k 
			}else{
				Dmat = Tx[[j]]
				wt1 = as.matrix(rep(1/nrow(Dmat),nrow(Dmat)))
				wt2 = list(as.matrix(rep(1,nrow(Dmat))))
				em = mstep(Dmat,wt1,wt2,...)
				for(l in setdiff(1:length(artem),mixind)){
					artem[[l]][[j]] = em[[l]][[1]]
					if(is.null(dim(artem[[l]][[j]]))){
						if(length(artem[[l]][[j]])==p) names(artem[[l]][[j]])<-colnames(clust.X[[1]][[1]])
					}else{
						if(ncol(artem[[l]][[j]])==p) colnames(artem[[l]][[j]])<-colnames(clust.X[[1]][[1]])
					}
				}
			}#if else		
		}# for j
		ret = list(emission = artem, leng=leng, state.clus=state.clus, ltr=ltr, nmix=nmix)
	ret
}
