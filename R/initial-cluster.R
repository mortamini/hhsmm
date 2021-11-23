#' initial clustering of the data set 
#'
#' Provides an initial clustering for a data of class \code{"hhsmmdata"} which 
#' determines the initial states and mixture components (if necessary) 
#' to be used for initial parameter and model estimation
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param train the train data set of class \code{"hhsmmdata"}, 
#' which can also contain missing data (NA or NaN)
#' @param nstate number of states 
#' @param nmix number of mixture components which is of one of the following forms:
#' \itemize{
#' \item a vector of positive (non-zero) integers of length \code{nstate}
#' \item a positive (non-zero) integer 
#' \item the text \code{"auto"}: the number of mixture components will be determined 
#' automatically based on the within cluster sum of squares 
#' }
#' @param ltr logical. if TRUE a left to right hidden hybrid Markov/semi-Markov model is assumed
#' @param final.absorb logical. if TRUE the final state of the sequence is assumed to be the absorbance state
#' @param verbose logical. if TRUE the outputs will be printed 
#' @param equispace logical. if TRUE the left to right clustering will be performed simply with equal time spaces. 
#' This option is suitable for speech recognition applications
#' @param regress logical. if TRUE the linear regression clustering will be performed
#' @param resp.ind the column indices of the response variables for the linear regression clustering approach. The 
#' default is 1, which means that the first column is the univariate response variable
#'
#' @return a list containing the following items:
#' \itemize{
#' \item \code{clust.X}{ a list of clustered observations for each sequence and state}
#' \item \code{mix.clus}{ a list of the clusters for the mixtures for each state}
#' \item \code{state.clus}{ the exact state clusters of each observation (available if \code{ltr}=FALSE)}
#' \item \code{nmix}{ the number of mixture components (a vector of positive (non-zero) integers of length \code{nstate})}
#' \item \code{ltr}{ logical. if TRUE a left to right hidden hybrid Markov/semi-Markov model is assumed}
#' \item \code{final.absorb}{ logical. if TRUE the final state of the sequence is assumed to be the absorbance state}
#' \item \code{miss}{ logical. if TRUE the \code{train$x} matrix contains missing 
#' data (NA or NaN)}
#' }
#'
#' @details In reliability applications, the hhsmm models are often left-to-right
#' and the modeling aims to predict the future states. In such cases, the
#' \code{ltr}=TRUE and \code{final.absorb}=TRUE should be set. 
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
#' test <-  simulate(model, nsim = c(7,3,3,8), seed = 1234, remission = rmixmvnorm)
#' clus = initial_cluster(train,nstate=3,nmix=c(2,2,2),ltr=FALSE,
#' final.absorb=FALSE,verbose=TRUE)
#'
#' @importFrom mice mice complete
#'
#' @export
#'
initial_cluster<-function(train,nstate,nmix,ltr=FALSE,equispace=FALSE,final.absorb=FALSE,verbose=FALSE,
	regress=FALSE,resp.ind=1){
		if(length(nmix)==1 & mode(nmix)=="numeric") nmix = rep(nmix,nstate)
		if(length(nmix)!=nstate & mode(nmix)=="numeric") stop("length of nmix must be 1 or equal the number of states.")
		if(class(train)!="hhsmmdata") stop("class of train data must be hhsmmdata !")
		if(equispace & !ltr) stop("equispace option is only applied for left to right model (ltr=TRUE)!")
		if(!ltr & final.absorb){
			final.absorb = FALSE
			warning("The final.absorb is for left to right case only! Changed to FALSE...")
		}
		Tx = list()
		data = as.matrix(train$x)
		miss = FALSE
		if(anyNA(data) | any(is.nan(data))){
			miss = TRUE
			allmiss = which(apply(data,1,function(t) all(is.na(t)|is.nan(t))))
			notallmiss = which(!apply(data,1,function(t) all(is.na(t)|is.nan(t))))
			for(ii in allmiss){
				neigh = notallmiss[order(abs(ii-notallmiss))[1:2]]
				data[ii,] = (data[neigh[1],]+data[neigh[2],])/2
			}
			if(ncol(data)>1) data = complete(mice(data,printFlag=FALSE))
		}
		num.units= length(train$N)
		for(j in 1:nstate){
			Tx[[j]] = matrix(0,nrow = 1,ncol=ncol(data))
		}# for j
		Ns = c(0,cumsum(train$N))
		xt = list()
		if(verbose) cat("Within sequence clustering ... \n")
		if(ltr){
			clusters = NULL
		} else {
			clusters = list()
		}
		if(equispace){
			for(m in 1:num.units){
				xt[[m]] = list()
				C=as.matrix(data[(Ns[m]+1):Ns[m+1],])
				if(final.absorb) D= as.matrix(C[-nrow(C),]) else D = C
				T = nrow(D)
				K = nstate-final.absorb
				clus = c(rep(1,T-K*trunc(T/K)),rep(1:K,each=trunc(T/K)))
				for(j in 1:(nstate-final.absorb)){
  					if(sum(clus==j)>0){
						xt[[m]][[j]]=matrix(D[clus==j,],
							nrow=sum(clus==j),ncol=ncol(D))
						colnames(xt[[m]][[j]]) <- colnames(train$x)
						Tx[[j]]= rbind(Tx[[j]],xt[[m]][[j]])
						colnames(Tx[[j]]) <- colnames(train$x)
					} else{
						xt[[m]][[j]] = NA
					}#if else
				}# for j
				if(final.absorb){
					Tx[[nstate]]=rbind(Tx[[nstate]],C[nrow(C),])
					colnames(Tx[[nstate]]) <- colnames(train$x)
					xt[[m]][[nstate]]=matrix(C[nrow(C),],
							nrow=1,ncol=ncol(D))	
					colnames(xt[[m]][[nstate]]) <- colnames(train$x)
				}
			}# for m 
		}else{
			if(ltr){
				for(m in 1:num.units){
					xt[[m]] = list()
					if(verbose) .progress(x=m,max=num.units)
					C=as.matrix(data[(Ns[m]+1):Ns[m+1],])
					if(final.absorb) D= as.matrix(C[-nrow(C),]) else D = C
					if(regress){
						clus = ltr_reg_clus(D,nstate-final.absorb,resp.ind=resp.ind)
					}else clus = ltr_clus(D,nstate-final.absorb)
					for(j in 1:(nstate-final.absorb)){
  						if(sum(clus==j)>0){
							xt[[m]][[j]]=matrix(D[clus==j,],
								nrow=sum(clus==j),ncol=ncol(D))
							colnames(xt[[m]][[j]]) <- colnames(train$x)
							Tx[[j]]= rbind(Tx[[j]],xt[[m]][[j]])
							colnames(Tx[[j]]) <- colnames(train$x)
						} else{
							xt[[m]][[j]] = NA
						}#if else
					}# for j
					if(final.absorb){
						Tx[[nstate]]=rbind(Tx[[nstate]],C[nrow(C),])
						colnames(Tx[[nstate]]) <- colnames(train$x)
						xt[[m]][[nstate]]=matrix(C[nrow(C),],
								nrow=1,ncol=ncol(D))	
						colnames(xt[[m]][[nstate]]) <- colnames(train$x)
					}
				}# for m 
			} else{
					if(regress){
						clus = .kregs(data,nstate,nstart=10,resp.ind=resp.ind)$cluster
					}else clus = kmeans(data,nstate,nstart=10)$cluster
					for(m in 1:num.units){
						clusters[[m]] = clus[(Ns[m]+1):Ns[m+1]]
						xt[[m]] = list()
						D = as.matrix(data[(Ns[m]+1):Ns[m+1],])
						for(j in 1:nstate){
							xt[[m]][[j]]=matrix(D[clusters[[m]]==j,],
								nrow=sum(clusters[[m]]==j),ncol=ncol(data))
							Tx[[j]]= rbind(Tx[[j]],xt[[m]][[j]])
							colnames(Tx[[j]]) <- colnames(train$x)
						}# for j
					}# for m 
			}# if else ltr
		}# if else equispace
		for(j in 1:nstate) Tx[[j]] = as.matrix(Tx[[j]][-1,])
		anmix = c()
		mix.clus = list()
		for(j in 1:nstate){
			if(verbose) cat("State ",j,"\n")
			if(verbose) cat("Between sequence clustering ... \n")
			if(length(nmix)==1){ if(nmix=="auto"){
				if(verbose) cat("Automatic determination of the number of mixture components ... \n")
				continue = TRUE
				DW = Inf
				oldW = (nrow(Tx[[j]])-1)*sum(apply(Tx[[j]],2,var))
				cntr = 0
				eps = 1e-2
				anmix[j] = 1
				while(continue & (cntr+1) < (nrow(Tx[[j]])*0.5) ){
					cntr = cntr + 1
					if(regress){
						tmpclus = .kregs(Tx[[j]],cntr+1,nstart=10,resp.ind=resp.ind)
					}else   tmpclus = kmeans(Tx[[j]],cntr+1,nstart=10)
					newW = sum(tmpclus$withinss)
					DW = c(DW,oldW - newW)
					oldW = newW
					if(cntr>2){
						DDW = -diff(DW)/DW[-1]
						DDDW = -diff(DDW)
						if(any(DDDW<=0) | cntr > 10){
							anmix[j] = which.max(DDW[-1]) + 2
							continue = FALSE
						}# if 
					}# if 
				}# while
				if(regress){
					mix.clus[[j]] = .kregs(Tx[[j]],anmix[j],nstart=10,resp.ind=resp.ind)$cluster
				}else mix.clus[[j]] = kmeans(Tx[[j]],anmix[j],nstart=10)$cluster	
			}# if 
			} else {
				if(regress){
					mix.clus[[j]] = .kregs(Tx[[j]],nmix[j],nstart=10,resp.ind=resp.ind)$cluster
				}else mix.clus[[j]] = kmeans(Tx[[j]],nmix[j],nstart=10)$cluster	
				anmix[j] = nmix[j]
			}# if else 
		}# for j
	out = list(clust.X=xt, mix.clus=mix.clus, state.clus = clusters, 
		nmix=anmix, ltr = ltr, final.absorb = final.absorb, miss = miss)
	return(out)
}
