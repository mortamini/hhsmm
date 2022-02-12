#' hhsmm model fit
#'
#' Fits a hidden hybrid Markov-semi-Markov model to a data of class \code{"hhsmmdata"} and using an initial
#' model created  by \code{\link{hhsmmspec}} or \code{\link{initialize_model}}
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat, \email{aftbayat@@gmail.com}
#'
#' @param x a data of class \code{"hhsmmdata"}, which can also contain missing values (NA or NaN)
#' @param model an initial model created  by \code{hhsmm.spec} or \code{initialize_model}
#' @param mstep the M step function for the EM algorithm, which also can be given in the model
#' @param M the maximum duration in each state
#' @param maxit the maximum number of iterations for the EM algorithm 
#' @param lock.transition logical. if TRUE the transition matrix will not be updated through the EM algorithm 
#' @param lock.d logical. if TRUE the sojourn probability matrix d will not be updated through the EM algorithm 
#' @param lock.init logical. if TRUE the initial probability vector will not be updated through the EM algorithm 
#' @param graphical logical. if TRUE a plot of the sojourn probabilities will be plotted through the EM algorithm
#' @param verbose logical. if TRUE the outputs will be printed
#' @param ... additional parameters for the dens.emission and mstep functions
#'
#' @return a list of class \code{"hhsmm"} containing the following items:
#' \itemize{
#' \item \code{loglike}{ the log-likelihood of the fitted model}
#' \item \code{AIC}{ the Akaike information criterion of the fitted model}
#' \item \code{BIC}{ the Bayesian information criterion of the fitted model}
#' \item \code{model}{ the fitted model}
#' \item \code{estep_variables}{ the E step (forward-backward) probabilities of the final iteration of the EM algorithm}
#' \item \code{M}{ the maximum duration in each state}
#' \item \code{J}{ the number of states}
#' \item \code{NN}{ the vector of sequence lengths}
#' \item \code{f}{ the emission probability density function}
#' \item \code{mstep}{ the M step function of the EM algorithm}
#' \item \code{yhat}{ the estimated sequence of states}
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
#' train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, 
#' remission = rmixmvnorm)
#' clus = initial_cluster(train, nstate = 3, nmix = c(2 ,2, 2),ltr = FALSE,
#' final.absorb = FALSE, verbose = TRUE)
#' initmodel1 = initialize_model(clus = clus, sojourn = "gamma", 
#' M = max(train$N), semi = semi)
#' fit1 = hhsmmfit(x = train, model = initmodel1, M = max(train$N),
#' maxit = 100, lock.transition = FALSE, lock.d = FALSE, lock.init = FALSE,
#' graphical = FALSE)
#'
#' @references
#' Guedon, Y. (2005). Hidden hybrid Markov/semi-Markov chains. 
#' \emph{Computational statistics and Data analysis}, 49(3), 663-688.
#'
#' OConnell, J., & Hojsgaard, S. (2011). Hidden semi Markov 
#' models for multiple observation sequences: The mhsmm package 
#' for R. \emph{Journal of Statistical Software}, 39(4), 1-22.
#'
#' @useDynLib hhsmm forward_backward
#' 
#' @importFrom stats approx chisq.test cov cov.wt density dnbinom dnorm dpois kmeans nlm pbeta pgamma plnorm pnbinom ppois pweibull qf qnorm quantile rgamma rgeom rlnorm rnbinom rnorm rpois runif sd ts uniroot var weighted.mean
#'
#' @importFrom utils head tail
#' 
#' @import Rcpp
#' 
#' @import Rdpack
#'
#' @export
#'
hhsmmfit <- function(x, model, mstep = NULL, M = NA, maxit = 100, 
	lock.transition = FALSE, lock.d = FALSE, lock.init = FALSE, 
	graphical = FALSE, verbose = TRUE, ...) 
{
  	sojourn.distribution=model$sojourn$type
  	tol=1e-4
  	ksmooth.thresh = 1e-20 #this is a threshold for which d(u) values to use 
	shiftthresh = 1e-20 #threshold for effective "0" when considering d(u)
 	J = nrow(model$transition)
  	model$J = J
	if(is.null(model$semi)) model$semi = rep(TRUE,J) 
   	if(is.null(mstep)) 
    if(is.null(model$mstep)){
	 stop("mstep not specified")
	} else{  mstep=model$mstep }     
  	.check.hhsmmspec(model)
	f=model$dens.emission
  	if(mode(x)=="numeric" | mode(x)=="integer") {
    		warning('x is a primitive vector.  Assuming single sequence.')
    		NN = NROW(x<-as.matrix(x))    
  	} else{
    		NN = x$N
    		x = as.matrix(x$x)
  	}
  	if (anyNA(M)) M = max(NN)
  	if(length(model$init)!=J) stop("length(model$init)!=J")
  	if(NROW(x)!=sum(NN)) stop("NROW(x)!=sum(NN)")
  	model <- .build_d(model,M)
  	new.model = model
  	ll = rep(NA,maxit)
  	rm(model)
  	for(it in 1:maxit) {
  		if (verbose) cat('iteration: ',it,"  ")
    		if(graphical)   plot.hhsmm(list(model=new.model,J=J))
		if(anyNA(x) | any(is.nan(x))){
			p = .densComputeMiss(x,new.model,...)
		}else{
			p = sapply(1:J,function(state) f(x,state,new.model,...))
		}# if else missing 
		p=p/max(p)
		p = matrix(p,ncol=J)
		p[is.na(p) | is.nan(p)] = 1e-300
		p[!is.finite(p)] = 1e+10
		if(any(apply(p,1,max)==0)) stop("Some values have 0 pdf for all states!  Check your model parameters")
		new.model$d[new.model$d==0]=1e-300
		estep_variables=.C("forward_backward", transition = as.double(new.model$transition), 
            init = as.double(new.model$init), p = as.double(p), 
            d = as.double(new.model$d), D = as.double(new.model$D), 
            timelength = as.integer(NN), J = as.integer(J), M = as.integer(rep(M, 
                J)), L1 = double(NROW(x) * J), N = double(NROW(x)), 
            eta = double(M * J), F1 = double(J * NROW(x)), si = double(J * 
                NROW(x)), gamma = double(J * NROW(x)), nsequences = as.integer(length(NN)), 
            totallength = NROW(x), G = double(J * NROW(x)), semi = as.double(new.model$semi),PACKAGE="hhsmm")
		if(!is.null(new.model$parms.emission$mix.p)){
			estep_mix_weights = list()
			for(j in 1:J){
				if(length(new.model$parms.emission$mix.p[[j]])>1){
					k = length(new.model$parms.emission$mix.p[[j]])
					estep_mix_weights[[j]] = matrix(nrow = nrow(x) , ncol = k)
					for(i in 1:k){
						tmp.model = new.model
						tmp.model$parms.emission$mix.p[[j]][-i]=0
						w = f(x,j,tmp.model,...)/f(x,j,new.model,...)
						w[is.nan(w)] = 1e-300
						estep_mix_weights[[j]][,i] = w
					}# for i
				}else{
					estep_mix_weights[[j]] = matrix(1, nrow = nrow(x) , ncol = 1)
				}
			}# for j
		}# if 
		estep_variables$gamma[is.nan(estep_variables$gamma) | is.na(estep_variables$gamma)] = 1e-100
		estep_variables$eta[is.nan(estep_variables$eta) | is.na(estep_variables$eta)] = 1e-100
		estep_variables$N[is.nan(estep_variables$N) | is.na(estep_variables$N)] = 1e-100
    		if(any(estep_variables$gamma<=0)) estep_variables$gamma = sapply(estep_variables$gamma,function(x) max(1e-100,x))
    		if(any(estep_variables$eta<=0)) estep_variables$eta = sapply(estep_variables$eta,function(x) max(1e-100,x))    
    		if(any(unlist(estep_variables$N)<=0))  estep_variables$N = sapply(estep_variables$N,function(x) max(1e-100,x))
    		old.model = new.model
    		state_wt <- matrix(estep_variables$gamma,ncol=J)
    		if(any(colSums(state_wt)==0)) stop("Error: at least one state has an
				 expected number of occurences equal to 0.\n 
					This may be caused by bad starting parameters 
					are insufficent sample size")
		if(anyNA(x) | any(is.nan(x))){
			if(!is.null(new.model$parms.emission$mix.p)){
				new.model$parms.emission = mstep(x,state_wt,estep_mix_weights,new.model$parms.emission,...)
			}else{
				new.model$parms.emission = mstep(x,state_wt,new.model$parms.emission,...)
			}			
		}else{
			if(!is.null(new.model$parms.emission$mix.p)){
				new.model$parms.emission = mstep(x,state_wt,estep_mix_weights,...)
			}else{
				new.model$parms.emission = mstep(x,state_wt,...)
			}
		}
		if(!all(!new.model$semi)){
    		  if(lock.d) {
      		new.model$d = old.model$d
     		 new.model$D = old.model$D
    		  }else {
      		if(sojourn.distribution=="nonparametric") {
        			new.model$d = apply(matrix(estep_variables$eta,ncol=J),2,function(x){if(sum(x)!=0){x/sum(x)}else{x}})
        		new.model$sojourn$d <- new.model$d
      		}else if(sojourn.distribution=="gamma") {
				estep_variables$eta[estep_variables$eta<=0]=1e-300
      			new.model$d = matrix(estep_variables$eta,ncol=J)
      		  	new.model$sojourn$shape = numeric(J)
      		  	new.model$sojourn$scale = numeric(J)
      		  	for(i in 1:J) { 
      		    		if(new.model$semi[i]){
      		    			tmp = .gammafit(1:M,wt=new.model$d[,i])
						new.model$sojourn$shape[i]=tmp$shape
						new.model$sojourn$scale[i]=tmp$scale
      			    		for(u in 1:M){
      		 	    			new.model$d[u,i] = pgamma(u,shape=new.model$sojourn$shape[i],scale=new.model$sojourn$scale[i]) - pgamma(u-1,shape=new.model$sojourn$shape[i],scale=new.model$sojourn$scale[i])
      		 	    			new.model$d[u,i] = new.model$d[u,i]/(pgamma(M,shape=new.model$sojourn$shape[i],scale=new.model$sojourn$scale[i]))
      			    		}# for u
					}# if semi
      			}# for i   
      		}else if(sojourn.distribution=="lnorm") {
				estep_variables$eta[estep_variables$eta<=0]=1e-300
      			new.model$d = matrix(estep_variables$eta,ncol=J)
      		  	new.model$sojourn$meanlog = numeric(J)
      		  	new.model$sojourn$sdlog = numeric(J)
      		  	for(i in 1:J) { 
      		    		if(new.model$semi[i]){
      		    			tmp = .lnormfit(log(1:M),wt=new.model$d[,i])
						new.model$sojourn$shape[i]=tmp$meanlog 
						new.model$sojourn$scale[i]=tmp$sdlog 
      			    		for(u in 1:M){
      		 	    			new.model$d[u,i] = plnorm(u,meanlog=new.model$sojourn$meanlog[i],sdlog=new.model$sojourn$sdlog[i]) - plnorm(u-1,meanlog=new.model$sojourn$meanlog[i],sdlog=new.model$sojourn$sdlog[i])
      		 	    			new.model$d[u,i] = new.model$d[u,i]/(plnorm(M,meanlog=new.model$sojourn$meanlog[i],sdlog=new.model$sojourn$sdlog[i]))
      			    		}# for u
					}# if semi
      			}# for i  
      		}else if(sojourn.distribution=="weibull") {
				estep_variables$eta[estep_variables$eta<=0]=1e-300
      			new.model$d = matrix(estep_variables$eta,ncol=J)
      		  	new.model$sojourn$shape = numeric(J)
      		  	new.model$sojourn$scale = numeric(J)
      		  	for(i in 1:J) { 
      		    		if(new.model$semi[i]){
      		    			tmp = .weibullfit(1:M,wt=new.model$d[,i])
						new.model$sojourn$shape[i]=tmp$shape
						new.model$sojourn$scale[i]=tmp$scale
      			    		for(u in 1:M){
      		 	    			new.model$d[u,i] = pweibull(u,shape=new.model$sojourn$shape[i],scale=new.model$sojourn$scale[i]) - pweibull(u-1,shape=new.model$sojourn$shape[i],scale=new.model$sojourn$scale[i])
      		 	    			new.model$d[u,i] = new.model$d[u,i]/(pweibull(M,shape=new.model$sojourn$shape[i],scale=new.model$sojourn$scale[i]))
      			    		}# for u
					}# if semi
      			}# for i         
      	}else if(sojourn.distribution=="ksmoothed-nonparametric") {
        		new.model$d = apply(matrix(estep_variables$eta+1e-100,ncol=J),2,function(x) x/sum(x))    
        		for(i in 1:J) {
          		new.model$d[,i] = approx(density(which(new.model$d[,i]>ksmooth.thresh),weights=new.model$d[which(new.model$d[,i]>ksmooth.thresh),i],from=1,n=M),xout=1:M)$y
          		new.model$d[is.na(new.model$d[,i]),i] = 0
          		new.model$d[,i] = (new.model$d[,i]+1e-300)/sum(new.model$d[,i])
        		}
        		new.model$sojourn$d <- new.model$d                
      	}else if(sojourn.distribution=="poisson") {
        		new.model$d = apply(matrix(estep_variables$eta,ncol=J),2,function(x){if(sum(x)!=0){x/sum(x)}else{x}})
        		new.model$sojourn$lambda = rep(0,J)
        		new.model$sojourn$shift = rep(0,J)            
        		for(i in 1:J) {
		 		if(new.model$semi[i]){
          			eta = new.model$d[,i]
          			maxshift =  match(TRUE,estep_variables$eta>shiftthresh)
					if(is.na(maxshift)) maxshift=min(NN)
          			Mtmp = tail(which(estep_variables$eta>shiftthresh),1)
					if(length(Mtmp)==0) Mtmp=max(NN)
          			new.model$sojourn$shift[i] = which.max(sapply(1:maxshift, function(shift) .dpois.hhsmm(x = maxshift:Mtmp,lambda=((maxshift:Mtmp)-shift)%*%estep_variables$eta[maxshift:Mtmp],shift=shift,log=TRUE)%*%estep_variables$eta[maxshift:Mtmp]))
          			new.model$sojourn$lambda[i] = ((new.model$sojourn$shift[i]:Mtmp)-new.model$sojourn$shift[i])%*%estep_variables$eta[new.model$sojourn$shift[i]:Mtmp]         
          			new.model$d[,i] = .dpois.hhsmm(1:M,new.model$sojourn$lambda[i],new.model$sojourn$shift[i])
				}
        		}
     	}else if(sojourn.distribution=="nbinom") {       
        		new.model$d = matrix(0,nrow=M,ncol=J)
        		new.model$sojourn$size = numeric(J)
        		new.model$sojourn$shift = integer(J)
        		new.model$sojourn$mu = numeric(J)                            
        		new.model$sojourn$prob = numeric(J)                            
        		eta = matrix(estep_variables$eta,ncol=J)
        		for(i in 1:J) { 
		 		if(new.model$semi[i]){
          			tmp = .nbinomfit(eta[,i])          
          			new.model$sojourn$shift[i] = tmp[1]
          			new.model$sojourn$size[i] =  tmp[2]
          			new.model$sojourn$mu[i] =  tmp[3]
          			new.model$sojourn$prob[i] =  tmp[4]
          			new.model$d[,i] =  .dnbinom.hhsmm(1:M,size=new.model$sojourn$size[i],mu=new.model$sojourn$mu[i],shift=new.model$sojourn$shift[i])
				}
       		}
    		} else if(sojourn.distribution=="logarithmic") {
        		new.model$d = apply(matrix(estep_variables$eta+1e-100,ncol=J),2,function(x) x/sum(x))
        		new.model$sojourn$shape = numeric(J)
        		for(i in 1:J) {
		 		if(new.model$semi[i]){           
          			new.model$sojourn$shape[i] = .logdistrfit(1:M,wt=new.model$d[,i])
          			new.model$d[,i] = .dlog(1:M,new.model$sojourn$shape[i])
				}
       		}
   		} else stop("Invalid sojourn distribution")         
      	new.model$D = apply(new.model$d,2,function(x) rev(cumsum(rev(x))))
    		}
	}
    if(lock.init | any(is.na(estep_variables$init) | is.nan(estep_variables$init))) {
      	new.model$init=old.model$init
    }else {
      	new.model$init=estep_variables$init
      	new.model$init[new.model$init<0]=0
    } 
    if(lock.transition | any(is.na(estep_variables$transition) | is.nan(estep_variables$transition))) {
      	new.model$transition = old.model$transition
    }else {
      	new.model$transition = matrix(estep_variables$transition,ncol=J)
      	new.model$transition[new.model$transition<0]=0
    }
    ll[it]=sum(log(unlist(estep_variables$N)))
    new.model$J = J
	if (verbose) cat("log-likelihood = ",ll[it],"\n")
    if(it>2) if(abs(ll[it]-ll[it-1])/abs(ll[it-1])<tol) break()
 } 
 if(all(!new.model$semi)){
	lock.d = TRUE 
	sojourn.distribution = ""
 }
 nep = length(new.model$parms.emission)
 np = 0
 for(ep in 1:nep) np = np + length(unlist(unlist(new.model$parms.emission[[ep]])))
 df = (!lock.init)*J+ (!lock.transition)*J^2-J + 
	(!lock.d) * (sojourn.distribution=="ksmoothed-nonparametric" | 
	sojourn.distribution=="nonparametric")*M*J + 
	(!lock.d) * (sojourn.distribution!="ksmoothed-nonparametric" & 
	sojourn.distribution!="nonparametric") * (length(new.model$sojourn)-1) +
	np
 AIC = 2*df - 2*ll[it]
 BIC = log(sum(NN))*df - 2*ll[it]
 if (verbose) cat("AIC = ",AIC,"\n")
 if (verbose) cat("BIC = ",BIC,"\n")
 class(new.model) <- "hhsmmspec"
 ret = list(loglik=ll[!is.na(ll)],AIC = AIC, BIC =BIC, 
             model=new.model,
             estep_variables=estep_variables,
             M=M,
             J=J,
             NN=NN,
             f=f,
             mstep=mstep,
             yhat=apply(matrix(estep_variables$gamma,ncol=J),1,which.max))
  class(ret) <- "hhsmm"
  ret
}
