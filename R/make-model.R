#' make a hhsmmspec model for a specified emission distribution
#'
#' Provides a hhsmmspec model by using the parameters  
#' obtained by \code{\link{initial_estimate}} for the emission distribution 
#' characterized by mstep and dens.emission
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param par the parameters obtained by \code{\link{initial_estimate}}
#' @param mstep the mstep function of the EM algorithm with an style simillar to that of \code{\link{mixmvnorm_mstep}}
#' @param dens.emission the density of the emission distribution with an style simillar to that of \code{\link{dmixmvnorm}}
#' @param semi logical and of one of the following forms:
#' \itemize{
#' \item a logical value: if TRUE all states are considered as semi-Markovian else Markovian
#' \item a logical vector of length nstate: the TRUE associated states are considered as semi-Markovian
#' and FALSE associated states are considered as Markovian
#' \item \code{NULL}{ if \code{ltr}=TRUE then \code{semi = c(rep(TRUE,nstate-1),FALSE)}, else 
#' \code{semi = rep(TRUE,nstate)}}
#' }
#' @param M maximum number of waiting times in each state
#' @param sojourn the sojourn time distribution which is one of the following cases:
#' \itemize{
#' \item \code{"nonparametric"}{ non-parametric sojourn distribution}
#' \item \code{"nbinom"}{ negative binomial sojourn distribution}
#' \item \code{"logarithmic"}{ logarithmic sojourn distribution}
#' \item \code{"poisson"}{ poisson sojourn distribution}
#' \item \code{"gamma"}{ gamma sojourn distribution}
#' \item \code{"weibull"}{ weibull sojourn distribution}
#' \item \code{"lnorm"}{ log-normal sojourn distribution}
#' \item \code{"auto"}{ automatic determination of the sojourn distribution using the chi-square test}
#' }
#'
#' @return a \code{\link{hhsmmspec}} model containing the following items:
#' \itemize{
#' \item \code{init}{ initial probabilities of states}
#' \item \code{transition}{ transition matrix}
#' \item \code{parms.emission}{ parameters of the mixture normal emission (\code{mu}, \code{sigma}, \code{mix.p})}
#' \item \code{sojourn}{ list of sojourn distribution parameters and its \code{type}}
#' \item \code{dens.emission}{ the emission probability density function}
#' \item \code{mstep}{ the M step function of the EM algorithm}
#' \item \code{semi}{ a logical vector of length nstate with the TRUE associated states are considered as semi-Markovian}
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
#' model = make_model(par,semi=NULL,M=max(train$N),sojourn="gamma")
#'
#' @export
#'
make_model<-function(par,mstep=mixmvnorm_mstep,dens.emission=dmixmvnorm,semi=NULL,M,sojourn){
	J = length(par$emission[[1]])
	nmix = par$nmix
	if(par$ltr){
		init <-c(1,rep(0,J-1))
	}else{
		init <-rep(1/J,J)
	} 
	if(is.null(semi)){
		if(par$ltr){
			semi = c(rep(TRUE,J-1),FALSE)
		}else{
			semi = rep(TRUE,J)
		}
	}
	B0 = par$emission
	d0 = matrix(0,nrow=M,ncol=J)
	lenn=list()
	for(j in 1:(J-par$ltr)){
		temp.lenn = c()
		if(par$ltr){
			temp.lenn = c(temp.lenn,par$leng[[j]])
		} else {
			clusters = par$state.clus
			num.units = length(clusters)
			for(m in 1:num.units){
				temp.lenn.m = c()
				runs = FALSE
				ii = 2
				while(ii<length(clusters[[m]])){
					add = FALSE
					if((clusters[[m]][ii] == clusters[[m]][ii-1]) & (clusters[[m]][ii]==j)){
						runs = TRUE
						cntr = 1 
						add = TRUE
					}
					while(runs & (ii<length(clusters[[m]]))){
						cntr = cntr+1
						ii = ii + 1
						runs = (clusters[[m]][ii] == clusters[[m]][ii-1])  	
					}
					if(add) temp.lenn.m = c(temp.lenn.m,cntr)
					ii = ii + 1
				}#while
				temp.lenn.m = sort(c(temp.lenn.m,rep(1,sum(clusters[[m]]== j)-sum(temp.lenn.m))))
				temp.lenn = c(temp.lenn,temp.lenn.m)				
			}# for m 
		}# if else ltr
		temp.lenn[temp.lenn==0]=1e-10
		lenn[[j]] = temp.lenn
	}# for j
  if(!all(!semi)){
	if(sojourn=="auto"){
		g.shape = g.scale = meanlog = sdlog = lambda = shift = 
	  	w.shape = w.scale = mu =  size = l.shape = numeric(J)
		pvalues = matrix(nrow = (J-par$ltr), ncol = 6)
		for(j in 1:(J-par$ltr)){	
			ps = table(lenn[[j]])/sum(table(lenn[[j]]))
			ind = as.numeric(names(ps))
			ind[ind==1e-10]=1
			d0j = rep(0,M)
			d0j[ind]=ps
			g.shape[j] = (mean(lenn[[j]])^2)/var(lenn[[j]])
		  	g.scale[j] = var(lenn[[j]])/mean(lenn[[j]])
			meanlog[j] = mean(log(lenn[[j]]))
			sdlog[j] = sd(log(lenn[[j]]))	
	  		shift[j] = min(lenn[[j]])
	  		lambda[j] = mean(lenn[[j]]-shift[j])
	  		mu[j] = mean(lenn[[j]]-shift[j])
			size[j] = round(mean((lenn[[j]]-shift[j])^2)/(var(lenn[[j]])-mean(lenn[[j]]-shift[j])))
		  	mean.time = mean(lenn[[j]])
  			fn <- function(p) mean.time + p/((1-p)*log(1-p))
  			l.shape[j] = uniroot(fn,c(1e-10,1-1e-10))$root
			mlequi <- function(alpha)	sum(lenn[[j]]^alpha*log(lenn[[j]]))/sum(lenn[[j]]^alpha)-1/alpha-mean(log(lenn[[j]]))
			lower = 1e-10
			upper = lower
			while(mlequi(upper)*mlequi(lower)>0){
				upper = upper*10
				if(is.na(mlequi(upper)) | is.nan(mlequi(upper))) upper = upper/15
			}
			w.shape[j]<-uniroot(function(t){sapply(t,function(a){mlequi(a)})},c(lower,upper))$root
			w.scale[j] <-(mean(lenn[[j]]^w.shape[j]))^(1/w.shape[j])
			d0jc = cumsum(d0j)
			breaks = unique(c(0,sapply(1:10,function(t) min(which(d0jc>=0.1*t)))))
			for(l in 1:length(breaks)){
				if(!is.finite(breaks[l])){
					breaks[l] = length(d0jc)
					breaks = breaks[1:l]
					break
				}
			}
			observed = expected.gamma = expected.lnorm = expected.poisson = expected.weibull =
			expected.nbinom = expected.log = c()
			M0 = breaks[length(breaks)]
			for(l in 1:(length(breaks)-1)){
				observed[l] = sum(d0j[breaks[l]:(breaks[l+1]-1)])
				expected.gamma[l] = pgamma(breaks[l+1],shape=g.shape[j],scale=g.scale[j])-pgamma(breaks[l],shape=g.shape[j],scale=g.scale[j])
				expected.gamma[l] = abs(expected.gamma[l]/pgamma(M0,shape=g.shape[j],scale=g.scale[j]))
				expected.weibull[l] = pweibull(breaks[l+1],shape=w.shape[j],scale=w.scale[j])-pweibull(breaks[l],shape=w.shape[j],scale=w.scale[j])
				expected.weibull[l] = abs(expected.weibull[l]/pweibull(M0,shape=w.shape[j],scale=w.scale[j]))
				expected.lnorm[l] = plnorm(breaks[l+1],meanlog=meanlog[j],sdlog=sdlog[j])-plnorm(breaks[l],meanlog=meanlog[j],sdlog=sdlog[j])
				expected.lnorm[l] = abs(expected.lnorm[l]/plnorm(M0,meanlog=meanlog[j],sdlog=sdlog[j]))
				expected.poisson[l] = ppois(breaks[l+1]-shift[j],lambda[j])-ppois(breaks[l]-shift[j],lambda[j])
				expected.poisson[l] = abs(expected.poisson[l]/ppois(M0-shift[j],lambda[j]))
				expected.nbinom[l] = pnbinom(breaks[l+1]-shift[j],size=size[j],mu=mu[j])-pnbinom(breaks[l]-shift[j],size=size[j],mu=mu[j])
				expected.nbinom[l] = abs(expected.nbinom[l]/pnbinom(M0-shift[j],size=size[j],mu=mu[j]))
				expected.log[l] = .plog(breaks[l+1]-shift[j],l.shape[j])-.plog(breaks[l]-shift[j],l.shape[j])
				expected.log[l] = abs(expected.log[l]/.plog(M0-shift[j],l.shape[j]))
			}
			pvalues[j,1] = chisq.test(x = observed,p = expected.gamma)$p.value
			pvalues[j,2] = chisq.test(x = observed,p = expected.lnorm)$p.value
			pvalues[j,3] = chisq.test(x = observed,p = expected.poisson)$p.value			
			pvalues[j,4] = chisq.test(x = observed,p = expected.weibull)$p.value		
			pvalues[j,5] = chisq.test(x = observed,p = expected.nbinom)$p.value
			pvalues[j,6] = chisq.test(x = observed,p = expected.log)$p.value								
		}# for j 
		pvals = apply(pvalues, 2 , min)
		if(max(pvals)>0.05){
			dist = which.max(pvals)
			sojourn = c("gamma","lnorm","poisson","weibull","nbinom","logarithmic")[dist]
		} else {
			sojourn = "nonparametric"
		}
	}# if auto 
	if(sojourn=="nonparametric"){
		for(j in 1:(J-par$ltr)){	
			ps = table(lenn[[j]])/sum(table(lenn[[j]]))
			ind = as.numeric(names(ps))
			ind[ind==1e-10]=1
			d0j = rep(0,M)
			d0j[ind]=ps
			d0[,j]=d0j
		}
		sojourn.list = list(d=d0,type=sojourn)
	}else if(sojourn=="gamma"){
		shape = scale = numeric(J)
		for(j in 1:(J-par$ltr)){
			shape[j] = (mean(lenn[[j]])^2)/var(lenn[[j]])
		  	scale[j] = var(lenn[[j]])/mean(lenn[[j]])
		}
		sojourn.list = list(shape = shape , scale = scale, type=sojourn)
	}else if (sojourn == "weibull"){
	  shape = scale = numeric(J)
	  for(j in 1:(J-par$ltr)){
		mlequi <- function(alpha) sum(lenn[[j]]^alpha*log(lenn[[j]]))/sum(lenn[[j]]^alpha)-1/alpha-mean(log(lenn[[j]]))
		lower = 1e-10
		upper = lower
		while(mlequi(upper)*mlequi(lower)>0){
			upper = upper*10
			if(is.na(mlequi(upper)) | is.nan(mlequi(upper))) upper = upper/15
		}
		shape[j]<-uniroot(function(t){sapply(t,function(a){mlequi(a)})},c(lower,upper))$root
		scale[j] <-(mean(lenn[[j]]^shape[j]))^(1/shape[j])
	  }
	  sojourn.list = list(shape = shape,scale = scale, type=sojourn)
   	}else if (sojourn == "lnorm"){
  		meanlog = sdlog = numeric(J)
   	  	for(j in 1:(J-par$ltr)){
			meanlog[j] = mean(log(lenn[[j]]))
			sdlog[j] = sd(log(lenn[[j]]))			
   	  	}
   	  	sojourn.list = list(meanlog = meanlog, sdlog = sdlog, type=sojourn) 
	}else if(sojourn=="poisson"){
	  	lambda = shift = numeric(J)
	  	for(j in 1:(J-par$ltr)){
	  		shift[j] = min(lenn[[j]])
	  		lambda[j] = mean(lenn[[j]]-shift[j])
	  	}
	  	sojourn.list = list(lambda = lambda, shift = shift, type=sojourn)
	}else if(sojourn=="nbinom"){
	  	mu = shift = size = numeric(J)
	  	for(j in 1:(J-par$ltr)){
	  		shift[j] = min(lenn[[j]])
	  		mu[j] = mean(lenn[[j]]-shift[j])
			size[j] = round(mean((lenn[[j]]-shift[j])^2)/(var(lenn[[j]])-mean(lenn[[j]]-shift[j])))
	  	}
	  	sojourn.list = list(shift = shift, mu = mu, size = size, type=sojourn)
	}else if(sojourn=="logarithmic"){
	  	shape = shift = numeric(J)
	  	for(j in 1:(J-par$ltr)){
		  	mean.time = mean(lenn[[j]])
  			fn <- function(p) mean.time + p/((1-p)*log(1-p))
  			shape[j] = uniroot(fn,c(1e-10,1-1e-10))$root
			shift[j] = min(lenn[[j]])
	  	}
	  	sojourn.list = list(shape = shape, shift = shift, type=sojourn)
	}else{
		stop(paste("sojourn type ",sojourn, " in unknown"))
	}#if
  } else{
	sojourn.list = NULL
  }# if else hmm
	P=matrix(0,nrow=J,ncol=J)
	if(par$ltr & all(semi[-J])){
		lenn.m = matrix(nrow=J-1,ncol=length(lenn[[1]]))
		for(k in 1:(J-1)) lenn.m[k,]=lenn[[k]]
		lennz = c()
		for(k in 1:(J-1)){
			lennz[k] = 0
			for(cl in 1:ncol(lenn.m)){
				lennz[k] = lennz[k] + (sum(lenn.m[,cl]!=0)==k)
			}# for cl
		}#for k 
		tmp.p = lennz/length(lenn[[k]])
		tmp.p[tmp.p == 0]=0.05
		tmp.p[J-1] = 1-sum(tmp.p[-(J-1)])
	}#if ltr no markov
	if(par$ltr & any(!semi[-J])){
		lenn.m = matrix(nrow=J-1,ncol=length(lenn[[1]]))
		for(k in 1:(J-1)) lenn.m[k,]=lenn[[k]]
		dg = rowSums(lenn.m)/sum(lenn.m)
		lennz = c()
		for(k in 1:(J-1)){
			lennz[k] = 0
			for(cl in 1:ncol(lenn.m)){
				lennz[k] = lennz[k] + (sum(lenn.m[,cl]!=0)==k)
			}# for cl
		}#for k 
		tmp.p = lennz/length(lenn[[k]])
		tmp.p[tmp.p == 0]=0.01
		tmp.p[J-1] = 1-sum(tmp.p[-(J-1)])
	}#if ltr any markov
	if(par$ltr){
		for(i in 1:J){
			if(i<J & !semi[i]){
				tmp.p2 =tmp.p
				tmp.p2[J-1] = 1-dg[i]-sum(tmp.p2[-(J-1)])
				tmp.p2[J] = dg[i]
			}else{
				tmp.p2 = tmp.p
				tmp.p2[J] = 0
			}
			for(j in J:1){
				if(i <= j){
					if(i <= (J-1)){
						P[i,j] = tmp.p2[J-j+i]
					}else{
						P[i,j] = 1
					}# ifelse i< J-1
				}# if i<j
				P[J,J] = 1
				if(semi[J-1]) P[J-1,J]=1			
			}# for j
		}# for i 
		P = P/rowSums(P)
	}else{
		state.clus = par$state.clus
		num.units = length(state.clus)
		P = matrix(0 , J , J)
		for(m in 1:num.units){
			tt <- table( c(state.clus[[m]][-length(state.clus[[m]])]), c(state.clus[[m]][-1]) )
			ttc = matrix(0 , J , J)
			ttc[as.numeric(rownames(tt)),as.numeric(colnames(tt))]<-tt
			ttc[ttc==0]=1e-10
    			P = P + ttc / rowSums(ttc)
		}# for m
		dg = diag(P)
		dg[!semi] = 0
		P = P - diag(dg)
		P = P / rowSums(P)
	}# ifelse ltr 
		model <- hhsmmspec(init=init, transition=P,
			parms.emission = B0,
		  	sojourn=sojourn.list,
			dens.emission = dens.emission,
			mstep = mstep, 
			semi=semi)
	return(model)
}
