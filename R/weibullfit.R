.weibullfit <- function(x,wt=NULL) {
	M = length(x)
	if(is.null(wt)) wt = rep(1,M)
	l.start = c()
	mlequi <- function(alpha) weighted.mean(x^alpha*log(x),wt)/weighted.mean(x^alpha,wt)-1/alpha-weighted.mean(log(x),wt)
	lower = 1e-10
	upper = lower
	while(mlequi(upper)*mlequi(lower)>0){
		upper = upper*10
		if(is.na(mlequi(upper))) upper = upper/15
	}
	l.start[1]<-uniroot(function(t){sapply(t,function(a){mlequi(a)})},c(lower,upper))$root
	l.start[2] <-(weighted.mean(x^l.start[1],wt))^(1/l.start[1])
  	par0 = if(any(is.nan(l.start) | l.start<=0)) c(mean(x)/var(x),mean(x)/gamma(1+mean(x)/var(x))) else l.start
  	logd = function(p){
		logd = matrix(nrow = M, ncol = 1)
    		for(u in 1:M){
      		logd[u,1] = log(pweibull(u,shape=exp(p[1]),scale=exp(p[2])) - pweibull(u-1,shape=exp(p[1]),scale=exp(p[2])))
      		logd[u,1] = logd[u,1]-pweibull(M,shape=exp(p[1]),scale=exp(p[2]),log.p=TRUE)
    		}
		logd
  	}
  	Quasi.loglike = function(p) -t(wt)%*%logd(p)
  	newpar = tryCatch({
		fit=nlm(Quasi.loglike,log(par0))
		c(fit$estimate,fit$code)
		},error=function(e){
		c(log(par0),1)
		})
  	if(newpar[3]<=2){
      	shape = exp(newpar[1])
      	scale = exp(newpar[2])
  	}else{
     	shape = par0[1]
      	scale = par0[2]
  	}	
  	return(list(shape=shape,scale=scale))
}
