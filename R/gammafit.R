.gammafit <- function(x,wt=NULL) {
	tol = 1e-08
	M = length(x)
	if(is.null(wt)) wt = rep(1,M)
	tmp = cov.wt(data.frame(x),wt=wt)
	xhat = tmp$center
	if(tmp$cov == Inf) tmp$cov = 1e300
	xs = sqrt(tmp$cov)
	s = log(xhat) - mean(weighted.mean(log(x),wt))
	if(is.nan(s) | is.na(s) | !is.finite(s)) s = log(mean(x))
	aold = (xhat/xs)^2
	a = Inf
	while(abs(a-aold)>tol & a > 1e-300) {
    		a = aold - (log(aold) - digamma(aold) - s)/((1/aold) - trigamma(aold))   
		if(is.nan(a) | is.na(a)){
			a = aold
			break
		}	
		if(a<=0){
			a = aold
			break
		}
    		aold=a
  	}
  	l.start = c(a, xhat/a)
  	par0 = if(any(is.nan(l.start) | l.start<=0)) c((xhat/xs)^2,xhat/(xhat/xs)^2) else l.start
  	logd = function(p){
		logd = matrix(nrow = M, ncol = 1)
    		for(u in 1:M){
      		logd[u,1] = log(pgamma(u,shape=exp(p[1]),scale=exp(p[2])) - pgamma(u-1,shape=exp(p[1]),scale=exp(p[2])))
      		logd[u,1] = logd[u,1]-pgamma(M,shape=exp(p[1]),scale=exp(p[2]),log.p=TRUE)
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
