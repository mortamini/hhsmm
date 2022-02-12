.lnormfit <- function(x,wt=NULL) {
	M = length(x)
	if(is.null(wt)) wt = rep(1,M)
	tmp = cov.wt(data.frame(x),wt=wt)
	xhat = tmp$center
	if(tmp$cov == Inf) tmp$cov = 1e300
	xs = sqrt(tmp$cov)
  	l.start = c(xhat, xs)
  	par0 = if(any(is.nan(l.start) | l.start<=0)) c(mean(x),sd(x)) else l.start
  	logd = function(p){
		logd = matrix(nrow = M, ncol = 1)
    		for(u in 1:M){
      		logd[u,1] = log(plnorm(u,meanlog=exp(p[1]),sdlog=exp(p[2])) - plnorm(u-1,meanlog=exp(p[1]),sdlog=exp(p[2])))
      		logd[u,1] = logd[u,1]-log(plnorm(M,meanlog=exp(p[1]),sdlog=exp(p[2])))
    		}
		logd
  	}
  	Quasi.loglike = function(p) -t(wt)%*%logd(p)
  	newpar = tryCatch({
		suppressWarnings(fit <- nlm(Quasi.loglike,log(par0)))
		c(fit$estimate,fit$code)
		},error=function(e){
		c(log(par0),1)
		})
  	if(newpar[3]<=2){
      	meanlog = exp(newpar[1])
      	sdlog = exp(newpar[2])
  	}else{
     	meanlog = par0[1]
      	sdlog = par0[2]
  	}	
  	return(list(meanlog = meanlog, sdlog = sdlog))
}
