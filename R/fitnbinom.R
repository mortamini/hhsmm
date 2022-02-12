.nbinomfit <- function(eta) {  
  shiftthresh=quantile(eta,0.5)
  maxshift =  match(TRUE,eta>shiftthresh)
  Mtmp = tail(which(eta>shiftthresh),1)
  fun1 <- function(shift) {
    m <- weighted.mean((maxshift:Mtmp)-shift,eta[maxshift:Mtmp])
    v <- as.numeric(cov.wt(data.frame((maxshift:Mtmp)-shift),wt=eta[maxshift:Mtmp])$cov)
	if(!is.finite(v)) v = 1e+10
    size <- if (v > m) m^2/(v - m) else 100
	size = trunc(size)
	if(size == 0) size =1
    densfun <- function(par) - sum(dnbinom((maxshift:Mtmp)-shift,size=par[1],mu=par[2],log=TRUE)*eta[maxshift:Mtmp])    
    out<- suppressWarnings(- nlm(densfun,c(size,m))$minimum)
  }
  shift = which.max(sapply(1:maxshift,fun1))
  m <- weighted.mean((maxshift:Mtmp)-shift,eta[maxshift:Mtmp])
  v <- as.numeric(cov.wt(data.frame((maxshift:Mtmp)-shift),wt=eta[maxshift:Mtmp])$cov)
  if(!is.finite(v)) v = 1e+10
  size <- if (v > m) m^2/(v - m) else 100
  size = trunc(size)
  if(size == 0) size =1
  size <- if (v > m) m^2/(v - m) else 100
  densfun <- function(par) -sum(dnbinom((maxshift:Mtmp)-shift,size=par[1],mu=par[2],log=TRUE)*eta[maxshift:Mtmp])        
  suppressWarnings(tmp <- nlm(densfun,c(size,m))$estimate)
  c(shift = shift,size=max(1,trunc(tmp[1])),mu=tmp[2],prob=max(1,trunc(tmp[1]))/(sum(tmp)))
}
