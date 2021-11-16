.dnbinom.hhsmm <- function(x,size,prob=NULL,shift,mu=NULL,log=FALSE) {
  if(shift<0) stop(".dnbinom.hhsmm: shift must be > 0")
  if(is.null(mu)){
  	dnbinom(x-shift,size,prob,log=log)
  } else {
    dnbinom(x-shift,size=size,mu=mu,log=log)
  }
}
