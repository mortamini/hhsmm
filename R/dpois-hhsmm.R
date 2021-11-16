.dpois.hhsmm <- function(x=NULL,lambda,shift,log=FALSE)   {
  if(shift<0) stop(".dpois.hhsmm: shift must be > 0")
  dpois(x-shift,lambda,log=log)
}
