.rnbinom.hhsmm <- function(n,size,prob,shift) {
  if(shift<0) stop(".dnbinom.hhsmm: shift must be > 0")
  rnbinom(n,size,prob) + shift 
}