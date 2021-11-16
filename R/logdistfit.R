.logdistrfit <- function(x,wt) {
  xbar = sum(wt*x)
  fn <- function(p) xbar + p/((1-p)*log(1-p))
  uniroot(fn,c(1e-10,1-1e-10))$root
}
