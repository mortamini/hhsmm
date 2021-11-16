.plog <- function(x,p) ifelse(x>=0,1+pbeta(p,x+1,1)*beta(x+1,1)/log(1-p),0)
