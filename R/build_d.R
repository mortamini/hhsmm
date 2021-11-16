.build_d <- function(model,M) {
  model$d = matrix(0,nrow=M,ncol=model$J) 
  if(!all(!model$semi)){
  sojourn.distribution=model$sojourn$type
  J = model$J
  if(sojourn.distribution=='nonparametric' | sojourn.distribution=="ksmoothed-nonparametric") {
    if(!is.null(model$sojourn$d)) {
      if(ncol(model$sojourn$d)!=J) stop("ncol(model$d)!=J")
      M = nrow(model$sojourn$d)
      model$d = model$sojourn$d
    }else stop("Sojourn distribution (model$sojourn$d) not specified.")
  }
  if(sojourn.distribution=="poisson") {
    if(is.null(model$sojourn$d)) {
      if(is.null(model$sojourn$lambda)) stop('Invalid waiting parameters supplied')
      if(is.null(model$sojourn$shift)) stop('No shift parameter provided for Poisson sojourn distribution (must be at least 1)')
      model$d = matrix(0,nrow=M,ncol=model$J)  	
      for(i in 1:J){
		 if(model$semi[i]) model$d[,i] = .dpois.hhsmm(1:M,model$sojourn$lambda[i],model$sojourn$shift[i])
	  }
    } else
      for(i in 1:J) model$d = model$sojourn$d
  }
 if(sojourn.distribution=="lnorm") {
    if(is.null(model$sojourn$d)) {
      if(is.null(model$sojourn$meanlog) | is.null(model$sojourn$sdlog)) stop('Invalid waiting parameters supplied')
      model$d = matrix(0,nrow=M,ncol=model$J)  	
      for(i in 1:J){
 		if(model$semi[i]){
        		for(u in 1:M){
          		model$d[u,i] = plnorm(u,meanlog=model$sojourn$meanlog[i],sdlog=model$sojourn$sdlog[i]) - plnorm(u-1,meanlog=model$sojourn$meanlog[i],sdlog=model$sojourn$sdlog[i])
          		model$d[u,i] = model$d[u,i] / plnorm(M,meanlog=model$sojourn$meanlog[i],sdlog=model$sojourn$sdlog[i])
        		}# for u
        }# if 
      }# for i
    }else
      for(i in 1:J) model$d = model$sojourn$d
  }
  if(sojourn.distribution=="gamma") {
    if(is.null(model$sojourn$shape) | is.null(model$sojourn$scale)) {
      if(is.null(model$sojourn$d))
        stop('Invalid waiting parameters supplied')
      else model$d = model$sojourn$d
    }else {    
      model$d = matrix(0,nrow=M,ncol=model$J)
      for(i in 1:J) {
 		if(model$semi[i]){
        		for(u in 1:M){
          		model$d[u,i] = pgamma(u,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]) - pgamma(u-1,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i])
          		model$d[u,i] = model$d[u,i] / (pgamma(M,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]))
        		}# for u
        }# if 
      }# for i
    }# if else null 
  }  # if gamma  
  if(sojourn.distribution=="weibull") {
    if(is.null(model$sojourn$shape) | is.null(model$sojourn$scale)) {
      if(is.null(model$sojourn$d))
        stop('Invalid waiting parameters supplied')
      else model$d = model$sojourn$d
    }else {    
      model$d = matrix(0,nrow=M,ncol=model$J)
      for(i in 1:J) {
 		if(model$semi[i]){
        		for(u in 1:M){
          		model$d[u,i] = pweibull(u,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]) - pweibull(u-1,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i])
          		model$d[u,i] = model$d[u,i] / (pweibull(M,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]))
        		}# for u
        }# if 
      }# for i
    }# if else null 
  }  # if weibull  
  if(sojourn.distribution=="logarithmic") {
    if(is.null(model$sojourn$shape)) {
      if(is.null(model$sojourn$d))
        stop('Invalid waiting parameters supplied')
      else model$d = model$sojourn$d
    } else {    
      model$d = matrix(0,nrow=M,ncol=model$J)
      	for(i in 1:J){
			if(model$semi[i]) model$d[,i] = .dlog(1:M,model$sojourn$shape[i])
	   	}
    }
  }  
  if (sojourn.distribution == "nbinom") {
    model$d = matrix(0,nrow = M, ncol = J)
    if(is.null(model$sojourn$mu))    for (i in 1:J){
			 								if(model$semi[i]) model$d[, i] = .dnbinom.hhsmm(1:M,size=model$sojourn$size[i],prob=model$sojourn$prob[i],shift=model$sojourn$shift[i])
									}
    if(is.null(model$sojourn$prob))    for (i in 1:J){
			 								if(model$semi[i]) model$d[, i] = .dnbinom.hhsmm(1:M,size=model$sojourn$size[i],mu=model$sojourn$mu[i],shift=model$sojourn$shift[i])
									}
  }
  model$d=head(model$d,M)
  }
  model$D = apply(model$d,2,function(x) rev(cumsum(rev(x))))
  model
}
