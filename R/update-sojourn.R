.udpate_sojourn <- function(model, sojourn, estep, M, NN) 
   {
  	ksmooth.thresh = 1e-20 #this is a threshold for which d(u) values to use 
	shiftthresh = 1e-20 #threshold for effective "0" when considering d(u)
	J <- model$J
	if (sojourn == "nonparametric") {
        	model$d = apply(matrix(estep$eta, ncol = J), 2, 
				function(x) {if (sum(x) != 0){ x/sum(x) } else {x}})
        	model$sojourn$d <- model$d
     }else if (sojourn == "gamma") {
		estep$eta[estep$eta <= 0] <- 1e-300
      	model$d = matrix(estep$eta, ncol = J)
      	model$sojourn$shape = numeric(J)
         model$sojourn$scale = numeric(J)
      	for (i in 1:J) { 
      		if(model$semi[i]){
      		    	tmp = .gammafit(1:M, wt = model$d[, i])
				model$sojourn$shape[i]=tmp$shape
				model$sojourn$scale[i]=tmp$scale
      			for(u in 1:M){
      		 	    	model$d[u,i] = pgamma(u,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]) - pgamma(u-1,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i])
      		 	    	model$d[u,i] = model$d[u,i]/(pgamma(M,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]))
      			}# for u
			}# if semi
      	}# for i   
      }else if(sojourn == "lnorm") {
		estep$eta[estep$eta <= 0] <- 1e-300
      	model$d = matrix(estep$eta, ncol = J)
      	model$sojourn$meanlog = numeric(J)
      	model$sojourn$sdlog = numeric(J)
      	for (i in 1:J) { 
      		if(model$semi[i]) {
      		    	tmp = .lnormfit(log(1:M), wt = model$d[, i])
				model$sojourn$shape[i] = tmp$meanlog 
				model$sojourn$scale[i] = tmp$sdlog 
      			for(u in 1:M){
      		 	    model$d[u,i] = plnorm(u,meanlog=model$sojourn$meanlog[i],sdlog=model$sojourn$sdlog[i]) - plnorm(u-1,meanlog=model$sojourn$meanlog[i],sdlog=model$sojourn$sdlog[i])
      		 	    model$d[u,i] = model$d[u,i]/(plnorm(M,meanlog=model$sojourn$meanlog[i],sdlog=model$sojourn$sdlog[i]))
      			}# for u
			}# if semi
      	}# for i  
	}else if(sojourn == "weibull") {
		estep$eta[estep$eta <= 0] <- 1e-300
      	model$d = matrix(estep$eta, ncol = J)
      	model$sojourn$shape = numeric(J)
      	model$sojourn$scale = numeric(J)
      	for (i in 1:J) { 
      		if(model$semi[i]){
      		    tmp = .weibullfit(1:M, wt = model$d[,i])
				model$sojourn$shape[i] = tmp$shape
				model$sojourn$scale[i] = tmp$scale
      			for(u in 1:M){
      		 	   model$d[u,i] = pweibull(u,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]) - pweibull(u-1,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i])
      		 	   model$d[u,i] = model$d[u,i]/(pweibull(M,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i]))
      			}# for u
			}# if semi
      	}# for i         
	}else if (sojourn == "ksmoothed-nonparametric") {
     	model$d = apply(matrix(estep$eta+ 1e-100, ncol = J), 2, 
				function(x) x / sum(x))    
        	for (i in 1:J) {
         	model$d[,i] = approx(density(which(model$d[,i] > ksmooth.thresh),weights = model$d[which(model$d[,i] > ksmooth.thresh),i],from=1,n=M),xout=1:M)$y
          	model$d[is.na(model$d[,i]),i] <- 0
          	model$d[,i] = (model$d[,i] + 1e-300) / sum(model$d[,i])
        	}
        	model$sojourn$d <- model$d                
    }else if(sojourn == "poisson") {
        	model$d = apply(matrix(estep$eta,ncol = J), 2,
			function(x){if(sum(x)!=0){x/sum(x)}else{x}})
        	model$sojourn$lambda = rep(0,J)
        	model$sojourn$shift = rep(0,J)            
        	for(i in 1:J) {
		 	if(model$semi[i]){
          		eta = model$d[,i]
          		maxshift =  match(TRUE,estep$eta>shiftthresh)
				if(is.na(maxshift)) maxshift=min(NN)
          		Mtmp = tail(which(estep$eta>shiftthresh),1)
				if(length(Mtmp)==0) Mtmp=max(NN)
          		model$sojourn$shift[i] = which.max(sapply(1:maxshift, function(shift) .dpois.hhsmm(x = maxshift:Mtmp,lambda=((maxshift:Mtmp)-shift)%*%estep$eta[maxshift:Mtmp],shift=shift,log=TRUE)%*%estep$eta[maxshift:Mtmp]))
          		model$sojourn$lambda[i] = ((model$sojourn$shift[i]:Mtmp)-model$sojourn$shift[i])%*%estep$eta[model$sojourn$shift[i]:Mtmp]         
          		model$d[,i] = .dpois.hhsmm(1:M,model$sojourn$lambda[i],model$sojourn$shift[i])
			}
        	}
    }else if (sojourn == "nbinom") {       
        	model$d = matrix(0,nrow = M, ncol = J)
        	model$sojourn$size = numeric(J)
        	model$sojourn$shift = integer(J)
        	model$sojourn$mu = numeric(J)                            
        	model$sojourn$prob = numeric(J)                            
        	eta = matrix(estep$eta, ncol=J)
        	for (i in 1:J) { 
		 	if(model$semi[i]){
          		tmp = .nbinomfit(eta[,i])          
          		model$sojourn$shift[i] = tmp[1]
          		model$sojourn$size[i] =  tmp[2]
          		model$sojourn$mu[i] =  tmp[3]
          		model$sojourn$prob[i] =  tmp[4]
          		model$d[,i] =  .dnbinom.hhsmm(1:M,size=model$sojourn$size[i],mu=model$sojourn$mu[i],shift=model$sojourn$shift[i])
			}
       	}
    	} else if(sojourn == "logarithmic") {
        	model$d = apply(matrix(estep$eta+1e-100,ncol=J),2,function(x) x/sum(x))
        	model$sojourn$shape = numeric(J)
        	for(i in 1:J) {
		 	if(model$semi[i]){           
          		model$sojourn$shape[i] = .logdistrfit(1:M,wt=model$d[,i])
          		model$d[,i] = .dlog(1:M,model$sojourn$shape[i])
			}
       	}
   	} else stop("Invalid sojourn distribution")
    model$D = apply(model$d,2,function(x) rev(cumsum(rev(x))))
	model
}
   