.ltr_reg_clus2<-function(Dat, resp.ind = 1){
	if(is.null(Dat)) return(NULL)
	if(is.null(dim(Dat))) Dat = t(as.matrix(Dat))
	n=nrow(Dat)
	y = as.matrix(Dat[,resp.ind])
	x = Dat[,-resp.ind]
	x = cbind(1,x)
	dx = ncol(x)
	dy = ncol(y)
	dif=c()
	if( n > max(5,(dy+1)) ){
		for(j in 5:(n-5)){
			x1 = matrix(x[1:j,],j,dx)
			x2 = matrix(x[(j+1):n],n-j,dx)
			y1 = matrix(y[1:j,],j,dy)
			y2 = matrix(y[(j+1):n],n-j,dy)
			coef1 = as.matrix(ginv(t(x1)%*%x1)%*%t(x1)%*%y1)			
			coef2 = as.matrix(ginv(t(x2)%*%x2)%*%t(x2)%*%y2)			
			yhat1 = x1%*%coef1
			yhat2 = x2%*%coef2
			res1 = y1-yhat1
			res2 = y2-yhat2
			sig1 = sum(res1^2)/j
			sig2 = sum(res2^2)/(n-j)
			Sig1 = sig1*ginv(t(x1)%*%x1)			
			Sig2 = sig2*ginv(t(x2)%*%x2)			
			Sp=(j*Sig1+(n-j)*Sig2)/n
			dif[j]=sum(t(coef1-coef2)%*%ginv(Sp)%*%(coef1-coef2))*(j*(n-j)/n)*(n-dx-1)/n/dx
		}
		jstar=which.max(dif)
		difstar = dif[jstar]
		if(difstar > qf(0.95,dx,n-1-dx)){
			clus=c(rep(1,jstar-1),rep(2,n-jstar+1))
		}else{
			clus=rep(1,n)
		}
	}else{
		clus=rep(1,n); difstar = 0
	}
	return(list(cluster=clus,mean.diff=difstar))
}
