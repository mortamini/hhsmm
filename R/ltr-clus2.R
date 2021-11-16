#' @importFrom MASS ginv
#'
.ltr_clus2<-function(Dat){
	if(is.null(Dat)) return(NULL)
	if(is.null(dim(Dat))) Dat = t(as.matrix(Dat))
	p=ncol(Dat)
	n=nrow(Dat)
	dif=c()
	if( n > (p+1) ){
		for(j in 2:(n-2)){
			mu1=apply(matrix(Dat[1:j,],j,p),2,mean)
			mu2=apply(matrix(Dat[(j+1):n],n-j,p),2,mean)
			Sig1=cov(matrix(Dat[1:j,],j,p))
			Sig2=cov(matrix(Dat[(j+1):n],n-j,p))
			Sp=((j-1)*Sig1+(n-j-1)*Sig2)/(n-2)
			dif[j]=(mu2-mu1)%*%ginv(Sp)%*%(mu2-mu1)*(j*(n-j)/n)*(n-p-1)/(n-2)/p
		}
		jstar=which.max(dif)
		difstar = dif[jstar]
		if(difstar > qf(0.95,p,n-1-p)){
			clus=c(rep(1,jstar-1),rep(2,n-jstar+1))
		}else{
			clus=rep(1,n)
		}
	}else{
		clus=rep(1,n); difstar = 0
	}
	return(list(cluster=clus,mean.diff=difstar))
}
