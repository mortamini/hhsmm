.kregs<-function(x, centers, iter.max = 10, nstart = 1, resp.ind = 1){
	n<-nrow(x)
	p<-ncol(x)
	if(centers>n) stop("the number of centers must be less than the number of observations!")	
	if(centers==n){
		out<-list(
		cluster = 1:n,
		centers = t(x[,resp.ind]),
		withinss = rep(0,n))
		return(out)
	}
	tmp.clus = matrix(nrow = n, ncol = nstart)
	sumdist = matrix(nrow = centers, ncol = nstart)
	for(re in 1:nstart){
		initclus = kmeans(x,centers)$cluster
		for(iter in 1:iter.max){
			coefj = list()
			for(j in 1:centers){
				xj = x[initclus==j,-resp.ind]
				yj = x[initclus==j,resp.ind]
				xj = cbind(1,xj)
				coefj[[j]] = as.matrix(ginv(t(xj)%*%xj)%*%t(xj)%*%yj)
			}#for j
			sumdist[,re] = rep(0,centers)
			for(i in 1:n){
				xi = x[i,-resp.ind]
				xi = c(1,xi)
				yi = x[i,resp.ind]
				yihats = lapply(coefj,function(bet) as.vector(t(bet)%*%xi))
				dists = lapply(yihats,function(yhat) (yi - yhat)^2)
				dists = matrix(unlist(dists),centers,length(resp.ind),byrow=TRUE)
				if(length(resp.ind)>1) dists = rowSums(dists)
				tmp.clus[i,re] = which.min(dists)-> iclus
				sumdist[iclus,re] = sumdist[iclus,re] + min(dists)
			}#for i
			if(all(tmp.clus[,re] == initclus)) break
			initclus = tmp.clus[,re]
		}#for iter 
	}#for re
	clus = tmp.clus[,best<-which.min(colSums(sumdist))]
	cent = sapply(1:centers,function(j) rowMeans(t(x[clus==j,resp.ind])))
	if(length(resp.ind)>1){
		colnames(cent)<-paste("Clust",1:centers)
	}else names(cent)<-paste("Clust",1:centers)
	withinss = sumdist[,best]
	names(withinss)<-paste("Clust",1:centers)
	out<-list(cluster = clus,centers=cent,withinss = withinss)
	return(out)
}