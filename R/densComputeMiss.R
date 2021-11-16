.densComputeMiss<-function(x,model,...){
			d = ncol(x)
			missed = apply(x,1,function(t) which(is.na(t)|is.nan(t)))
			p = matrix(0,nrow(x),J)
			means = list()
			J = model$J
			for(j in 1:J){
				means[[j]] = list()
				if(!is.null(model$parms.emission$mix.p)){
					k = length(model$parms.emission$mix.p[[j]])
					for(i in 1:k){
						means[[j]][[i]]=sapply(1:length(missed), function(ii){
							l = missed[[ii]]
							if(length(l)==0){NA}else{
								if(length(l) == d){
									model$parms.emission$mu[[j]][[i]][l]
								}else{
									model$parms.emission$sigma[[j]][[i]][l,-l]%*%
										ginv(model$parms.emission$sigma[[j]][[i]][-l,-l])%*%
										(x[ii,-l]-model$parms.emission$mu[[j]][[i]][-l])+
										model$parms.emission$mu[[j]][[i]][l]	
								}
							}
						})
						xr = x
						for(ii in 1:nrow(xr)) xr[ii,is.na(xr[ii,])|is.nan(xr[ii,])] = means[[j]][[i]][[ii]]
						tmp.model = model
						tmp.model$parms.emission$mix.p[[j]][i]=1
						tmp.model$parms.emission$mix.p[[j]][-i]=0
						p[,j] = p[,j] + model$parms.emission$mix.p[[j]][i]*model$dens.emission(xr,j,tmp.model,...)
					}#for i
				}else{
					means[[j]]=sapply(1:length(missed), function(ii){
							l = missed[[ii]]
							if(length(l)==0){NA}else{
								if(length(l) == d){
									model$parms.emission$mu[[j]][l]
								}else{
									model$parms.emission$sigma[[j]][l,-l]%*%
										ginv(model$parms.emission$sigma[[j]][-l,-l])%*%
										(x[ii,-l]-model$parms.emission$mu[[j]][-l])+
										model$parms.emission$mu[[j]][l]	
								}
							}
						})
						xr = x
						for(ii in 1:nrow(xr)) xr[ii,is.na(xr[ii,])|is.nan(xr[ii,])] = means[[j]][[ii]]
						p[,j] = model$dens.emission(xr,j,model,...)
				}#if else mix.p
			}# for j
		p
}