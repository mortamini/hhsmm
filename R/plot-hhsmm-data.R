#' @importFrom grDevices dev.new
#'
#' @export
#'
plot.hhsmmdata <- function (x, ...) 
{
	opar <- par(mfrow=c(1,1),no.readonly = TRUE)
	N = x$N
	Ns = c(0,cumsum(N))
	xx = as.matrix(x$x)
	d = ncol(xx)
	if(anyNA(xx) | any(is.nan(xx))){
		allmiss = which(apply(xx,1,function(t) all(is.na(t)|is.nan(t))))
		notallmiss = which(!apply(xx,1,function(t) all(is.na(t)|is.nan(t))))
		for(ii in allmiss){
			neigh = notallmiss[order(abs(ii-notallmiss))[1:2]]
			xx[ii,] = (xx[neigh[1],]+xx[neigh[2],])/2
		}
		if(ncol(xx)>1) xx = complete(mice(xx,printFlag=FALSE))
	}
	for(j in 1:d){
		for(i in 1:length(N)){
			if(d * length(N) <= 9){
				q1 = trunc(sqrt(d * length(N)))
				tmp = d * length(N) / q1
				q2 = trunc(tmp) 
				if(tmp > q2) q2 = q2 + 1
				if(i ==1 & j==1){
					dev.new(width = 70*q2, height = 35*q1, unit = "px")
 					par(mfrow = c(q1,q2))
					on.exit(par(opar))
				}
			} else{
				dev.new(width=10, height=5, unit="in")
			}
			xxx = xx[(Ns[i]+1):Ns[i+1],j]
			sc = 1
			if(!is.null(x$s)) sc = 1.2
			plot(ts(xxx), xlab = "Time",
			ylab = "Observations", 
			main = paste("Sequence ", i," variable ",j),
			ylim = c(min(xxx)/sc , max(xxx)), ...)
    			if (!is.null(x$s)) 
        		.add.states(x$s[(Ns[i]+1):Ns[i+1]], ht = axTicks(2)[1], time.scale = 1)
		}
	}
}