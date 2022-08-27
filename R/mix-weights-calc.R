.mix_weights_clac <- function(model, x, ...)
  {
	f = model$dens.emission
	estep_mix_weights = list()
	for (j in 1:model$J) {
		if (length(model$parms.emission$mix.p[[j]]) > 1) {
			k = length(model$parms.emission$mix.p[[j]])
			estep_mix_weights[[j]] = matrix(nrow = nrow(x), ncol = k)
			for(i in 1:k){
				tmp.model = model
				tmp.model$parms.emission$mix.p[[j]][-i] <- 0
				w = f(x, j, tmp.model, ...) / f(x, j, model, ...)
				w[is.nan(w)] <- 1e-300
				estep_mix_weights[[j]][,i] <- w
			}# for i
		}else{
			estep_mix_weights[[j]] = matrix(1, nrow = nrow(x), ncol = 1)
		}# if 
	}# for j 
	estep_mix_weights
}# end of function 
