.dfcalc <- function(model,par,M)
  {
	nep = length(model$parms.emission)
	np = 0
 	for(ep in 1:nep){
		newpar = unlist(unlist(model$parms.emission[[ep]]))
		if (is.matrix(newpar)) {
				if(isSymmetric(newpar)) 
					len = sqrt(length(newpar)) * (sqrt(length(newpar))+1) / 2
		} else {
			len = length(newpar)
		}
		np = np + len
	}
 	df = (!par$lock.init) * model$J + (!par$lock.transition) * 
		(model$J ^ 2 - model$J)  + (!par$lock.d) + np
	if(!is.null(model$sojourn$type)){
		df = df + 
		(model$sojourn$type == "ksmoothed-nonparametric" | 
		model$sojourn$type == "nonparametric") * M * (model$J - 1) + 
		(!par$lock.d) * (model$sojourn$type != "ksmoothed-nonparametric" & 
		model$sojourn$type != "nonparametric") * (length(model$sojourn) - 1)
	} 
	df
}
