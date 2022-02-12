#' Computing maximum homogeneity of two state sequences 
#'
#' A function to compute the maximum homogeneity of two state sequences.
#' 
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param state.seq1 first state sequence
#' @param state.seq2 second state sequence
#'
#' @return a vector of a length equal to the maximum number of states giving 
#' the maximum homogeneity ratios 
#' 
#' @examples
#' state.seq1 = c(3, 3, 3, 1, 1, 2, 2, 2, 2)
#' state.seq2 = c(2, 2, 2, 3, 3, 1, 1, 1, 1)
#' homogeneity(state.seq1, state.seq2)
#' 
#' @export
homogeneity <- function(state.seq1, state.seq2)
{
	state.seq1 <- as.factor(state.seq1)
	state.seq2 <- as.factor(state.seq2)
	if(nlevels(state.seq2) > nlevels(state.seq1)){
		tmp <- state.seq1
		state.seq1 <- state.seq2
		state.seq2 <- tmp
	}
	if (length(state.seq2) != length(state.seq1)){
		stop("Different sequence lengths !")
	}
	K <- max(nlevels(state.seq1), nlevels(state.seq2))
	ind1 = ind2 = list()
	for(i in 1:K){
		ind1[[i]] <- which(state.seq1 == i)
		ind2[[i]] <- which(state.seq2 == i)
	}
	comp = matrix(0, K, K)
	for (i in 1:K) {
		for (j in 1:K) {
			comp[i, j] = sum(ind2[[i]] %in% ind1[[j]])
		}
	}
	compn <- comp
	maxh = indm = numeric(K)
	for (j in 1:(K-1)) {
		maxhin <- which(compn == max(compn), arr.ind = TRUE)[1,]
		indm[j] <- maxhin[2]
		maxh[indm[j]] <- max(compn)
		compn[maxhin[1], ] <- compn[, maxhin[2]] <- rep(-1, K)
	}
	maxh[setdiff(1:K,indm)] <- max(compn)
	homogeneity <- maxh / sapply(1:K, function(i) length(ind1[[i]]))
	return(homogeneity)
}