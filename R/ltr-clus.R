#' left to right clustering
#'
#' A left to right initial clustering method using the 
#' mean differences and Hotelling's T-squared test
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param Dat a data matrix
#' @param k number of clusters
#'
#' @return a (left to right) clustering vector 
#'
#' @examples
#' data(CMAPSS)
#' clus = ltr_clus(CMAPSS$train$x[1:CMAPSS$train$N[1], ], 3)
#'
#' @export
#'
ltr_clus <- function(Dat, k) {
	n <- nrow(Dat)
	CDat <- list()
	CDat[[1]] <- Dat 
	clus <- rep(1,n)
	n.clus <- 1
	diff <- c()
	change <- TRUE
	while (n.clus < k & change) {
		tmp.clus <- c()
		nc <- n.clus
		ncc <- 0
		for (j in 1:n.clus) {
			tmp <- .ltr_clus2(CDat[[j]])
			if (length(unique(tmp$cluster)) == 2) {
				ncc <- ncc + 1 
				clus1 <- CDat[[j]][tmp$cluster == 1, ]
				clus2 <- CDat[[j]][tmp$cluster == 2, ]
				if (is.null(dim(clus1))) clus1 <- t(as.matrix(clus1))
				if (is.null(dim(clus2))) clus2 <- t(as.matrix(clus2))
				CDat[[j]] <- clus1
				CDat[[n.clus+ncc]] <- clus2
				nc <- nc + 1
				tmp.clus <- c(tmp.clus, tmp$cluster + (j - 1) * 2)
				diff <- c(diff, tmp$mean.diff)
			} else {
				tmp.clus <- c(tmp.clus, tmp$cluster + j - 1)
			}
			clus <- tmp.clus
			clus <- as.numeric(as.factor(clus))
		}
		if (nc == n.clus) change <- FALSE
		n.clus <- nc
	} 	
	if (n.clus > k) {
		dc <- n.clus - k
		merge <- order(diff)[1:dc] + 1
		if (min(merge) > 1 & max(merge) < n.clus) {
			tar <- which.min(c(diff[min(merge) - 1], diff[max(merge)]))
			tar <- c(min(merge) - 1, max(merge))[tar]
		}else if (min(merge) == 1) {
			tar <- max(merge) + 1
		}else{
			tar <- min(merge) - 1 
		}
		clus[clus %in% merge] <- tar
		clus <- as.numeric(as.factor(clus))
	}
	return(clus)
}
