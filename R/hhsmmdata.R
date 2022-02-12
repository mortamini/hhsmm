#' convert to hhsmm data 
#'
#' Converts a matrix of data and its associated vector of sequence lengths
#' to a data list of class \code{"hhsmmdata"} 
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param x a matrix of data
#' @param N a vector of sequence lengths. If NULL then \code{N = nrow(x)}
#'
#' @return a data list of class \code{"hhsmmdata"} containing \code{x} and \code{N}
#'
#'
#' @examples
#' x = sapply(c(1, 2), function(i) rnorm(100, i, i/2))
#' N = c(10, 15, 50, 25)
#' data = hhsmmdata(x, N)
#'
#' @export
#'
hhsmmdata <- function(x, N = NULL){
	if (is.null(N)) N = nrow(x)
	if (nrow(x) != sum(N)) stop("nrow of x != sum(N) !")
	data <- list(x = x, N = N)
	class(data) <- "hhsmmdata"
	data
}