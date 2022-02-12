#' Splitting the data sets to train and test 
#'
#' A function to split the train data of class \code{"hhsmmdata"}
#' to train and test subsets with an option to right trim the sequences
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param train the train data of class \code{"hhsmmdata"}
#' @param train.ratio a number in (0,1] which determines the ratio of the train subset. It can be 
#' equal to 1, if we need the test set to be equal to the train set and we only 
#' need to right trim the sequences
#' @param trim logical. if TRUE the sequences will be right trimmed with random lengths
#' @param trim.ratio a vector of trim ratios with a length equal to that of \code{train$N},
#' or a single trim ratio for all sequences. If it is \code{NULL}, then random trim 
#' ratios will be used
#'
#' @return a list containing:
#' \itemize{
#' \item\code{train}{ the randomly selected subset of train data of class \code{"hhsmmdata"}}
#' \item\code{test}{ the randomly selected subset of test data of class \code{"hhsmmdata"}}
#' \item\code{trimmed}{ right trimmed test subset, if \code{trim}=TRUE, with trim ratios equal to \code{trim.ratio} }
#' \item\code{trimmed.count}{ the number of right trimmed individuals in each sequence of the test subset, if \code{trim}=TRUE }
#' }
#' 
#' @details This function splits the sample to train and test samples and 
#' trims the test sample from right, in order to provide a sample for examination 
#' of the prediction tools. 
#' In reliability applications, the hhsmm models are often left-to-right
#' and the modeling aims to predict the future states. In such cases, the 
#' test sets are right trimmed and the prediction aims to predict the 
#' residual useful lifetime (RUL) of a new sequence. 
#'
#' @examples
#'\donttest{
#' data(CMAPSS)
#' tt = train_test_split(CMAPSS$train, train.ratio = 0.7, trim = TRUE)
#'}
#' 
#' @export
#'
train_test_split <- function(train, train.ratio = 0.7, trim = FALSE, trim.ratio = NULL){
	if (class(train)!="hhsmmdata") stop("train must be of class hhsmmdata")
	if (train.ratio <= 0 | train.ratio > 1) stop("train.ratio must be in (0,1]")
	x = train$x
	N = train$N
	n = length(N)
	if (is.null(trim.ratio)){ 
		trim.ratio = runif(n) / 2 + 0.5
	} else {
		if (length(trim.ratio) != n & length(trim.ratio) != 1) {
			stop("length of train.ratio must be equal to the number of sequences or 1 !")
		} else if (length(trim.ratio) == 1) {
			trim.ratio = rep(trim.ratio, n)
		}
	}	
	if (!is.null(train$s)) s = train$s else s = NULL
	if (train.ratio < 1) {
		ntest = trunc((1 - train.ratio) * n)
		sam = sample(1:n, ntest)
		trim.ratio = trim.ratio[sam]
		Ns = cumsum(c(0, N))
		xtest = xtrain = rep(0, ncol(x))
		Ntest = Ntrain = c()
		if (!is.null(s)) stest = strain = c()
		for (i in 1:n) {
			if (i %in% sam) {
				xtest = rbind(xtest, as.matrix(x[(Ns[i] + 1):Ns[i + 1], ]))
				Ntest = c(Ntest, N[i])
				if (!is.null(s)) stest = c(stest, s[(Ns[i] + 1):Ns[i + 1]])
			} else {
				xtrain = rbind(xtrain, as.matrix(x[(Ns[i] + 1):Ns[i + 1], ]))
				Ntrain = c(Ntrain, N[i])
				if (!is.null(s)) strain = c(strain, s[(Ns[i] + 1):Ns[i + 1]])
			}
		}
		xtest = as.matrix(xtest[-1, ])
		xtrain = as.matrix(xtrain[-1, ])
	} else {
		xtest = xtrain = train$x
		if (!is.null(s)) strain = stest = train$s
		Ntrain = Ntest = train$N
		ntest = n 
	}
	if (trim) {
		if (!is.null(s)) strimmed = c()
		xtrimmed = rep(0, ncol(x))	
		Ntrim = trunc(Ntest * trim.ratio)
		Nts = cumsum(c(0, Ntest))
		for (i in 1:ntest) {
			xtrimmed = rbind(xtrimmed, as.matrix(xtest[(Nts[i] + 1):(Nts[i] + Ntrim[i]), ]))
			if (!is.null(s)) strimmed = c(strimmed, s[(Nts[i] + 1):(Nts[i] + Ntrim[i])])
		}
		xtrimmed = as.matrix(xtrimmed[-1, ])
	}
	if (!is.null(s)) {
		train = list(x = xtrain, N = Ntrain , s = strain)
		test = list(x = xtest, N = Ntest , s = stest)
		if (trim) trimmed = list(x = xtrimmed, N = Ntrim , s = strimmed)
			else trimmed = test 
	} else {
		train = list(x = xtrain, N = Ntrain)
		test = list(x = xtest, N = Ntest)
		if(trim) trimmed = list(x = xtrimmed, N = Ntrim)
			else trimmed = test 
	}
	class (train) <- class (test) <- class (trimmed) <- "hhsmmdata"
	list(train = train , test = test, trimmed = trimmed, trimmed.count = Ntest - Ntrim)
}
