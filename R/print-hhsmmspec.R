#' @export
#'
print.hhsmmspec <- function(x, ...){
  cat("Hidden Hybrid markov-semi-Markov Model specification:\n")
  cat("Markov states: ",which(!x$semi),"\n")
  cat(sprintf("J (number of states): \n%i \n", x$J))
  cat("init:\n")
  print(x$init)
  cat ("transition matrix:\n")
  print(x$transition)
  cat("emission distribution:\n")
  print(x$parms.emission)
  if(!all(!x$semi)){
	cat("sojourn distribution:\n")
  	print(x$sojourn)
  }
  return(invisible(x))
}
