#'
#' @useDynLib hhsmm simulate_markov
#'
.simulate_markov <- function(init,transition,N) {
  if(!all.equal(rowSums(transition),rep(1,nrow(transition)))) stop("Rows of the transition matrix must sum to one")
  if(!all.equal(sum(init),1)) stop("initial probabilities must sum to one")
  a0 =  t(apply(transition,1,cumsum))
  st= cumsum(init)
  state = integer(sum(N))
  .C("simulate_markov",as.double(st),as.double(a0),as.integer(nrow(transition)),
		state=state,as.integer(N),as.integer(length(N)),PACKAGE="hhsmm")$state
}
