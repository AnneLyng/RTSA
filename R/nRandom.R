#' nRandom
#'
#' Calculate the RIS using diversity. Implementation of method from TSA software.
#'
#' @param alpha type 1 error
#' @param beta type 2 error
#' @param pI proportion of events in intervention group
#' @param pC proportion of events in control group
#' @param diversity estimate of diversity. Must be between 0 and 1.
#'
#' @return RIS for random-effects meta-analysis
#' @export
#'
#' @importFrom stats qnorm
#'
#' @examples
#' log.RR = log(0.9)
#' p0 = 0.1
#' pI = exp(log(p0)+log.RR/2)
#' pC = exp(log(p0)-log.RR/2)
#' nRandom(alpha = 0.05, beta = 0.2, pI = pI, pC = pC, diversity = 0.7)
#'
nRandom <- function(alpha, beta, pI, pC, diversity = NULL){

  p <- (pC + pI)/2
  v <- p*(1-p)
  theta <- pC-pI
  NF <- 2*(qnorm(1-alpha/2)+qnorm(1-beta))^2*2*v/theta^2

  NR <- 1/(1-diversity)*NF
  return(NR)
}
