#' nFixed
#'
#' Calculate the number of total participants when there is no presence
#'  of heterogeneity to achieve a specific power.
#'
#' @param alpha type 1 error
#' @param beta type 2 error
#' @param pI probability of event in intervention group (for binary outcomes)
#' @param pC probability of event in control group (for binary outcomes)
#' @param iE intervention effect (for continuous outcomes)
#' @param sdE standard deviation of intervention effect (for continuous outcomes)
#' @param binary TRUE/FALSE. If continuous outcome set to FALSE. Default is TRUE.
#'
#' @return The total required number of  participants to achieve a specific power
#' in a fixed-effect model
#' @export
#'
#' @importFrom stats qnorm
#'
#' @examples
#' log.RR = log(0.9)
#' p0 = 0.1
#' pI = exp(log(p0)+log.RR/2)
#' pC = exp(log(p0)-log.RR/2)
#' nFixed(alpha = 0.05, beta = 0.2, pI = pI, pC = pC, binary = TRUE)
#'
nFixed <- function(alpha, beta, pI = NULL, pC = NULL, iE = NULL,
                   sdE = NULL, binary = TRUE){
  if(binary == TRUE){
    p <- (pC + pI)/2
    v <- p*(1-p)
    theta <- pC-pI
    return(2*(qnorm(1-alpha/2)+qnorm(1-beta))^2*2*v/theta^2)
  } else {
    return(2*(qnorm(1-alpha/2)+qnorm(1-beta))^2*2*sdE^2/iE^2)
  }
}
