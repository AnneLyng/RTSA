#' RTSA
#'
#' @param data A data set containing eI, eC, nI, nC or mI, mC, sdI, sdC, nI, nC
#' @param outcome Outcome of interest, RR only possibility now
#' @param mc minimum clinical relevant outcome
#' @param alpha The type I error
#' @param beta The type II error
#'
#' @return A TSA object (list for now)
#'
#' @export
#'
#' @examples
#' data(perioOxy)
#' RTSA(data = perioOxy, outcome = "RR", mc = 0.9)
RTSA <- function(data, outcome = "RR", mc, alpha = 0.05, beta = 0.2, ...){
  # calculate the meta-analysis
  mp = metaPrepare(outcome = outcome, eI = data$eI, eC = data$eC,
                         nI = data$nI, nC = data$nC)

  # Calculate the cumulative number of participants
  count <- cumsum(data$nI+data$nC)

  # Calculate the RIS
  if(outcome == "RR"){
  p0 = sum(data$eC+data$eI)/sum(data$nC+data$nI)
  pI = exp(log(p0)+log(mc))
  pC = exp(log(p0)-log(mc))
  RIS = nRandom(alpha = alpha, beta = beta, pI = pI, pC = pC, diversity = 0)
  }

  # Set the timings of the studies relative to the RIS
  timing <- c(count/RIS)
  if(max(timing) < 1){
    timing = c(timing,1)
  } else {
    warning("More information than needed. Setting timing to 1, where sample size is larger than RIS")
    timing[timing>1] = 1
  }

  RTSAout = TSA(timing = timing, synth = mp, anaTimes = 2:length(timing[timing <= 1]),
                alpha = alpha, ...)
  class(RTSAout) <- c("list", "RTSA")
  return(RTSAout)
}
