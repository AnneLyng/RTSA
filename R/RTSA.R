#' RTSA
#'
#' @param data A data set containing eI, eC, nI, nC or mI, mC, sdI, sdC, nI, nC
#' @param outcome Outcome of interest, RR only possibility now
#' @param mc minimum clinical relevant outcome
#' @param alpha The type I error
#' @param beta The type II error
#' @param fixed Default is FALSE. Sets heterogeneity to 0 if TRUE.
#' @param ... other arguments
#'
#' @return A TSA object (list for now)
#'
#' @export
#'
#' @examples
#' data(perioOxy)
#' RTSA(data = perioOxy, outcome = "RR", mc = 0.9)
RTSA <- function(data, outcome = "RR", mc, alpha = 0.05, beta = 0.2, fixed = FALSE, ...){
  # calculate the meta-analysis
  mp = RTSA:::metaPrepare(outcome = outcome, eI = data$eI, eC = data$eC,
                         nI = data$nI, nC = data$nC)

  syn = RTSA:::synthesize(mp)

  # Calculate the cumulative number of participants
  subjects <- cumsum(data$nI+data$nC)

  logit <- function(x) log(x/(1-x))
  invlogit <- function(x) 1/(1+exp(-x))

  # Calculate the RIS
  if(outcome == "RR"){
  p0 = sum(data$eC+data$eI)/sum(data$nC+data$nI)
  pI = exp(log(p0)+log(mc)/2)
  pC = exp(log(p0)-log(mc)/2)
  RIS = nRandom(alpha = alpha, beta = beta, pI = pI, pC = pC,
                diversity = ifelse(fixed == TRUE | syn$U[1] == 0, 0, syn$U[4]))
  } else if(outcome == "OR"){
    p0 = sum(data$eC+data$eI)/sum(data$nC+data$nI)
    pI = invlogit(logit(p0)+log(mc)/2)
    pC = invlogit(logit(p0)-log(mc)/2)
    RIS = nRandom(alpha = alpha, beta = beta, pI = pI, pC = pC,
                  diversity = ifelse(fixed == TRUE | syn$U[1] == 0, 0, syn$U[4]))
  }

  # Set the timings of the studies relative to the RIS and the number of subjects
  timing <- c(subjects/RIS)
  if(max(timing) < 1){
    orgTiming = timing
    timing = c(timing,1)
  } else {
    warning("More information than needed. Setting timing to 1, first time number of subjects is larger than RIS. Boundaries will only be calculated for timings less or equal to 1. \n Note that the final results are based on the full number of subjects.")
    orgTiming = timing
    timing[timing>1] = 1
  }

  RTSAout = TSA(timing = timing, synth = mp, ana_time = 1:length(timing[timing <= 1]),
                alpha = alpha, mc = mc, ...)
  RTSAout$orgTiming = orgTiming
  RTSAout$RIS = RIS
  class(RTSAout) <- c("list", "RTSA")
  return(RTSAout)
}
