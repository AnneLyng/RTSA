#' RTSA
#'
#' @param data A data set containing eI, eC, nI, nC or mI, mC, sdI, sdC, nI, nC
#' @param side I or II sided hypothesis test
#' @param outcome Outcome of interest, RR only possibility now
#' @param mc minimum clinical relevant outcome
#' @param alpha The type I error
#' @param beta The type II error
#' @param futility Futility boundaries added to design. Default is "none".
#' @param fixed Default is FALSE. Sets heterogeneity to 0 if TRUE.
#' @param method Weighting method. Defaults to Mantel-Haenzel (MH).
#' @param vartype TBD.
#' @param sign TBD.
#' @param fixedStudy TBD.
#' @param hksj TBD.
#' @param tau.ci.method TBD.
#' @param ... other arguments
#'
#' @return A TSA object (list for now)
#'
#' @export
#'
#' @examples
#' data(perioOxy)
#' RTSA(data = perioOxy, outcome = "RR", mc = 0.9, side = 2)
RTSA <- function(data, side, outcome = "RR", mc, alpha = 0.05, beta = 0.2,
                 futility = "none", fixed = FALSE, anaTimes = NULL,
                 method = "MH", vartype = "equal", sign = NULL,
                 fixedStudy = FALSE,
                 hksj = FALSE,
                 tau.ci.method = "BJ",...){
  # calculate the meta-analysis
  mp = metaPrepare(outcome = outcome, data = data, method = method,
                   vartype = vartype, alpha = alpha)

  syn = synthesize(mp, sign = sign, hksj = hksj, tau.ci.method = tau.ci.method,
                   fixedStudy = fixedStudy)

  # Calculate the cumulative number of participants
  subjects <- cumsum(data$nI+data$nC)

  logit <- function(x) log(x/(1-x))
  invlogit <- function(x) 1/(1+exp(-x))

  # Calculate the RIS
  if(outcome == "RR"){
  p0 = sum(data$eC+data$eI)/sum(data$nC+data$nI)
  pI = exp(log(p0)+log(mc)/2)
  pC = exp(log(p0)-log(mc)/2)
  } else if(outcome == "OR"){
    p0 = sum(data$eC+data$eI)/sum(data$nC+data$nI)
    pI = invlogit(logit(p0)+log(mc)/2)
    pC = invlogit(logit(p0)-log(mc)/2)
  }

  if(side == 1){
    RIS = nRandom(alpha = alpha*2, beta = beta, pI = pI, pC = pC,
                diversity = ifelse(fixed == TRUE | syn$U[1] == 0, 0, syn$U[4]))
  } else {
    RIS = RTSA:::nRandom(alpha = alpha, beta = beta, pI = pI, pC = pC,
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

  if(is.null(anaTimes)){
    anaTimes <- 1:length(timing[timing <= 1])
  }

  RTSAout = TSA(timing = timing, side = side, synth = mp, ana_time = anaTimes,
                alpha = alpha, beta = beta, futility = futility, mc = mc, sign = sign,
                fixedStudy = fixedStudy,
                hksj = hksj,
                tau.ci.method = tau.ci.method)
  RTSAout$Pax <- list(RIS = RIS, p0 = p0, pI = pI, pC = pC)
  RTSAout$orgTiming = orgTiming
  RTSAout$adjRIS = RIS*RTSAout$root
  RTSAout$root = RTSAout$root
  RTSAout$AIS = sum(data$nC+data$nI)
  class(RTSAout) <- c("list", "RTSA")
  return(RTSAout)
}
