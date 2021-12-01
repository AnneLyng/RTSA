#' TSA
#'
#' @param timing Vector containing values from 0 to 1 in chronological order.
#' @param anaTimes Analysis times.
#' @param synth metaPrepare object.
#' @param side 1 or 2 sided test. Only 2 sided available for now.
#' @param alpha Type 1 error rate.
#' @param stopTime Stopping time. Defaults to second last analysis time.
#' @param confInt If confidence intervals should be included.
#' @param count If no timing is provided, provide cumulative count of participants.
#' @param RIS If no timing is provided, provide RIS.
#' @param hakn If HKSJ adjustment should be used.
#'
#' @return A list consisting of:
#' \item{boundout}{Timing and boundaries}
#' \item{zout}{Meta-analysis synthesize object at either the time of stopping
#' (stopTime input) or the meta-analysis synthesize object containing all studies.}
#' \item{anaTimes}{Analysis times}
#' \item{stopTime}{Stopping time. Either the time of stopping
#' (stopTime input) or the analysis time which containing all studies.}
#' \item{naiveCI}{1-alpha% CI for fixed-effect and random-effects model at stopping
#' time}
#' \item{adjCI}{TSA adjusted confidence intervals. Also known as repeated confidence
#'  intervals}
#'
TSA = function(timing,
                anaTimes,
                synth,
                side = 2,
                alpha,
                stopTime = NULL,
                confInt = TRUE,
                count,
                RIS,
                hakn) {
  # start with preparing the timing data
  if (is.null(timing)) {
    timing = c(count / RIS, 1)
    timing = timing[timing < 1]
  }

  if(is.null(anaTimes)){
    anaTimes = 1:length(timing)
  }

  timing = timing[anaTimes]

  IFincrementThreshold <- 0.01
  IFtotalThreshold <- 0.05

  timingincr <- timing - c(0, timing[-length(timing)])
  trials <- cbind(timing, timingincr)

  anaTimes = anaTimes[which(trials[, 2] > IFincrementThreshold)]
  trials <- trials[trials[, 2] > IFincrementThreshold, ]

  trials[, 2] <- trials[, 1] - c(0, trials[, 1][-length(trials[, 1])])

  # calculate the boundaries
  boundout = boundary(informationFractions = trials[, 1],
                            side = side,
                            alpha = alpha)

  # calculate the cum. z-score (do we want this per study?)
  zout = lapply(anaTimes[anaTimes <= dim(synth$data)[1]],
                function(x) {
                  synout = synthesize(
                    metaPrepare(
                      data = synth$data[1:x,],
                      outcome = synth$outcome,
                      method = synth$method
                    )
                  )
                  return(synout)
                })

  names(zout) = anaTimes[anaTimes <= dim(synth$data)[1]]

  zvalues = sapply(names(zout), function(x){c(zout[[x]]$peF[4], zout[[x]]$peR[4])})

  if (confInt == TRUE) {
    if(is.null(stopTime)){ stopTime =
      as.character(max(anaTimes[anaTimes <= dim(synth$data)[1]]))}
    naiveCI = list(CIfixed = zout[[stopTime]]$peF[c(2, 3)],
                   CIrandom = zout[[stopTime]]$peR[c(2, 3)])
    adjCI = list(
      CIfixed = exp(
        log(zout[[stopTime]]$peF[1]) +
          c(-1, 1) * boundout$alpha.boundaries.upper[which(stopTime == anaTimes)] *
          sqrt(zout[[stopTime]]$peF[7])
      ),
      CIrandom = exp(
        log(zout[[stopTime]]$peR[1]) +
          c(-1, 1) * boundout$alpha.boundaries.upper[which(stopTime == anaTimes)] *
          sqrt(zout[[stopTime]]$peR[6])
      )
    )
  }

  RTSAout =     list(
    side = side,
    boundout = boundout,
    zout = zout[stopTime],
    zvalues = zvalues,
    anaTimes = anaTimes,
    stopTime = stopTime,
    naiveCI = naiveCI,
    adjCI = adjCI
  )

  class(RTSAout) <- c("list", "RTSA")

  return(
    RTSAout
  )

}

#' RTSA
#'
#' @param data A data set containing eI, eC, nI, nC or mI, mC, sdI, sdC, nI, nC
#' @param outcome Outcome of interst, RR only possibility now
#' @param mc minimum clinical relevant outcome
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
  RIS = RTSA::nRandom(alpha = alpha, beta = beta, pI = pI, pC = pC, diversity = 0)
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
