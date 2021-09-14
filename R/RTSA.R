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
#' @export
#'
#' @examples
#' count <- cumsum(perioOxy$nI+perioOxy$nC)
#' timing <- c(count/5673,1)
#' mp = metaPrepare(outcome = "RR", eI = perioOxy$eI, nI = perioOxy$nI,
#' eC = perioOxy$eC, nC = perioOxy$nC, method = "IV")
#' TSA(timing = timing, synth = mp, anaTimes = c(4,5,7,8), side = 2, alpha = 0.05)
#'
TSA = function(timing,
                anaTimes,
                synth,
                side,
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
  boundout = RTSA::boundary(informationFractions = trials[, 1],
                            side = side,
                            alpha = alpha)

  # calculate the cum. z-score (do we want this per study?)
  zout = lapply(anaTimes[anaTimes <= length(synth$eI)],
                function(x) {
                  synout = synthesize(
                    metaPrepare(
                      outcome = synth$outcome,
                      eI = synth$eI[1:x],
                      eC = synth$eC[1:x],
                      nI = synth$nI[1:x],
                      nC = synth$nC[1:x],
                      method = synth$method
                    )
                  )
                  return(synout)
                })

  names(zout) = anaTimes[anaTimes <= length(synth$eI)]

  zvalues = sapply(names(zout), function(x){c(zout[[x]]$peF[4], zout[[x]]$peR[4])})

  if (is.null(stopTime) & confInt == TRUE) {
    stopTime = as.character(max(anaTimes[anaTimes < length(synth$eI)]))
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
#' @param data A data set containing eI, eC, nI, nC
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
RTSA <- function(data, outcome = "RR", mc){
  # calculate the meta-analysis
  mp = RTSA::metaPrepare(outcome = outcome, eI = data$eI, eC = data$eC,
                         nI = data$nI, nC = data$nC)

  # Calculate the cumulative number of participants
  count <- cumsum(data$nI+data$nC)

  # Calculate the RIS
  if(outcome == "RR"){
  p0 = sum(data$eC+data$eI)/sum(data$nC+data$nI)
  pI = exp(log(p0)+log(mc))
  pC = exp(log(p0)-log(mc))
  RIS = RTSA::nRandom(alpha = 0.05, beta = 0.2, pI = pI, pC = pC, diversity = 0)
  }

  # Set the timings of the studies relative to the RIS
  timing <- c(count/RIS)
  if(max(timing) < 1){
    timing = c(timing,1)
  }

  RTSA::TSA(timing = timing, synth = mp, anaTimes = 2:length(timing[timing < 1]),
      side = 2, alpha = 0.05)
}
