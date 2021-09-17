#' minTrial
#'
#' Calculate minimum number of trials for wanted power in a meta-analysis with
#' heterogeneity
#'
#' @param metric Metric of interest, options include only risk ratio (RR) for now
#' @param value Value of minimum relevant difference
#' @param tau2 Heterogeneity estimate
#' @param alpha Type 1 error
#' @param beta Type 2 error
#' @param p0 If verbose i set to TRUE and metric is set to RR (risk ration),
#' define common probability of event p0
#' @param verbose TRUE or FALSE. If the number of required participants per
#' trial is wanted as well, a matrix with the number of trials and the number
#' of required participants per trial is printed.
#'
#' @return Either a number (minimum required trials) or the minimum required
#' required trials together with a matrix of required participants per trial given
#' different number of trials.
#' @export
#'
#' @importFrom stats uniroot
#'
#' @examples
#' minTrial(metric = "RR", value = 0.7, tau2 = 0.05)
#'
minTrial = function(metric, value, tau2, alpha = 0.05, beta = 0.2, p0 = NULL, verbose = FALSE){

  if(metric == "RR"){
    log.RR = log(value)

    ntrial.func = function(log.RR, trial, alpha, beta, tau2){
      log.RR^2*trial/((qnorm(1-alpha/2)+qnorm(1-beta))^2)-tau2
    }

    minTrial = ceiling(uniroot(
      function(trial) ntrial.func(log.RR, trial, alpha, beta, tau2),
      interval = c(0, 1000)
    )$root)

  if(verbose == TRUE){
    log.p0 = log(p0)
    pI <- round(exp(log.p0+log.RR/2),4)
    pC <- round(exp(log.p0-log.RR/2),4)
    var.k = 1/(pC)+1/(pI)-2

    out.mat = matrix(NA, ncol = 4, nrow = 2)
    out.mat[1,] = c(minTrial, minTrial+1, minTrial+2, minTrial+3)
    out.mat[2,] = ceiling(2*var.k/(log.RR^2*c(minTrial, minTrial+1, minTrial+2, minTrial+3)/((qnorm(1-alpha/2)+qnorm(1-beta))^2)-tau2))
  }
  }

  if(verbose == TRUE){
    return(list(minTrial = minTrial, nPax = out.mat))
  } else {
    return(minTrial)
  }
}

