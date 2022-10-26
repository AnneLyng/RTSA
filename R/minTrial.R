#' minTrial
#'
#' Calculate minimum number of trials for wanted power in a meta-analysis with
#' heterogeneity
#'
#' @param outcome Metric of interest, options include only risk ratio (RR) for now
#' @param delta Value of minimum relevant difference
#' @param tau2 Heterogeneity estimate
#' @param alpha Type 1 error
#' @param beta Type 2 error
#' @param p0 If verbose i set to TRUE and outcome is set to RR (risk ration),
#' define common probability of event p0
#' @param var_delta Variance of the estimated effect
#' @param var_random Estimated variance from the random-effects meta-analysis
#'
#' @return Either a number (minimum required trials) or the minimum required
#' required trials together with a matrix of required participants per trial given
#' different number of trials.
#' @export
#'
#' @importFrom stats uniroot qnorm
#'
#' @examples
#' minTrial(outcome = "RR", delta = 0.7, tau2 = 0.05)
#'
minTrial = function(outcome,
                    delta,
                    tau2,
                    alpha = 0.05,
                    beta = 0.2,
                    p0 = NULL,
                    var_delta = NULL,
                    var_random = NULL) {
  if (outcome %in% c("RR", "OR")) {
    delta = log(delta)
  }

  if(!is.null(var_random)){
    ntrial.func = function(delta, trial, alpha, beta, tau2, var_random) {
      trial / ((qnorm(1 - alpha / 2) + qnorm(1 - beta)) ^ 2/delta^2-1/var_random) - tau2
    }

    minTrial <- ceiling(uniroot(
      function(trial)
        ntrial.func(delta, trial, alpha, beta, tau2, var_random),
      interval = c(0, 1000)
    )$root)
  } else {
    ntrial.func = function(delta, trial, alpha, beta, tau2) {
      trial / ((qnorm(1 - alpha / 2) + qnorm(1 - beta)) ^ 2/delta^2) - tau2
    }

    minTrial <- ceiling(uniroot(
      function(trial)
        ntrial.func(delta, trial, alpha, beta, tau2),
      interval = c(0, 1000)
    )$root)
  }

  if (outcome == "RR") {
    pI <- exp(log(p0) + delta / 2)
    pC <- exp(log(p0) - delta / 2)
    var_delta <- 1 / pC + 1 / pI - 2
  } else if (outcome == "OR") {
    logit <- function(x)
      log(x / (1 - x))
    invlogit <- function(x)
      1 / (1 + exp(-x))

    pI <- invlogit(logit(p0) + delta / 2)
    pC <- invlogit(logit(p0) - delta / 2)
    var_delta <- 1 / pI + 1 / pC + 1 / (1 - pI) + 1 / (1 - pC)
  }

  out.mat = matrix(NA, ncol = 4, nrow = 3)
  out.mat[1, ] <- c(minTrial, minTrial + 1, minTrial + 2, minTrial + 3)
  if(is.null(var_random)){
  out.mat[2, ] <- ceiling(2 * var_delta / (
    delta ^ 2 * c(minTrial, minTrial + 1, minTrial + 2, minTrial + 3) / ((qnorm(1 -
                                                                                   alpha / 2) + qnorm(1 - beta)) ^ 2) - tau2
  ))} else {
    out.mat[2, ] <- ceiling(2 * var_delta / (
      c(minTrial, minTrial + 1, minTrial + 2, minTrial + 3) / ((qnorm(1 -alpha / 2) + qnorm(1 - beta)) ^ 2/delta ^ 2 -1/var_random) - tau2
    ))
  }
  out.mat[3, ] <- out.mat[1, ] * out.mat[2, ]

  return(list(minTrial = minTrial, nPax = out.mat))
}

