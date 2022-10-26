#' ris: Calculate required information size and required number of trials.
#'
#' RTSA's sample size calculation.
#'
#' @param outcome Choose between: "cont", "RR", "OR" or "RD".
#' @param delta Minimum clinical relevant effect. For "OR" or "RR" set to natural scale, not log scale.
#' @param side Test type. Set to 1 or 2 depending on the test being 1- or 2-sided.
#' @param alpha Level of type-I-error. Default value is 0.05.
#' @param beta Level of type-II-error. Default value is 0.2.
#' @param random Should sample size be based on a random effects model? Defaults to FALSE.
#' @param sd_delta Standard deviation of estimated effect. Only needed when outcome = "cont".
#' @param p0 Probability of event in population.
#' @param data Data set containing studies used for meta-analysis.
#' @param tau2 The value of the heterogeneity. Use when estimating the sample size under a random effects model. If data is provided, the estimated heterogeneity is used instead.
#' @param I2 Inconsistency.
#' @param D2 Diversity.
#' @param var_random Estimated variance from the random-effects meta-analysis
#' @param type Prospective or retrospective sample size calculation.
#' @param ... additional arguments
#'
#' @return A list containing required information size and more.
#' @export
#' @aliases print.ris
#'
#' @examples
#' data("perioOxy")
#' ris(outcome = "RR", delta = 0.7, p0 = 0.1, data = perioOxy, random = FALSE)
ris <-
  function(outcome,
           delta,
           side = 2,
           alpha = 0.05,
           beta = 0.2,
           random = FALSE,
           sd_delta = NULL,
           p0 = NULL,
           data = NULL,
           tau2 = NULL,
           I2 = NULL,
           D2 = NULL,
           var_random = NULL,
           type = "prospective",
           ...) {
    if (outcome == "cont" & is.null(sd_delta)) {
      stop("For continuous outcomes provide delta and sd_delta")
    }

    if (outcome %in% c("RR", "OR") & is.null(p0)) {
      stop("For binary outcomes provide p0.")
    }

    if(type == "retrospective" & (is.null(data)|is.null(var_random))){
      stop("Data and specification of the variance of the pooled effect size (random-effects) is required for retrospective sample size calculation")
    }

    if (outcome == "RR" | outcome == "OR") {
      # prepare helper functions
      logit <- function(x)
        log(x / (1 - x))
      invlogit <- function(x)
        1 / (1 + exp(-x))

      # calculate common prob of event if NULL
      if (is.null(p0)) {
        warning(
          "Common probability of event (p0) is calculated from data. This affects the required number of participants calculation (power calculation). See vignette 'Calculating required sample size and required number of trials' for more information."
        )
        p0 <- sum(data$eC + data$eI) / sum(data$nC + data$nI)
      }
      if (outcome == "RR") {
        pI <- exp(log(p0) + log(delta) / 2)
        pC <- exp(log(p0) - log(delta) / 2)
      } else if (outcome == "OR") {
        pI = invlogit(logit(p0) + log(delta) / 2)
        pC = invlogit(logit(p0) - log(delta) / 2)
      }
      p <- (pC + pI) / 2
      var_delta <- p * (1 - p)
      delta_nf <- pC - pI
    }

    if(outcome == "cont"){
      var_delta <- sd_delta^2
      delta_nf <- delta
    }

    NF <-
      2 * (qnorm(1 - alpha / side) + qnorm(1 - beta)) ^ 2 * 2 * var_delta / delta_nf ^
      2

    NF <- ceiling(NF) + ceiling(NF) %% 2

    args <- mget(names(formals()),sys.frame(sys.nframe()))

    outlist <-
      list(
        var_delta = var_delta,
        NF = NF,
        args = args
      )

    if (random == TRUE & !is.null(data)) {
      ma <-
        RTSA::metaanalysis(
          data = data,
          outcome = outcome,
          alpha = alpha,
          beta = beta,
          side = side
        )

      if(outcome %in% c("RR", "OR")){
        NR <- minTrial(outcome = outcome, delta = delta, alpha = alpha,
                       beta = beta, p0 = p0, tau2 = ma$synthesize$U[1],
                       var_random = var_random)
        NR_bc <- minTrial(outcome = outcome, delta = delta, alpha = alpha,
                          beta = beta, p0 = p0,
                          tau2 = ma$synthesize$ci.tau$random[1,2],
                          var_random = var_random)
        NR_bc <- append(NR_bc, list(tau2 = ma$synthesize$ci.tau$random[1,2]))
        NR_wc <- minTrial(outcome = outcome, delta = delta, alpha = alpha,
                          beta = beta, p0 = p0,
                          tau2 = ma$synthesize$ci.tau$random[1,3],
                          var_random = var_random)
        NR_wc <- append(NR_wc, list(tau2 = ma$synthesize$ci.tau$random[1,3]))
      } else{
        NR <- minTrial(outcome = outcome, delta = delta, alpha = alpha,
                       beta = beta, var_delta = var_delta, tau2 = ma$synthesize$U[1],
                       var_random = var_random)
        NR_bc <- minTrial(outcome = outcome, delta = delta, alpha = alpha,
                          beta = beta, var_delta = var_delta,
                          tau2 = ma$synthesize$ci.tau$random[1,2],
                          var_random = var_random)
        NR_wc <- try(minTrial(outcome = outcome, delta = delta, alpha = alpha,
                              beta = beta, var_delta = var_delta,
                              tau2 = ma$synthesize$ci.tau$random[1,3],
                              var_random = var_random), TRUE)
        if(class(NR_wc) == "try-error"){
          NR_wc = "Too much heterogeneity to estimate an upper number of trials."
        }
      }


      NR_div <- 1 / (1 - ma$synthesize$U[4]) * NF
      NR_inc <- 1 / (1 - ma$synthesize$U[3]) * NF

      NR_div <- ceiling(NR_div) + ceiling(NR_div) %% 2
      NR_inc <- ceiling(NR_inc) + ceiling(NR_inc) %% 2

      outlist <-
        append(outlist, list(
          tau2 = ma$synthesize$U[1],
          NR_div = NR_div,
          NR_inc = NR_inc,
          NR = NR,
          NR_bc = NR_bc,
          NR_wc = NR_wc
        ))

    } else if(random == TRUE & is.null(data)){

      if(is.null(tau2) & is.null(I2) & is.null(D2)){
        stop("No value for heterogeneity (tau2, I2 and/or D2) is provided")
      }

      if(is.null(tau2)) tau2 <- 0
      NR_div = NULL
      NR_inc = NULL

      if(outcome %in% c("RR", "OR")){
        NR <- minTrial(outcome = outcome, delta = delta, alpha = alpha,
                       beta = beta, p0 = p0, tau2 = tau2, var_random = var_random)
      } else{
        NR <- minTrial(outcome = outcome, delta = delta, alpha = alpha,
                       beta = beta, var_delta = var_delta, tau2 = tau2, var_random = var_random)
      }
      if(!is.null(I2)){
        NR_inc <- 1 / (1 - I2) * NF
        NR_inc <- ceiling(NR_inc) + ceiling(NR_inc) %% 2
      }
      if(!is.null(D2)){
        NR_div <- 1 / (1 - D2) * NF
        NR_div <- ceiling(NR_div) + ceiling(NR_div) %% 2
      }

      outlist <-
        append(outlist, list(
          tau2 = tau2,
          NR = NR,
          NR_div = NR_div,
          NR_inc = NR_inc
        ))
    }

    class(outlist) <- "ris"

    return(outlist)
  }


# FUNCTION | print ris ----
#' @method print ris
#' @export
print.ris <- function(x, ...) {
  if (x$args$type == "prospective") {
    cat(
      "The sample size calculation assumes a", paste0(x$args$side, "-sided test,"), "equal group sizes, a type-I-error of",
      x$args$alpha,
      "and a type-II-error of",
      x$args$beta,
      "\n\n"
    )
    cat("Fixed-effect required information size:\n")
    cat(paste(x$NF, "participants in total. \n"))
    if (x$args$random == TRUE) {
      cat("\n")
      cat("Random-effects required information size:\n")
      if (x$tau2 != 0) {
        cat(
          paste(
            x$NR$nPax[3, 4],
            "participants in total split over",
            x$NR$nPax[1, 4],
            "trials. Adjusted by tau^2.\n"
          )
        )
      }
      if (!is.null(x$NR_div)) {
        cat(paste(
          x$NR_div,
          "participants in total. Adjusted by diversity (D^2).\n"
        ))
      }
      if (!is.null(x$NR_inc)) {
        cat(paste(
          x$NR_inc,
          "participants in total. Adjusted by inconsistency (I^2).\n"
        ))
      }
    }
    cat("\n")
    cat(
      "For more information about the sample size calculation see vignette 'Calculating required sample size and required number of trials'."
    )
  } else {
    cat(
      "This is a retrospective meta-analysis. \nThe sample size calculation assumes equal group sizes, a type-I-error of",
      x$args$alpha,
      "and a type-II-error of",
      x$args$beta,
      "\n\n"
    )
    cat("Fixed-effect required information size:\n")
    cat(paste(x$NF, "participants in total are additionally required. \n"))
    if (x$args$random == TRUE) {
      cat("\n")
      cat("Random-effects required information size:\n")
      if (x$tau2 != 0) {
        cat(
          paste(
            x$NR$nPax[3, 4],
            "participants in total are additionally required. These can be split over",
            x$NR$nPax[1, 4],
            "trials. Adjusted by tau^2.\n"
          )
        )
      }
      if (!is.null(x$NR_div)) {
        cat(paste(
          x$NR_div,
          "participants in total are additionally required. Adjusted by diversity (D^2).\n"
        ))
      }
      if (!is.null(x$NR_inc)) {
        cat(paste(
          x$NR_inc,
          "participants in total are additionally required. Adjusted by inconsistency (I^2).\n"
        ))
      }
    }
    cat("\n")
    cat(
      "For more information about the sample size calculation see vignette 'Calculating required sample size and required number of trials'."
    )
  }
}
