#' metaPreapre
#'
#' Preparation for synthesization of results for binary data. Calculates study specific effect sizes, standard errors, confidence intervals and weights to use for meta-analysis. Handling of continuous outcomes will be added soon.
#'
#' @param outcome binary outcome of interest (option include: RD, OR, RR)
#' @param eI number of events in intervention group
#' @param nI number of pax in intervention group
#' @param eC number of events in control group
#' @param nC number of pax in control group
#' @param mI effect estimate in intervention group (continuous outcome)
#' @param mC effect estimate in control group (continuous outcome)
#' @param sdI standard error of estimate in intervention group (continuous outcome)
#' @param sdC standard error of estimate in control group (continuous outcome)
#' @param vartype Variance type for continuous outcomes. Default is "equal", other choice is "non-equal"
#' @param method Method for calculating weights. Options include: MH (Mantel-Haenzel), Inverse variance weighting (IV) or GLM
#'
#' @return A list is returned with the following items:
#' \item{w}{weights}
#' \item{te}{trial specific effect size}
#' \item{lower}{lower CI for trial specific effect size}
#' \item{upper}{upper CI for trial specific effect size}
#' \item{pe}{pooled effect size (for MH method svpe is also included which is the sum of variance)}
#' \item{sig}{sig: sigma (standard deviation of effect size per trial - consider renaming)}
#' \item{outcome}{statistic used for effect size (OR, RR og RD)}
#' \item{method}{method for calculating weights. GLM, Inverse variance weighting (IV) or Mantel-Haenzel (MH)}
#' \item{eI}{number of events in the intervention group}
#' \item{eC}{number of events in the control group}
#' \item{nI}{number of pax in the intervention group}
#' \item{nC}{number of pax in the control group}
#'
#' @export
#'
#' @examples
#' data(perioOxy)
#' with(perioOxy, metaPrepare(outcome = "RR", eI = eI, nI = nI, eC = eC, nC = nC,
#' method = "MH"))
metaPrepare <- function(outcome = "RR", eI = NULL, nI = NULL, eC = NULL, nC = NULL,
                         mI = NULL, mC = NULL, sdI = NULL, sdC = NULL, vartype = "equal",
                         method = "MH") {
  if(outcome %in% c("OR", "RR")){

    # if there is total zero events
    if(sum(eI == 0 & eC == 0) > 0){
      stop("Zero total events in data.\nOdds Ratio or Relative Risk cannot be computed.
 Please remove zero total events trials from data or change to Risk Difference")
    }

    # if one of the event counts is zero
    if(sum(eI == 0 | eC == 0) >0){
      zc <- which(eI == 0 | eC == 0)
      eI[zc] <- eI[zc] + 0.5
      nI[zc] <- nI[zc] + 1
      eC[zc] <- eC[zc] + 0.5
      nC[zc] <- nC[zc] + 1
    }

    # calculate event probability
    pI <- eI/nI
    pC <- eC/nC

    # set the number of trails
    K <- length(eI)

    # prepare returned object
    out <- list()
    class(out) <- "synthPrepped"

    # calculate effects sizes and beloning std.
    if(outcome == "RD") { # Risk difference
      te <- pI-pC
      sig <- sqrt(pI*(1-pI)/nI+pC*(1-pC)/nC)
    } else if(outcome == "OR") { # Odds Ratio
      te <- (pI/(1-pI))/(pC/(1-pC))
      sig <- sqrt(1/eI+1/(nI-eI)+1/eC+1/(nC-eC))
    } else if(outcome == "RR") { # Risk Ratio
      te <- pI/pC
      sig <- sqrt(1/eI-1/nI+1/eC-1/nC)
    }

    # calculate the weights
    if(method == "IV") { # inverse variance
      w <- 1/(sig^2)
      w <- w/sum(w)
      pe <- sum(te*w)

      w <- w*100
    } else if(method == "MH" | method == "GLM"){ # Mantel-Haenszel or GLM
      A <- eI # make the complicated variance (emerson)
      B <- nI - eI
      C <- eC
      D <- nC - eC
      N <- A + B + C + D

      if(outcome == "OR"){
        w <- (nI-eI)*eC/(nC+nI)

        T1 <- (A+D)/N
        T2 <- (B+C)/N
        T3 <- A*D/N
        T4 <- B*C/N

        vpe <- 0.5*((T1*T3)/(sum(T3)^2)+(T1*T4+T2*T3)/(sum(T3)*sum(T4))+T2*T4/(sum(T4)^2))
        svpe <- sum(vpe)

      } else if(outcome == "RR"){
        w <- (nI)*eC/(nI+nC)

        D1 <- ((A+B)*(C+D)*(A+C)-(A*C*N))/(N^2)
        R <- (A*(C+D))/(N)
        S <- (C*(A+B))/N

        svpe <- sum(D1)/(sum(R)*sum(S))

      } else if(outcome == "RD"){
        w <- (nI)*(nC)/(nI+nC)
      }
      pe <- sum(te*w)/sum(w)
    }

    lower <- exp(log(te) - qnorm(1-0.05/2)*sig)
    upper <- exp(log(te) + qnorm(1-0.05/2)*sig)

    # return results
    if(method == "MH"){
      out <- list(w = w, te = te, lower = lower, upper = upper, pe = c(pe, svpe),
                  sig = sig, outcome = outcome, method = method, eI = eI,
                  eC = eC, nC = nC, nI = nI)
    } else {
      out <- list(w = w, te = te, lower = lower, upper = upper, pe = pe,
                  sig = sig, outcome = outcome, method = method, eI = eI,
                  eC = eC, nC = nC, nI = nI)
    }
  } else if(outcome == "cont"){
    if(is.null(vartype)){
      vartype = "equal"
      message("No choice about type of variance was indicated. Equal variance is chosen.")
    }

    te = mI - mC

    if(vartype == "equal"){
      spooled <- sqrt(((nI-1)*sdI^2+(nC-1)*sdC^2)/(nI+nC-2))
      vte <- (nI+nC)/(nI*nC)*spooled^2
      sete <- sqrt(vte)
      df <- (nI-1)+(nC-1)
      lower <- te - qt(0.975, df = df)*sete
      upper <- te + qt(0.975, df = df)*sete
    } else {
      vte <- sdI^2/nI + sdC^2/nC
      sete <- sqrt(vte)
      #df2 <- (sdI^2/nI+sdC^2/nC)^2/((sdI^4)/(nI^2*(nI-1))+(sdC^4)/(nC^2*(nC-1)))
      lower <- te - qnorm(0.975)*sete
      upper <- te + qnorm(0.975)*sete
    }

    w <- 1/vte
    pe <- sum(w*te)/sum(w)
    sig <- sqrt(vte)
    w <- w/sum(w)

    out <- list(w = w, te = te, lower = lower, upper = upper, pe = pe,
                sig = sig, outcome = outcome, method = method, mI = mI,
                mC = mC, nC = nC, nI = nI)
  }


  class(out) <- "synthPrepped"
  return(out)
}
