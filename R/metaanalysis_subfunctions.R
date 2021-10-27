# metaPrepare - func ----
#' @importFrom stats glm coef vcov pnorm qnorm pchisq pt qt
metaPrepare <- function(data = NULL,
                        eI = NULL, nI = NULL, eC = NULL, nC = NULL,
                        mI = NULL, mC = NULL, sdI = NULL, sdC = NULL,
                        outcome = "RR",method = "MH",vartype = "unequal") {
  #Maybe add a test if no data is inputted.

  #Import from 'data' depending on the chosen outcome.
  #NB. If data i chosen this will overwrite any single defined variables.
  if (!is.null(data) & outcome %in% c("OR", "RR", "RD")) {
    eI = data$eI
    nI = data$nI
    eC = data$eC
    nC = data$nC
  }else if (!is.null(data) & outcome == "cont") {
    mI = data$mI
    mC = data$mC
    sdI = data$sdI
    sdC = data$sdC
    nI = data$nI
    nC = data$nC
  }

  if(is.null(data) & outcome %in% c("OR", "RR", "RD")){
    data = data.frame(eI, nI, eC, nC)
  }else if(is.null(data) & outcome == "cont"){
    data = data.frame(mI, mC, sdI, sdC, nI, nC)
  }


  if(outcome == "cont"){
    method = "IV"
  }

  #Prepare dichotomous outcomes.
  if(outcome %in% c("OR", "RR", "RD")){

    # Stop if any trial has zero total events.
    if(sum(eI == 0 & eC == 0) > 0 & outcome %in% c("RR", "OR")){
      stop("One or more trials with zero total events.
      Odds Ratio (OR) or Relative Risk (RR) cannot be computed.
      Please remove zero total events trials from data or change to Risk Difference (RD)")
    }

    # Adding 0.5 if one of the event counts is zero
    if(outcome %in% c("OR", "RR"))
    if(sum(eI == 0 | eC == 0) > 0){
      zc <- which(eI == 0 | eC == 0)
      eI[zc] <- eI[zc] + 0.5
      nI[zc] <- nI[zc] + 1
      eC[zc] <- eC[zc] + 0.5
      nC[zc] <- nC[zc] + 1
    }

    # Calculate event probability
    pI <- eI/nI
    pC <- eC/nC

    # Set the number of trials
    K <- length(eI)

    # Prepare returned object
    out <- list()
    class(out) <- "synthPrepped"

    # Calculate effects sizes and beloning std.
    if(outcome == "RD") { # Risk difference
      te <- pI-pC
      sig <- sqrt(pI*(1-pI)/nI+pC*(1-pC)/nC)
    }else if(outcome == "OR") { # Odds Ratio
      te <- (pI/(1-pI))/(pC/(1-pC))
      sig <- sqrt(1/eI+1/(nI-eI)+1/eC+1/(nC-eC))
    }else if(outcome == "RR") { # Risk Ratio
      te <- pI/pC
      sig <- sqrt(1/eI-1/nI+1/eC-1/nC)
    }

    # Calculate the weights
    if(method == "IV") { # Inverse Variance
      w <- 1/(sig^2)
      w <- w/sum(w)
      pe <- sum(te*w)
      w <- w*100
    }else if(method == "MH" | method == "GLM"){ # Mantel-Haenszel or GLM
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
      }else if(outcome == "RR"){
        w <- (nI)*eC/(nI+nC)
        D1 <- ((A+B)*(C+D)*(A+C)-(A*C*N))/(N^2)
        R <- (A*(C+D))/(N)
        S <- (C*(A+B))/N
        svpe <- sum(D1)/(sum(R)*sum(S))
      }else if(outcome == "RD"){
        w <- (nI)*(nC)/(nI+nC)
      }

      pe <- sum(te*w)/sum(w)
    }

    lower <- exp(log(te) - qnorm(1-0.05/2)*sig)
    upper <- exp(log(te) + qnorm(1-0.05/2)*sig)

    # return results
    if(method == "MH"){
      out <- list(w = w, te = te, lower = lower, upper = upper, pe = c(pe, svpe),
                  sig = sig, outcome = outcome, method = method, data = data)
    }else{
      out <- list(w = w, te = te, lower = lower, upper = upper, pe = pe,
                  sig = sig, outcome = outcome, method = method, data = data)
    }

  }else if(outcome == "cont"){
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
      lower <- te - qnorm(0.975)*sete
      upper <- te + qnorm(0.975)*sete
    }

    w <- 1/vte
    pe <- sum(w*te)/sum(w)
    sig <- sqrt(vte)
    w <- w/sum(w)

    out <- list(w = w, te = te, lower = lower, upper = upper, pe = pe,
                sig = sig, outcome = outcome, method = method, data = data)
  }

  class(out) <- "synthPrepped"
  return(out)
}

# synthesize - func ----
#' @importFrom stats qnorm qt binomial
synthesize <- function(y,
                       sign = NULL,
                       fixedStudy = TRUE,
                       hksj = FALSE) {

  # Denne test kan vel fjernes når brugeren bare kører metaanalysis?
  # if (class(y) != "synthPrepped") {
  #   warning('sig object is not of synthPrepped-class, results may be inaccurate')
  # }

  w <- y$w   # collect objects
  sig <- y$sig
  te <- y$te
  pe <- y$pe
  eI <- y$eI
  eC <- y$eC
  nI <- y$nI
  nC <- y$nC
  df <- length(w) - 1
  if (length(w) == 1) df <- 1 #check up on this

  #For GLM
  if (y$method == "GLM") {
    bi <- nI - eI
    di <- nC - eC

    grpOut <- cbind(xi = c(rbind(eI, eC)), mi = c(rbind(bi, di)))
    k <- length(eI)
    eff <- rep(rep(1, k), 2)
    trial <- factor(rep(seq_len(k), each = 2))
    group <- rep(c(1, 0), times = k)
    eff <- eff * group

    if (fixedStudy == FALSE) {
      const <- rep(rep(1, k), 2)
      glmFit <- lme4::glmer(grpOut ~ -1 + eff + const + (1 | trial),
                            nAGQ = 7,
                            family = "binomial")
      es <- glmFit@beta[1]
      sigma2 <- lme4::VarCorr(glmFit)[[1]][1]
      tau2 <- 0
    } else {
      glmFit <- glm(grpOut ~ -1 + eff + trial, family = binomial)
      es <- coef(glmFit)[1]
    }

    se <- sqrt(vcov(glmFit)[1, 1])
    zval <- es / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)

    lci <- exp(es - qnorm(0.05 / 2, lower.tail = FALSE) * se)
    uci <- exp(es + qnorm(0.05 / 2, lower.tail = FALSE) * se)
    fpe <- exp(es)

    sv <- rep(c(.5, -.5), times = k)
    if (fixedStudy == FALSE) {
      glmRandomFit <- lme4::glmer(grpOut ~ -1 + eff + const + (1 | trial) +
                                    (sv - 1 | trial), family = binomial)
    } else {
      glmRandomFit <-
        lme4::glmer(grpOut ~ -1 + eff + trial + (sv - 1 | trial),
                    family = binomial, nAGQ = 7)
    }

    esR <- glmRandomFit@beta[1]
    seR <- sqrt(vcov(glmRandomFit)[1, 1])
    if (fixedStudy == FALSE) {
      tau2 <- lme4::VarCorr(glmRandomFit)[[2]][1]
      sigma2 <- lme4::VarCorr(glmRandomFit)[[1]][1]
    } else {
      tau2 <- lme4::VarCorr(glmRandomFit)[[1]][1]
    }

    zvalR <- esR / seR
    pvalR <- 2 * pnorm(abs(zvalR), lower.tail = FALSE)
    lciR <- exp(esR - qnorm(0.05 / 2, lower.tail = FALSE) * seR)
    uciR <- exp(esR + qnorm(0.05 / 2, lower.tail = FALSE) * seR)
    peR <- exp(esR)

    return(list(
      peF = c(fpe, lci, uci, zval, pval, se, es),
      peR = c(peR, lciR, uciR, zvalR, pvalR),
      U = tau2
    ))

  } else {
    # NOT GLM

    if (y$method != "MH") {
      w <- 1 / (sig ^ 2) # weight inverse variance
      rw <- w / sum(w) # relative weight
      vw <- 1 / sum(w) # variance of pooled effect
    }

    if (y$outcome == "cont") {
      peF <- sum(w * te) / sum(w)
      lci <- peF - 1.96 * sqrt(vw)
      uci <- peF + 1.96 * sqrt(vw)
      if (is.null(sign)) {
        zval <- peF / sqrt(vw)
      } else {
        zval <- sign * peF / sqrt(vw)
      }
      pval <- (1 - pnorm(abs(zval))) * 2

    }else if (y$outcome != "cont" & y$method == "MH") {
      # method the same for OR and RR
      vw <- pe[2]
      lpeF <- log(sum(te * w) / sum(w))
      lci <- exp(lpeF - 1.96 * sqrt(vw))
      uci <- exp(lpeF + 1.96 * sqrt(vw))
      peF <- exp(lpeF)
      if (is.null(sign)) {
        zval <- lpeF / sqrt(vw)
      } else {
        zval <- sign * lpeF / sqrt(vw)
      }

      pval <- (1 - pnorm(abs(zval))) * 2

    }else{
      lpeF <- sum(log(te) * rw) # fixed effect log pooled estimate
      uci <- exp(lpeF + 1.96 * sqrt(vw))
      lci <- exp(lpeF - 1.96 * sqrt(vw))
      peF <- exp(lpeF)
      if(is.null(sign)){
        zval <- lpeF / sqrt(vw)
      }else{
        zval <- sign * lpeF / sqrt(vw)
      }
      pval <- (1 - pnorm(abs(zval))) * 2

    }

    if(y$outcome != "cont"){
      w <- 1 / (sig ^ 2)
      Q <- sum(w * log(te) ^ 2) - (sum(w * log(te))) ^ 2 / sum(w)
      U <- sum(w) - sum(w ^ 2) / sum(w)
      tau2 <-
        ifelse(Q > df, (Q - df) / U, 0) # DerSimonian-Laird estimate
      pQ <- pchisq(Q, df, lower.tail = FALSE)
      if(!is.na(Q) & Q / df <= 0){
        H <- 0
      }else if(is.na(Q)){
        H <- NA
      }else{
        H <- sqrt(Q / df)
      }
      I2 <- ifelse((Q - df) / Q >= 0, (Q - df) / Q, 0)

      #Random effect
      if (tau2 != 0) {
        # cal. random effect weights and belonging pooled estimate
        wR <- 1 / (sig ^ 2 + tau2)
        vwR <- 1 / sum(wR) # variance of pooled effect (random)
        rwR <- wR * vwR
        teR <- sum(te * wR) / sum(wR)
        leR <- sum(log(te) * rwR)
        peR <- exp(leR)
        if (hksj == TRUE) {
          vwR <- 1 / df * sum(wR * (log(te) - leR) ^ 2 / sum(wR))
          if (is.null(sign)) {
            zvalR <- leR / sqrt(vwR)
          } else {
            zvalR <- sign * leR / sqrt(vwR)
          }
          pvalR <- 2 * pt(abs(zvalR), df = df, lower.tail = FALSE)
          lciR <- exp(leR - qt(1 - 0.05 / 2, df = df) * sqrt(vwR))
          uciR <- exp(leR + qt(1 - 0.05 / 2, df = df) * sqrt(vwR))
        } else {
          if (is.null(sign)) {
            zvalR <- leR / sqrt(vwR)
          } else {
            zvalR <- sign * leR / sqrt(vwR)
          }
          pvalR <- (1 - pnorm(abs(zvalR))) * 2
          lciR <- exp(leR - 1.96 * sqrt(vwR))
          uciR <- exp(leR + 1.96 * sqrt(vwR))
        }

        vw <- 1 / sum(w)
        D2 <- 1 - vw / vwR
        synth <-
          list(
            fw = round(w / sum(w) * 100, 1),
            peF = c(peF, lci, uci, zval, pval, lpeF, vw),
            rwR = rwR * 100,
            peR = c(peR, lciR, uciR, zvalR, pvalR, vwR),
            Q = c(Q, df, pQ),
            U = c(tau2, H, I2, D2)
          )
        class(synth) <- "synthesized"
        return(synth)
      } else {
        synth <-
          list(
            peF = c(peF, lci, uci, zval, pval, lpeF, vw),
            Q = c(Q, df, pQ),
            U = c(tau2, H, I2)
          )
        class(synth) <- "synthesized"
        return(synth)
      }
    } else {
      Q <- sum(w * te ^ 2) - (sum(w * te)) ^ 2 / sum(w)
      U <- sum(w) - sum(w ^ 2) / sum(w)
      tau2 <- ifelse(Q > df, (Q - df) / U, 0)
      wR <- 1 / (sig ^ 2 + tau2)
      vwR <- 1 / sum(wR) # variance of pooled effect (random)
      rwR <- wR * vwR
      peRest <- sum(te * wR) / sum(wR)
      lciR <- peRest - 1.96 * sqrt(vwR)
      uciR <- peRest + 1.96 * sqrt(vwR)
      zvalR <- peRest / sqrt(vwR)
      pvalR <- 2 * (1 - pnorm(abs(zvalR)))
      pQ <- pchisq(Q, df, lower.tail = FALSE)
      H <- sqrt(Q / df)
      I2 <- (Q - df) / Q
      D2 <- (1 - vw / vwR)

      synth <-
        list(
          fw = round(w / sum(w) * 100, 1),
          peF = c(peF, vw, lci, uci, zval, pval),
          w = w,
          rwR = round(rwR * 100, 1),
          peR = c(peRest, vwR, lciR, uciR, zvalR, pvalR),
          Q = c(Q, df, pQ),
          U = c(tau2, H, I2, D2)
        )
      class(synth) <- "synthesized"
      return(synth)
    }
  }
}
