# DOCUMENTATION ----

#' metaanalysis
#'
#' @param data Data.frame. If data is provided as a data set. Dataset must then containing arguments for
#'  meta-analysis. Either `study`, `eI`, `eC`, `nC` or `nI` for discrete data, or, `study`, `mI`, `mC`, `sdI` and `sdC` for cont.
#'  See details. If a data frame is not provided, the columns must be set in the function.
#' @param outcome Outcome metric for the studies. Choose between: cont, RR, RD or OR
#' @param vartype Variance type for continuous outcomes. Choices are "equal" or "non-equal". Defaults to "equal".
#' @param method Method for calculating weights. Options include: MH (Mantel-Haenzel), IV (Inverse variance weighting)
#'  or GLM
#' @param fixedStudy TRUE or FALSE. For using GLM methods for the meta-analysis, should the study effect be
#' fixed-effect. Defaults to TRUE.
#' @param hksj TRUE or FALSE. Should the Hartung-Knapp-Sidik-Jonkman adjustment be used to the random-effects meta-analysis. Defaults to FALSE.
#' @param sign TBA.
#' @param study XXX.
#' @param mI,mC,sdI,sdC See details.
#' @param eI,nI,eC,nC See details.
#' @param ... Additional variables. See Details.
#'
#' @return A list object containing two data frames.
#' 1. studyResults contains information about the individual studies
#' 2. metaResults contains information about the meta-analysis (traditional, non-sequential)
#' @export
#' @aliases print.metaanalysis
#'
#' @examples
#' data(perioOxy)
#' metaanalysis(outcome = "RR", data = perioOxy, study = perioOxy$trial)
#'

# FUNCTION | metaanalysis ----
#' @importFrom metafor rma.uni confint.rma.uni
#' @importFrom stats glm coef vcov pnorm qnorm pchisq pt qt binomial

metaanalysis <- function(data = NULL,
                         study = NULL,
                         mI = NULL, mC = NULL,
                         sdI = NULL, sdC = NULL,
                         eI = NULL, nI = NULL,
                         eC = NULL, nC = NULL,
                         outcome = "RR",
                         vartype = "equal",
                         method = "MH",
                         fixedStudy = TRUE,
                         hksj = FALSE,
                         sign = NULL,
                         alpha = 0.05,
                         tau.ci.method = "BJ",
                         ...) {

  # Check inputs ----
  # check | outcome
  if (!(outcome %in% c("cont", "RD", "RR", "OR"))) {
    stop("`outcome` must be 'cont', 'RD', 'RR' or 'OR'.")
  }
  # check | data; mI; mC; sdI; sdC; eI; nI; eC; nC.
  if(is.null(data)){
    if(outcome == "cont"){
      if(list(NULL) %in% list(mI,mC,sdI,sdC,nI,nC)){
        stop("'mI','mC','sdI,'sdC','nI','nC' must all be provided
           if `data` is not provided.")
      }
      data <- data.frame(mI=mI, mC=mC, sdI=sdI, sdC=sdC, nI=nI, nC=nC)
    }else{
      if(list(NULL) %in% list(eI,nI,eC,nC)){
        stop("'eI','nI','eC','nC' must all be provided
           if `data` is not provided.")
      }
      data <- data.frame(eI=eI, nI=nI, eC=eC, nC=nC)
    }

  }else{
    if(outcome == "cont" &
       !all(c("mI","mC","sdI","sdC","nI","nC") %in%
            colnames(data))){
      stop("`data` must have the following columns:
           'mI','mC','sdI,'sdC','nI', and 'nC'.")
    }else if(!all(c("eI","nI","eC","nC") %in%
                  colnames(data))){
      stop("`data` must have the following columns:
           'eI','nI,'eC', and 'nC'.")
    }
  }
  # check | study
  if(!any(colnames(data) == "study")){
    if(!is.null(study)){
      data <- cbind(study,data)
      missing_vec <- NULL
    }else{
      study <- c(1:nrow(data))
      missing_vec <-1
      data <- cbind(study,data)
    }
  }else{
    missing_vec <- NULL
  }
  # check | vartype
  if (!(vartype %in% c("equal","non-equal"))) {
    stop("`vartype` must be either 'equal' or 'non-equal'")
  }
  # check | method
  if (!(method %in% c("MH","IV","GLM"))) {
    stop("`method` must be either 'MH', 'GLM', or 'IV'")
  }
  # check | fixedStudy
  if (!(fixedStudy %in% c(TRUE,FALSE))) {
    stop("`fixedStudy` must be either TRUE or FALSE")
  }
  # check | hksj
  if (!(hksj %in% c(TRUE,FALSE))) {
    stop("`hksj` must be either TRUE or FALSE")
  }
  # check | alpha
  if (class(alpha) != "numeric" & alpha < 1) {
    stop("`alpha` must numeric and below 1")
  }
  # check | alpha
  if (!(tau.ci.method %in% c("BJ"))) {
    stop("`tau.ci.method` must be 'BJ'")
  }

  # Define - metaprepate ----
  metaPrepare <- function(data,outcome,method,vartype,alpha) {

    #Prepare dichotomous outcomes.
    if(outcome %in% c("OR", "RR", "RD")){

      # Remove studies with zero total events.
      if(sum(data$eI == 0 & data$eC == 0) > 0){
        nonevent <- which(data$eI == 0 & data$eC == 0)
        data <- data[-nonevent,]
      } else {
        nonevent <- NULL
      }

      # Adding 0.5 if one of the event counts is zero
      if(sum(data$eI == 0 | data$eC == 0) > 0){
        zc <- data$which(data$eI == 0 | data$eC == 0)
        data$eI[zc] <- data$eI[zc] + 0.5
        data$nI[zc] <- data$nI[zc] + 1
        data$eC[zc] <- data$eC[zc] + 0.5
        data$nC[zc] <- data$nC[zc] + 1
      }

      # Calculate event probability
      pI <- data$eI/data$nI
      pC <- data$eC/data$nC

      # Set the number of trials
      K <- length(data$eI)

      # Prepare returned object
      out <- list()

      # Calculate effects sizes and beloning std.
      if(outcome == "RD") { # Risk difference
        te <- pI-pC
        sig <- sqrt(pI*(1-pI)/data$nI+pC*(1-pC)/data$nC)
      }else if(outcome == "OR") { # Odds Ratio
        te <- (pI/(1-pI))/(pC/(1-pC))
        sig <- sqrt(1/data$eI+1/(data$nI-data$eI)+1/data$eC+1/(data$nC-data$eC))
      }else if(outcome == "RR") { # Risk Ratio
        te <- pI/pC
        sig <- sqrt(1/data$eI-1/data$nI+1/data$eC-1/data$nC)
      }

      # Calculate the weights
      if(method == "IV") { # Inverse Variance
        w <- 1/(sig^2)
        w <- w/sum(w)
        pe <- sum(te*w)
        w <- w*100
      }else if(method == "MH" | method == "GLM"){ # Mantel-Haenszel or GLM
        A <- data$eI # make the complicated variance (emerson)
        B <- data$nI - data$eI
        C <- data$eC
        D <- data$nC - data$eC
        N <- A + B + C + D

        if(outcome == "OR"){
          w <- (data$nI-data$eI)*data$eC/N
          T1 <- (A+D)/N
          T2 <- (B+C)/N
          T3 <- A*D/N
          T4 <- B*C/N
          vpe <- 0.5*((T1*T3)/(sum(T3)^2)+(T1*T4+T2*T3)/(sum(T3)*sum(T4))+T2*T4/(sum(T4)^2))
          svpe <- sum(vpe)
        }else if(outcome == "RR"){
          w <- (data$nI)*data$eC/N
          D1 <- ((A+B)*(C+D)*(A+C)-(A*C*N))/(N^2)
          R <- (A*(C+D))/(N)
          S <- (C*(A+B))/N
          svpe <- sum(D1)/(sum(R)*sum(S))
        }else if(outcome == "RD"){
          w <- (data$nI)*(data$nC)/N
          svpe <- sum(((A*B*data$nI)^3+(C*D*data$nC)^3)/
                        (data$nI*data$nC*N)^2)/(sum(data$nI*data$nC)/N)^2

        }
        pe <- sum(te*w)/sum(w)
      }

      # Calculate confidence limits
      if (outcome %in% c("RR", "OR")) {
        lower <- exp(log(te) - qnorm(1 - alpha / 2) * sig)
        upper <- exp(log(te) + qnorm(1 - alpha / 2) * sig)
      } else {
        lower <- te - qnorm(1 - alpha / 2) * sig
        upper <- te + qnorm(1 - alpha / 2) * sig
      }

      # Return results
      if(method == "MH"){
        out <- list(w = w, te = te, lower = lower, upper = upper,
                    pe = c(pe, svpe), sig = sig, outcome = outcome,
                    method = method, data = data, nonevent = nonevent)
      }else{
        out <- list(w = w, te = te, lower = lower, upper = upper, pe = pe,
                    sig = sig, outcome = outcome, method = method, data = data,
                    nonevent = nonevent)
      }

    }else if(outcome == "cont"){

      te <- data$mI - data$mC

      if(vartype == "equal"){
        spooled <- sqrt(((data$nI-1)*data$sdI^2+(data$nC-1)*data$sdC^2)/(data$nI+data$nC-2))
        vte <- (data$nI+data$nC)/(data$nI*data$nC)*spooled^2
        sete <- sqrt(vte)
        df <- (data$nI-1)+(data$nC-1)
        lower <- te - qt(0.975, df = df)*sete
        upper <- te + qt(0.975, df = df)*sete
      } else if(vartype == "non-equal"){
        vte <- sdI^2/data$nI + sdC^2/data$nC
        sete <- sqrt(vte)
        lower <- te - qnorm(0.975)*sete
        upper <- te + qnorm(0.975)*sete
      }

      w <- 1/vte
      pe <- sum(w*te)/sum(w)
      sig <- sqrt(vte)
      w <- w/sum(w)

      out <- list(w = w, te = te, lower = lower, upper = upper, pe = pe,
                  sig = sig, outcome = outcome, method = method, data = data,
                  nonevent = nonevent)
    }

    class(out) <- "synthPrepped"
    return(out)
  }

  # Use - metaprepare ----
  mp <- metaPrepare(data, outcome = outcome, method = method,
                      vartype = vartype, alpha = alpha)

  # Define - synthesize ----
  synthesize <- function(y,sign,fixedStudy,hksj,tau.ci.method) {

    w <- y$w   # collect objects
    sig <- y$sig
    te <- y$te
    pe <- y$pe
    # eI <- y$eI
    # eC <- y$eC
    # nI <- y$nI
    # nC <- y$nC
    df <- length(w) - 1
    data <- y$data
    ci.tau <- ""

    if (length(w) == 1) df <- 1 #check up on this

    #### For GLM
    if (y$method == "GLM") {
      bi <- data$nI - data$eI
      di <- data$nC - data$eC

      grpOut <- cbind(xi = c(rbind(data$eI, data$eC)), mi = c(rbind(bi, di)))
      k <- length(data$eI)
      eff <- rep(rep(1, k), 2)
      trial <- factor(rep(seq_len(k), each = 2))
      group <- rep(c(1, 0), times = k)
      eff <- eff * group

      if (fixedStudy == FALSE) {
        const <- rep(rep(1, k), 2)
        glmFit <- lme4::glmer(grpOut ~ -1 + eff + const + (1 | trial),
                              nAGQ = 7,
                              family = binomial)
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
      ### NOT GLM

      if (y$method != "MH")
        w <- 1 / (sig ^ 2) # weight inverse variance
        rw <- w / sum(w) # relative weight
        vw <- 1 / sum(w) # variance of pooled effect

      if (y$method == "cont") {
        peF <- sum(w * te) / sum(w)
        lci <- peF - 1.96 * sqrt(vw)
        uci <- peF + 1.96 * sqrt(vw)
        if (is.null(sign)) {
          zval <- peF / sqrt(vw)
        } else {
          zval <- sign * peF / sqrt(vw)
        }
        pval <- (1 - pnorm(abs(zval))) * 2

      } else if (y$method == "MH") {
        if(y$outcome == "RD"){
          vw <- pe[2]
          peF <- sum(te * w) / sum(w)
          lci <- peF - 1.96 * sqrt(vw)
          uci <- peF + 1.96 * sqrt(vw)
          if (is.null(sign)) {
            zval <- peF / sqrt(vw)
          } else {
            zval <- sign * peF / sqrt(vw)
          }
        } else {       # method the same for OR and RR
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
        }

        pval <- (1 - pnorm(abs(zval))) * 2

      }else{
        # if IV
        if(y$outcome == "RD"){
          peF <- sum(te * rw) # fixed effect log pooled estimate
          uci <- peF + 1.96 * sqrt(vw)
          lci <- peF - 1.96 * sqrt(vw)
          if(is.null(sign)){
            zval <- peF / sqrt(vw)
          }else{
            zval <- sign * peF / sqrt(vw)
          }
        } else {
          lpeF <- sum(log(te) * rw) # fixed effect log pooled estimate
          uci <- exp(lpeF + 1.96 * sqrt(vw))
          lci <- exp(lpeF - 1.96 * sqrt(vw))
          peF <- exp(lpeF)
          if(is.null(sign)){
            zval <- lpeF / sqrt(vw)
          }else{
            zval <- sign * lpeF / sqrt(vw)
          }
        }
        pval <- (1 - pnorm(abs(zval))) * 2

      }

      if(y$method != "cont"){
        w <- 1 / (sig ^ 2)
        if(y$outcome == "RD"){
          Q <- sum(w * te ^ 2) - (sum(w * te)) ^ 2 / sum(w)
        } else {
          Q <- sum(w * log(te) ^ 2) - (sum(w * log(te))) ^ 2 / sum(w) }
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
          peR <- sum(te * wR) / sum(wR)
          if(y$outcome != "RD") {
            leR <- sum(log(te) * rwR)
            peR <- exp(leR)
          }
          if (hksj == TRUE) {
            if (y$outcome == "RD") {
              vwR <- 1 / df * sum(wR * (te - peR) ^ 2 / sum(wR))
              if (is.null(sign)) {
                zvalR <- peR / sqrt(vwR)
              } else {
                zvalR <- sign * peR / sqrt(vwR)
              }
              pvalR <- 2 * pt(abs(zvalR), df = df, lower.tail = FALSE)
              lciR <- peR - qt(1 - 0.05 / 2, df = df) * sqrt(vwR)
              uciR <- peR + qt(1 - 0.05 / 2, df = df) * sqrt(vwR)
            } else {
              vwR <- 1 / df * sum(wR * (log(te) - leR) ^ 2 / sum(wR))
              if (is.null(sign)) {
                zvalR <- leR / sqrt(vwR)
              } else {
                zvalR <- sign * leR / sqrt(vwR)
              }
              pvalR <- 2 * pt(abs(zvalR), df = df, lower.tail = FALSE)
              lciR <- exp(leR - qt(1 - 0.05 / 2, df = df) * sqrt(vwR))
              uciR <- exp(leR + qt(1 - 0.05 / 2, df = df) * sqrt(vwR))
            }
          } else {
            if(y$outcome == "RD"){
              if (is.null(sign)) {
                zvalR <- peR / sqrt(vwR)
              } else {
                zvalR <- sign * peR / sqrt(vwR)
              }
              pvalR <- (1 - pnorm(abs(zvalR))) * 2
              lciR <- peR - 1.96 * sqrt(vwR)
              uciR <- peR + 1.96 * sqrt(vwR)
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
          }

          vw <- 1 / sum(w)
          D2 <- 1 - vw / vwR

          if(y$outcome %in% c("RR", "OR") & tau.ci.method == "BJ"){
            ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
              yi = log(te), sei = sig, method = "GENQ", weights = 1/sig^2))
          }

          if(y$outcome %in% c("RR", "OR") & tau.ci.method == "QP"){
            ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
              yi = log(te), sei = sig, method = "DL"))
          }

          synth <-
            list(
              fw = round(rw * 100, 1),
              peF = c(peF, lci, uci, zval, pval, log(peF), vw),
              rwR = rwR * 100,
              peR = c(peR, lciR, uciR, zvalR, pvalR, vwR),
              Q = c(Q, df, pQ),
              U = c(tau2, H, I2, D2),
              ci.tau = ci.tau
            )
          class(synth) <- "synthesized"

          return(synth)
        } else {
          synth <-
            list(
              peF = c(peF, lci, uci, zval, pval, log(peF), vw),
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

  # Use - synthesize ----
  sy <- synthesize(y = mp, sign = sign, fixedStudy = fixedStudy, hksj = hksj,
                   tau.ci.method = tau.ci.method)

  # Create - output ----

  # Check for non-event studies
  if(!is.null(mp$nonevent)){
    if(class(data$study[mp$nonevent]) == "integer"){
      nonevent <- paste(paste("Study",data$study[mp$nonevent]),collapse=", ")
    }else{
      nonevent <- paste(data$study[mp$nonevent],collapse=", ")
    }
  }else{ nonevent <- NULL }

  #Apply weights
  if(!is.null(sy$rwR)) {
    weightRandom = sy$rwR
  } else {
    weightRandom <- rep(NA, length(mp$w))
  }

  studyResults <- data.frame(
    study = mp$data$study,
    "ES" = mp$te,
    "stdError" = sqrt(mp$sig),
    "lowerCI" = mp$lower,
    "upperCI" = mp$upper,
    weightFixed = mp$w,
    weightRandom = weightRandom
  )
  colnames(studyResults)[2] <- outcome

  if(!is.null(sy$peR)){
    metaResults <- data.frame(
      type = c("Fixed", "Random"),
      "ES" = c(sy$peF[1], sy$peR[1]),
      "stdError" = c(sqrt(sy$peF[7]), sqrt(sy$peR[6])),
      "lowerCI" = c(sy$peF[2], sy$peR[2]),
      "upperCI" = c(sy$peF[3], sy$peR[3]),
      "pValue" = c(sy$peF[5], sy$peR[5]))
  }else{
    metaResults = data.frame(
      type = c("Fixed"),
      "ES" = c(sy$peF[1]),
      "stdError" = c(sqrt(sy$peF[7])),
      "lowerCI" = c(sy$peF[2]),
      "upperCI" = c(sy$peF[3]),
      "pValue" = c(sy$peF[5]))
  }
  colnames(metaResults)[2] <- outcome

  out = list(studyResults = studyResults, metaResults = metaResults,
             metaPrepare = mp, synthesize = sy, nonevent = nonevent,
             missing_vec = missing_vec)
  class(out) <- "metaanalysis"
  return(out)
}

# FUNCTION | print metaanalysis ----
#' @method print metaanalysis
#' @export
print.metaanalysis <- function(x,...){
  cat("Individual trial results: \n \n")
  y <- x$studyResults
  print(y)
  cat("\nNon-sequential metaanalysis results: \n \n")
  y <- x$metaResults
  print(y)
  invisible(x)
  if(x$metaPrepare$method == "GLM"){
    message("\n NB. Only fixed-effect is analysed since method is GLM")
  }
  #ZERO TRIAL
  if(!is.null(x$missing_vec)){
    message("\n NB. Please provide a study vector to name studies.")
  }
  if(!is.null(x$nonevent)){
    message(paste("\n NB.",x$nonevent,"was excluded from the analysis due to zero events, consider changing outcome to RD"))
  }

}

# FUNCTION | plot metaanalysis ----
#' @method plot metaanalysis
#' @export
plot.metaanalysis <- function(x, type="both", ...){
#TODO: Size depending on weight, prioritise Random

  if(type=="both") cat("hej")

  plot <- merge(x$metaPrepare$data,x$studyResults)
  results <- x$metaResults
  colnames(results)[1] <- "study"
  results[,colnames(plot)[!(colnames(plot) %in% colnames(results))]] <- NA
  results$study[results$study == "Fixed"] <- "Fixed-effect"
  results$study[results$study == "Random"] <- "Random-effect"

  plot <- rbind(plot, results[names(plot)])
  plot <- cbind(nrow(plot):1,plot)
  colnames(plot)[1] <- "yaxis"

  outcome <- colnames(plot)[8]
  colnames(plot)[8] <- "outcome"

  plot$out_ci <- paste0(round(plot$outcome,2)," (",round(plot$lowerCI,2),"-",
                        round(plot$upperCI,2),")")
  shapes <- grepl("Fixed-ef|Random-ef",plot$study)*2+21
  sizes <- (plot$weightRandom/100+0.5)*4
  sizes[is.na(sizes)] <- 3

  colors <- rep("black",nrow(plot))
  colors[grepl("Fixed-ef",plot$study)] <- "#0053a3"
  colors[grepl("Random-ef",plot$study)] <- "#b30000"

  nchar_n <- max(nchar(c(plot$nI,plot$nC)),na.rm=T)

  left_theme <- theme(plot.margin = margin(),
    axis.title.x.bottom = element_text(color="white"),
    axis.text.x.bottom = element_text(color="white"),
    axis.text.x.top = element_text(hjust=1),
    axis.title.x.top = element_text(hjust=0.5),
    axis.ticks.x = element_blank(), axis.line.x = element_blank(),
    axis.title.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(hjust=0, color="black"),
    axis.line.y = element_blank())

  left_1 <- ggplot(data = plot, aes(y = yaxis)) +
    annotate("text",x = 0, y = plot$yaxis, label = plot$eC, hjust = 1) +
    annotate("text",x = nchar_n/2, y = plot$yaxis, label = plot$nC,hjust = 1) +
    scale_x_continuous(sec.axis = dup_axis(),limits=c(-1,NA),
      name="Control", breaks = c(0,nchar_n/2), labels=c("eC","nC")) +
    scale_y_continuous(breaks=plot$yaxis,labels=plot$study) +
    theme_classic() + left_theme +
    theme(axis.text.y = element_text(size=10))

  left_2 <- ggplot(data = plot, aes(y = yaxis)) +
    annotate("text",x = 0, y = plot$yaxis, label = plot$eI, hjust = 1) +
    annotate("text",x = nchar_n/2, y = plot$yaxis, label = plot$nI,hjust = 1) +
    scale_x_continuous(sec.axis = dup_axis(),limits=c(-1,NA),
      name="Experimental", breaks = c(0,nchar_n/2), labels=c("eI","nI")) +
    scale_y_continuous(breaks=plot$yaxis,labels=plot$study) +
    theme_classic() + left_theme + theme(axis.text.y = element_blank())

  middle <-
   ggplot(plot,aes(x=outcome,xmin=lowerCI,xmax=upperCI,y=yaxis)) +
    geom_vline(xintercept = 1, color="gray", linetype=3) +
    geom_segment(aes(x=plot$outcome[grepl("Fixed-effect",plot$study)],
                     xend=plot$outcome[grepl("Fixed-effect",plot$study)],
                     y=Inf,yend=plot$yaxis[grepl("Fixed-effect",plot$study)]),
                     color = "#7cbfff") +
    geom_segment(aes(x=plot$outcome[grepl("Random-effect",plot$study)],
                     xend=plot$outcome[grepl("Random-effect",plot$study)],
                     y=Inf,yend=plot$yaxis[grepl("Random-effect",plot$study)]),
                 color = "#ff8c8c") +
     geom_segment(aes(x=min(plot$lowerCI),xend=max(plot$upperCI),y=-Inf,yend=-Inf)) +
    geom_point(shape = shapes, color = colors, fill = colors,size=sizes) +
    geom_errorbar(color = colors, width=0) +
    theme_classic() +
    scale_colour_identity() +
    scale_x_continuous(trans="log10",sec.axis = dup_axis(), name=paste(outcome, "(95%CI)")) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          plot.margin = margin(), axis.line = element_blank(),
          axis.ticks.x.top = element_blank(), axis.text.x.top = element_text(color="white"), axis.line.y.left = element_blank(), axis.ticks.y.left = element_blank())

  right_1 <- ggplot(data = plot, aes(y = yaxis)) +
    annotate("text",x = 0, y = plot$yaxis, label = plot$out_ci, hjust = 1) +
    annotate("text",x = 3, y = plot$yaxis, label = round(plot$weightFixed,2), hjust = 1) +
    annotate("text",x = 6, y = plot$yaxis, label = round(plot$weightRandom,2), hjust = 1) +
    scale_x_continuous(sec.axis = dup_axis(), limits=c(-6,NA),
                       name=" ", breaks = c(0,3,6), labels=c(paste(outcome, "(95%CI)"),"F-Weight (%)","R-Weight (%)")) +
    scale_y_continuous(breaks=plot$yaxis,labels=plot$study) +
    theme_classic() + left_theme + theme(axis.text.y = element_blank())
      # coord_fixed(.50)

  grid.arrange(left_1,left_2,middle,right_1,ncol=4,widths=c(1.5,1,4,3))

}
