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
#' @param alpha Type-I-error.
#' @param beta Type-II-error. Not used unless we want to make a sample size calculation.
#' @param tau.ci.method Methods for computation of CI for tau.
#' @param delta Minimum clinically relevant value
#' @param z_thres Threshold value for z-score
#' @param power_calc Sample size calculation for fixed-effect and random-effects meta-analysis.
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
                         method = "IV",
                         fixedStudy = TRUE,
                         hksj = FALSE,
                         sign = NULL,
                         alpha = 0.05,
                         beta = 0.2,
                         tau.ci.method = "BJ",
                         delta = NULL,
                         z_thres = NULL,
                         power_calc = FALSE,
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
           if `data` is not provided. If outcome is continuous, set outcome = 'cont'")
      }
      data <- data.frame(eI=eI, nI=nI, eC=eC, nC=nC)
    }

  }else{
    if(outcome == "cont" &
       !all(c("mI","mC","sdI","sdC","nI","nC") %in%
            colnames(data))){
      stop("`data` must have the following columns:
           'mI','mC','sdI,'sdC','nI', and 'nC'.")
    }else if(outcome != "cont" & !all(c("eI","nI","eC","nC") %in%
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
  if (!is.numeric(alpha) & alpha < 1) {
    stop("`alpha` must numeric and below 1")
  }
  # check | alpha
  if (!(tau.ci.method %in% c("BJ", "QP"))) {
    stop("`tau.ci.method` must be 'BJ' or 'QP'")
  }

  # Use - metaprepare ----
  mp <- metaPrepare(data, outcome = outcome, method = method,
                      vartype = vartype, alpha = alpha)

  # Use - synthesize ----
  sy <- synthesize(y = mp, sign = sign, fixedStudy = fixedStudy, hksj = hksj,
                   tau.ci.method = tau.ci.method)

  # Create - output ----

  # Check for non-event studies
  if(!is.null(mp$nonevent)){
    if(is.integer(data$study[mp$nonevent])){
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

  # if study name is not unique - add number and give warning:
  if(length(unique(mp$data$study)) != length(mp$data$study)){
    warning("Study names are not unique - change names of studies.")
    dup <- duplicated(mp$data$study)
    mp$data$study[dup] <- paste0(mp$data$study[dup], 1:length(sum(dup)))
  }

  studyResults <- data.frame(
    study = mp$data$study,
    "ES" = mp$te,
    "stdError" = sqrt(mp$sig),
    "lowerCI" = mp$lower,
    "upperCI" = mp$upper,
    weightFixed = sy$fw,
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

    heteResults <- list(
      "Q" = sy$Q[1],
      "Q_df" = sy$Q[2],
      "Q_pval" = sy$Q[3],
      "tau2" = sy$U[1],
      "I^2" = sy$U[3],
      "D^2" = sy$U[4],
      "CI_heterogen" = sy$ci.tau
    )
  }else{
    metaResults = data.frame(
      type = c("Fixed"),
      "ES" = c(sy$peF[1]),
      "stdError" = c(sqrt(sy$peF[7])),
      "lowerCI" = c(sy$peF[2]),
      "upperCI" = c(sy$peF[3]),
      "pValue" = c(sy$peF[5]))

    heteResults = NULL
  }

  colnames(metaResults)[2] <- outcome

  out <- list(studyResults = studyResults, metaResults = metaResults,
             heteResults = heteResults,
             metaPrepare = mp, synthesize = sy, nonevent = nonevent,
             missing_vec = missing_vec)

  if(power_calc){
    if(outcome %in% c("RR", "OR")){
      if(!is.null(sy$peR)){out_ris <- ris(outcome = outcome, delta = delta, data = data,
                                          alpha = alpha, beta = beta, random = TRUE,
                                          var_random = sy$peR[6], type = "retrospective",
                                          ...)}
      out_ris$NF <- out_ris$NF - (sum(data$nI)+sum(data$nC))
      out_ris$NR_inc <- out_ris$NR_inc - (sum(data$nI)+sum(data$nC))
      out_ris$NR_div <- out_ris$NR_div - (sum(data$nI)+sum(data$nC))
      out <- append(out, list(out_ris = out_ris))
    } else {
      # TODO
      ris(outcome = outcome, delta = 0.7, data = perioOxy, alpha = alpha,
          beta = beta, random = TRUE)
    }
  }

  if(!is.null(delta)){
    if(is.null(z_thres)) z_thres <- 1.96
    if(outcome %in% c("RR", "OR"))  delta <- log(delta)
    z_score_fe <- delta/sqrt(sy$peF[7])
    pwr_fe <- 1-pnorm(z_thres - z_score_fe)+pnorm(-z_thres - z_score_fe)

    if(!is.null(sy$peR)){z_score_re <- delta/sqrt(sy$peR[6])
    pwr_re <- 1-pnorm(z_thres - z_score_re)+pnorm(-z_thres - z_score_re)
    out <- append(out, list(pwr_re = pwr_re))
    }

    out <- append(out, list(pwr_fe = pwr_fe))
  }

  class(out) <- "metaanalysis"
  return(out)
}

# Define - metaprepare ----
metaPrepare <- function(data, outcome, method, vartype, alpha, nonevent = NULL) {

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
      zc <- which(data$eI == 0 | data$eC == 0)
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


# Define - synthesize ----
synthesize <- function(y,sign,fixedStudy,hksj,tau.ci.method) {

  w <- y$w   # collect objects
  sig <- y$sig
  te <- y$te
  pe <- y$pe
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

    if (y$method != "MH") w <- 1 / (sig ^ 2) # weight inverse variance
    rw <- w / sum(w) # relative weight
    vw <- 1 / sum(w) # variance of pooled effect

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

    } else if (y$method == "MH") {
      rw <- w / sum(w)
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

    }else {
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

    if(y$outcome != "cont"){
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
            fw = round(rw/sum(rw) * 100, 4),
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
            fw = round(rw/sum(rw) * 100, 4),
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

      if(tau.ci.method == "BJ"){
        ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
          yi = te, sei = sig, method = "GENQ", weights = 1/sig^2))
      }

      if(tau.ci.method == "QP"){
        ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
          yi = te, sei = sig, method = "DL"))
      }

      synth <-
        list(
          fw = round(w / sum(w) * 100, 4),
          peF = c(peF, lci, uci, zval, pval, NA, vw),
          rwR = round(rwR * 100, 4),
          peR = c(peRest, lciR, uciR, zvalR, pvalR, vwR),
          Q = c(Q, df, pQ),
          U = c(tau2, H, I2, D2),
          ci.tau = ci.tau
        )
      class(synth) <- "synthesized"
      return(synth)
    }
  }
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
  if(!is.null(x$missing_vec)){
    message("\n NB. Please provide a study vector to name studies.")
  }
  #Total-zero trials
  if(!is.null(x$nonevent)){
    message(paste("\n NB.",x$nonevent,"was excluded from the analysis due to zero events, consider changing outcome to RD"))
  }

}

# FUNCTION | plot metaanalysis ----
#' @method plot metaanalysis
#' @importFrom ggplot2 annotate labs scale_colour_identity geom_errorbar
#' @export
plot.metaanalysis <- function(x, type="both", ...){

  # FOR TESTING
  x <- readRDS("C:/Oel/Artikler/CTU/CTU_RTSA/cord.Rdata")
  x <- metaanalysis(x)
  type = "both"
  library(ggplot2)
  x <- get(load("C:/Oel/Artikler/CTU/CTU_RTSA/RTSA/data/perioOxy.RData"))
  x <- metaanalysis(x)
  type = "random"
  # FOR TESTING

  # Create dataframe for plot
  fplot <- merge(x$metaPrepare$data,x$studyResults)
  results <- x$metaResults
  colnames(results)[1] <- "study"
  results[,colnames(fplot)[!(colnames(fplot) %in% colnames(results))]] <- NA
  results$study[results$study == "Fixed"] <- "Fixed-effect"
  results$study[results$study == "Random"] <- "Random-effects"
  fplot <- rbind(fplot, results[names(fplot)])

  if(type=="fixed"){fplot <- fplot[!grepl("Random-ef",fplot$study),]}
  if(type=="random"){fplot <- fplot[!grepl("Fixed-ef",fplot$study),]}

  fplot <- cbind(nrow(fplot):1,fplot)
  colnames(fplot)[1] <- "yaxis"

  fplot$eI[grepl("Fixed-ef|Random-ef",fplot$study)] <- sum(fplot$eI,na.rm=T)
  fplot$nI[grepl("Fixed-ef|Random-ef",fplot$study)] <- sum(fplot$nI,na.rm=T)
  fplot$eC[grepl("Fixed-ef|Random-ef",fplot$study)] <- sum(fplot$eC,na.rm=T)
  fplot$nC[grepl("Fixed-ef|Random-ef",fplot$study)] <- sum(fplot$nC,na.rm=T)

  outcome <- colnames(x$studyResults)[2]
  colnames(fplot)[which(colnames(fplot) == outcome)] <- "outcome"

  fplot$out_ci <- paste0(sprintf(fplot$outcome, fmt = '%#.2f')," (",
                         sprintf(fplot$lowerCI, fmt = '%#.2f'),"-",
                         sprintf(fplot$upperCI, fmt = '%#.2f'),")")
  fplot$controls <- paste0(fplot$eC,"/",fplot$nC)
  fplot$controls[fplot$controls == "NA/NA"] <- ""
  fplot$experimental <- paste0(fplot$eI,"/",fplot$nI)
  fplot$experimental[fplot$experimental == "NA/NA"] <- ""
  fplot$wF <- round(fplot$weightFixed)
  fplot$wF[fplot$wF < 10 & !is.na(fplot$wF)] <- round(fplot$weightFixed[fplot$wF <10 & !is.na(fplot$wF)],1)
  fplot$wR <- round(fplot$weightRandom)
  fplot$wR[fplot$wR < 10 & !is.na(fplot$wR)] <- round(fplot$weightRandom[fplot$wR <10 & !is.na(fplot$wR)],1)

  # Ensure fixed effect model, if Tau = 0
  if(x$synthesize$U[1] == 0) type <- "fixed"; message <- "No heterogeneity, only fixed effect estimates presented"


  #Define shapes and colors
  shapes <- grepl("Fixed-ef|Random-ef",fplot$study)*2+21
  sizes <- (fplot$weightRandom/100+0.5)*4
  sizes[is.na(sizes)] <- 3

  colors <- rep("black",nrow(fplot))
  colors[grepl("Fixed-ef",fplot$study)] <- "#0053a3"
  colors[grepl("Random-ef",fplot$study)] <- "#b30000"

  # define x-axis
  if(is.null(xlims)) xlims <- c(min(fplot$lowerCI),max(fplot$upperCI))
  if(is.na(xlims[1])) xlims[1] <- min(fplot$lowerCI)
  if(is.na(xlims[2])) xlims[2] <- max(fplot$upperCI)

  xl <- list(`l0` = 10^(log10(xlims[1])-0.9),
             `l1` = 10^(log10(xlims[1])-0.6),
             `l2` = 10^(log10(xlims[1])-0.1),
             `r1` = 10^(log10(xlims[2])+0.9),
             `r2` = 10^(log10(xlims[2])+1.3),
             `r3` = 10^(log10(xlims[2])+1.6))
  if(type=="random"){ xl$r3 <- xl$r2; xl$r2 <- NULL }
  if(type=="fixed") xl$r3 <- NULL

  xl2 <- c(unlist(unname(xl[-1])),mean(c(xl$r2,xl$r3),na.rm=T))
  xl2 <- xl2[order(xl2)]

  #Axis text
  xlabs <- c("Experimental","Control",paste0("",outcome,"(95% CI)"),
             if(type!="fixed"){"Random"},"Weights (%)\n",
             if(type!="random"){"Fixed"})

  #Lower left text
  heterogen <- paste0("Heterogeneity: Tau² = ",
                      sprintf(x$synthesize$U[1], fmt = '%#.2f'), " (",
                      sprintf(x$synthesize$ci.tau$random["tau^2","ci.lb"], fmt = '%#.2f'),"-",
                      sprintf(x$synthesize$ci.tau$random["tau^2","ci.ub"], fmt = '%#.2f'),")",
                      ", Q = ", sprintf(x$synthesize$Q[1], fmt = '%#.1f'),
                      ", df = ", round(x$synthesize$Q[2]),
                      ", I² = ", round(x$synthesize$U[3]*100),
                      "%\n")
  if(type!="random"){
    heterogen <- paste0(heterogen,"Overall effect: Fixed-effect, z = ",
                        sprintf(x$synthesize$peF[4], fmt = '%#.2f'),
                        " (p = ",
                        sprintf(x$metaResults[x$metaResults == "Fixed","pValue"], fmt = '%#.4f'),
                        ")")
  }
  if(type!="fixed"){
    heterogen <- paste0(heterogen,"; Random-effects, z = ",
                        sprintf(x$synthesize$peR[4], fmt = '%#.2f'),
                        " (p = ",
                        sprintf(x$metaResults[x$metaResults == "Random","pValue"], fmt = '%#.4f'),
                        ")")
  }

  #The plot
  ggplot(fplot,aes(x=outcome,xmin=lowerCI,xmax=upperCI,y=yaxis)) +
    geom_vline(xintercept = 1, color="gray", linetype=3) +
    geom_segment(aes(x=xlims[1],xend=xlims[2],y=-Inf,yend=-Inf)) +
    annotate("text",x=xl$l1,y=fplot$yaxis,label=fplot$experimental,hjust=1) +
    annotate("text",x=xl$l2,y=fplot$yaxis,label=fplot$control,hjust=1) +

    #Confidence intervals
    {if(type!="fixed")
      geom_segment(aes(x=fplot$outcome[grepl("Random-effect",fplot$study)],
                       xend=fplot$outcome[grepl("Random-effect",fplot$study)],
                       y=Inf,yend=fplot$yaxis[grepl("Random-effect",fplot$study)]),
                   color = "#ff8c8c")} +
    {if(type!="random") geom_segment(aes(x=fplot$outcome[grepl("Fixed-effect",fplot$study)],
                                         xend=fplot$outcome[grepl("Fixed-effect",fplot$study)],
                                         y=Inf,yend=fplot$yaxis[grepl("Fixed-effect",fplot$study)]),
                                     color = "#7cbfff")} +

    annotate("text",x=xl$r1,y=fplot$yaxis,label=fplot$out_ci,hjust=1) +
    {if(type!="fixed") annotate("text",x=xl$r3,y=fplot$yaxis,label=fplot$wR,hjust=1) } +
    {if(type!="random") annotate("text",x=xl$r2,y=fplot$yaxis,label=fplot$wF,hjust=1) } +
    geom_segment(aes(x=0,xend=xl$l2, # dots before summary study side
                     y=max(fplot$yaxis[grepl("Fixed-ef|Random-ef",fplot$study)])+0.5,
                     yend=max(fplot$yaxis[grepl("Fixed-ef|Random-ef",fplot$study)])+0.5),
                 linetype=3) +
    geom_segment(aes(x=xlims[2],xend=Inf, # dots after
                     y=max(fplot$yaxis[grepl("Fixed-ef|Random-ef",fplot$study)])+0.5,
                     yend=max(fplot$yaxis[grepl("Fixed-ef|Random-ef",fplot$study)])+0.5),
                 linetype=3) +
    labs(tag=heterogen)+
    geom_point(shape = shapes, color = colors, fill = colors,size=sizes) +
    geom_errorbar(color = colors, width=0) +
    theme_classic() +
    scale_colour_identity() +
    scale_x_continuous(trans="log10",
                       sec.axis = sec_axis(~.,breaks = xl2,label=xlabs),
                       limits=c(xl$l0,if(type!="fixed"){xl$r3}else{xl$r2}), name="o\no",
                       breaks=c(round(1/((xlims[2]-1)/2),nchar(round(xlims[2]))),1,round((xlims[2]-1)/2,nchar(round(xlims[2]))))) + # ANNE: Lower bound cannot handle upper over 100
    scale_y_continuous(breaks=fplot$yaxis, label=fplot$study) +
    theme(axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.ticks.x.top = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y.left = element_blank(),
          axis.title.x = element_text(color="white"),
          axis.text.x.top = element_text(hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10,face="bold"),
          plot.tag.position = c(0,0), plot.tag = element_text(hjust=0, vjust=0, size=9))


  if(nchar(message) > 0){message(message)}

}
