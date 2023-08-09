#' Fixed-effect or random-effects meta-analysis
#' 
#' @description
#' Computes a fixed-effect or random-effects meta-analysis including heterogeneity statistics. If \code{mc} is specified, a retrospective sample and trial size is calculated.
#' 
#' 
#'
#' @param outcome Outcome metric for the studies. Choose between: MD (mean difference), RR (relative risk), RD (risk difference) or OR (odds ratio).
#' @param data A data.frame containing the study results. The data set must containing a specific set of columns. These are respectively `eI` (events in intervention group), `eC` (events in control group), `nC` (participants intervention group) or `nI` (participants control group) for discrete data, or, `mI` (mean intervention group), `mC` (mean control group), `sdI` (standard error intervention group), `sdC` (standard error control group),`nC` (participants intervention group) and `nI` (participants control group)  for continuous outcomes. Preferable also a `study` column as an indicator of study.
#' @param study Optional vector of study IDs. If no study indicator is provided in `data`, a vector of study indicators e.g. names.
#' @param alpha The level of type I error as a percentage, the default is 0.05 corresponding to 5\%.
#' @param beta The level of type II error as a percentage, the default is 0.1 corresponding to 10\%. Not used unless a sample and trial size calculation is wanted. 
#' @param side Whether a 1- or 2-sided hypothesis test is used. Options are 1 or 2. Default is 2.
#' @param weights Method for calculating weights. Options are "MH" (Mantel-Haenzel and only optional for binary data) or "IV" (Inverse variance weighting). Default is "IV".
#' @param re_method Methods are "DL" for DerSimonian-Laird or "DL_HKSJ" for DerSimonian-Laird with Hartung-Knapp-Sidik-Jonkman adjustment. Default is "DL_HKSJ".
#' @param tau_ci_method Methods for computation of confidence interval for heterogeneity estimate tau. Calls rma.uni from the metafor package. Options are "BJ" and "QP". Default is "BJ"
#' @param cont_vartype Variance type for continuous outcomes. Choices are "equal" (homogeneity of treatment group variances) or "non-equal" (heterogeneity of treatment group variances). Default is "equal".
#' @param mc Minimum clinically relevant value. Used for sample and trial size calculation.
#' @param RRR Relative risk reduction. Used for binary outcomes with outcome metric RR. Argument mc can be used instead. Must be a value between 0 and 1.
#' @param conf_level Confidence interval coverage
#' @param sd_mc The expected standard deviation. Used for sample and trial size calculation for mean differences.
#' @param zero_adj Zero adjustment for null events in binary data. Options for now is 0.5. Default is 0.5.
#' @param ... Additional variables.
#'
#' @return A \code{metaanalysis} object which is a list with 6 or 7 elements.
#' \item{study_results}{A data.frame containing study results which is information about the individual studies}
#' \item{meta_results}{A data.frame containing the results of the meta-analysis such as the pooled estimate, its standard error, confidence interval and p-value}
#' \item{hete_results}{A list containing statistics about hetergeneity.}
#' \item{metaPrepare}{A list containing the elements used for calculating the study results.}
#' \item{synthesize}{A list containing the elements used for calculating the meta-analysis results.}
#' \item{settings}{A list containing the arguments used in the \code{metaanalysis} call.}
#' \item{ris}{(Only when \code{mc} has been specified or meta-analysis is created as part of \code{RTSA}). List of sample and trial size calculation. See documentation for \code{ris}.}
#' 
#' @export
#' @aliases print.metaanalysis
#'
#' @examples
#' ### Basic uses
#' # Use perioOxy data from package and run meta-analysis with default settings
#' data(perioOxy)
#' metaanalysis(outcome = "RR", data = perioOxy, study = perioOxy$trial)
#'
#' # Run same meta-analysis but with odds ratio as outcome metric, Mantel-Haenzel 
#' # weights and DerSimionian-Laird for the variance estimate 
#' metaanalysis(outcome = "OR", data = perioOxy, study = perioOxy$trial,
#'  weights = "MH", re_method = "DL")
#'  
#' # Run meta-analysis with mean difference as outcome metric
#' data(eds)
#' metaanalysis(outcome = "MD", data = eds)
#' 
#' ### Retrospective sample size calculation
#' # minimal clinical relevant difference set to an odds ratio of 0.7.
#' ma <- metaanalysis(outcome = "OR", data = perioOxy, mc = 0.7)
#' ma$ris

# FUNCTION | metaanalysis ----
#' @importFrom metafor rma.uni confint.rma.uni
#' @importFrom stats glm coef vcov pnorm qnorm pchisq pt qt binomial

metaanalysis <- function(outcome,
                         data,
                         side = 2,
                         alpha = 0.05,
                         beta = 0.1,
                         weights = "IV",
                         re_method = "DL_HKSJ",
                         tau_ci_method = "BJ",
                         cont_vartype = "equal",
                         mc = NULL,
                         RRR = NULL,
                         sd_mc = NULL,
                         study = NULL,
                         conf_level =  0.95,
                         zero_adj = 0.5,
                         ...) {
  # Check inputs ----
  # check | outcome
  if (!(outcome %in% c("MD", "RD", "RR", "OR"))) {
    stop("`outcome` must be 'MD', 'RD', 'RR' or 'OR'.")
  }
  # check | data; mI; mC; sdI; sdC; eI; nI; eC; nC.
  if (is.null(data)) {
    stop("A data.frame containing either: eI, nI, eC, nC")
  } else{
    if (outcome == "MD" &
        !all(c("mI", "mC", "sdI", "sdC", "nI", "nC") %in%
             colnames(data))) {
      stop("`data` must have the following columns:
           'mI','mC','sdI,'sdC','nI', and 'nC'.")
    } else if (outcome != "MD" & !all(c("eI", "nI", "eC", "nC") %in%
                                      colnames(data))) {
      stop("`data` must have the following columns:
           'eI','nI,'eC', and 'nC'.")
    }
  }
  # check | study
  if (!any(colnames(data) == "study")) {
    if (!is.null(study)) {
      data <- cbind(study, data)
      missing_vec <- NULL
    } else{
      study <- c(1:nrow(data))
      missing_vec <- 1
      data <- cbind(study, data)
    }
  } else{
    missing_vec <- NULL
  }
  
  # check | mc or RRR
  if(is.null(mc) & !is.null(RRR) & outcome == "RR"){
    mc <- 1 - RRR
  }
  
  # check | side
  if (!(side %in% c(1, 2))) {
    stop("`side` must be either 1 or 2")
  }
  # check | cont_vartype
  if (!(cont_vartype %in% c("equal", "non-equal"))) {
    stop("`cont_vartype` must be either 'equal' or 'non-equal'")
  }
  # check | weights
  if (!(weights %in% c("MH", "IV"))) {
    stop("`weights` must be either 'MH', or 'IV'")
  }
  # check | re_method
  if (!(re_method %in% c("DL", "DL_HKSJ"))) {
    stop("`re_method` must be either DL or DL_HKSJ")
  }
  # check | alpha
  if (!is.numeric(alpha) | alpha > 1 | alpha < 0) {
    stop("`alpha` must numeric and between 0 and 1")
  }
  # check | beta
  if (!is.numeric(beta) | beta > 1 | beta < 0) {
    stop("`beta` must numeric and between 0 and 1")
  }
  # check | ci_method for tau
  if (!(tau_ci_method %in% c("BJ", "QP"))) {
    stop("`tau_ci_method` must be 'BJ' or 'QP'")
  }
  
  argg <- c(as.list(environment()), list(...))

  # Use - metaprepare ----
  # calculate estimates for each study
  mp <- metaPrepare(
    data,
    outcome = outcome,
    weights = weights,
    cont_vartype = cont_vartype,
    alpha = alpha,
    conf_level = conf_level,
    zero_adj = zero_adj
  )
  
  # Use - synthesize ----
  # calculate actual meta-analysis
  sy <- synthesize(
    y = mp,
    re_method = re_method,
    conf_level = conf_level,
    tau_ci_method = tau_ci_method
  )

  # Create - output ----
  # Check for non-event studies (total zero)
  if (!is.null(mp$nonevent)) {
    if (is.integer(data$study[mp$nonevent])) {
      nonevent <-
        paste(paste("Study", data$study[mp$nonevent]), collapse = ", ")
    } else{
      nonevent <- paste(data$study[mp$nonevent], collapse = ", ")
    }
  } else{
    nonevent <- NULL
  }

  #Apply weights
  if (!is.null(sy$rwR)) {
    w_random = sy$rwR
  } else {
    w_random <- rep(NA, length(mp$w))
  }

  # if study name is not unique - add number and give warning:
  if (length(unique(mp$data$study)) != length(mp$data$study)) {
    warning("Study names are not unique - change names of studies.")
    dup <- duplicated(mp$data$study)
    mp$data$study[dup] <-
      paste0(mp$data$study[dup], 1:length(sum(dup)))
  }
  
  # collect study results
  study_results <- data.frame(
    study = mp$data$study,
    "ES" = mp$te,
    "se" = sqrt(mp$sig),
    "lower" = mp$lower,
    "upper" = mp$upper,
    w_fixed = sy$fw,
    w_random = w_random
  )
  colnames(study_results)[2] <- outcome
  colnames(study_results)[c(4,5)] <- paste0(colnames(study_results)[c(4,5)], ".", 100*(1-alpha), "CI")

  # collect meta-analysis results
  if (!is.null(sy$peR)) { # if any heterogeneity collect both fe and re
    meta_results <- data.frame(
      type = c("Fixed", "Random"),
      "ES" = c(sy$peF[1], sy$peR[1]),
      "se" = c(sqrt(sy$peF[7]), sqrt(sy$peR[6])),
      "lower" = c(sy$peF[2], sy$peR[2]),
      "upper" = c(sy$peF[3], sy$peR[3]),
      "pValue" = c(sy$peF[5], sy$peR[5])
    )

    hete_results <- list(
      hete_est = data.frame(
        "Q" = sy$Q[1],
        "Q_df" = sy$Q[2],
        "Q_pval" = sy$Q[3],
        "tau2" = sy$U[1],
        "I^2" = sy$U[3],
        "D^2" = sy$U[4]
      ),
      "CI_heterogen" = sy$ci.tau
    )
  } else { # collect only fixed-effect
    meta_results = data.frame(
      type = c("Fixed"),
      "ES" = c(sy$peF[1]),
      "se" = c(sqrt(sy$peF[7])),
      "lower" = c(sy$peF[2]),
      "upper" = c(sy$peF[3]),
      "pValue" = c(sy$peF[5])
    )

    hete_results <- list(
      hete_est = data.frame(
        "Q" = sy$Q[1],
        "Q_df" = sy$Q[2],
        "Q_pval" = sy$Q[3],
        "tau2" = sy$U[1],
        "I^2" = sy$U[3],
        "D^2" = sy$U[4]
      ),
      "CI_heterogen" = sy$ci.tau
    )
  }

  if(outcome %in% c("RR", "OR", "RD")){
    colnames(study_results)[3] <- paste0("se(log(", outcome, "))")
    colnames(meta_results)[3] <- paste0("se(log(", outcome, "))")
  }
  
  colnames(meta_results)[c(4,5)] <- paste0(colnames(meta_results)[c(4,5)], ".", 100*(1-alpha), "CI")
  colnames(meta_results)[2] <- outcome

  out <-
    list(
      study_results = study_results,
      meta_results = meta_results,
      hete_results = hete_results,
      metaPrepare = mp,
      synthesize = sy,
      settings = c(argg, list("nonevent" = nonevent, "missing_vec" = missing_vec))
    )
  
  if(outcome == "MD" & is.null(sd_mc)){
    if(sy$U[1] > 0){
      sd_mc = sqrt(sy$peR[6])
    } else {
      sd_mc = sqrt(sy$peF[7])
    }
  }

  # calculate retrospective sample and trial size if mc is provided
  if (!is.null(mc)) {
    if (outcome %in% c("RR", "OR")) {
      if (!is.null(sy$peR)) {
        out_ris <- ris(
          outcome = outcome,
          mc = mc,
          ma = list("metaPrepare" = mp, "synthesize" = sy),
          alpha = alpha,
          beta = beta,
          fixed = FALSE,
          side = side,
          type = "retrospective",
          ...
        )
      } else {
        out_ris <- ris(
          outcome = outcome,
          mc = mc,
          ma = list("metaPrepare" = mp, "synthesize" = sy),
          alpha = alpha,
          beta = beta,
          fixed = TRUE,
          side = side,
          type = "retrospective",
          ...
        )
      }
    } else {
      out_ris <- ris(
        outcome = outcome,
        mc = mc,
        alpha = alpha,
        beta = beta,
        fixed = FALSE,
        sd_mc = sd_mc,
        side = side,
        ma = list("metaPrepare" = mp, "synthesize" = sy), ...
      )
    }
    out <- append(out, list(ris = out_ris))
  }

  class(out) <- "metaanalysis"
  return(out)
}

# Define - metaprepare ----
metaPrepare <-
  function(data,
           outcome,
           weights,
           cont_vartype = "equal",
           alpha,
           conf_level,
           nonevent = NULL,
           zero_adj = 0.5) {
    
    # check | data
    if(sum(unlist(apply(is.na(data),2,which)))>0){
      stop("Missing data in data.frame. The function can not handle missing data.")
    }
    
    # store the original data in case of zero events handling
    org_data <- data

    #Prepare dichotomous outcomes.
    if (outcome %in% c("OR", "RR", "RD")) {
      # Remove studies with zero total events.
      if (sum(data$eI == 0 & data$eC == 0) > 0) {
        nonevent <- which(data$eI == 0 & data$eC == 0)
        data <- data[-nonevent, ]
      } else {
        nonevent <- NULL
      }

      # if some has all events
      if (sum(data$eI == data$nI & data$eC == data$nC) > 0 & outcome %in% c("OR", "RR")) {
        allevent <- which(data$eI == data$nI & data$eC == data$nC)
        data <- data[-allevent, ]
      } else {
        allevent <- NULL
      }
      
      # Adding 0.5 if one of the event counts is zero
      if (sum(data$eI == 0 | data$eC == 0) > 0) {
        zc <- which(data$eI == 0 | data$eC == 0)
        data$eI[zc] <- data$eI[zc] + zero_adj
        data$nI[zc] <- data$nI[zc] + 2*zero_adj
        data$eC[zc] <- data$eC[zc] + zero_adj
        data$nC[zc] <- data$nC[zc] + 2*zero_adj
      }

      # Calculate event probability
      pI <- data$eI / data$nI
      pC <- data$eC / data$nC

      # Set the number of trials
      K <- length(data$eI)

      # Prepare returned object
      out <- list()

      # Calculate effects sizes and belonging std.
      if (outcome == "RD") {
        # Risk difference
        te <- pI - pC
        sig <- sqrt(pI * (1 - pI) / data$nI + pC * (1 - pC) / data$nC)
      } else if (outcome == "OR") {
        # Odds Ratio
        te <- (pI / (1 - pI)) / (pC / (1 - pC))
        sig <-
          sqrt(1 / data$eI + 1 / (data$nI - data$eI) + 1 / data$eC + 1 / (data$nC -
                                                                            data$eC))
      } else if (outcome == "RR") {
        # Risk Ratio
        te <- pI / pC
        sig <- sqrt(1 / data$eI - 1 / data$nI + 1 / data$eC - 1 / data$nC)
      }

      # Calculate the weights
      if (weights == "IV") {
        # Inverse Variance
        w <- 1 / (sig ^ 2)
        w <- w / sum(w)
        pe <- sum(te * w)
        w <- w * 100
      } else if (weights == "MH") {
        # Mantel-Haenszel
        A <- data$eI
        B <- data$nI - data$eI
        C <- data$eC
        D <- data$nC - data$eC
        N <- A + B + C + D

        if (outcome == "OR") {
          w <- (data$nI - data$eI) * data$eC / N
          T1 <- (A + D) / N
          T2 <- (B + C) / N
          T3 <- A * D / N
          T4 <- B * C / N
          vpe <-
            0.5 * ((T1 * T3) / (sum(T3) ^ 2) + (T1 * T4 + T2 * T3) / (sum(T3) * sum(T4)) +
                     T2 * T4 / (sum(T4) ^ 2))
          svpe <- sum(vpe)
        } else if (outcome == "RR") {
          w <- (data$nI) * data$eC / N
          D1 <- ((A + B) * (C + D) * (A + C) - (A * C * N)) / (N ^ 2)
          R <- (A * (C + D)) / (N)
          S <- (C * (A + B)) / N
          svpe <- sum(D1) / (sum(R) * sum(S))
        } else if (outcome == "RD") {
          w <- (data$nI) * (data$nC) / N
          svpe <- sum(((A * B * data$nI) ^ 3 + (C * D * data$nC) ^ 3) /
                        (data$nI * data$nC * N) ^ 2) / (sum(data$nI * data$nC) /
                                                          N) ^ 2

        }
        pe <- sum(te * w) / sum(w)
      }

      # Calculate confidence limits
      if (outcome %in% c("RR", "OR")) {
        lower <- exp(log(te) - qnorm((1-conf_level)/2, lower.tail = FALSE) * sig)
        upper <- exp(log(te) + qnorm((1-conf_level)/2, lower.tail = FALSE) * sig)
      } else {
        lower <- te - qnorm((1-conf_level)/2, lower.tail = FALSE) * sig
        upper <- te + qnorm((1-conf_level)/2, lower.tail = FALSE) * sig
      }

      # Return results
      if (weights == "MH") {
        out <- list(
          w = w,
          te = te,
          lower = lower,
          upper = upper,
          pe = c(pe, svpe),
          sig = sig,
          outcome = outcome,
          weights = weights,
          data = data,
          org_data = org_data,
          nonevent = nonevent
        )
      } else{
        out <- list(
          w = w,
          te = te,
          lower = lower,
          upper = upper,
          pe = pe,
          sig = sig,
          outcome = outcome,
          weights = weights,
          data = data,
          org_data = org_data,
          nonevent = nonevent
        )
      }

    } else if (outcome == "MD") {
      # calculate effect sizes
      te <- data$mI - data$mC

            # calculate std.
      if (cont_vartype == "equal") {
        spooled <-
          sqrt(((data$nI - 1) * data$sdI ^ 2 + (data$nC - 1) * data$sdC ^ 2) / (data$nI +
                                                                                  data$nC - 2))
        vte <- (data$nI + data$nC) / (data$nI * data$nC) * spooled ^ 2
        sete <- sqrt(vte)
        df <- (data$nI - 1) + (data$nC - 1)
        lower <- te - qt((1-conf_level)/2, df = df, lower.tail = F) * sete
        upper <- te + qt((1-conf_level)/2, df = df, lower.tail = F) * sete
      } else if (cont_vartype == "non-equal") {
        vte <- data$sdI ^ 2 / data$nI + data$sdC ^ 2 / data$nC
        sete <- sqrt(vte)
        lower <- te - qnorm((1-conf_level)/2, lower.tail = F) * sete
        upper <- te + qnorm((1-conf_level)/2, lower.tail = F) * sete
      }

      # only weight option for MD is "IV"
      w <- 1 / vte
      pe <- sum(w * te) / sum(w)
      sig <- sqrt(vte)
      w <- w / sum(w)

      out <-
        list(
          w = w,
          te = te,
          lower = lower,
          upper = upper,
          pe = pe,
          sig = sig,
          outcome = outcome,
          data = data,
          org_data = org_data,
          weights = weights,
          nonevent = nonevent
        )
    }

    return(out)
  }

# Define - synthesize ----
synthesize <- function(y, re_method, tau_ci_method, conf_level) {
  # collect objects
  w <- y$w   # collect objects
  sig <- y$sig
  te <- y$te
  pe <- y$pe
  df <- length(w) - 1
  ci.tau <- ""

  if (y$weights != "MH") w <- 1 / (sig ^ 2) # weight inverse variance
  rw <- w / sum(w) # relative weight
  vw <- 1 / sum(w) # variance of pooled effect

  if (y$weights == "MH") {
    rw <- w / sum(w)
    if (y$outcome == "RD") {
      vw <- pe[2]
      peF <- sum(te * w) / sum(w)
      lci <- peF - qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vw)
      uci <- peF + qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vw)
      zval <- peF / sqrt(vw)
    } else {
      # weights the same for OR and RR
      vw <- pe[2]
      lpeF <- log(sum(te * w) / sum(w))
      lci <- exp(lpeF - qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vw))
      uci <- exp(lpeF + qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vw))
      peF <- exp(lpeF)
      zval <- lpeF / sqrt(vw)
    }

    pval <- (1 - pnorm(abs(zval))) * 2

  } else {
    # if IV
    if (y$outcome %in% c("RD", "MD")) {
      peF <- sum(te * rw) # fixed effect log pooled estimate
      uci <- peF + qnorm((1-conf_level)/2, lower.tail = F) * sqrt(vw)
      lci <- peF - qnorm((1-conf_level)/2, lower.tail = F) * sqrt(vw)
      zval <- peF / sqrt(vw)
    } else {
      lpeF <- sum(log(te) * rw) # fixed effect log pooled estimate
      uci <- exp(lpeF + qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vw))
      lci <- exp(lpeF - qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vw))
      peF <- exp(lpeF)
      zval <- lpeF / sqrt(vw)
    }
    pval <- (1 - pnorm(abs(zval))) * 2
  }
  
  if(df == 0){
    Q <- 0; U <- 0; tau2 <- 0; H <- 0; I2 <- 0; D2 <- 0; pQ <- 1; rwR <- NULL
    
    if(!(y$outcome %in% c("RR", "OR"))){
      lpeF <- NA
    } else {
      lpeF <- log(peF)
    }
    
    synth <-
      list(
        fw = round(rw / sum(rw) * 100, 4),
        peF = c(peF, lci, uci, zval, pval, lpeF, vw),
        Q = c(Q, df, pQ),
        U = c(tau2, H, I2, D2)
      )
    class(synth) <- "synthesized"
    return(synth)
    
  } else {
  # Heterogeneity estimates and statistics
  if (y$outcome != "MD") { # for binary data
    w <- 1 / (sig ^ 2)
    if (y$outcome == "RD") {
      Q <- sum(w * te ^ 2) - (sum(w * te)) ^ 2 / sum(w)
    } else {
      Q <- sum(w * log(te) ^ 2) - (sum(w * log(te))) ^ 2 / sum(w)
    }
    U <- sum(w) - sum(w ^ 2) / sum(w)
    tau2 <-
      ifelse(Q > df, (Q - df) / U, 0) # DerSimonian-Laird estimate
    pQ <- pchisq(Q, df, lower.tail = FALSE)
    if (!is.na(Q) & Q / df <= 0) {
      H <- 0
    } else if (is.na(Q)) {
      H <- NA
    } else{
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
      if (y$outcome != "RD") {
        leR <- sum(log(te) * rwR)
        peR <- exp(leR)
      }
      if (re_method == "DL_HKSJ") {
        if (y$outcome == "RD") {
          vwR <- 1 / df * sum(wR * (te - peR) ^ 2 / sum(wR))
          zvalR <- peR / sqrt(vwR)
          pvalR <- 2 * pt(abs(zvalR), df = df, lower.tail = FALSE)
          lciR <- peR - qt((1 - conf_level) / 2, df = df, lower.tail = FALSE) * sqrt(vwR)
          uciR <- peR + qt((1 - conf_level) / 2, df = df, lower.tail = FALSE) * sqrt(vwR)
        } else {
          vwR <- 1 / df * sum(wR * (log(te) - leR) ^ 2 / sum(wR))
          zvalR <- leR / sqrt(vwR)
          pvalR <- 2 * pt(abs(zvalR), df = df, lower.tail = FALSE)
          lciR <- exp(leR - qt((1 - conf_level) / 2, df = df, lower.tail = FALSE) * sqrt(vwR))
          uciR <- exp(leR + qt((1 - conf_level) / 2, df = df, lower.tail = FALSE) * sqrt(vwR))
        }
      } else {
        if (y$outcome == "RD") {
          zvalR <- peR / sqrt(vwR)
          pvalR <- (1 - pnorm(abs(zvalR))) * 2
          lciR <- peR - qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vwR)
          uciR <- peR + qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vwR)
        } else {
          zvalR <- leR / sqrt(vwR)
          pvalR <- (1 - pnorm(abs(zvalR))) * 2
          lciR <- exp(leR - qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vwR))
          uciR <- exp(leR + qnorm((1-conf_level)/2, lower.tail = FALSE) * sqrt(vwR))
        }
      }

      vw <- 1 / sum(w)
      D2 <- 1 - vw / vwR

      if (y$outcome %in% c("RR", "OR") & tau_ci_method == "BJ") {
        ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
          yi = log(te),
          sei = sig,
          method = "GENQ",
          weights = 1 / sig ^ 2
        ))
      }

      if (y$outcome %in% c("RR", "OR") & tau_ci_method == "QP") {
        ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
          yi = log(te),
          sei = sig,
          method = "DL"
        ))
      }

      if (y$outcome == c("RD") & tau_ci_method == "BJ") {
        ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
          yi = te,
          sei = sig,
          method = "GENQ",
          weights = 1 / sig ^ 2
        ))
      }

      if (y$outcome == c("RD") & tau_ci_method == "QP") {
        ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
          yi = te,
          sei = sig,
          method = "DL"
        ))
      }

      if(y$outcome != "RD"){
      synth <-
        list(
          fw = round(rw / sum(rw) * 100, 4),
          peF = c(peF, lci, uci, zval, pval, log(peF), vw),
          rwR = rwR * 100,
          peR = c(peR, lciR, uciR, zvalR, pvalR, vwR),
          Q = c(Q, df, pQ),
          U = c(tau2, H, I2, D2),
          ci.tau = ci.tau
        )
      class(synth) <- "synthesized" } else {
        synth <-
          list(
            fw = round(rw / sum(rw) * 100, 4),
            peF = c(peF, lci, uci, zval, pval, NA, vw),
            rwR = rwR * 100,
            peR = c(peR, lciR, uciR, zvalR, pvalR, vwR),
            Q = c(Q, df, pQ),
            U = c(tau2, H, I2, D2),
            ci.tau = ci.tau
          )
      }

      return(synth)
    } else {
      D2 <- 0
      synth <-
        list(
          fw = round(rw / sum(rw) * 100, 4),
          peF = c(peF, lci, uci, zval, pval, log(peF), vw),
          Q = c(Q, df, pQ),
          U = c(tau2, H, I2, D2)
        )
      class(synth) <- "synthesized"
      return(synth)
    }
  } else {
    # heterogeneity estimate and statistics for MD
    Q <- sum(w * te ^ 2) - (sum(w * te)) ^ 2 / sum(w)
    U <- sum(w) - sum(w ^ 2) / sum(w)
    tau2 <- ifelse(Q > df, (Q - df) / U, 0)
    wR <- 1 / (sig ^ 2 + tau2)
    vwR <- 1 / sum(wR) # variance of pooled effect (random)
    rwR <- wR * vwR
    peRest <- sum(te * wR) / sum(wR)
    lciR <- peRest - qnorm((1-conf_level)/2, lower.tail = F) * sqrt(vwR)
    uciR <- peRest + qnorm((1-conf_level)/2, lower.tail = F) * sqrt(vwR)
    zvalR <- peRest / sqrt(vwR)
    pvalR <- 2 * (1 - pnorm(abs(zvalR)))
    pQ <- pchisq(Q, df, lower.tail = FALSE)
    H <- sqrt(Q / df)
    I2 <- ifelse((Q - df) / Q >= 0, (Q - df) / Q, 0)
    D2 <- ifelse(is.na(1 - vw / vwR), 0, (1 - vw / vwR))
    
    if (tau_ci_method == "BJ") {
      ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
        yi = te,
        sei = sig,
        method = "GENQ",
        weights = 1 / sig ^ 2
      ))
    }

    if (tau_ci_method == "QP") {
      ci.tau <- metafor::confint.rma.uni(metafor::rma.uni(
        yi = te,
        sei = sig,
        method = "DL"
      ))
    }
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

#' @method print metaanalysis
#' @export
print.metaanalysis <- function(x, ...) {
  cat("Individual trial results: \n \n")
  y <- x$study_results
  y[,2:5] <- round(y[,2:5],3)
  y[,6:7] <- round(y[,6:7],2)
  print(y)
  cat("\nNon-sequential meta-analysis results: \n \n")
  y <- x$meta_results
  y[,2:5] <- round(y[,2:5], 3)
  y[,6] <- round(y[,6], 4)
  print(y)
  invisible(x)
  if (!is.null(x$settings$missing_vec)) {
    message("\n NB. Please provide a study vector to name studies.")
  }
  #Total-zero trials
  if (!is.null(x$settings$nonevent)) {
    message(
      paste(
        "\n NB.",
        x$settings$nonevent,
        "was excluded from the analysis due to zero total events, consider changing outcome to RD."
      )
    )
  }

}

#' Forestplot for metaanalysis object.
#' 
#' @param x metaanalysis object from the RTSA package.
#' @param type Define whether or not both fixed-effect and random-effects meta-analysis results should be printed on the plot. Options are: "fixed", "random" or "both". Default is "both".
#' @param xlims Set default limits on the outcome scale. Default is NULL.
#' @param ... Additional arguments
#'
#' @method plot metaanalysis
#' @importFrom ggplot2 annotate labs scale_colour_identity geom_errorbar aes_string geom_polygon arrow unit
#' @importFrom stats na.omit complete.cases
#' @export
#' 
#' @examples
#' # Example with OR 
#' ma <- metaanalysis(data = coronary, outcome = "OR")
#' plot(ma)
#' 
#' # Example with RR 
#' ma <- metaanalysis(data = perioOxy, outcome = "RR")
#' plot(ma)
#' 
#' # Example with MD
#' ma <- metaanalysis(data = eds, outcome = "MD")
#' plot(ma, type = "random")
plot.metaanalysis <- function(x, type = "both", xlims = NULL, ...) {

  plot_message <- NULL

  # Create dataframe for plot
  fplot <- merge(x$metaPrepare$org_data, x$study_results, sort = FALSE)
  CInames <- colnames(fplot[grepl("CI", colnames(fplot))])
  results <- x$meta_results
  colnames(results)[1] <- "study"
  results[, colnames(fplot)[!(colnames(fplot) %in% colnames(results))]] <-
    NA
  results$study[results$study == "Fixed"] <- "Fixed-effect"
  results$study[results$study == "Random"] <- "Random-effects"
  fplot <- rbind(fplot, results[names(fplot)])

  # Ensure fixed effect model, if Tau = 0
  if (x$synthesize$U[1] == 0) {
    type <- "fixed"
    plot_message <-
      "No heterogeneity estimated, a fixed-effect meta-analysis is presented."
  }

  if (type == "fixed") {
    fplot <- fplot[!grepl("Random-ef", fplot$study), ]
  }
  if (type == "random") {
    fplot <- fplot[!grepl("Fixed-ef", fplot$study), ]
  }

  fplot <- cbind(nrow(fplot):1, fplot)
  colnames(fplot)[1] <- "yaxis"

  if(x$settings$outcome %in% c("RR","OR","RD")){
  fplot$nI[grepl("Fixed-ef|Random-ef", fplot$study)] <-
    sum(fplot$nI, na.rm = T)
  fplot$nC[grepl("Fixed-ef|Random-ef", fplot$study)] <-
    sum(fplot$nC, na.rm = T)
  fplot$eI[grepl("Fixed-ef|Random-ef", fplot$study)] <-
    sum(fplot$eI, na.rm = T)
  fplot$eC[grepl("Fixed-ef|Random-ef", fplot$study)] <-
    sum(fplot$eC, na.rm = T)
  } 

  outcome <- colnames(x$study_results)[2]
  colnames(fplot)[which(colnames(fplot) == outcome)] <- "outcome"

  fplot$out_ci <- paste0(
    sprintf(fplot$outcome, fmt = '%#.2f'),
    " (",
    sprintf(fplot[,CInames[1]], fmt = '%#.2f'),
    "; ",
    sprintf(fplot[,CInames[2]], fmt = '%#.2f'),
    ")"
  )
  
  if(x$settings$outcome %in% c("OR", "RR", "RD")){
  fplot$controls <- paste0(fplot$eC, "/", fplot$nC)
  fplot$controls[fplot$controls == "NA/NA"] <- ""
  fplot$experimental <- paste0(fplot$eI, "/", fplot$nI)
  fplot$experimental[fplot$experimental == "NA/NA"] <- ""
  } else {
    fplot$controls <- paste0(round(fplot$mC,1), "/",round(fplot$sdC,1), "/", fplot$nC)
    fplot$controls[fplot$controls == "NA/NA/NA"] <- ""
    fplot$experimental <- paste0(round(fplot$mI,2), "/",round(fplot$sdI,2), "/", fplot$nI)
    fplot$experimental[fplot$experimental == "NA/NA/NA"] <- ""
  }
  fplot$wF <- round(fplot$w_fixed)
  fplot$wF[fplot$wF < 10 &
             !is.na(fplot$wF)] <-
    round(fplot$w_fixed[fplot$wF < 10 & !is.na(fplot$wF)], 1)
  fplot$wR <- round(fplot$w_random)
  fplot$wR[fplot$wR < 10 &
             !is.na(fplot$wR)] <-
    round(fplot$w_random[fplot$wR < 10 & !is.na(fplot$wR)], 1)

  #Define shapes and colors
  shapes <- 21
  if(type == "fixed"){ sizes <- pmax(fplot$w_fixed/4,1)
  } else {
    sizes <-pmax(fplot$w_random/4,1)
  }
  colors <- rep("black", nrow(fplot))
  
  # define x-axis
  if (is.null(xlims))
    xlims <- c(min(fplot[,CInames[1]], na.rm = T), max(fplot[,CInames[2]], na.rm = T))
  if(x$settings$outcome %in% c("OR","RR")){
    if(xlims[1] < 1 & xlims[2] < 1){
      xlims[2] <- 1.05
    } else if(xlims[1] > 1 & xlims[2] > 1){
      xlims[1] <- 0.95
    }
  } else {
    if(xlims[1] < 0 & xlims[2] < 0){
      xlims[2] <- 0
    } else if(xlims[1] > 0 & xlims[2] > 0){
      xlims[1] <- 0
    }
  }
  
  midpoint <- ifelse(x$settings$outcome %in% c("OR","RR"), 1, 0)
  
  n <- length(na.omit(fplot[,CInames[2]]))
  xmax2 <- sort(na.omit(fplot[,CInames[2]]))[n-1]
  xmin2 <- sort(na.omit(fplot[,CInames[1]]), decreasing = T)[n-1]
  arrow_data <- data.frame(x = NA, y = NA, xend = NA, yend = NA)
  
  if(midpoint != xlims[1] & (midpoint - xlims[1] > (midpoint - xmin2)*2)){
    xlims[1] <- xmin2 - abs((xlims[1]-xmin2)/3)
    arrow_data <- rbind(arrow_data, c(xmax2,
                                      fplot$yaxis[fplot[,CInames[1]] < xlims[1]],
                        xlims[1],
                        fplot$yaxis[fplot[,CInames[1]] < xlims[1]]))
    fplot[,CInames[1]][fplot[,CInames[1]] < xlims[1]] <- xlims[1]
  }
  if(midpoint != xlims[2] & (xlims[2]-midpoint) > (xmax2-midpoint)*2){
    xlims[2] <- xmax2 + abs((xlims[2]-xmax2)/3)
    arrow_data <- rbind(arrow_data, c(xmin2,
                                      fplot$yaxis[fplot[,CInames[2]] > xlims[2]],
                                      xlims[2],
                                      fplot$yaxis[fplot[,CInames[2]] > xlims[2]]))
    fplot[,CInames[2]][fplot[,CInames[2]] > xlims[2]] <- xlims[2]
  }
  arrow_data <- arrow_data[complete.cases(arrow_data),]
  
  # set the ratio between the forest plot and the other columns in plot
  # we want the forest plot to be at least 30% of the plot
  # start by calculating the wanted range of the forest plot
  if(x$settings$outcome %in% c("OR", "RR")){
    ci_range <- xlims[2]-xlims[1]
     xl <- list(
       `l0` = exp(log(xlims[1])-ci_range*3/5),
       `l1` = exp(log(xlims[1])-ci_range*2/5),
       `l2` = exp(log(xlims[1])-ci_range*1/5),
       `r1` = exp(log(xlims[2])+ci_range*2/5),
       `r2` = exp(log(xlims[2])+ci_range*3/5),
       `r3` = exp(log(xlims[2])+ci_range*4/5)
  )
  } else {
    int <- abs(xlims[1]-xlims[2])
    xl <- list(
      `l0` = xlims[1] - int,
      `l1` = xlims[1] - int/2,
      `l2` = xlims[1] - int/10,
      `r1` = xlims[2] + int,
      `r2` = xlims[2] + int*1.3,
      `r3` = xlims[2] + int*1.6
    )
  }
  
  if (type == "random") {
    xl$r3 <- xl$r2
    xl$r2 <- NULL
  }
  if (type == "fixed")
    xl$r3 <- NULL

  xl2 <- c(unlist(unname(xl[-1])), mean(c(xl$r2, xl$r3), na.rm = T))
  xl2 <- xl2[order(xl2)]

  #Axis text
  xlabs <- c(
    ifelse(x$settings$outcome == "MD", "Experimental\n mean/sd/n","Experimental\n events/n"),
    ifelse(x$settings$outcome == "MD", "Control\n mean/sd/n","Control\n events/n"),
    paste0("", outcome, " (95% CI)"),
    if (type != "random") {
      "Fixed"
    },
    "Weights (%)\n",
    if (type != "fixed") {
      "Random"
    }
  )

  #Lower left text
  heterogen <- paste0(
    "Heterogeneity: tau^2 = ",
    sprintf(x$synthesize$U[1], fmt = '%#.2f'),
    " (",
    sprintf(x$synthesize$ci.tau$random["tau^2", "ci.lb"], fmt = '%#.2f'),
    "; ",
    sprintf(x$synthesize$ci.tau$random["tau^2", "ci.ub"], fmt = '%#.2f'),
    ")",
    ", Q = ",
    sprintf(x$synthesize$Q[1], fmt = '%#.1f'),
    ", df = ",
    round(x$synthesize$Q[2]),
    ", I^2 = ",
    round(x$synthesize$U[3] * 100),
    "%\n"
  )
  if (type != "random") {
    heterogen <- paste0(
      heterogen,
      "Test (Fixed-effect): z = ",
      sprintf(x$synthesize$peF[4], fmt = '%#.2f'),
      " (p = ",
      sprintf(x$meta_results[x$meta_results == "Fixed", "pValue"], fmt = '%#.4f'),
      ")"
    )
  }
  if (type != "fixed") {
    heterogen <- paste0(
      heterogen,
      "; Random-effects, z = ",
      sprintf(x$synthesize$peR[4], fmt = '%#.2f'),
      " (p = ",
      sprintf(x$meta_results[x$meta_results == "Random", "pValue"], fmt = '%#.4f'),
      ")"
    )
  }

  lfp <- dim(fplot)[1]
  if(type == "both"){
    rm <- c(lfp-1,lfp)
    poly_data <- data.frame(x_fixed = c(fplot[rm[1],CInames[1]],fplot$outcome[rm[1]],
                                        fplot[rm[1],CInames[2]],fplot$outcome[rm[1]]),
                            y_fixed = c(fplot$yaxis[rm[1]],fplot$yaxis[rm[1]]+0.25,
                                        fplot$yaxis[rm[1]],fplot$yaxis[rm[1]]-0.25),
                            x_random = c(fplot[rm[2],CInames[1]],fplot$outcome[rm[2]],
                                        fplot[rm[2],CInames[2]],fplot$outcome[rm[2]]),
                            y_random = c(fplot$yaxis[rm[2]],fplot$yaxis[rm[2]]+0.25,
                                        fplot$yaxis[rm[2]],fplot$yaxis[rm[2]]-0.25))
  } else {
    rm <- c(lfp)
    if(type == "fixed"){
      poly_data <- data.frame(x_fixed = c(fplot[rm,CInames[1]],fplot$outcome[rm],
                                          fplot[rm,CInames[2]],fplot$outcome[rm]),
                              y_fixed = c(fplot$yaxis[rm],fplot$yaxis[rm]+0.25,
                                          fplot$yaxis[rm],fplot$yaxis[rm]-0.25))
    } else {
      poly_data <- data.frame(x_random = c(fplot[rm,CInames[1]],fplot$outcome[rm],
                                           fplot[rm,CInames[2]],fplot$outcome[rm]),
                              y_random = c(fplot$yaxis[rm],fplot$yaxis[rm]+0.25,
                                           fplot$yaxis[rm],fplot$yaxis[rm]-0.25))
    }
  }
  
  if(x$settings$outcome %in% c("RD", "MD")){
    if(abs(xlims[2]-xlims[1]) > 5){
      break_xlim <- c(round(seq(from = xlims[1], to = xlims[2], length.out = 5),0),0)
    } else {
      break_xlim <- c(round(seq(from = xlims[1], to = xlims[2], length.out = 5),2),0)
    }
  } else {
    break_xlim <- c(round(exp(seq(from = log(xlims[1]), to = log(1), length.out = 3)),2),
                    round(exp(seq(from = log(1), to = log(xlims[2]), length.out = 3)),2)[-1])
    if(sum(diff(log(break_xlim)) < abs(log(ci_range))/5) > 0){
      break_xlim <- c(break_xlim[1],break_xlim[-1][which(diff(log(break_xlim)) > abs(log(ci_range))/5)])
    }
  }


  #The plot
  p <-
    ggplot(fplot, aes(
      x = c(outcome[-rm], rep(NA, length(rm))),
      xmin = c(fplot[-rm,CInames[1]], rep(NA, length(rm))),
      xmax = c(fplot[-rm,CInames[2]], rep(NA, length(rm))),
      y = yaxis
    )) +
    geom_vline(xintercept = ifelse(x$settings$outcome %in% c("RR", "OR"), 1, 0),
               color = "gray",
               linetype = 1) +
    geom_segment(aes(
      x = xlims[1],
      xend = xlims[2],
      y = -Inf,
      yend = -Inf
    )) +
    annotate(
      "text",
      x = xl$l1,
      y = fplot$yaxis,
      label = fplot$experimental,
      hjust = 1
    ) +
    annotate(
      "text",
      x = xl$l2,
      y = fplot$yaxis,
      label = fplot$control,
      hjust = 1
    )+
    annotate(
      "text",
      x = xl$r1,
      y = fplot$yaxis,
      label = fplot$out_ci,
      hjust = 1
    )

  # Vertical estimate lines and polygon
  if(type != "fixed"){
    p <- p + geom_segment(aes(
      x = outcome[grepl("Random-effect", study)],
      xend = outcome[grepl("Random-effect", study)],
      y = Inf,
      yend = yaxis[grepl("Random-effect", study)]
    ),
    color = "blue", linetype = 2, linewidth = 0.25) + annotate(
      "text",
      x = xl$r3,
      y = fplot$yaxis,
      label = fplot$wR,
      hjust = 1
    ) + geom_polygon(data = poly_data,inherit.aes = FALSE,
                     mapping = aes(x = x_random, y = y_random), fill = "blue",
                     col = "black", alpha = 0.5)
  }
   if (type != "random"){
     p <- p +
        geom_segment(aes(
          x = outcome[grepl("Fixed-effect", study)],
          xend = outcome[grepl("Fixed-effect", study)],
          y = Inf,
          yend = yaxis[grepl("Fixed-effect", study)]
        ),
        color = "red", linetype = 2, linewidth = 0.25) + annotate(
          "text",
          x = xl$r2,
          y = fplot$yaxis,
          label = fplot$wF,
          hjust = 1
        ) + geom_polygon(data = poly_data,inherit.aes = FALSE,
                         mapping = aes(x = x_fixed, y = y_fixed), fill = "red",
                         col = "black", alpha = 0.5)
    }
  
    # add dots between study results and meta-analysis results
    p <- p + geom_segment(aes(
      x = ifelse(x$settings$outcome %in% c("RR", "OR"),0,-Inf),
      xend = xl$l2,
      y = max(fplot$yaxis[grepl("Fixed-ef|Random-ef", fplot$study)]) +
        0.5,
      yend = max(fplot$yaxis[grepl("Fixed-ef|Random-ef", fplot$study)]) +
        0.5
    ), linetype = 3) +
    geom_segment(aes(
      x = xlims[2],
      xend = Inf,
      y = max(fplot$yaxis[grepl("Fixed-ef|Random-ef", fplot$study)]) +
        0.5,
      yend = max(fplot$yaxis[grepl("Fixed-ef|Random-ef", fplot$study)]) +
        0.5
    ), linetype = 3) +
    labs(tag = heterogen) +
    geom_point(
      shape = shapes,
      color = colors,
      fill = colors,
      size = sizes
    ) +
    geom_errorbar(color = colors, width = 0) +
      geom_segment(data = arrow_data, inherit.aes = F,
                           aes(x = x, y = y, xend = xend,
                                          yend = yend),
                   arrow = arrow(length = unit(0.13, "cm"))) +
    theme_classic() +
    scale_colour_identity() + {if(x$settings$outcome %in% c("OR", "RR"))
      scale_x_continuous(
        trans = "log10",
        sec.axis = sec_axis( ~ ., breaks = xl2, labels = xlabs),
        limits = c(xl$l0, if (type != "fixed") {
          xl$r3
        } else{
          xl$r2
        }),
        name = "o\no",
        breaks = break_xlim
      )
      } +  {if(x$settings$outcome %in% c("RD", "MD"))
        scale_x_continuous(
          sec.axis = sec_axis( ~ ., breaks = xl2, labels = xlabs),
          limits = c(xl$l0, if (type != "fixed") {
            xl$r3
          } else{
            xl$r2
          }),
          name = "o\no",
          breaks = break_xlim
        )
      } +
    scale_y_continuous(breaks = fplot$yaxis, labels = fplot$study) +
    theme(
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.line.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.title.x = element_text(color = "white"),
      axis.text.x.top = element_text(hjust = 1, color = "black"),
      axis.text.y = element_text(
        color = "black",
        size = 10,
        face = "bold"
      ),
      plot.tag.position = c(0, 0),
      plot.tag = element_text(
        hjust = 0,
        vjust = 0,
        size = 9
      )
    )

  return(suppressWarnings(print(p)))

  if (!is.null(plot_message)) {
    message(plot_message)
  }

}

if(getRversion() >= "2.15.1"){
  utils::globalVariables(c("yaxis", "study", "x_random", "y_random", "x_fixed",
                           "y_fixed", "y", "xend", "yend"))
}
