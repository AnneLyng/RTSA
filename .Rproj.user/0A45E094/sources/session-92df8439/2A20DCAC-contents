#' RTSA
#'
#' @param type Type of RTSA. Options are "design" or "analysis".
#' @param data A data.frame containing the study results. The data set must containing a specific set of columns. These are respectively `eI` (events in intervention group), `eC` (events in control group), `nC` (participants intervention group) or `nI` (participants control group) for discrete data, or, `mI` (mean intervention group), `mC` (mean control group), `sdI` (standard error intervention group), `sdC` (standard error control group),`nC` (participants intervention group) and `nI` (participants control group)  for continuous outcomes. Preferable also a `study` column as an indicator of study.
#' @param study An optional vector of study names and perhaps year of study. Defaults to NULL.
#' @param ana_time An optional vector of analysis times. Used if the sequential analysis is not done for all studies included in the meta-analysis.
#' @param timing Expected timings of interim analyses when type = "design". Defaults to NULL.
#' @param side Whether a 1- or 2-sided hypothesis test is used. Options are 1 or 2.
#' @param outcome Outcome metric. Options are: RR (risk ratio/relative risk), OR (odds ratio), RD (risk difference) and MD (mean difference).
#' @param mc Minimal clinical relevant outcome value
#' @param sd_mc The expected standard deviation. Used for sample size calculation for mean differences.
#' @param p0 The expected probability of event in the control group. Used for sample size calculation for binary outcomes.
#' @param alpha The level of type I error
#' @param beta The level of type II error
#' @param futility Futility boundaries added to design. Options are: none, non-binding and binding. Default is "none".
#' @param fixed Should only a fixed-effect meta-analysis be computed. Default is FALSE.
#' @param tau2 Heterogeneity estimate. Used for sample and trial size calculation. Defaults to NULL.
#' @param I2 Inconsistency estimate. Used for sample and trial size calculation. Defaults to NULL.
#' @param D2 Diversity estimate. Used for sample and trial size calculation. Defaults to NULL.
#' @param weights Weighting method options include IV (inverse-variance) and MH (Mantel-Haenszel). Defaults to IV.
#' @param cont_vartype For mean difference outcomes, do we expect the variance in the different groups to be "equal" or "non-equal".
#' @param re_method Method for calculating the estimate of heterogeneity, tau^2, and the random-effects meta-analysis variance. Options are "DL" for DerSimonian-Laird and "HKSJ" for the Harting-Knapp-Sidik-Jonkman adjustment of the DerSimonian-Laird estimator.
#' @param tau_ci_method Method for calculating confidence intervals for the estimated heterogeneity tau^2. Options are "QP" for Q-profiling and "BJ" for Biggelstaff ....
#' @param es_alpha The spending function for alpha-spending. Options are: esOF (Lan & DeMets version of O'Brien-Fleming), esPoc (Lan & DeMets version of Pocock), HSDC (Hwang Sihi and DeCani) and rho (rho family).
#' @param es_beta The spending function for beta-spending. For options see es_alpha.
#' @param gamma Parameter for the HSDC error spending function.
#' @param rho Parameter for the rho family error spending function.
#' @param design RTSA object where type is design.
#' @param zero_adj Zero adjustment. Options for now is 0.5.
#' @param conf_int Stopping time confidence interval. Options for now is sw (stage-wise).
#' @param conf_level Confidence level on stopping time confidence interval.
#' @param ... other arguments
#'
#' @returns A RTSA object, a list of five elements:
#' \item{settings}{A list containing all of the settings used in the \code{RTSA} call. See Arguments.}
#' \item{ris}{List containing sample and trial size calculations. See documentation for \code{ris} function.}
#' \item{bounds}{List of stopping boundaries, timing of trials and more. See documentation for \code{boundaries} function.}
#' \item{results}{List of 3 to 7 elements. \code{DARIS} diversity adjusted required information size. \code{DARIS_F} fixed-effect meta-analysis required sample size. \code{AIS} Achieved information size. \code{results_df} a data.frame of inference, see documentation for \code{inference} function.  \code{seq_inf} a list of conditional inference, see documentation for \code{inference} function. \code{metaanalysis} A metaanalysis object, see documentation for \code{metaanalysis} function. \code{design_df} a data.frame containing the stopping boundaries and timings from the design.}
#' \item{warnings}{List of warnings}
#'
#'
#'
#' @export
#' @aliases print.RTSA
#'
#' @examples
#' ### Retrospective sequential meta-analysis:
#' # A RRR of 20% is expected which gives mc = 1 - RRR = 0.8. 
#' # No futility boundaries
#' data(perioOxy)
#' RTSA(type = "analysis", data = perioOxy, outcome = "RR", mc = 0.8, side = 2,
#'  alpha = 0.05, beta = 0.2, es_alpha = "esOF")
#'
#' # Set binding futility boundaries
#' # And use Lan and DeMets' version of Pocock stopping boundaries 
#' RTSA(type = "analysis", data = perioOxy, outcome = "RR", mc = 0.8, side = 2,
#'  alpha = 0.05, beta = 0.2, es_alpha = "esOF", futility = "binding",
#'  es_beta = "esPoc")
#'  
#' # Set non-binding futility boundaries
#' RTSA(type = "analysis", data = perioOxy, outcome = "RR", mc = 0.8, side = 2,
#'  alpha = 0.05, beta = 0.2, es_alpha = "esOF", futility = "non-binding",
#'  es_beta = "esPoc")
#'  
#' ### Design a prospective sequential meta-analysis
#' # For continuous data without expected heterogeneity
#' RTSA(type = "design", outcome = "MD", mc = 5, sd_mc = 10, side = 1, 
#' timing = c(0.33, 0.66, 1), fixed = TRUE,
#' alpha = 0.025, beta = 0.1, es_alpha = "esOF", futility = "non-binding", 
#' es_beta = "esPoc")
#' 
#' # For binary outcome
#' RTSA(type = "design", outcome = "RR", mc = 0.75, side = 1, 
#' timing = c(0.33, 0.66, 1), p0 = 0.1, D2 = 0.1, 
#' alpha = 0.025, beta = 0.2, es_alpha = "esOF", futility = "non-binding", 
#' es_beta = "esOF")
#' 
#' # extract sample size calculation
#' out_rtsa <-  RTSA(type = "design", outcome = "RR", mc = 0.75, side = 1, 
#' timing = c(0.33, 0.66, 1), p0 = 0.1, D2 = 0.1, 
#' alpha = 0.025, beta = 0.2, es_alpha = "esOF", futility = "non-binding", 
#' es_beta = "esOF")
#' out_rtsa$ris
#' 
#' # plot the design
#' plot(out_rtsa)
#' 
#' # update the design with data as it accumulates (here toy-data)
#' fake_data <- data.frame(eI = c(10,10), eC = c(13, 11), nI = c(750, 750),
#' nC = c(750,750))
#' RTSA(type = "analysis", design = out_rtsa, data = fake_data)
#' 
#' # plot the analysis
#' an_rtsa <- RTSA(type = "analysis", design = out_rtsa, data = fake_data)
#' plot(an_rtsa)


RTSA <-
  function(type = "design",
           data = NULL,
           study = NULL,
           ana_time = NULL,
           timing = NULL,
           side = NULL,
           outcome = NULL,
           mc,
           sd_mc = NULL,
           p0 = NULL,
           alpha = NULL,
           beta = NULL,
           zero_adj = 0.5,
           futility = "none",
           fixed = FALSE,
           tau2 = NULL,
           I2 = NULL,
           D2 = NULL,
           weights = "IV",
           cont_vartype = "equal",
           re_method = "DL",
           tau_ci_method = "BJ",
           es_alpha = NULL,
           es_beta = NULL,
           gamma = NULL,
           rho = NULL,
           design = NULL,
           conf_int = "sw",
           conf_level = 0.95,
           ...) {
    # Check inputs ----
    # check | type
    if(type == "analysis" & is.null(data)){
      stop("Provide `data` for analysis.")
    }

    # check | design arguments equal to analysis arguments if design is present
    if(type == "analysis" & !is.null(design)){
      if(!is.null(es_alpha)) {
        if (es_alpha != design$settings$es_alpha) {
          warning(
            paste0(
              "`es_alpha` it not equal to the design setting of es_alpha. The design setting overrules the error spending function.",
              " es_alpha is set to ",
              design$setting$es_alpha
            )
          )
        }
      }
      if(!is.null(es_beta)){
        if(es_beta != design$settings$es_beta){
        warning(paste0("`es_beta` it not equal to the design setting of es_beta. The design setting overrules the error spending function.",
                       " es_beta is set to ", design$setting$es_beta))
        }
      }
      if(!is.null(alpha)){
        if(alpha != design$settings$alpha){
        warning(paste0("`alpha` it not equal to the design setting of alpha. The design setting overrules.",
                       " alpha is set to ", design$setting$alpha))
        }
      }
      if(!is.null(beta)){
        if(beta != design$settings$beta){
        warning(paste0("`beta` it not equal to the design setting of beta. The design setting overrules.",
                       " beta is set to ", design$setting$beta))
        }
      }
      if(!is.null(outcome)){
        if(outcome != design$settings$outcome){
        warning(paste0("`outcome` it not equal to the design setting of outcome. The design setting overrules.",
                       " outcome is set to ", design$setting$outcome))
        }
      }
      if(!is.null(side)){
        if(side != design$settings$side){
        warning(paste0("`side` it not equal to the design setting of side. The design setting overrules.",
                       " side is set to ", design$setting$side))
        }
      }
      if(futility != "none" & futility != design$settings$futility){
        warning(paste0("`futility` it not equal to the design setting of futility. The design setting overrules.",
                       " futility is set to ", design$setting$futility))
      }
    }

    # check | timing
    if(type == "design" & is.null(timing)){
      stop("Provide expected `timing` for design.")
    }

    # if design is present
    if(!is.null(design)){
      outcome <- design$settings$outcome
      alpha <- design$settings$alpha
      beta <- design$settings$beta
      es_alpha <- design$settings$es_alpha
      es_beta <- design$settings$es_beta
      futility <- design$settings$futility
      mc <- design$settings$mc
    }
    
    # check | mc
    if(outcome %in% c("OR", "RR") & mc < 0){
      stop("For outcomes such as risk ratios (RR) or odds ratio (OR), the minimum clinical value (mc) can not be less than 0.")
    }

    # check | outcome
    if(!(outcome %in% c("MD", "RD", "RR", "OR"))) {
      stop("`outcome` must be 'MD', 'RD', 'RR' or 'OR'.")
    }
    # check | data; mI; mC; sdI; sdC; eI; nI; eC; nC.
    if (!is.null(data)) {
      if (outcome == "MD" &
          !all(c("mI", "mC", "sdI", "sdC", "nI", "nC") %in%
               colnames(data))) {
        stop("`data` must have the following columns:
           'mI','mC','sdI,'sdC','nI', and 'nC'.")
      } else if (outcome != "MD" &
                 !all(c("eI", "nI", "eC", "nC") %in%
                      colnames(data))) {
        stop("`data` must have the following columns:
           'eI','nI,'eC', and 'nC'.")
      }
    }
    # check | study
    if (!any(colnames(data) == "study") & !is.null(data)) {
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
    # check | cont_vartype
    if (!(cont_vartype %in% c("equal", "non-equal"))) {
      stop("`cont_vartype` must be either 'equal' or 'non-equal'")
    }
    # check | method
    if (!(weights %in% c("MH", "IV"))) {
      stop("`weights` must be either 'MH', or 'IV'")
    }
    # check | re_method
    if (!(re_method %in% c("DL", "DL_HKSJ"))) {
      stop("`re_method` must be either DL or DL_HKSJ")
    }
    # check | alpha
    if (is.null(alpha) | !is.numeric(alpha) | alpha > 1) {
      stop("`alpha` must be provided, numeric and below 1")
    }

    # check | beta
    if (is.null(beta) | !is.numeric(beta) | beta > 1) {
      stop("`beta` must be provided, numeric and below 1")
    }

    # check | confidence interval for tau
    if (!(tau_ci_method %in% c("BJ", "QP"))) {
      stop("`tau_ci_method` must be 'BJ' or 'QP'")
    }

    # check | futility
    if(!(futility %in% c("none", "non-binding", "binding"))){
      stop("`futility` must be 'none', 'non-binding' or 'binding'")
    }

    # check | alpha spending
    if(!(es_alpha %in% c("esPoc", "esOF", "HSDC", "rho"))){
      stop("`es_alpha` must be 'esOF', 'esPOC', 'HSDC' or 'rho'")
    }

    # check | beta spending
    if(!is.null(es_beta)){
      if(!(es_beta %in% c("esPoc", "esOF", "HSDC", "rho"))){
      stop("`es_beta` must be 'esOF', 'esPOC', 'HSDC' or 'rho' if not NULL.")
      }
    }

    if (futility != "none" & is.null(es_beta)) {
      stop(paste(
        "Set spending function for beta spending. Alpha spending is set to:",
        es_alpha
      ))
    }
    
    # check | timing
    if(type == "design"){
      dif <-
        timing - c(0, timing[-length(timing)])
      if(sum(dif < 0.01) > 0){
        stop("Timings must be such that studies provide more than 1% of information.")
      }
    }

    argg <- c(as.list(environment()), list(...))

    # helper functions
    logit <- function(x)
      log(x / (1 - x))
    invlogit <- function(x)
      1 / (1 + exp(-x))

    war_order = NULL
    war_p0 = NULL
    war_tim = NULL
    war_design = NULL
    war_ana = NULL
    war_design = NULL
    war_het_design = NULL

    if (!is.null(data$order) & type == "analysis") {
      data <- data[sort(data$order), ]
    } else if(is.null(data$order) & type == "analysis"){
      war_order <-
        c(
          "The order of the Trial Sequential Analysis will be based on the order of the studies in the data-set. Please add a 'order' column in the data-set to specify the order."
        )
    }

    # calculate the meta-analysis
    if (!is.null(data)) {
      ma <- metaanalysis(outcome = outcome, data = data,
                   weights = weights, cont_vartype = cont_vartype,
                   alpha = alpha, zero_adj = zero_adj,
                   conf_level = conf_level,
                   re_method = re_method, tau_ci_method = tau_ci_method)

      mp <- ma$metaPrepare
      sy <- ma$synthesize
      hete_results <- ma$hete_results

      # Calculate the cumulative number of participants
      subjects <- cumsum(data$nI + data$nC)

      # Calculate the RIS
      if (outcome %in% c("RR", "OR")) {
        if (is.null(p0) & is.null(design)) {
          p0 = sum(data$eC + data$eI) / sum(data$nC + data$nI)
          war_p0 <- NULL
        } else {
          p0 <- ifelse(!is.null(design),design$settings$p0,p0)
          war_p0 <- c(
            paste0(
              "Prob. of event in the control group is set to ",
              ifelse(!is.null(design),design$settings$p0,p0),
              ". The observed prob. of event is ",
              round(sum(data$eC + data$eI) / sum(data$nC + data$nI), 4),
              ". The power of the sequential might be affected."
            )
          )
        }
      }
      
      if (outcome %in% c("MD")) {
        if (is.null(sd_mc) & is.null(design)) {
          if(sy$U[1] > 0){
            sd_mc = sqrt(sy$peR[6])
          } else {
            sd_mc = sqrt(sy$peF[7])
          }
          war_p0 <- paste0("Standard deviation is set to", round(sd_mc,4))
        } else {
          sd_mc <- ifelse(!is.null(design),design$settings$sd_mc,sd_mc)
          war_p0 <- paste0("Standard deviation is set to", round(sd_mc,4))
        }
      }

      mc <- ifelse(!is.null(design),design$settings$mc,mc)

      if (outcome == "RR") {
        pI = exp(log(p0) + log(mc)) # pI = exp(log(p0)+log(mc)/2)
      } else if (outcome == "OR") {
        pI = invlogit(logit(p0) + log(mc))
      }

      if(is.null(design)){
      if (side == 1 & (outcome == "RR" | outcome == "OR")) {
        outris = ris(
          outcome = outcome,
          mc = mc,
          side = 1,
          alpha = alpha,
          beta = beta,
          fixed = fixed,
          p0 = p0,
          tau2 = sy$U[1],
          D2 = sy$U[4]
        )

      } else if (outcome == "RR" | outcome == "OR") {
        outris = ris(
          outcome = outcome,
          mc = mc,
          side = 2,
          alpha = alpha,
          beta = beta,
          fixed = fixed,
          p0 = p0,
          tau2 = sy$U[1],
          D2 = sy$U[4]
        )
      } else if (side == 1 & outcome == "MD") {
        outris = ris(
          outcome = outcome,
          mc = mc,
          side = side,
          alpha = alpha,
          beta = beta,
          sd_mc = sd_mc,
          fixed = fixed,
          tau2 = sy$U[1],
          D2 = sy$U[4]
        )
      } else {
        outris = ris(
          outcome = outcome,
          mc = mc,
          side = side,
          alpha = alpha,
          beta = beta,
          fixed = FALSE,
          sd_mc = sd_mc,
          tau2 = sy$U[1],
          D2 = sy$U[4]
        )
      }

      if(sy$U[1] == 0 | fixed){
        RIS = outris$NF
      } else {
        RIS = outris$NR_D2
      }
      } else {
        RIS = ceiling(design$results$RIS)
      }
    } else {
      # calculate ris without data
      if(fixed == FALSE){
        war_het_design <- "Note that you have adjusted for heterogeneity in the prospective sample size calculation. If the estimate of the heterogeneity is different than the expected size, the analysis is most likely not valid."
      }

      if(outcome == "MD"){
      outris = ris(
          outcome = outcome,
          mc = mc,
          sd_mc = sd_mc,
          side = side,
          alpha = alpha,
          beta = beta,
          fixed = fixed,
          tau2 = tau2,
          I2 = I2,
          D2 = D2
        )
      } else if(outcome %in% c("RR", "OR")){
        if (outcome == "RR") {
          pI = exp(log(p0) + log(mc)) # pI = exp(log(p0)+log(mc)/2)
        } else if (outcome == "OR") {
          pI = invlogit(logit(p0) + log(mc))
        }

        outris = ris(
          outcome = outcome,
          mc = mc,
          p0 = p0,
          side = side,
          alpha = alpha,
          beta = beta,
          fixed = fixed,
          tau2 = tau2,
          I2 = I2,
          D2 = D2
        )
      }

      if(fixed == TRUE){
        RIS = outris$NF
      } else {
        if(!is.null(tau2)){
          RIS = outris$NR$nPax[3,1]
        } else if(!is.null(I2)){
          RIS = outris$NR_I2
        } else {
          RIS = outris$NR_D2
        }
      }
    }

    # calculate boundaries for design
    war_ana <- NULL
    war_design <- NULL
    
    if (type == "design") {
      bounds <-
        boundaries(
          timing = timing,
          alpha = alpha,
          beta = beta,
          side = side,
          futility = futility,
          es_alpha = es_alpha,
          es_beta = es_beta,
          type = type
        )
      war_design <- c(
        "The RTSA function is used for design. Boundaries are computed but sequential inference will not be calculated. Use the metaanalysis() function if interested in meta-analysis results."
      )
    } else if(type == "analysis" & is.null(design)){
      
      # calculate timings
      timing <- c(subjects / RIS)
      
      # set analysis time
      if (is.null(ana_time)) {
        ana_time <- 1:length(timing)
      } 
      org_ana_time <- ana_time
      orgTiming = timing
      
      if (max(timing) > 1) {
        war_tim <- c("There might be more information than needed. ")
        timing <- timing[timing <= 1]
      } else {
        war_tim <- NULL
      }
      
      if(!is.null(design)){
        timing <- c(timing, design$bounds$root) 
      }
    
    if(length(timing) < max(ana_time)){
      ana_time <- ana_time[ana_time <= length(timing)]
      if(length(ana_time) == 0) { 
        timing = 1
        ana_time <- 1
      } else {
        timing = timing[ana_time]
      }
    } else {
      timing = timing[ana_time]
    }
    
    trials <- cbind(timing, NA, 1:length(timing))
    
    time_tf = 0.01
    war_tim2 <- NULL
    
    if(length(ana_time) == 1){
      trials[2] <- trials[1] - 0
      if (trials[2] < time_tf) {
        war_tim2 <- c(
          "Interim analyses will not be made with less than 1% increase in information.\nThe number of studies and number of interim analysis will not be the same."
        )
        stop("Not enough information in data to make analysis")
      }
      
    } else  {
      trials[, 2] <-
        trials[, 1] - c(0, trials[, 1][-length(trials[, 1])])
      
      # if new study adds less than 1% of RIS, the analysis is not performed.
      if (length(which(trials[, 2] > time_tf)) != length(trials[, 1])) {
        war_tim2 <- c(
          "Interim analyses will not be made with less than 1% increase in information.\nThe number of studies and number of interim analysis will not be the same."
        )
      }
      trials <-
        matrix(trials[trials[, 2] > time_tf, ], ncol = 3, byrow = F)
      if(length(ana_time) > 0) ana_time = trials[, 3]
      
      trials[, 2] <-
        trials[, 1] - c(0, trials[, 1][-length(trials[, 1])]) 
    }
      
        war_ana <- c(
          "design was not provided for the analysis. The analysis will be retrospective and the results validity is affected."
        )
        
        if(length(ana_time) == 1){
          timing <- trials[1] } else {
            timing <- trials[, 1] }

        bounds <-
          boundaries(
            timing = timing,
            alpha = alpha,
            beta = beta,
            side = side,
            futility = futility,
            es_alpha = es_alpha,
            es_beta = es_beta,
            type = "design"
          )
        design_R <- bounds$root
    } 
    
    if(type == "analysis"){
    if (!is.null(design)) {
          war_ana <- NULL
          design_R <- design$bounds$root
    }
      
      # recalculate timings
      timing <- c(subjects / RIS)
      orgTiming <- timing
      
      # set analysis time
      if (is.null(ana_time)) {
        ana_time <- 1:length(timing)
        org_ana_time <- NULL
      } 
      if(is.null(org_ana_time)){
        org_ana_time <- ana_time
      }
      orgTiming = timing
      
      timing <- timing[org_ana_time]
      
      if(max(timing) < design_R){
        timing <- c(timing, design_R)
        ana_time <- org_ana_time[org_ana_time <= length(timing)]
      } else if(max(timing) > design_R){
        timing <- c(timing[timing < design_R], design_R) 
        ana_time <- org_ana_time[org_ana_time <= length(timing)]
      }
      
      if(length(ana_time) == 0){
        ana_time <- 1
      }

      trials <- cbind(timing, NA, 1:length(timing))
      trials[, 2] <-
        trials[, 1] - c(0, trials[, 1][-length(trials[, 1])])

      # if new study adds less than 1% of RIS, the analysis is not performed.
      time_tf = 0.01
      war_tim2 <- NULL
      if (length(which(trials[, 2] > time_tf)) != length(trials[, 1])) {
        war_tim2 <- c(
          "Interim analyses will not be made with less than 1% increase in information.\nThe number of studies and number of interim analysis will not be the same."
        )
      }
      
      trials <-
        matrix(trials[trials[, 2] > time_tf, ], ncol = 3, byrow = F)
      if(length(ana_time) > 0) ana_time = trials[, 3]

      trials[, 2] <-
        trials[, 1] - c(0, trials[, 1][-length(trials[, 1])]) 

      if (max(trials[, 1]) < 1) {
        timing <- c(trials[, 1], design_R) 
      } else {
        timing <- trials[, 1]
      }

      if(!is.null(design)){
        
        bounds <- boundaries(
          timing = timing,
          alpha = design$settings$alpha,
          beta = design$settings$beta,
          side = design$settings$side,
          futility = design$settings$futility,
          es_alpha = design$settings$es_alpha,
          es_beta = design$settings$es_beta,
          type = type,
          design_R = design_R
        )
      } else {
        
        bounds <- boundaries(
          timing = timing,
          alpha = alpha,
          beta = beta,
          side = side,
          futility = futility,
          es_alpha = es_alpha,
          es_beta = es_beta,
          type = type,
          design_R = design_R
        )
      }

      if((outcome %in% c("OR", "RR") & mc < 1) | (outcome %in% c("MD", "RD") & mc < 0)){
        direction <- -1
      } else {direction <- 1}
      
      # inference
      inf <- inference(
        design = bounds,
        timing = timing,
        ana_time = ana_time,
        orgTiming = orgTiming,
        ma = ma,
        direction = direction,
        fixed = ifelse(!is.null(design), design$settings$fixed, FALSE),
        conf_int = conf_int,
        conf_level = conf_level
        )
    }

    RTSAout <- list()
    RTSAout$settings = argg

    if (outcome %in% c("RR", "OR", "RD"))
      RTSAout$settings$Pax <-
      list(
        RIS = RIS,
        p0 = p0,
        pI = pI,
        pC = p0
      )

    # store sample and trial size calculation
    if(!is.null(design)){
      RTSAout$design_ris <- design$ris
    } else {
      RTSAout$ris <- outris
    }

    if(!is.null(data)) RTSAout$settings$Pax$subjects <- subjects
    RTSAout$bounds = bounds
    #if(!is.null(data)) RTSAout$orgTiming = orgTiming
    RTSAout$results$RIS = RIS
    RTSAout$results$DARIS = RIS * RTSAout$bounds$root
    RTSAout$results$DARIS_F = RTSAout$ris$NF * RTSAout$bounds$root

    #if(!is.null(data)) RTSAout$results$timing <- c(subjects / RTSAout$results$DARIS)
    if(!is.null(data)) RTSAout$results$AIS = sum(data$nC + data$nI)
    #if(!is.null(data)) RTSAout$het = heteResults
    RTSAout$warnings <- list(
      war_order = war_order,
      war_p0 = war_p0,
      war_tim = war_tim,
      war_design = war_design,
      war_ana = war_ana,
      war_design = war_design,
      war_het_design = war_het_design
    )

    if(type == "analysis" & !is.null(design)){
      RTSAout$settings$side <- design$settings$side
      RTSAout$settings$outcome <- design$settings$outcome
      RTSAout$settings$mc <- design$settings$mc
      RTSAout$settings$sd_mc <- design$settings$sd_mc
      RTSAout$settings$p0 <- design$settings$p0
      RTSAout$settings$alpha <- design$settings$alpha
      RTSAout$settings$beta <- design$settings$beta
      RTSAout$settings$futility <- design$settings$futility
      RTSAout$settings$fixed <- design$settings$fixed
      RTSAout$settings$es_alpha <- design$settings$es_alpha
      RTSAout$settings$es_beta <- design$settings$es_beta
    }

    if (type == "design"){
      if (argg$side == 1) {
        RTSAout$results$design_df <- data.frame(
          "sma_timing" = bounds$inf_frac*bounds$root,
          "upper" = bounds$alpha_ubound
        )
        if (argg$futility != "none") {
          RTSAout$results$design_df <- cbind(RTSAout$results$design_df,
                                              data.frame("fut_lower" = bounds$beta_lbound))
        }
      } else {
        RTSAout$results$design_df <- data.frame(
          "sma_timing" = bounds$inf_frac*bounds$root,
          "upper" = bounds$alpha_ubound,
          "lower" = bounds$alpha_lbound
        )
        if (argg$futility != "none") {
          RTSAout$results$design_df <- cbind(
            RTSAout$results$design_df,
            data.frame(
              "fut_upper" = bounds$beta_ubound,
              "fut_lower" = bounds$beta_lbound
            )
          )
        }
      }
    }

    if (type == "analysis") {
      RTSAout$results$results_df = inf$results_df
      RTSAout$results$seq_inf = inf$seq_inf
      RTSAout$results$metaanalysis <- ma
      RTSAout$results$design_df <- design$results$design_df
    }
    class(RTSAout) <- c("RTSA", "list")
    return(RTSAout)
  }

# FUNCTION | print RTSA ----
#' @method print RTSA
#' @importFrom scales percent
#'
#' @export
print.RTSA <- function(x, ...) {
  if(x$settings$type == "design"){
    cat("Design with ")
  }
  cat("Trial Sequential Analysis was computed with the following settings: \n \n")
  cat(
    paste0(
      "Boundaries for a ",
      x$settings$side,
      "-sided design with a type-I-error of ",
      x$settings$alpha,
      ", and type-II-error of ",
      x$settings$beta,
      ".\n"
    )
  )
  cat(
    paste0(
      "Futility is set to: ",
      x$settings$futility,
      ". Alpha-spending function: ",
      x$settings$es_alpha,
      ".\n",
      "Beta-spending function: ",
      x$settings$es_beta,
      ".\n\n"
    )
  )

  cat("Timing,", ifelse(x$settings$type == "design",
                        "and boundaries:\n","boundaries, and test statistic:\n"))
  if(x$settings$type == "design"){
    y <- x$results$design_df
    print(round(y, 3), row.names = FALSE)
    cat("sma_timing is the ratio of the required sample for a sequential meta-analysis to a non-sequential meta-analysis sample size.")
    cat("\n\nSample size calculation for standard meta-analysis:\n")
    cat(print(x$ris))
    cat("\n")
    cat("\nSample size calculation for sequential meta-analysis:\n")
    cat("Fixed-effect:", paste0(ceiling(x$results$DARIS_F)," participants."))
    if(x$settings$fixed == FALSE) cat("\nRandom-effects:", paste0(ceiling(x$results$DARIS)," participants."))
  }
  if(x$settings$type == "analysis"){
    y <- x$results$results_df[,c("sma_timing", "upper", "lower", "fut_upper", "fut_lower")]
  y <-
    cbind(
      y,
      data.frame(
        "z_fixed" = x$results$results_df$z_fixed,
        "z_random" = x$results$results_df$z_random
      )
    )
  print(round(y, 3), row.names = FALSE)
  cat("sma_timing is the ratio of the required sample for a sequential meta-analysis to a non-sequential meta-analysis sample size.")
  cat("\n\nTiming, outcomes, and confidences intervals for fixed-effect and random-effects models:\n")
  z <- data.frame(
    "sma_timing" = x$results$results_df$sma_timing,
    "outcome1" = x$results$results_df$outcome_fixed,
    "LCI1" = x$results$results_df$TSAadjCIfixed_lower,
    "UCI1" = x$results$results_df$TSAadjCIfixed_upper,
    "outcome2" = x$results$results_df$outcome_random,
    "LCI2" = x$results$results_df$TSAadjCIrandom_lower,
    "UCI2" = x$results$results_df$TSAadjCIrandom_upper
  )
  colnames(z)[2:7] <- c(
    paste0(x$settings$outcome, "_fixed"),
    paste0("TSA_", x$settings$conf_level, "lci_fixed"),
    paste0("TSA_", x$settings$conf_level, "uci_fixed"),
    paste0(x$settings$outcome, "_random"),
    paste0("TSA_", x$settings$conf_level, "lci_random"),
    paste0("TSA_", x$settings$conf_level, "uci_random")
  )
  print(round(z[,c(1:4)], 3), row.names = FALSE)
  print(round(z[,c(1,5:7)], 3), row.names = FALSE)
  cat("lci is the lower limit of the confidence interval. uci is the upper limit of the confidence interval.")

  cat("\n\nMeta-analysis results:\n")
  tmp_ca <- x$settings$alpha / x$settings$side

  #LABELS
  df <- x$results$results_df
  f_tmp_outcome <- df$outcome_fixed[!is.na(df$outcome_fixed)]
  f_tmp_outcome <- f_tmp_outcome[length(f_tmp_outcome)]
  f_tmp_lcl <- df$TSAadjCIfixed_lower
  f_tmp_lcl1 <-
    df$TSAadjCIfixed_lower[!is.na(df$TSAadjCIfixed_lower)]
  f_tmp_lcl1 <- f_tmp_lcl1[length(f_tmp_lcl1)]
  f_tmp_ucl <- df$TSAadjCIfixed_upper
  f_tmp_ucl1 <-
    df$TSAadjCIfixed_upper[!is.na(df$TSAadjCIfixed_upper)]
  f_tmp_ucl1 <- f_tmp_ucl1[length(f_tmp_ucl1)]
  f_tmp_pvalue <- df$pvalues_fixed[!is.na(df$pvalues_fixed)]
  f_tmp_pvalue <- f_tmp_pvalue[length(f_tmp_pvalue)]

  r_tmp_outcome <- df$outcome_random[!is.na(df$outcome_random)]
  r_tmp_outcome <- r_tmp_outcome[length(r_tmp_outcome)]
  r_tmp_lcl <- df$TSAadjCIrandom_lower
  r_tmp_lcl1 <-
    df$TSAadjCIrandom_lower[!is.na(df$TSAadjCIrandom_lower)]
  r_tmp_lcl1 <- r_tmp_lcl1[length(r_tmp_lcl1)]
  r_tmp_ucl <- df$TSAadjCIrandom_upper
  r_tmp_ucl1 <-
    df$TSAadjCIrandom_upper[!is.na(df$TSAadjCIrandom_upper)]
  r_tmp_ucl1 <- r_tmp_ucl1[length(r_tmp_ucl1)]
  r_tmp_pvalue <- df$pvalues_random[!is.na(df$pvalues_random)]
  r_tmp_pvalue <- r_tmp_pvalue[length(r_tmp_pvalue)]


  results_fixed <- paste0(
    "Fixed Pooled effect (",
    x$settings$outcome,
    "): ",
    format(round(f_tmp_outcome, 2), nsmall = 2),
    " (95% TSA-adjusted CI: ",
    format(round(f_tmp_lcl1[length(f_tmp_lcl1)], 2), nsmall = 2),
    "-",
    format(round(f_tmp_ucl1[length(f_tmp_ucl1)], 2), nsmall = 2),
    "); naive p-value: ",
    format(round(f_tmp_pvalue, 4), nsmall = 4)
  )

  results_random <- paste0(
    "Random Pooled effect (",
    x$settings$outcome,
    "): ",
    format(round(r_tmp_outcome, 2), nsmall = 2),
    " (95% TSA-adjusted CI: ",
    format(round(r_tmp_lcl1[length(r_tmp_lcl1)], 2), nsmall = 2),
    "-",
    format(round(r_tmp_ucl1[length(r_tmp_ucl1)], 2), nsmall = 2),
    "); naive p-value: ",
    format(round(r_tmp_pvalue, 4), nsmall = 4)
  )

  if(!is.null(x$results$seq_inf)){
    if(x$results$seq_inf$lower > x$results$seq_inf$upper){
      temp <- x$results$seq_inf$lower
      x$results$seq_inf$lower <- x$results$seq_inf$upper
      x$results$seq_inf$upper <- temp
    }

    results_sw <- paste0(
      "Median unbiased pooled effect (",
      x$settings$outcome,
      "): ",
      format(round(x$results$seq_inf$median_unbiased, 2), nsmall = 2),
      " (95% SW-adjusted CI: ",
      format(round(x$results$seq_inf$lower, 2), nsmall = 2),
      "-",
      format(round(x$results$seq_inf$upper, 2), nsmall = 2),
      "); SW p-value: ",
      format(round(x$results$seq_inf$p.value, 4), nsmall = 4)
    )
  } else {results_sw <- ""}

  #CREATE LABELS
  settings <- paste0(
    "Pc: ",
    percent(x$Pax$p0, 0.1),
    "; ",
    "MC: ",
    percent(x$settings$mc, 0.1),
    "; ",
    "Alpha: ",
    percent(x$settings$alpha, 0.1),
    "; ",
    "Beta: ",
    percent(x$settings$beta),
    "\n",
    "Methods: Random effects, ",
    x$settings$re_method,
    "; ",
    "Weight, ",
    x$settings$weights,
    "; ",
    "Alpha spending, ",
    x$settings$es_alpha,
    "; ",
    "Beta spending, ",
    x$settings$es_beta
  )

  results <- paste0(
    results_fixed,
    "\n",
    results_random,"\n",
    results_sw,
    "\nHeterogeneity results:\n",
    "tau^2: ",
    format(round(x$results$metaanalysis$hete_results$hete_est$tau2, 2), nsmall = 2),
    "; ",
    "I^2: ",
    percent(x$results$metaanalysis$hete_results$hete_est$I.2, 0.1),
    "; ",
    "D^2: ",
    percent(x$results$metaanalysis$hete_results$hete_est$D.2, 0.1),
    "; ",
    "Heterogeneity p-value: ",
    format(round(x$results$metaanalysis$hete_results$hete_est$Q_pval, 4), nsmall = 4)
  )
  cat(results)}

  if(!is.null(unlist(x$warnings))) cat("\n\nPlease note the following warnings:\n")
  if (!is.null(x$warnings$war_order))
    cat("-", x$warnings$war_order, "\n\n")
  if (!is.null(x$warnings$war_p0))
    cat("-", x$warnings$war_p0, "\n\n")
  if (!is.null(x$warnings$war_tim))
    cat("-", x$warnings$war_tim, "\n\n")
  if (!is.null(x$warnings$war_design))
    cat("-", x$warnings$war_design, "\n\n")
  if (!is.null(x$warnings$war_ana))
    cat("-", x$warnings$war_ana)
  if (!is.null(x$warnings$war_het_design))
    cat("-", x$warnings$war_het_design)
}
