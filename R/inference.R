#' Inference from RTSA
#' @description
#' Naive inference, TSA-adjusted confidence intervals and stage-wise ordered inference if the meta-analysis crosses a alpha-spending boundary or a binding beta-spending boundary.
#' 
#'
#' @param design design
#' @param timing timing
#' @param ana_time analysis times
#' @param ma meta-analysis item
#' @param fixed Whether the analysis is fixed or random.
#' @param orgTiming The original timing of all included studies
#' @param direction whether the upper or lower boundary was crossed.
#' @param conf_int type of confidence interval
#' @param conf_level confidence interval level
#' @param tol tolerance
#'
#' @returns A data.frame of cumulative meta-analysis results including stopping boundaries and a list of conditional sequential inference to be parsed to RTSA
#' \item{results_df}{A data.frame containing information about: Cumulative test values, cumulative outcomes, timing of trials,
#' stopping boundaries (alpha_upper, alpha_lower, beta_upper, beta_lower), naive confidence intervals, TSA-adjusted confidence intervals,
#' cumulative p-values and standard deviations.}
#' \item{seq_inf}{If the meta-analysis crosses a alpha-spending boundary or a binding beta-spending boundary inference conditional on stopping is provided. A median unbiased estimate, lower and upper confidence interval, and p-value is provided based on stage-wise ordering.}
#'
#' @export
#'
#' @examples
#' ma <- metaanalysis(data = perioOxy, outcome = "RR", mc = 0.8)
#' sts <- ma$ris$NR_D2$NR_D2_full
#' timing <- cumsum(perioOxy$nI + perioOxy$nC)/sts
#' bound_oxy <- boundaries(timing = timing, alpha = 0.05, beta = 0.2, side = 2,
#'                        futility = "none", es_alpha = "esOF")
#' inference(timing = bound_oxy$inf_frac, design = bound_oxy, ma = ma, direction = -1, fixed = FALSE,
#' ana_time = 1:length(timing), orgTiming = timing)
#' 
inference <- function(design,
                      timing,
                      ana_time,
                      ma,
                      fixed,
                      orgTiming,
                      conf_int = "sw",
                      conf_level = 0.95,
                      direction,
                      tol = 1e-15){
  
  # extract methods and data from meta-analysis
  mp <- ma$metaPrepare

  # calculate the cum. z-score, standard errors, pooled estimates etc.
  zout = lapply(1:length(mp$w),
                function(x) {
                  synout = synthesize(
                    metaPrepare(
                      data = mp$data[1:x,],
                      outcome = ma$settings$outcome,
                      weights = mp$weights,
                      alpha = design$alpha,
                      conf_level= conf_level,
                      cont_vartype = ma$settings$cont_vartype
                    ), re_method = ma$settings$re_method,
                    conf_level= conf_level,
                    tau_ci_method = ma$settings$tau_ci_method
                  )
                  return(synout)
                })

  # set names of the cumulative meta-analyses
  names(zout) = 1:dim(mp$data)[1]
  if("1" %in% names(zout) &
     length(grep("peR", names(unlist(zout)))) > 0){
    zout[["1"]]$peR <- c(0, 0, 0, zout[["1"]]$peF[4])
  }

  # collect the outcomes in matrices
  zvalues = sapply(names(zout), function(x){c(zout[[x]]$peF[4],
                                              ifelse(length(grep("peR", names(unlist(zout[[x]]))))!=0,zout[[x]]$peR[4], NA))})

  sd_values = sapply(names(zout), function(x){c(zout[[x]]$peF[7],
                                              ifelse(length(grep("peR", names(unlist(zout[[x]]))))!=0,zout[[x]]$peR[6], NA))})

  outcome_values = sapply(names(zout), function(x){c(zout[[x]]$peF[1],
                                                     ifelse(length(grep("peR", names(unlist(zout[[x]]))))!=0,zout[[x]]$peR[1], NA))})

  p_values = sapply(names(zout), function(x){c(zout[[x]]$peF[5],
                                                     ifelse(length(grep("peR", names(unlist(zout[[x]]))))!=0,zout[[x]]$peR[5], NA))})

  outcome_values[2,1] <- outcome_values[1,1]
  outcome_values[2,is.na(outcome_values[2,])] <- outcome_values[1,is.na(outcome_values[2,])]
  zvalues[2,is.na(zvalues[2,])] <- zvalues[1,is.na(zvalues[2,])]
  p_values[2,is.na(p_values[2,])] <- p_values[1,is.na(p_values[2,])]

  # if we have more planned analyses than actual analyses
  if(length(timing) > length(ana_time)){
    extra <- TRUE
  } else {
    extra <- FALSE
  }
  ll <- -Inf

  if(length(ana_time) > 0){
  if(dim(zvalues)[2] < length(ana_time)){
    ana_time <- 1:dim(zvalues)[2]
  } else if(max(ana_time) > dim(zvalues)[2]){
    ana_time <- ana_time[ana_time <= dim(zvalues)[2]]
    extra <- TRUE
  }
  }

  # set the boundaries that we compare the test statistic with
  if(extra | length(ana_time) < length(design$alpha_ubound)){
    ul <- design$alpha_ubound[-length(design$alpha_ubound)]
    if(design$side == 2){ ll <- design$alpha_lbound[-length(design$alpha_lbound)] 
    if(design$futility == "binding"){
      bll <- design$beta_lbound[-length(design$beta_lbound)]
      bul <- design$beta_ubound[-length(design$beta_ubound)]
    } else {
      bll <- -Inf
      bul <- Inf }
    }
    if(design$side == 1 & design$futility == "binding"){ ll <- design$beta_lbound[-length(design$beta_lbound)]
    } else if(design$side == 1 & design$futility != "binding"){ ll <- -Inf }
  } else {
    ul <- design$alpha_ubound
    if(design$side == 2){ ll <- design$alpha_lbound 
    if(design$futility == "binding"){
      bll <- design$beta_lbound
      bul <- design$beta_ubound
    } else {
      bll <- -Inf
      bul <- Inf }
    }
    if(design$side == 1 & design$futility == "binding"){ ll <- design$beta_lbound
    } else if(design$side == 1 & design$futility != "binding"){ ll <- -Inf }
  }
  
  # if no random effects
  if(identical(zvalues[1,], zvalues[2,])){
    fixed = TRUE
  }
  
  test_values <- direction*zvalues[ifelse(fixed,1,2),ana_time]
  
  # if no set stop time is made, see if the analysis crossed a boundary
  if(design$side == 2){
    if(sum(test_values > ul) > 0){
      stop_time <- min(which((test_values > ul) == TRUE))
      stop_direction <- "ul"
    }  else if(sum(test_values < ll) > 0){
      stop_time <- min(which((test_values < ll) == TRUE))
      stop_direction <- "ll"
    }  else if(design$futility == "binding" & sum(test_values > bll & test_values < bul, na.rm = T) > 0){
      stop_time <- min(which((test_values > bll & test_values < bul) == TRUE))
      stop_direction <- "fut"
    } else {
      stop_time <- NULL
    }
  } else {
    if(sum(test_values > ul) > 0){
      stop_time <- min(which((test_values > ul) == TRUE))
      stop_direction <- "ul"
    }  else if(sum(test_values < ll) > 0){
      stop_time <- min(which((test_values < ll) == TRUE))
      stop_direction <- "fut"
    } else {
      stop_time <- NULL
    } 
  }
  
  # calculate naive p-values and CI
  naiveCI = list(CIfixed = sapply(ana_time, function(x) zout[[x]]$peF[c(2, 3)]),
                   CIrandom = sapply(ana_time, function(x) {
                     if(zout[[x]]$U[1] > 0 & !fixed) zout[[x]]$peR[c(2, 3)] else matrix(nrow = 2)}))
    if(ma$settings$outcome %in% c("RR","OR")){
      TSAadjCI = list(
        CIfixed = sapply(ana_time, function(x) {exp(
          log(zout[[x]]$peF[1]) +
            c(-1, 1) * design$alpha_ubound[which(x == ana_time)] *
            sqrt(zout[[x]]$peF[7])
        ) } ) )
      if(!fixed){
        TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed, CIrandom = sapply(ana_time, function (x) {
          if(zout[[x]]$U[1] > 0) {exp(
            log(zout[[x]]$peR[1]) +
              c(-1, 1) * design$alpha_ubound[which(x == ana_time)] *
              sqrt(zout[[x]]$peR[6])) } else {matrix(nrow = 2)}}))
      } else {
        TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed,CIrandom = TSAadjCI$CIfixed)
      }
    } else {
      TSAadjCI = list(
        CIfixed = sapply(ana_time, function(x) {zout[[x]]$peF[1] +
            c(-1, 1) * design$alpha_ubound[which(x == ana_time)] *
            sqrt(zout[[x]]$peF[7]) }))
      if(!fixed){
        TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed,CIrandom = sapply(ana_time, function (x) {
          if(zout[[x]]$U[1] > 0){
            zout[[x]]$peR[1] +
              c(-1, 1) * design$alpha_ubound[which(x == ana_time)] *
              sqrt(zout[[x]]$peR[6]) } else {matrix(nrow = 2)}} ))} else {
                TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed,CIrandom = TSAadjCI$CIfixed)
              }
    }
  
  # calculate minimum clinical relevant value "boundary"
  if(fixed){
    stnd_dv <- sapply(ana_time, function(x) sqrt(zout[[x]]$peF[7])) 
    } else {
      stnd_dv <- sapply(ana_time, function(x) sqrt(zout[[x]]$peR[6]))
    }
  
  # if the analysis crossed a boundary
  if (!is.null(stop_time)) {
    stnd_dv <- ifelse(fixed, sqrt(zout[[stop_time]]$peF[7]), sqrt(zout[[stop_time]]$peR[6]))
    info_ana <- sd_inf(design$inf_frac * design$root)
    if (conf_int == "sw") {
      # if the type of inference is stage-wise (only option for now)
      # zscore <- direction * zvalues[ifelse(fixed, 1, 2), stop_time]
      if (design$side == 2) {
        zscore <- direction*zvalues[ifelse(fixed, 1, 2), stop_time]
        zb <- ul[1:stop_time]; za <- ll[1:stop_time]; zd <- bul[1:stop_time]; zc <- bll[1:stop_time]
        zd <- zd[!is.na(zd)]; zc <- zc[!is.na(zc)]
        if(length(zd) == 1 & -Inf %in% zd) zd <- NULL
        if(length(zc) == 1 & -Inf %in% zc) zc <- NULL
        if(stop_direction == "ul") zb <- c(zb[-stop_time], zscore)
        if(stop_direction == "ll") za <- c(za[-stop_time], zscore)
        if(stop_direction == "fut") zd <- c(zd[-stop_time], zscore)
        
        if(stop_direction == "fut"){
          zb[stop_time] <- zd[stop_time]  
          
          sw_p <- pmin(sw_pvalue(side = design$side, info = info_ana, za = za, zb = zb,
                            zc = zc, zd = zd),1)
          
          lowci <- uniroot(sw_cilower,
                           upper = max(ul),
                           lower = min(ll),
                           conf_level = conf_level, info = info_ana,
                           za = za, zb = zb, zc = zc, zd = zd)$root
          sw.lower <- direction *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
          
          upci <- uniroot(sw_ciupper,
                          upper = 10,
                          lower = 0,
                          conf_level = conf_level, info = info_ana,
                          za = za, zb = zb, zc = zc, zd = zd)$root
          sw.upper <- direction *  upci * stnd_dv * info_ana$sd_proc[stop_time] 
        } else {
          
          sw_p <- sw_pvalue(side = design$side, info = info_ana, za = za, zb = zb,
                            zc = zc, zd = zd)
          
          lowci <- uniroot(sw_cilower,
                           upper = max(ul),
                           lower = min(ll),
                           conf_level = conf_level, info = info_ana,
                           za = za, zb = zb, zc = zc, zd = zd)$root
          sw.lower <- direction *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
          
          upci <- uniroot(sw_ciupper,
                          upper = 10,
                          lower = 0,
                          conf_level = conf_level, info = info_ana,
                          za = za, zb = zb, zc = zc, zd = zd)$root
          sw.upper <- direction *  upci * stnd_dv * info_ana$sd_proc[stop_time] 
        }
        
        } else {
          zscore <- direction*zvalues[ifelse(fixed, 1, 2), stop_time]
          zb <- ul[1:stop_time]; za <- ll[1:stop_time]; zc <- NULL; zd <- NULL
          if(stop_direction == "ul") zb <- c(zb[-stop_time], zscore)
          if(stop_direction == "fut") za <- c(za[-stop_time], zscore)
          if(length(za) == 1 & -Inf %in% za) za <- NULL
          
          if(stop_direction == "fut"){
            zb[stop_time] <- za[stop_time]
            
            sw_p <- sw_pvalue(side = design$side, info = info_ana, za = za, zb = zb,
                              zc = zc, zd = zd)
            
            lowci <- uniroot(sw_cilower,
                             upper = max(ul),
                             lower = ifelse(min(ll) == -Inf, -20, min(ll)),
                             conf_level = conf_level, info = info_ana,
                             za = za, zb = zb, zc = zc, zd = zd)$root
            sw.lower <- direction *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
            
            upci <- uniroot(sw_ciupper,
                            upper = 10,
                            lower = 0,
                            conf_level = conf_level, info = info_ana,
                            za = za, zb = zb, zc = zc, zd = zd)$root
            sw.upper <- direction *  upci * stnd_dv * info_ana$sd_proc[stop_time]   
          } else {
            sw_p <- sw_pvalue(side = design$side, info = info_ana, za = za, zb = zb,
                              zc = zc, zd = zd)
            
          lowci <- uniroot(sw_cilower,
                           upper = max(ul),
                           lower = ifelse(min(ll) == -Inf, -20, min(ll)),
                           conf_level = conf_level, info = info_ana,
                           za = za, zb = zb, zc = zc, zd = zd)$root
          sw.lower <- direction *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
          
          upci <- uniroot(sw_ciupper,
                          upper = 10,
                          lower = 0,
                          conf_level = conf_level, info = info_ana,
                          za = za, zb = zb, zc = zc, zd = zd)$root
          sw.upper <- direction *  upci * stnd_dv * info_ana$sd_proc[stop_time] 
          }
        }
    }
      
      median_unbiased_z <-
        uniroot(sw_ciupper, upper = 10, lower = -10, conf_level = 0,info = info_ana,
                za = za, zb = zb, zc = zc, zd = zd)$root
      
      median_unbiased <- direction * median_unbiased_z * stnd_dv * info_ana$sd_proc[stop_time]
      
    
    if (ma$settings$outcome %in% c("RR", "OR")) {
      sw.lower <- exp(sw.lower)
      sw.upper <- exp(sw.upper)
      median_unbiased <- exp(median_unbiased)
    }
    
    
    seq_inf <-
      list(
        median_unbiased = median_unbiased,
        lower = sw.lower,
        upper = sw.upper,
        p.value = sw_p,
        lower_z = lowci,
        upper_z = upci,
        median_unbiased_z = median_unbiased_z,
        stnd_dv = stnd_dv,
        info = info_ana
      )
    
  } else {
    seq_inf <- NULL
  }
  
  n_out <- ifelse(is.null(dim(zvalues)[2]), 1, dim(zvalues)[2])
  if(length(ana_time) > 0) n_out2 <- max(ana_time)
  n_row <- ifelse((is.null(stop_time) | length(ana_time) < length(timing)), n_out+1, n_out)
  design$inf_frac <- sort(unique(c(orgTiming,timing)))
  n_out <- which(design$inf_frac %in% orgTiming)
  
  results <- as.data.frame(matrix(ncol = 21, nrow = n_row))
  results[n_out,1:2] <- t(zvalues)
  results[n_out,3:4] <- t(outcome_values)
  results[ana_time,10:11] <- t(naiveCI$CIfixed)
  results[ana_time,12:13] <- t(naiveCI$CIrandom)
  results[ana_time,14:15] <- t(TSAadjCI$CIfixed)
  results[ana_time,16:17] <- t(TSAadjCI$CIrandom)
  results[n_out,18:19] <- t(p_values)
  results[n_out,20:21] <- t(sd_values)

  if(length(ana_time) < length(design$inf_frac) & length(ana_time) > 0){ana_time <- c(ana_time, max(ana_time)+1)}

  if(length(timing) > 1){
    indi_seq <- which(design$inf_frac %in% timing) } else {
      indi_seq <- which(timing == design$inf_frac)
    }
  
  results[,5] <- design$inf_frac
  results[indi_seq,6] <- design$alpha_ubound
  results[indi_seq,7] <- design$alpha_lbound
  results[indi_seq,8] <- design$beta_ubound
  results[indi_seq,9] <- design$beta_lbound


  colnames(results) <- c("z_fixed", "z_random", "outcome_fixed", "outcome_random",
                         "sma_timing", "upper", "lower", "fut_upper",
                         "fut_lower", "naiveCIfixed_lower", "naiveCIfixed_upper",
                         "naiveCIrandom_lower", "naiveCIrandom_upper", "TSAadjCIfixed_lower",
                         "TSAadjCIfixed_upper", "TSAadjCIrandom_lower",
                         "TSAadjCIrandom_upper", "pvalues_fixed", "pvalues_random",
                         "sdvalues_fixed", "sdvalues_random")

  inf_out =     list(results_df = results,
                    seq_inf = seq_inf)

  return(
    inf_out
  )
}

# define - p-value function for stage-wise ---
sw_pvalue <- function(side, info, za, zb, zc, zd){
  return(side *
    ma_power(zb = zb, za = za, zc = zc, zd = zd, info = info,
      delta = 0)[[2]])
} 

# define - lower ci limit for stage-wise ---
sw_cilower <- function(x, conf_level, za, zb, zc, zd, info) {
    return((1-conf_level) / 2 - ma_power( zb = zb, za = za, zc = zc, zd = zd,
      info = info, delta = x)[[2]])
  }

# define - upper ci limit for stage-wise ---
sw_ciupper <- function(x, conf_level, za, zb, zc, zd, info) {
    return((1 - (1-conf_level)/ 2) - ma_power(zb = zb, za = za, zd = zd,  zc = zc,
      info = info, delta = x)[[2]])
}
