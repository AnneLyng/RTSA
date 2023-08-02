#' Inference calculations for sequential meta-analysis
#' 
#' @description
#' Calculates point-estimates, p-values and confidence intervals. Computes naive inference and TSA-adjusted confidence intervals. If the meta-analysis crosses a alpha-spending boundary, a binding beta-spending boundary or reached the sequential RIS, stage-wise ordered inference is also calculated. This function is not supposed to be used individually for Trial Sequential Analysis (TSA). RTSA() is recommended for TSA. 
#' 
#' @param timing The timing of the studies relative to the sequential RIS. A vector consisting of values equal to the proportion of study participants out of the sequential RIS. 
#' @param bounds The boundaries for the analysis as calculated by the boundaries() function in RTSA.
#' @param ana_times The analysis times presented as a vector. Describes at which studies the meta-analyses were performed. If one expects that the meta-analysis was updated per study a vector from 1 to the number of studies included can be used.
#' @param ma A metaanalysis object from the metaanalysis function.
#' @param fixed Whether the analysis is for fixed-effect or random-effects meta-analysis. Options are TRUE (meta-analysis is fixed-effect) or FALSE (meta-analysis is random-effects). 
#' @param org_timing The timing of all included studies as a proportion of RIS and not sequential RIS. 
#' @param inf_type For now only option is "sw" (stage-wise). Type of inference used for point estimates, confidence intervals and p-values.
#' @param conf_level The confidence interval level. Defaults to 0.95 which is 95\%.
#' @param final_analysis Whether or not the this analysis is considered the final analysis.
#' @param tol The tolerance level. Set to 1e+09.
#' 
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
#' inference(timing = bound_oxy$inf_frac, bounds = bound_oxy, ma = ma,fixed = FALSE,
#' ana_times = 1:length(timing), org_timing = timing)
#' 
inference <- function(bounds,
                      timing,
                      ana_times,
                      ma,
                      fixed,
                      org_timing,
                      inf_type = "sw",
                      conf_level = 0.95,
                      final_analysis = FALSE,
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
                      alpha = bounds$alpha,
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
  if(length(timing) > length(ana_times)){
    extra <- TRUE
  } else {
    extra <- FALSE
  }
  ll <- -Inf

  if(length(ana_times) > 0){
  if(dim(zvalues)[2] < length(ana_times)){
    ana_times <- 1:dim(zvalues)[2]
  } else if(max(ana_times) > dim(zvalues)[2]){
    ana_times <- ana_times[ana_times <= dim(zvalues)[2]]
    extra <- TRUE
  }
  }

  # set the boundaries that we compare the test statistic with
  if(extra | length(ana_times) < length(bounds$alpha_ubound)){
    ul <- bounds$alpha_ubound[-length(bounds$alpha_ubound)]
    if(bounds$side == 2){ ll <- bounds$alpha_lbound[-length(bounds$alpha_lbound)] 
    if(bounds$futility == "binding"){
      bll <- bounds$beta_lbound[-length(bounds$beta_lbound)]
      bul <- bounds$beta_ubound[-length(bounds$beta_ubound)]
    } else {
      bll <- -Inf
      bul <- Inf }
    }
    if(bounds$side == 1 & bounds$futility == "binding"){ ll <- bounds$beta_lbound[-length(bounds$beta_lbound)]
    } else if(bounds$side == 1 & bounds$futility != "binding"){ ll <- -Inf }
  } else {
    ul <- bounds$alpha_ubound
    if(bounds$side == 2){ ll <- bounds$alpha_lbound 
    if(bounds$futility == "binding"){
      bll <- bounds$beta_lbound
      bul <- bounds$beta_ubound
    } else {
      bll <- -Inf
      bul <- Inf }
    }
    if(bounds$side == 1 & bounds$futility == "binding"){ ll <- bounds$beta_lbound
    } else if(bounds$side == 1 & bounds$futility != "binding"){ ll <- -Inf }
  }
  
  # if no random effects
  if(identical(zvalues[1,], zvalues[2,])){
    fixed = TRUE
  }
  
  test_values <- zvalues[ifelse(fixed,1,2),ana_times]
  
  # if no set stop time is made, see if the analysis crossed a boundary
  if(bounds$side == 2){
    if(sum(test_values > ul) > 0){
      stop_time <- min(which((test_values > ul) == TRUE))
      stop_direction <- "ul"
      stop_sign <- 1
    }  else if(sum(test_values < ll) > 0){
      stop_time <- min(which((test_values < ll) == TRUE))
      stop_direction <- "ll"
      stop_sign <- -1
    }  else if(bounds$futility == "binding" & sum(test_values > bll & test_values < bul, na.rm = T) > 0){
      stop_time <- min(which((test_values > bll & test_values < bul) == TRUE))
      stop_direction <- "fut"
      stop_sign <- 1
    } else {
      stop_time <- NULL; stop_sign <- NULL; stop_direction <- NULL
    }
  } else {
    if(sum(test_values > ul) > 0){
      stop_time <- min(which((test_values > ul) == TRUE))
      stop_direction <- "ul"
      stop_sign <- 1
    }  else if(sum(test_values < ll) > 0){
      stop_time <- min(which((test_values < ll) == TRUE))
      stop_direction <- "fut"
      stop_sign <- 1
    } else {
      stop_time <- NULL; stop_sign <- NULL; stop_direction <- NULL
    } 
  }
  
  # calculate naive p-values and CI
  naiveCI = list(CIfixed = sapply(ana_times, function(x) zout[[x]]$peF[c(2, 3)]),
                   CIrandom = sapply(ana_times, function(x) {
                     if(zout[[x]]$U[1] > 0 & !fixed) zout[[x]]$peR[c(2, 3)] else matrix(nrow = 2)}))
    if(ma$settings$outcome %in% c("RR","OR")){
      TSAadjCI = list(
        CIfixed = sapply(ana_times, function(x) {exp(
          log(zout[[x]]$peF[1]) +
            c(-1, 1) * bounds$alpha_ubound[which(x == ana_times)] *
            sqrt(zout[[x]]$peF[7])
        ) } ) )
      if(!fixed){
        TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed, CIrandom = sapply(ana_times, function (x) {
          if(zout[[x]]$U[1] > 0) {exp(
            log(zout[[x]]$peR[1]) +
              c(-1, 1) * bounds$alpha_ubound[which(x == ana_times)] *
              sqrt(zout[[x]]$peR[6])) } else {matrix(nrow = 2)}}))
      } else {
        TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed,CIrandom = TSAadjCI$CIfixed)
      }
    } else {
      TSAadjCI = list(
        CIfixed = sapply(ana_times, function(x) {zout[[x]]$peF[1] +
            c(-1, 1) * bounds$alpha_ubound[which(x == ana_times)] *
            sqrt(zout[[x]]$peF[7]) }))
      if(!fixed){
        TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed,CIrandom = sapply(ana_times, function (x) {
          if(zout[[x]]$U[1] > 0){
            zout[[x]]$peR[1] +
              c(-1, 1) * bounds$alpha_ubound[which(x == ana_times)] *
              sqrt(zout[[x]]$peR[6]) } else {matrix(nrow = 2)}} ))} else {
                TSAadjCI <- list(CIfixed = TSAadjCI$CIfixed,CIrandom = TSAadjCI$CIfixed)
              }
    }
  
  stnd_dv_func <- function(x){
    if(zout[[x]]$U[1] == 0){
      return(sqrt(zout[[x]]$peF[7]))
    } else {
      return(sqrt(zout[[x]]$peR[6]))
    }
  }
  
  # calculate minimum clinical relevant value "boundary"
  if(fixed){
    stnd_dv <- sapply(ana_times, function(x) sqrt(zout[[x]]$peF[7])) 
    } else {
      stnd_dv <- sapply(ana_times, stnd_dv_func)
    }
  
  # if the analysis crossed a boundary
  if (!is.null(stop_time)) {
    stnd_dv <- ifelse(fixed | zout[[stop_time]]$U[1] == 0, sqrt(zout[[stop_time]]$peF[7]),
                      sqrt(zout[[stop_time]]$peR[6]))
    info_ana <- sd_inf(bounds$inf_frac * bounds$root)
    if (inf_type == "sw") {
      # if the type of inference is stage-wise (only option for now)
      # zscore <- direction * zvalues[ifelse(fixed, 1, 2), stop_time]
      if (bounds$side == 2) {
        zscore <- stop_sign*zvalues[ifelse(fixed, 1, 2), stop_time]
        zb <- ul[1:stop_time]; za <- ll[1:stop_time]; zd <- bul[1:stop_time]; zc <- bll[1:stop_time]
        zd <- zd[!is.na(zd)]; zc <- zc[!is.na(zc)]
        if(length(zd) == 1 & -Inf %in% zd) zd <- NULL
        if(length(zc) == 1 & -Inf %in% zc) zc <- NULL
        if(stop_direction == "ul") zb <- c(zb[-stop_time], zscore)
        if(stop_direction == "ll") za <- c(za[-stop_time], zscore)
        if(stop_direction == "fut") zd <- c(zd[-stop_time], zscore)
        
        if(stop_direction == "fut"){
          zb[stop_time] <- zd[stop_time]  
          
          sw_p <- pmin(sw_pvalue(side = bounds$side, info = info_ana, za = za, zb = zb,
                            zc = zc, zd = zd),1)
          
          lowci <- uniroot(sw_cilower,
                           upper = max(ul),
                           lower = min(ll),
                           conf_level = conf_level, info = info_ana,
                           za = za, zb = zb, zc = zc, zd = zd)$root
          sw.lower <- stop_sign *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
          
          upci <- uniroot(sw_ciupper,
                          upper = 10,
                          lower = 0,
                          conf_level = conf_level, info = info_ana,
                          za = za, zb = zb, zc = zc, zd = zd)$root
          sw.upper <- stop_sign *  upci * stnd_dv * info_ana$sd_proc[stop_time] 
        } else {
          
          sw_p <- sw_pvalue(side = bounds$side, info = info_ana, za = za, zb = zb,
                            zc = zc, zd = zd)
          
          lowci <- uniroot(sw_cilower,
                           upper = 20,
                           lower = -20,
                           conf_level = conf_level, info = info_ana,
                           za = za, zb = zb, zc = zc, zd = zd)$root
          sw.lower <- stop_sign *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
          
          browser()
          
          upci <- uniroot(sw_ciupper,
                          upper = 30,
                          lower = 0,
                          conf_level = conf_level, info = info_ana,
                          za = za, zb = zb, zc = zc, zd = zd)$root
          sw.upper <- stop_sign *  upci * stnd_dv * info_ana$sd_proc[stop_time] 
        }
        
        } else {
          zscore <- stop_sign*zvalues[ifelse(fixed, 1, 2), stop_time]
          zb <- ul[1:stop_time]; za <- ll[1:stop_time]; zc <- NULL; zd <- NULL
          if(stop_direction == "ul") zb <- c(zb[-stop_time], zscore)
          if(stop_direction == "fut") za <- c(za[-stop_time], zscore)
          if(length(za) == 1 & -Inf %in% za) za <- NULL
          
          if(stop_direction == "fut"){
            zb[stop_time] <- za[stop_time]
            
            sw_p <- sw_pvalue(side = bounds$side, info = info_ana, za = za, zb = zb,
                              zc = zc, zd = zd)
            
            lowci <- uniroot(sw_cilower,
                             upper = max(ul),
                             lower = ifelse(min(ll) == -Inf, -20, min(ll)),
                             conf_level = conf_level, info = info_ana,
                             za = za, zb = zb, zc = zc, zd = zd)$root
            sw.lower <- stop_sign *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
            
            upci <- uniroot(sw_ciupper,
                            upper = 10,
                            lower = 0,
                            conf_level = conf_level, info = info_ana,
                            za = za, zb = zb, zc = zc, zd = zd)$root
            sw.upper <- stop_sign *  upci * stnd_dv * info_ana$sd_proc[stop_time]   
          } else {
            sw_p <- sw_pvalue(side = bounds$side, info = info_ana, za = za, zb = zb,
                              zc = zc, zd = zd)
            
          lowci <- uniroot(sw_cilower,
                           upper = max(ul),
                           lower = ifelse(min(ll) == -Inf, -20, min(ll)),
                           conf_level = conf_level, info = info_ana,
                           za = za, zb = zb, zc = zc, zd = zd)$root
          sw.lower <- stop_sign *  lowci * stnd_dv * info_ana$sd_proc[stop_time]    
          
          upci <- uniroot(sw_ciupper,
                          upper = 15,
                          lower = 0,
                          conf_level = conf_level, info = info_ana,
                          za = za, zb = zb, zc = zc, zd = zd)$root
          sw.upper <- stop_sign *  upci * stnd_dv * info_ana$sd_proc[stop_time] 
          }
        }
    }
      
      median_unbiased_z <-
        uniroot(sw_ciupper, upper = 20, lower = -20, conf_level = 0,info = info_ana,
                za = za, zb = zb, zc = zc, zd = zd)$root
      
      median_unbiased <- stop_sign * median_unbiased_z * stnd_dv * info_ana$sd_proc[stop_time]
      
    
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
        info = info_ana, 
        stop_time = stop_time,
        stop_direction = stop_direction,
        stop_sign = stop_sign
      )
    
  } else {
    seq_inf <- list(
      median_unbiased = NULL,
      stop_time = stop_time,
      stop_direction = stop_direction,
      stop_sign = stop_sign
    )
  }
  
  n_out <- ifelse(is.null(dim(zvalues)[2]), 1, dim(zvalues)[2])
  if(length(ana_times) > 0) n_out2 <- max(ana_times)
  n_row <- ifelse((is.null(stop_time) | length(ana_times) < length(timing)), n_out+1, n_out)
  bounds$inf_frac <- sort(unique(c(org_timing,timing)))
  n_out <- which(bounds$inf_frac %in% org_timing)
  
  results <- as.data.frame(matrix(ncol = 21, nrow = n_row))
  results[n_out,1:2] <- t(zvalues)
  results[n_out,3:4] <- t(outcome_values)
  results[ana_times,10:11] <- t(naiveCI$CIfixed)
  results[ana_times,12:13] <- t(naiveCI$CIrandom)
  results[ana_times,14:15] <- t(TSAadjCI$CIfixed)
  results[ana_times,16:17] <- t(TSAadjCI$CIrandom)
  results[n_out,18:19] <- t(p_values)
  results[n_out,20:21] <- t(sd_values)

  if(length(ana_times) < length(bounds$inf_frac) & length(ana_times) > 0){ana_times <- c(ana_times, max(ana_times)+1)}

  if(length(timing) > 1){
    indi_seq <- which(bounds$inf_frac %in% timing) } else {
      indi_seq <- which(timing == bounds$inf_frac)
    }
  
  if(length(bounds$inf_frac) < n_row){
    results[indi_seq,5] <- bounds$inf_frac 
  } else { results[,5] <- bounds$inf_frac }
  results[indi_seq,6] <- bounds$alpha_ubound
  results[indi_seq,7] <- bounds$alpha_lbound
  results[indi_seq,8] <- bounds$beta_ubound
  results[indi_seq,9] <- bounds$beta_lbound


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
