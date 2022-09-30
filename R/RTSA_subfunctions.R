# nFixed ----
# Calculate RIS with fixed effect.
#' @importFrom stats qnorm
nFixed <- function(alpha, beta, pI = NULL, pC = NULL, iE = NULL,
                   sdE = NULL, binary = TRUE){
  if(binary == TRUE){
    p <- (pC + pI)/2
    v <- p*(1-p)
    theta <- pC-pI
    return(2*(qnorm(1-alpha/2)+qnorm(1-beta))^2*2*v/theta^2)
  } else {
    return(2*(qnorm(1-alpha/2)+qnorm(1-beta))^2*2*sdE^2/iE^2)
  }
}


# nRandom ----
# Calculate RIS with random effect using diversity
nRandom <- function(alpha, beta, pI, pC, diversity = NULL){

  p <- (pC + pI)/2
  v <- p*(1-p)
  theta <- pC-pI
  NF <- 2*(qnorm(1-alpha/2)+qnorm(1-beta))^2*2*v/theta^2

  NR <- 1/(1-diversity)*NF
  return(NR)
}

# ma_power ----
# Calculate power and type-1-error rate
ma_power <- function(upper_bound, theta, info, za = rep(-20, length(upper_bound))){
  nints <- numeric(length(upper_bound))
  p_cross <- 1-pnorm(upper_bound[1], mean = theta*info$sd_proc[1], sd = 1) # the first one

  if(length(upper_bound) > 1){
    lnn <- 5000
    nints[1] <- round((abs(za[1] - upper_bound[1])/0.05*info$sd_incr[1])) + 1 # how many splits of the integral
    last <- first(za = za[1], zb = upper_bound[1], h = 0.05, stdv = info$sd_incr[1],
                         nints = nints, delta = theta, lnn = lnn, bs = TRUE)
    p_cross <- c(p_cross,qpos(xq = upper_bound[2]*info$sd_proc[2], last = last, i = 2, nint = nints[1], zam1 = za[1],
                                     zbm1 = upper_bound[1], stdv = info, bs = FALSE, delta = theta))
    if(length(upper_bound) > 2){
    for(i in 2:(length(upper_bound)-1)){
      nints[i] <- round(abs(upper_bound[i]-za[i])/0.05*info$sd_incr[i])+1
      last <- other(za[1:(i+1)], upper_bound[1:(i+1)], i, stdv = info, h = 0.05, last = last, nints, delta = theta, bs = TRUE)
      p_cross <- c(p_cross, qpos(xq = upper_bound[i+1]*info$sd_proc[i+1], last = last, i = (i+1), nint = nints[i], zam1 = za[i],
                                        zbm1 = upper_bound[i], stdv = info, bs = FALSE, delta = theta))
    }
    }
  }
  return(list(p_cross, sum(p_cross)))
}

# boundary ----
# Calculate boundaries for sequential meta-analysis.
boundary <- function(inf_frac, side, alpha, beta,
                     zninf = -20, tol = 1e-13, delta = NULL, bs = FALSE,
                     tol_alpha = 1e-16){
  # set nuissance variables
  # TODO: see if next 6 times can be moved to other function or to more
  # meaningfull section

  if(is.null(delta)) delta = abs(qnorm(alpha/side)+qnorm(beta))

  nn <- length(inf_frac); lnn <- 5000; h <- 0.05
  maxnn <- max(c(length(inf_frac), 50))

  # create variables (to be filled)
  za <- numeric(length = nn); zb <- numeric(length = nn)
  ya <- numeric(length = maxnn); yb <- numeric(length = maxnn)
  nints <- numeric(length = maxnn)

  # calculate the alpha spending
  alpha_spend <- obf_as(alpha = alpha, side = side, ti = inf_frac)
  # change in information
  info <- sd_inf(ti = inf_frac)

  # set criteria for first spending (not under 0, not over alpha)
  alpha_spend$as_incr[1] <- pmin(alpha, alpha_spend$as_incr[1])
  alpha_spend$as_incr[1] <- pmax(0, alpha_spend$as_incr[1])

  if(alpha_spend$as_incr[1] < tol_alpha){
    zb[1] <- -zninf
  } else if(alpha_spend$as_incr[1] == alpha){
    zb[1] <- 0 # TODO: should this not be alpha?
  } else {
    zb[1] <- qnorm(alpha_spend$as_incr[1]/side, lower.tail = FALSE)
  }

  yb[1] <- zb[1]*info$sd_incr[1]

  if(side == 1){
    za[1] <- zninf
    ya[1] <- za[1]* info$sd_incr[1]
  } else {
    za[1] <- -zb[1]
    ya[1] <- -yb[1]
  }

  # why this definition? (what is h and why this formula)
  nints[1] <- round((abs(zb[1] - za[1])/h*info$sd_incr[1])) + 1

  for(i in 2:nn){
    if(i == 2){
      last <- first(za = za[1], zb = zb[1], h = h, stdv = info$sd_incr[1],
                    nints = nints, delta = 0, lnn = lnn, bs = bs)
    }
    if(alpha_spend$as_incr[i] <= 0 || alpha_spend$as_incr[i] >= 1){
      alpha_spend$as_incr[i] <- min(c(1, alpha_spend$as_incr[i]))
      alpha_spend$as_incr[i] <- max(c(0, alpha_spend$as_incr[i]))
    }
    if(alpha_spend$as_incr[i] < tol_alpha){
      zb[i] <- -zninf
      yb[i] <- zb[i]*info$sd_incr[i]
    } else if(alpha_spend$as_incr[i] == 1){
      zb[i] <- 0
      yb[i] <- zb[i]*info$sd_incr[i]
    } else {
      zb[i] <- searchfunc(last = last, nints = nints, i = i,
                          as = alpha_spend$as_incr[i]/side,
                          stdv = info, za = za, zb = zb, bs = bs, tol = tol, delta = 0) # delta = delta
      yb[i] <- zb[i]*info$sd_proc[i]
      if(side == 1){
        ya[i] = zninf*info$sd_proc[i] # should this be incr?
        za[i] = zninf
      } else {
        ya[i] = -yb[i]
        za[i] = -zb[i]
      }
    }
    nints[i] <- round((zb[i]-za[i])/h*info$sd_incr[i])+1

    if(i != nn){
      last <- other(za = za, zb = zb, i = i, stdv = info, h = h,
                    last = last, nints = nints, delta = 0, bs = bs) # delta = delta
    }
  }

  alpha_ubound <- zb
  if(side == 2){
    alpha_lbound <- -zb
  } else {
    alpha_lbound <-NULL
  }
  return(list(inf_frac = inf_frac,
              alpha_ubound = alpha_ubound,
              alpha_lbound = alpha_lbound,
              alpha = alpha,
              alpha_spend = alpha_spend,
              delta = delta))
}

# fcap ----
# Quantify where boundaries will lie
fcab <- function(last, nint, zam1, zbm1, h, x,
                 stdv, delta, i, bs){

  nlim <- 5000
  f <- numeric(length = nlim)

  for(fcon in 1:(nint+1)){
    grid2 <- (zam1 + h*(fcon-1))*stdv$sd_proc[i-1]
    if(bs == TRUE){
      f[fcon] <- last[fcon] *stdv$sd_proc[i]/stdv$sd_incr[i] *
        dnorm((x-grid2)/stdv$sd_incr[i], mean = delta*stdv$sd_incr[i])
    } else {
      f[fcon] <- last[fcon] *stdv$sd_proc[i]/stdv$sd_incr[i] *
        dnorm((grid2-x)/stdv$sd_incr[i], mean = delta*stdv$sd_incr[i]) # check if should be delta*stdv
    }
  }
  return(trap(f = f, n = nint, h = h))
}

#other ----
other <- function(za, zb, i, stdv, h, last, nints, delta, bs){
  nlim <- 5000
  fn <- numeric(length = nlim)
  hh <- (zb[i] - za[i])/nints[i]
  hlast <- (zb[i-1] - za[i-1])/nints[i-1]

  for(j in 1:(nints[i]+1)){
    grid1 <- (za[i] + hh*(j-1))*stdv$sd_proc[i] # renamed to grid1 from grid
    fn[j] <- fcab(last = last, nint = nints[i-1], zam1 = za[i-1], zbm1 = zb[i-1],
                  h = hlast, x = grid1, stdv = stdv, delta = delta, i = i,
                  bs = bs)
  }
  last <- fn
  return(last)
}

# trap ----
trap <- function(f, n, h){
  sum1 = f[1] # rename from sum to sum1

  for(j in 2:n){
    sum1 <- sum1 + 2*f[j]
  }
  sum1 <- sum1 + f[n + 1]

  return(h/2*sum1)
}

# qpos ----
# reverse integrals
#' @importFrom stats pnorm
qpos <- function(xq, last, nint, i, zam1, zbm1, stdv, bs = bs, delta){
  nwork <- 5000
  fun1 <- numeric(length = nwork)
  hlast <- (zbm1 - zam1)/nint

  for(j in 1:(nint+1)){
    grid3 <- (zam1 + hlast * (j-1))*stdv$sd_proc[i-1]
    if(bs == FALSE & delta != 0) {
      fun1[j] <- last[j]*(pnorm((grid3-xq)/stdv$sd_incr[i], mean= -delta*stdv$sd_incr[i], sd = 1,
                                lower.tail = TRUE))
    }
    else if(bs == TRUE){
      fun1[j] <- last[j]*(pnorm((xq-grid3)/stdv$sd_incr[i], mean= delta*stdv$sd_incr[i], sd = 1,
                                lower.tail = TRUE))
    } else {
      fun1[j] <- last[j]*(pnorm((grid3-xq)/stdv$sd_incr[i], mean= delta*stdv$sd_incr[i], sd = 1, lower.tail = TRUE))
    }
  }
  return(trap(f = fun1, n = nint, h = hlast))
}

# first ----
# focus on the first trial
first <- function(za, zb, h, stdv, nints, delta, lnn, bs = bs){
  hh <- (zb - za)/nints[1]

  last <- numeric(length = lnn)

  for(j in 1:(nints[1])+1){
    #last[j] <- gfunc(x = (ya+hh*(j-1))/stdv, delta = delta)/stdv #original
    last[j] <- dnorm(x = (za+hh*(j-1)))
    if(bs == TRUE) {last[j] <- dnorm(x = (za+hh*(j-1)), mean = delta*stdv)}
  }
  return(last)
}

# obf_as ----
# Alpha spending boundaries
# original: alphas
#' @importFrom stats qnorm pnorm
obf_as <- function(alpha, side, ti){
  as_cum <- numeric(length(ti))
  as_incr <- numeric(length(ti))

  for(i in 1:length(ti)){
    as_cum[i] <- 2*(1 - pnorm(qnorm(1-alpha/side/2)/sqrt(ti[i]), mean = 0, sd = 1))
    as_cum[i] <- side*as_cum[i]
    if(i == 1){
      as_incr[i] <- as_cum[i]
    } else{
      as_incr[i] <- as_cum[i] - as_cum[i-1]
    }
  }
  return(list(as_cum = as_cum, as_incr = as_incr))
}

#betas_Obf ----
#betas_Obf <- function(beta, side, inf_frac){
#  bs_cum <- numeric(length(inf_frac))
#  bs_incr <- numeric(length(inf_frac))

#  for(i in 1:length(inf_frac)){
#    pn_betaSpend[i] <- 2*(1 - pnorm(qnorm(1-beta/2)/sqrt(inf_frac[i]), mean = 0, sd = 1))
#    if(i == 1){
#      pn_betaSpendDelta[i] <- pn_betaSpend[i]
#    } else {
#      pn_betaSpendDelta[i] <- pn_betaSpend[i] - pn_betaSpend[i-1]
#    }
#    if(pn_betaSpendDelta[i] < tol){
#      pn_betaSpendDelta[i] <- 0 }
#  }
#  return(list(as_cum = pn_betaSpend, as_incr = pn_betaSpendDelta))
#}

# sd_inf ----
# Ratio of how much information is available
# original: sdfunc
sd_inf <- function(ti){
  sd_incr <- sqrt(ti - c(0, ti[-length(ti)]))
  sd_proc <- sqrt(ti)

  return(list(sd_incr = sd_incr, sd_proc = sd_proc))
}

#searchfunc ----
# Readjusting RIS based on the number of interim-analyses
# TODO: Talk with CG, should it be implemented
searchfunc <- function(last, nints, i, as, stdv, za, zb, tol, bs, delta){
  maxnn <- max(c(length(nints), 50)); upper <- zb[i - 1]*stdv$sd_proc[i]
  if(bs == TRUE) upper <- za[i - 1]*stdv$sd_proc[i]
  del <- 10
  qout <- qpos(xq = upper, last = last, i = i, nint = nints[i-1], zam1 = za[i-1],
               zbm1 = zb[i-1], stdv = stdv, bs = bs, delta = delta)

  cond <- TRUE

  while(cond){ # what does this do?
    if(abs(qout-as) <= tol){
      cond <- FALSE
      break
    }
    if(qout > as + tol){
      del <- del/10
      for(k in 1:maxnn){
        if(bs == TRUE) { upper <- upper - 2*del }
        upper <- upper + del
        qout <- qpos(xq = upper, last = last, i = i, nint = nints[i-1], zam1 = za[i-1],
                     zbm1 = zb[i-1], stdv = stdv, bs = bs, delta = delta)
        if( qout <= as + tol){
          break
        }
      }
    }
    if(qout < as - tol){
      del <- del/10
      for(k in 1:maxnn){
        if(bs == TRUE) { upper <- upper + 2*del }
        upper <- upper - del
        qout <- qpos(xq = upper, last = last, i = i, nint = nints[i-1], zam1 = za[i-1],
                     zbm1 = zb[i-1], stdv = stdv, bs = bs, delta = delta)
        if( qout >= as - tol){
          break
        }
      }
    }
  }
  return(upper/stdv$sd_proc[i])
}

#getInnerWedge ----
# TODO: To be implemented
#' @importFrom stats qnorm
getInnerWedge <- function(inf_frac, beta, delta = NULL, side, fakeIFY = 0,
                          zninf = -20, tol = tol, outer_boundaries, rm_bs = NULL){
  maxnn <- 50; lnn <- 5000; h <- 0.01
  nints <- numeric(length = maxnn)

  if(is.null(delta)) delta = abs(qnorm(outer_boundaries$alpha/side)+qnorm(beta))
  outer_bound <- outer_boundaries$alpha_ubound

  #if(any(inf_frac > 1)){
  #  en <- inf_frac[inf_frac <= 1]
  #} else {
  #  en <- inf_frac
  #}

  nn <- length(inf_frac)
  za <- numeric(length = nn); zb <- numeric(length = nn)
  stdv <- numeric(length = nn)
  ya <- numeric(length = maxnn); yb <- numeric(length = maxnn)
  last <- numeric(length = lnn)

  #beta_spend = RTSA:::betas_Obf(nn, beta, inf_frac = inf_frac/inf_frac[nn])
  beta_spend = obf_as(beta, side, inf_frac/inf_frac[nn])
  if(!is.null(rm_bs)){
    beta_spend = obf_as(beta, side, c(rep(0, max(rm_bs)),inf_frac[-c(1:rm_bs)]/inf_frac[nn]))
  }

  if(side == 2){
    beta_spend$as_cum <- beta_spend$as_cum/2
    beta_spend$as_incr <- beta_spend$as_incr/2
  }

  info = sd_inf(ti = inf_frac)

  if(beta_spend$as_incr[1] <= 0 || beta_spend$as_incr[1] >= beta){
    beta_spend$as_incr[1] <- pmin(beta, beta_spend$as_incr[1])
    beta_spend$as_incr[1] <- pmax(0, beta_spend$as_incr[1])
  }

  if(beta_spend$as_incr[1] == 0){
    za[1] <- zninf
    ya[1] <- za[1]*info$sd_incr[1]
  } else if(beta_spend$as_incr[1] == beta){
    za[1] <- 0
    ya[1] <- za[1]*info$sd_incr[1]
  } else{
    za[1] <- qnorm(beta_spend$as_incr[1]/side, mean = info$sd_proc[1]*delta, sd = 1)
    ya[1] <- za[1] * info$sd_incr[1]
  }

  zb <- outer_bound
  yb <- zb*info$sd_incr

  #if(side == 1){
  #  zb[1] <- outer_bound[1]
  #  yb[1] <- zb[1]* info$sd_incr[1]
  #} else {
  #  zb[1] <- -za[1]
  #  yb[1] <- -ya[1]
  #}

  nints[1] <- round((abs(za[1] - zb[1])/h*info$sd_incr[1])) + 1 # how many splits of the integral
  for(i in 2:nn){
    if(i == 2){
      last <- first(za = za[1], zb = zb[1], h = h, stdv = info$sd_incr[1],
                    nints = nints, delta = delta, lnn = lnn, bs = TRUE)
    }
    if(beta_spend$as_incr[i] <= 0 || beta_spend$as_incr[i] >= 1){
      beta_spend$as_incr[i] <- min(c(1, beta_spend$as_incr[i]))
      beta_spend$as_incr[i] <- max(c(0, beta_spend$as_incr[i]))
    }
    if(beta_spend$as_incr[i] < tol){
      za[i] <- zninf
      ya[i] <- za[i]*info$sd_incr[i]
    } else if(beta_spend$as_incr[i] == beta){
      za[i] <- 0
      ya[i] <- za[i]*info$sd_incr[i]
    } else {
      za[i] <- searchfunc(last = last, nints = nints,
                          i = i, as = beta_spend$as_incr[i],
                          stdv = info, za = za, zb = zb, tol = tol, bs =TRUE,
                          delta = delta)
      ya[i] <- za[i]*info$sd_proc[i]
      #if(side == 1){
      #  zb[i] = outer_bound[i]
      #  yb[i] = zb[i]*info$sd_proc[i]
      #} else {
      #  yb[i] = -ya[i]
      #  zb[i] = -za[i]
      #}
    }
    nints[i] <- round(abs(zb[i]-za[i])/h*info$sd_incr[i])+1
    if(i != nn){
      last <- other(za, zb, i, stdv = info, h, last, nints, delta = delta, bs = TRUE)
    }
  }
  testDrift <- fakeIFY + abs(ya[length(inf_frac)])
  ret1 <- inf_frac
  ret2 <- info$sd_proc*delta
  return(list(ret1 = ret1, ret2 = ret2, as_incr = beta_spend$as_incr,
              as_cum = beta_spend$as_cum, za = za, zb = zb,
              ya = ya, yb = yb, drift = testDrift, info = info))
}



# TSA ----
# Gathers all of the information in a nested list for TSA
TSA = function(timing,
               ana_time,
               synth,
               side,
               alpha,
               beta,
               futility,
               mc,
               stopTime = NULL,
               confInt = TRUE,
               subjects,
               RIS,
               hakn,
               sign,
               fixedStudy,
               hksj,
               tau.ci.method) {

  # if analysis times are not specified, analysis is performed at all timings.
  if(is.null(ana_time)){
    ana_time = 1:length(timing)
  }

  # reduce analysis to analysis times
  timing = timing[ana_time]

  timing_incr <- timing - c(0, timing[-length(timing)])
  trials <- cbind(timing, timing_incr)

  # if new study adds less than 1% of RIS, the analysis is not performed.
  time_tf = 0.01
  trials <- trials[trials[, 2] > time_tf, ]
  # ana_time = ana_time[which(trials[, 2] > time_tf)] # TODO: check if redundant

  trials[, 2] <- trials[, 1] - c(0, trials[, 1][-length(trials[, 1])])

  # calculate the boundaries
  boundout = RTSA:::boundary(inf_frac = trials[, 1],
                      side = side,
                      alpha = alpha, beta = beta)

  # calculate the power of the analysis
  if(side == 1 & is.null(futility)){
    right_power <- function(x, inner, outer, theta){
      info = sd_inf(trials[, 1]*x)
      pwr <- ma_power(upper_bound = outer, theta = theta,
                             info = info, za = inner)[[2]]
      pwr - (1-beta)
    }

    lb <- rep(-20,length(trials[, 1]))

    root <- uniroot(right_power, lower = 0.9, upper = 1.2, tol = 1e-9,
                    inner = lb, outer = boundout$alpha_ubound, theta =
                      boundout$delta)$root

    info <- sd_inf(trials[, 1]*root)
    pwr <- ma_power(upper_bound = boundout$alpha_ubound, theta = boundout$delta,
                    info = info, za = lb)
    t1e <- ma_power(upper_bound = boundout$alpha_ubound, theta = 0,
                    info = info, za = lb)
  } else if(side == 1 & futility == "non-binding"){
    lb <- getInnerWedge(inf_frac = trials[, 1], beta = beta,
                               side = 1, fakeIFY = 0,zninf = -20, tol = 1e-13,
                               outer_boundaries = boundout,
                               delta = NULL)

    getWarp <- function(x, outer){
      a <- getInnerWedge(inf_frac = trials[, 1]*x, beta = beta,
                                side = 1, zninf = -20, tol = 1e-13,
                                outer_boundaries = outer,
                                delta = NULL)$za[length(trials[, 1])]
      outer$alpha_ubound[length(trials[, 1])] - a
    }

    root <- uniroot(getWarp, lower = 0.9, upper = 1.2, outer = boundout,
                    tol = 1e-9)$root

    lb <- getInnerWedge(inf_frac = trials[,1]*root, beta = beta,
                               side = 1, zninf = -20, tol = 1e-13,
                               outer_boundaries = boundout,
                               delta = NULL)

    right_power <- function(x, inner, outer, theta){
      info <- sd_inf(trials[,1]*x)
      pwr <- ma_power(upper_bound = outer$alpha_ubound, theta = theta,
                             info = info, za = inner)[[2]]
      pwr - (1-beta)
    }

    root <- uniroot(right_power, lower = 0.9, upper = 1.2, tol = 1e-9,
                    inner = lb$za, outer = boundout, theta = boundout$delta)$root

    info <- sd_inf(trials[,1]*root)
    pwr <- ma_power(upper_bound = boundout$alpha_ubound, theta = boundout$delta,
                    info = info, za = lb$za)
    t1e <- ma_power(upper_bound = boundout$alpha_ubound, theta = 0,
                    info = info, za = lb$za)
    boundout$alpha_lbound <- lb$za
  }


  # calculate the cum. z-score (do we want this per study?)
  zout = lapply(ana_time[ana_time <= dim(synth$data)[1]],
                function(x) {
                  synout = synthesize(
                    metaPrepare(
                      data = synth$data[1:x,],
                      outcome = synth$outcome,
                      method = synth$method,
                      alpha = alpha
                    ), sign = sign,
                    fixedStudy = fixedStudy,
                    hksj = hksj,
                    tau.ci.method = tau.ci.method
                  )
                  return(synout)
                })

  names(zout) = ana_time[ana_time <= dim(synth$data)[1]]
  if("1" %in% names(zout) &
     length(grep("peR", names(unlist(zout)))) > 0){
    zout[["1"]]$peR <- c(0, 0, 0, zout[["1"]]$peF[4])
  }

  zvalues = sapply(names(zout), function(x){c(zout[[x]]$peF[4], zout[[x]]$peR[4])})

  if (confInt == TRUE) {
    if(is.null(stopTime)){ stopTime =
      as.character(max(ana_time[ana_time <= dim(synth$data)[1]]))}
    naiveCI = list(CIfixed = zout[[stopTime]]$peF[c(2, 3)],
                   CIrandom = zout[[stopTime]]$peR[c(2, 3)])
    adjCI = list(
      CIfixed = exp(
        log(zout[[stopTime]]$peF[1]) +
          c(-1, 1) * boundout$alpha_ubound[which(stopTime == ana_time)] *
          sqrt(zout[[stopTime]]$peF[7])
      ),
      if(zout[[stopTime]]$U[1] > 0){ CIrandom = exp(
        log(zout[[stopTime]]$peR[1]) +
          c(-1, 1) * boundout$alpha_ubound[which(stopTime == ana_time)] *
          sqrt(zout[[stopTime]]$peR[6])
      ) }
    )
  }

  RTSAout =     list(
    alpha = alpha,
    side = side,
    boundout = boundout,
    zout = zout[stopTime],
    zvalues = zvalues,
    ana_time = ana_time,
    stopTime = stopTime,
    naiveCI = naiveCI,
    adjCI = adjCI,
    pwr = pwr,
    t1e = t1e,
    futility = futility,
    root = root
  )

  class(RTSAout) <- c("list", "RTSA")

  return(
    RTSAout
  )

}


