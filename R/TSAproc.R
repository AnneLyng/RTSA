#' boundary
#'
#' Calculate boundaries for sequential meta-analysis. The functions are based
#' on the TSA software
#'
#' @param informationFractions Vector of information levels out of the required
#'  information. Must be in cronological order and between 0 and 1.
#' @param side 2 or 1, for 2-sided or 1-sided testing. For now is only 2-sided
#'  possible
#' @param alpha alpha level e.g. 0.05
#' @param zninf First boundary on z-scale when the first alpha spending is too
#' small based on the variable tol. Defaults to -8.
#' @param tol The tolerance level on probabilities when defining the stopping
#' boundaries. Defaults to 1e-13.
#'
#' @return A list of two elements
#' \item{InformationFractions}{The information fractions. On a scale from 0 to 1.}
#' \item{alpha.boundaries.upper}{The alpha spending boundaries. For 2-sided testing,
#'  the values are the upper boundaries.}
#' \item{alpha.boundaries.lower}{The lower boundaries. Only computed if side = 2.}
#'
#' @export
#'
#' @examples
#' timing = c(0.19, 0.5, 1)
#' boundary(informationFractions = timing, side = 2, alpha = 0.05)
#'
boundary <- function(informationFractions, side, alpha,
                     zninf = -8, tol = 1e-7){
  # set variables
  nn <- length(informationFractions); lnn <- 5000; h <- 0.05
  maxnn <- max(c(length(informationFractions), 50))

  # create variables (to be filled)
  za <- numeric(length = nn); zb <- numeric(length = nn)

  ya <- numeric(length = maxnn); yb <- numeric(length = maxnn)
  nints <- numeric(length = maxnn)

  outalpha <- alphas(alpha = alpha, side = side, ti = informationFractions, tol = tol)
  outsd <- sdfunc(ti = informationFractions)

  if(outalpha$alphaValuesDelta[1] <= 0 || outalpha$alphaValuesDelta[1] > alpha){
    outalpha$alphaValuesDelta[1] <- pmin(alpha, outalpha$alphaValuesDelta[1])
    outalpha$alphaValuesDelta[1] <- pmax(0, outalpha$alphaValuesDelta[1])
  }

  if(outalpha$alphaValuesDelta[1] == 0){
    zb[1] <- -zninf
    yb[1] <- zb[1]*outsd$sdincr[1]
  } else if(outalpha$alphaValuesDelta[1] == alpha){ # should this not be alpha?
    zb[1] <- 0
    yb[1] <- zb[1]*outsd$sdincr[1]
  } else{
    zb[1] <- qnorm(1-outalpha$alphaValuesDelta[1]/side)
    yb[1] <- zb[1] * outsd$sdincr[1]
  }

  if(side == 1){
    za[1] <- zninf
    ya[1] <- za[1]* outsd$sdincr[1]
  } else {
    za[1] <- -zb[1]
    ya[1] <- -yb[1]
  }

  nints[1] <- round((abs(yb[1] - ya[1])/h*outsd$sdincr[1])) + 1

  for(i in 2:nn){
    if(i == 2){
      last <- first(ya = ya[1], yb = yb[1], h = h, stdv = outsd$sdincr[1],
                    nints = nints, delta = 0, lnn = lnn)
    }
    if(outalpha$alphaValuesDelta[i] <= 0 || outalpha$alphaValuesDelta[i] >= 1){
      outalpha$alphaValuesDelta[i] <- min(c(1, outalpha$alphaValuesDelta[i]))
      outalpha$alphaValuesDelta[i] <- max(c(0, outalpha$alphaValuesDelta[i]))
    }
    if(outalpha$alphaValuesDelta[i] < tol){
      zb[i] <- -zninf
      yb[i] <- zb[i]*outsd$sdincr[i]
    } else if(outalpha$alphaValuesDelta[i] == 1){
      zb[i] <- 0
      yb[i] <- zb[i]*outsd$sdincr[i]
    } else {
      yb <- searchfunc(last = last, nints = nints, i = i,
                       valSF = outalpha$alphaValuesDelta[i]/side,
                       stdv = outsd$sdincr[i], ya = ya, yb = yb)
      zb[i] <- yb[i]/outsd$sdproc[i]
      if(side == 1){
        ya[i] = zninf*outsd$sdproc[i] # should this be incr?
        za[i] = zninf
      } else {
        ya[i] = -yb[i]
        za[i] = -zb[i]
      }
    }
    nints[i] <- round((yb[i]-ya[i])/h*outsd$sdincr[i])+1

    if(i != nn){
      last <- other(ya = ya, yb = yb, i = i, stdv = outsd$sdincr[i], h = h,
                    last = last, nints = nints)
    }
  }

  alpha.boundaries.upper <- zb
  if(side == 2){
    alpha.boundaries.lower <- -zb
  } else {
    alpha.boundaries.lower <-NULL
  }
  return(list(informationFractions = informationFractions,
              alpha.boundaries.upper = alpha.boundaries.upper,
              alpha.boundaries.lower = alpha.boundaries.lower))
}
#' @export

#' fcab
#'
#' Test
#'
fcab <- function(last, nint, yam1, ybm1, h, x,
                 stdv, delta){
  nlim <- 5000
  f <- numeric(length = nlim)

  for(fcon in 1:(nint+1)){
    grid2 <- yam1 + h*(fcon-1)
    f[fcon] <- last[fcon] * gfunc((x-grid2)/stdv, delta = delta)/stdv
  }
  return(trap(f = f, n = nint, h = h))
}

other <- function(ya, yb, i, stdv, h, last, nints){
  nlim <- 5000
  fn <- numeric(length = nlim)
  hh <- (yb[i] - ya[i])/nints[i]
  hlast <- (yb[i-1] - ya[i-1])/nints[i-1]

  for(j in 1:(nints[i]+1)){
    grid1 <- ya[i] + hh*(j-1) # renamed to grid1 from grid
    fn[j] <- fcab(last = last, nint = nints[i-1], yam1 = ya[i-1], ybm1 = yb[i-1],
                  h = hlast, x = grid1, stdv = stdv, delta = 0)
  }
  last <- fn
  return(last)
}

trap <- function(f, n, h){
  sum1 = f[1] # rename from sum to sum1

  for(j in 2:n){
    sum1 <- sum1 + 2*f[j]
  }
  sum1 <- sum1 + f[n + 1]

  return(h/2*sum1)
}

qpos <- function(xq, last, nint, yam1, ybm1, stdv){
  nwork <- 5000
  fun1 <- numeric(length = nwork)
  hlast <- (ybm1 - yam1)/nint

  for(j in 1:(nint+1)){
    grid3 <- yam1 + hlast * (j-1)
    fun1[j] <- last[j]*(1 - pnorm((xq - grid3)/stdv, mean= 0, sd = 1))
  }
  return(trap(f = fun1, n = nint, h = hlast))
}

gfunc <- function(x, delta){
  exp(-0.5*(x - delta)^2)/sqrt(2*pi)
}

first <- function(ya, yb, h, stdv, nints, delta, lnn){
  hh <- (yb - ya)/nints[1]

  last <- numeric(length = lnn)

  for(j in 1:(nints[1])+1){
    last[j] <- gfunc(x = (ya+hh*(j-1))/stdv, delta = delta)/stdv
  }
  return(last)
}


alphas <- function(alpha, side, ti, tol){
  alphaSpendCum <- numeric(length(ti))
  alphaSpendDelta <- numeric(length(ti))

  for(i in 1:length(ti)){
    alphaSpendCum[i] <- 2*(1 - pnorm(qnorm(1-alpha/side/2)/sqrt(ti[i]), mean = 0, sd = 1))
    alphaSpendCum[i] <- side*alphaSpendCum[i]
    if(i == 1){
      alphaSpendDelta[i] <- alphaSpendCum[i]
    } else{
      alphaSpendDelta[i] <- alphaSpendCum[i] - alphaSpendCum[i-1]
    }
    if(alphaSpendDelta[i] < tol){
      alphaSpendDelta[i] <- 0 }
  }
  return(list(alphaValuesCum = alphaSpendCum, alphaValuesDelta = alphaSpendDelta))
}


betas_Obf <- function( nn, beta, informationFractions, tol = 1e-13){
  pn_betaSpend <- numeric(length(informationFractions))
  pn_betaSpendDelta <- numeric(length(informationFractions))

  for(i in 1:length(informationFractions)){
    pn_betaSpend[i] <- 2*(1 - pnorm(qnorm(1-beta/2)/sqrt(informationFractions[i]), mean = 0, sd = 1))
    if(i == 1){
      pn_betaSpendDelta[i] <- pn_betaSpend[i]
    } else {
      pn_betaSpendDelta[i] <- pn_betaSpend[i] - pn_betaSpend[i-1]
    }
    if(pn_betaSpendDelta[i] < tol){
      pn_betaSpendDelta[i] <- 0 }
  }
  return(list(betaValuesCumulated = pn_betaSpend, betaValuesDelta = pn_betaSpendDelta))
}

betas_An <- function( nn, beta, informationFractions, tol = tol){
  pn_betaSpend <- numeric(length(informationFractions))
  pn_betaSpendDelta <- numeric(length(informationFractions))

  for(i in 1:length(informationFractions)){
    pn_betaSpend[i] <- 2*(1 - pnorm(qnorm(1-beta/2)/sqrt(informationFractions[i]), mean = 0, sd = 1))
    if(i == 1){
      pn_betaSpendDelta[i] <- pn_betaSpend[i]
    } else {
      pn_betaSpendDelta[i] <- pn_betaSpend[i] - pn_betaSpend[i-1]
    }
    if(pn_betaSpendDelta[i] < tol){
      pn_betaSpendDelta[i] <- 0 }
  }
  return(list(betaValuesCumulated = pn_betaSpend, betaValuesDelta = pn_betaSpendDelta))
}


sdfunc <- function(ti){
  sdincr <- numeric(length = length(ti))
  sdproc <- numeric(length = length(ti))

  for(ticon in 1:length(ti)){
    if(ticon == 1){
      sdincr[1] <- sqrt(ti[1])
      sdproc[1] <- sdincr[1]
    } else {
      sdincr[ticon] <- sqrt(ti[ticon]-ti[ticon-1])
      sdproc[ticon] <- sqrt(ti[ticon])
    }
  }
  return(list(sdincr = sdincr, sdproc = sdproc))
}

searchfunc <- function(last, nints, i, valSF, stdv, ya, yb){
  maxnn <- max(c(length(nints), 50)); upper <- yb[i - 1]
  del <- 10
  eps = 1e-7
  qout <- qpos(xq = upper, last = last, nint = nints[i-1], yam1 = ya[i-1],
               ybm1 = yb[i-1], stdv = stdv)

  cond <- TRUE

  while(cond){
    if(abs(qout-valSF) <= eps){
      yb[i] = upper
      cond <- FALSE
      break
    }
    if(qout > valSF + eps){
      del <- del/10
      for(k in 1:maxnn){
        upper <- upper + del
        qout <- qpos(xq = upper, last = last, nint = nints[i-1], yam1 = ya[i-1],
                     ybm1 = yb[i-1], stdv = stdv)
        if( qout <= valSF + eps){
          break
        }
      }
    }
    if(qout < valSF - eps){
      del <- del/10
      for(k in 1:maxnn){
        upper <- upper - del
        qout <- qpos(xq = upper, last = last, nint = nints[i-1], yam1 = ya[i-1],
                     ybm1 = yb[i-1], stdv = stdv)
        if( qout >= valSF - eps){
          break
        }
      }
    }
  }
  return(yb)
}

getInnerWedge <- function(informationFractions, beta, bsInf, delta, side, fakeIFY,
                          zninf = -8, tol = tol){
  maxnn <- 50; lnn <- 5000; h <- 0.05
  nints <- numeric(length = maxnn)

  if(any(informationFractions > 1)){
    en <- informationFractions[informationFractions <= 1]
  } else {
    en <- informationFractions
  }

  nn <- length(informationFractions)
  za <- numeric(length = nn); zb <- numeric(length = nn)
  stdv <- numeric(length = nn)
  #betaValuesCumulated <- numeric(length = maxnn); betaValuesDelta <- numeric(length = maxnn)
  ya <- numeric(length = maxnn); yb <- numeric(length = maxnn)
  last <- numeric(length = lnn)

  outbeta = betas_Obf( nn, beta, informationFractions = informationFractions)
  outsd = sdfunc(ti = informationFractions)

  if(outbeta$betaValuesDelta[1] <= 0 || outbeta$betaValuesDelta[1] >= beta){
    outbeta$betaValuesDelta[1] <- pmin(beta, outbeta$betaValuesDelta[1])
    outbeta$betaValuesDelta[1] <- pmax(0, outbeta$betaValuesDelta[1])
  }

  if(outbeta$betaValuesDelta[1] == 0){
    za[1] <- zninf
    ya[1] <- za[1]*outsd$sdincr[1]
  } else if(outbeta$betaValuesDelta[1] == beta){
    za[1] <- 0
    ya[1] <- za[1]*outsd$sdincr[1]
  } else{
    zb[1] <- zninf
    yb[1] <- zb[1]*outsd$sdincr[1]
    za[1] <- -zb[1]
    ya[1] <- -yb[1]
    nints[1] <- round((abs(yb[1] - ya[1])/h*outsd$sdincr[1])) + 1
    za[1] <- qnorm(outbeta$betaValuesDelta[1])
    ya[1] <- za[1] * outsd$sdincr[1]
  }

  if(side == 1){
    zb[1] <- -zninf
    yb[1] <- zb[1]* outsd$sdincr[1]
  } else {
    zb[1] <- -za[1]
    yb[1] <- -yb[1]
  }

  nints[1] <- round((abs(yb[1] - ya[1])/h*outsd$sdincr[1])) + 1
  for(i in 2:nn){
    if(i == 2){
      last <- first(ya = ya[1], yb = yb[1], h = h, stdv = outsd$sdincr[1],
                    nints =nints, delta = 0, lnn = lnn)
    }
    if(outbeta$betaValuesDelta[i] <= 0 || outbeta$betaValuesDelta[i] >= 1){
      outbeta$betaValuesDelta[i] <- min(c(1, outbeta$betaValuesDelta[i]))
      outbeta$betaValuesDelta[i] <- max(c(0, outbeta$betaValuesDelta[i]))
    }
    if(outbeta$betaValuesDelta[i] < tol){
      za[i] <- zninf
      ya[i] <- za[i]*outsd$sdincr[i]
    } else if(outbeta$betaValuesDelta[i] == 1){
      za[i] <- 0
      ya[i] <- za[i]*outsd$sdincr[i]
    } else {
      yb <- searchfunc(last = last, nints = nints,
                       i = i, valSF = outbeta$betaValuesDelta[i],
                       stdv = outsd$sdincr[i], ya = ya, yb = yb)
      zb[i] <- yb[i]/outsd$sdproc[i]
      if(side == 1){
        ya[i] = zninf*outsd$sdproc[i]
        za[i] = zninf
      } else {
        ya[i] = -yb[i]
        za[i] = -zb[i]
      }
    }
    nints[i] <- round(abs(yb[i]-ya[i])/h*outsd$sdincr[i])+1
    if(i != nn){
      last <- other(ya, yb, i, outsd$sdincr[i], h, last, nints)
    }
  }
  testDrift <- fakeIFY + abs(ya[length(informationFractions)])
  ret1 <- informationFractions
  ret2 <- za + outsd$sdproc*testDrift
  return(list(ret1 = ret1, ret2 = ret2, betaValuesDelta = outbeta$betaValuesDelta, za = za, zb = zb,
              ya = ya, yb = yb, drift = testDrift))
}
