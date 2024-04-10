#' @useDynLib RTSA
#' @importFrom Rcpp sourceCpp
NULL

# Define - root power function ----
right_power <- function(x, zb, za, zc = NULL, zd = NULL, delta,
                        beta, timing){
  info <- sd_inf(timing*x)
  pwr <- ma_power(zb = zb,  za = za,
                  delta = delta,
                  info = info, zc = zc, zd = zd)[[2]]
  pwr - (1-beta)
}

# Define - letting beta and alpha boundaries meet ----
inf_warp <- function(x, alpha_boundaries, rm_bs = 0, side = 1, delta = NULL,
                     beta, timing, es_beta, design_R = NULL){
  a <- beta_boundary(inf_frac = timing, beta = beta,
                     side = side, alpha_boundaries = alpha_boundaries,
                     delta = delta, rm_bs = rm_bs,
                     es_beta = es_beta,
                     design_R = design_R, warp_root = x)$za[length(timing)]
  alpha_boundaries$alpha_ubound[length(timing)] - a
}

# Define - ma_power (meta-analysis power) ----
# Calculate power and type-1-error rate
ma_power <- function(zb, za = rep(-20, length(zb)), delta, info, zc = NULL, zd = NULL){
  # zb is upper alpha boundary
  # za is lower alpha boundary
  # zd is upper beta boundary
  # zc is lower beta boundary

  # probability of crossing the first boundary
  p_cross <- 1-pnorm(zb[1], mean = delta*info$sd_proc[1], sd = 1)

  # calculate the probability of crossing the following boundaries
  if(is.null(zc)){ # if no fut. boundaries (for two-sided designs)
    if(length(zb) > 1){
      zj_wj <- z_n_w(r = 18, info = info, za = za, zb = zb, i = 1,
                   delta = delta)
      zj <- zj_wj$zj
      wj <- zj_wj$wj

      last <- init_int(wj = wj, zj = zj, delta = delta, stdv = info$sd_incr)

      p_cross <- c(p_cross,prob(xq = zb[2]*info$sd_proc[2], last = last, zj = zj, k = 2,
                              stdv = matrix(unlist(info), ncol = 2),
                              bs = FALSE, delta = delta))
      if(length(zb) > 2){
        for(i in 2:(length(zb)-1)){

          zj_wj_up <- z_n_w(r = 18, info = info, za = za, zb = zb, i = i,
                        delta = delta)

          last <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                        zj = zj, last = last, zj_up = zj_wj_up$zj,
                        wj_up = zj_wj_up$wj, delta = delta, bs = FALSE)

          zj <- zj_wj_up$zj

      p_cross <- c(p_cross, prob(xq = zb[i+1]*info$sd_proc[i+1], last = last, k = (i+1),
                                 zj = zj, stdv = matrix(unlist(info), ncol = 2),
                                 bs = FALSE, delta = delta))
    }
    }
  }
  } else {
    fk <- min(which(!is.na(zc)))
    if(length(zb) > 1){
      if(fk == 1){
        zj_wj_u <- z_n_w(r = 18, info = info, za = zd, zb = zb, i = 1,
                              delta = delta)
        zj_u <- zj_wj_u$zj
        wj_u <- zj_wj_u$wj

        last_u <- init_int(wj = wj_u, zj = zj_u, delta = delta, stdv = info$sd_incr)

        u_prob <- prob(xq = zb[2]*info$sd_proc[2], last = last_u, zj = zj_u, k = 2,
                    stdv = matrix(unlist(info), ncol = 2),
                    bs = FALSE, delta = delta)

        zj_wj_l <- z_n_w(r = 18, info = info, za = za, zb = zc, i = 1,
                                delta = delta)
        zj_l <- zj_wj_l$zj
        wj_l <- zj_wj_l$wj

        last_l <- init_int(wj = wj_l, zj = zj_l, delta = delta, stdv = info$sd_incr)

        l_prob <- prob(xq = zb[2]*info$sd_proc[2], last = last_l, zj = zj_l, k = 2,
                              stdv = matrix(unlist(info), ncol = 2),
                              bs = FALSE, delta = delta)

        p_cross <- c(p_cross,u_prob+l_prob)
      } else {
        zj_wj <- z_n_w(r = 18, info = info, za = za, zb = zb, i = 1,
                              delta = delta)
        zj <- zj_wj$zj
        wj <- zj_wj$wj

        # first integral
        last <- init_int(wj = wj, zj = zj, delta = delta, stdv = info$sd_incr)

        p_cross <- c(p_cross,prob(xq = zb[2]*info$sd_proc[2], last = last, zj = zj, k = 2,
                                         stdv = matrix(unlist(info), ncol = 2),
                                         bs = FALSE, delta = delta))
      }
      if(length(zb) > 2){
        for(i in 2:(length(zb)-1)){
          if(fk <= i){
          zj_wj_up_u <- z_n_w(r = 18, info = info, za = zd, zb = zb, i = i,
                                   delta = delta)

          zj_wj_up_l <- z_n_w(r = 18, info = info, za = za, zb = zc, i = i,
                                     delta = delta)

          if(fk < i ){
          last_u <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                                   zj = zj_u, last = last_u, zj_up = zj_wj_up_u$zj,
                                   wj_up = zj_wj_up_u$wj, delta = delta, bs = FALSE)
          last_l <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                                     zj = zj_l, last = last_l, zj_up = zj_wj_up_l$zj,
                                     wj_up = zj_wj_up_l$wj, delta = delta, bs = FALSE)
          } else {
            last_u <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                                       zj = zj, last = last, zj_up = zj_wj_up_u$zj,
                                       wj_up = zj_wj_up_u$wj, delta = delta, bs = FALSE)
            last_l <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                                       zj = zj, last = last, zj_up = zj_wj_up_l$zj,
                                       wj_up = zj_wj_up_l$wj, delta = delta, bs = FALSE)
          }

          zj_u <- zj_wj_up_u$zj
          zj_l <- zj_wj_up_l$zj

          l_prob <- prob(xq = zb[i+1]*info$sd_proc[i+1], last = last_l, zj = zj_l, k = (i+1),
                                stdv = matrix(unlist(info), ncol = 2),
                                bs = FALSE, delta = delta)
          u_prob <- prob(xq = zb[i+1]*info$sd_proc[i+1], last = last_u, zj = zj_u, k = (i+1),
                                stdv = matrix(unlist(info), ncol = 2),
                                bs = FALSE, delta = delta)

          p_cross <- c(p_cross, l_prob+u_prob)

          } else {
            zj_wj_up <- z_n_w(r = 18, info = info, za = za, zb = zb, i = i,
                                     delta = delta)

            last <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                                     zj = zj, last = last, zj_up = zj_wj_up$zj,
                                     wj_up = zj_wj_up$wj, delta = delta, bs = FALSE)

            zj <- zj_wj_up$zj

            p_cross <- c(p_cross, prob(xq = zb[i+1]*info$sd_proc[i+1], last = last, k = (i+1),
                                              zj = zj, stdv = matrix(unlist(info), ncol = 2),
                                              bs = FALSE, delta = delta))
          }
      }
    }
    }
  }
  return(list(p_cross, sum(p_cross))) 
}

# alpha_boundary ----
# Calculate boundaries for sequential meta-analysis based on type-1-error spending.
alpha_boundary <- function(inf_frac, side, alpha, beta,
                     zninf = -20, tol = 1e-09, delta = NULL, bs = FALSE,
                     es_alpha = NULL,
                     gamma = NULL, rho = NULL, r = 18,
                     type = "design", design_R = NULL){
  # set delta if missing
  if(is.null(delta)) delta = abs(qnorm(alpha/side)+qnorm(beta))

  nn <- length(inf_frac); lnn <- 5000; h <- 0.05

  # create variables (to be filled)
  za <- numeric(length = nn); zb <- numeric(length = nn)
  ya <- numeric(length = nn); yb <- numeric(length = nn)

  org_inf_frac <- NULL
  if(type == "analysis"){
    org_inf_frac <- inf_frac
    nn <- length(inf_frac)
  }

  if(sum(inf_frac > 1) > 0){
    alpha_timing <- c(inf_frac/max(inf_frac))
    nn <- length(alpha_timing)
  } else {
    alpha_timing <- inf_frac
    nn <- length(alpha_timing)
  }

  # calculate the alpha spending
  alpha_spend <- switch(which(es_alpha == c("esOF", "esPoc", "HSDC", "rho")),
                        esOF(alpha/side, alpha_timing), esPoc(alpha/side,alpha_timing),
                        HSDC(alpha/side, alpha_timing, gamma),
                        rho(alpha/side, alpha_timing, rho))

  # change in information
  info <- sd_inf(timing= inf_frac)
  if(type == "analysis"){
    info <- sd_inf(timing= inf_frac*design_R)
  }
  
  # set criteria for first spending (not under 0, not over alpha)
  alpha_spend$as_incr[1] <- pmin(alpha, alpha_spend$as_incr[1])
  alpha_spend$as_incr[1] <- pmax(0, alpha_spend$as_incr[1])

  if(alpha_spend$as_incr[1] < tol){
    zb[1] <- -zninf
  } else {
    zb[1] <- qnorm(alpha_spend$as_incr[1], lower.tail = FALSE)
  }

  yb[1] <- zb[1]*info$sd_incr[1]

  if(side == 1){
    za[1] <- zninf
    ya[1] <- za[1]* info$sd_incr[1]
  } else {
    za[1] <- -zb[1]
    ya[1] <- -yb[1]
  }

  zj_wj <- z_n_w(r = r, info = info, za = za, zb = zb, i = 1,
                 delta = 0)
  zj <- zj_wj$zj
  wj <- zj_wj$wj

  if(nn > 1){
  for(i in 2:nn){
    if(i == 2){
      last <- init_int(wj = wj, zj = zj, delta = 0, stdv = info$sd_incr)
    }
    if(alpha_spend$as_incr[i] <= 0 || alpha_spend$as_incr[i] >= 1){
      alpha_spend$as_incr[i] <- min(c(1, alpha_spend$as_incr[i]))
      alpha_spend$as_incr[i] <- max(c(0, alpha_spend$as_incr[i]))
    }
    if(alpha_spend$as_incr[i] < tol){
      zb[i] <- -zninf
      yb[i] <- zb[i]*info$sd_proc[i]
    } else if(alpha_spend$as_incr[i] == 1){
      zb[i] <- 0
      yb[i] <- zb[i]*info$sd_proc[i]
    } else {
      zb[i] <- searchfunc(last = last, i = i,
                          as = alpha_spend$as_incr[i],
                          stdv = matrix(unlist(info), ncol = 2), za = za,
                          zb = zb, bs = bs, tol = tol, delta = 0, zj = zj) # delta = delta
      yb[i] <- zb[i]*info$sd_proc[i]
    }
    if(side == 1){
      ya[i] = zninf*info$sd_proc[i]
      za[i] = zninf
    } else {
      ya[i] = -yb[i]
      za[i] = -zb[i]
    }
    #nints[i] <- round((zb[i]-za[i])/h*info$sd_incr[i])+1

    if(i != nn){
      zj_wj_up <- z_n_w(r = r, info = info, za = za, zb = zb, i = i,
                        delta = 0)

      last <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                                zj = zj, last = last, zj_up = zj_wj_up$zj,
                                wj_up = zj_wj_up$wj, delta = 0, bs = bs)

      zj <- zj_wj_up$zj
      wj <- zj_wj_up$wj
    }
  }
  }

  alpha_ubound <- zb
  if(side == 2){
    alpha_lbound <- -zb
  } else {
    alpha_lbound <-NULL
  }
  return(list(inf_frac = inf_frac,
              org_inf_frac = org_inf_frac,
              alpha_ubound = alpha_ubound,
              alpha_lbound = alpha_lbound,
              alpha = alpha,
              alpha_spend = alpha_spend,
              delta = delta,
              design_R = design_R,
              info = info))
}

# esOF ----
# Error spending function - Lan and DeMets version of OBrien-Fleming
#' @importFrom stats qnorm pnorm
esOF <- function(alpha, timing){
  as_cum <- numeric(length(timing))
  as_incr <- numeric(length(timing))

  for(i in 1:length(timing)){
    as_cum[i] <- 2*(1 - pnorm(qnorm(1-alpha/2)/sqrt(timing[i]), mean = 0, sd = 1)) # removed /side
    if(i == 1){
      as_incr[i] <- as_cum[i]
    } else{
      as_incr[i] <- as_cum[i] - as_cum[i-1]
    }
  }
  return(list(as_cum = as_cum, as_incr = as_incr))
}

# esPoc ----
# Error spending function - Lan and DeMets version of Pocock
esPoc <- function(alpha, timing){
  as_cum <- numeric(length(timing))
  as_incr <- numeric(length(timing))

  for(i in 1:length(timing)){
    as_cum[i] <- alpha*log(1+(exp(1)-1)*timing[i])
    if(i == 1){
      as_incr[i] <- as_cum[i]
    } else{
      as_incr[i] <- as_cum[i] - as_cum[i-1]
    }
  }
  return(list(as_cum = as_cum, as_incr = as_incr))
}

# HSDC ----
# Error spending function - Hwang Sihi and DeCani
HSDC <- function(alpha, timing, gamma){
  as_cum <- numeric(length(timing))
  as_incr <- numeric(length(timing))

  for(i in 1:length(timing)){
    if(gamma != 0){
      as_cum[i] <- alpha*(1-(exp(-gamma*timing[i])))/(1-(exp(-gamma)))
    } else {
      as_cum[i] <- alpha*timing
    }
    if(i == 1){
      as_incr[i] <- as_cum[i]
    } else{
      as_incr[i] <- as_cum[i] - as_cum[i-1]
    }
  }
  return(list(as_cum = as_cum, as_incr = as_incr))
}

# rho ----
# Error spending function - rho family error spending
rho <- function(alpha, timing, rho){
  as_cum <- numeric(length(timing))
  as_incr <- numeric(length(timing))

  for(i in 1:length(timing)){
      as_cum[i] <- alpha*timing[i]^rho
    if(i == 1){
      as_incr[i] <- as_cum[i]
    } else{
      as_incr[i] <- as_cum[i] - as_cum[i-1]
    }
  }
  return(list(as_cum = as_cum, as_incr = as_incr))
}

# sd_inf ----
# Ratio of how much information is available
sd_inf <- function(timing){
  sd_incr <- sqrt(timing - c(0, timing[-length(timing)]))
  sd_proc <- sqrt(timing)

  return(list(sd_incr = sd_incr, sd_proc = sd_proc))
}


#searchfunc ----
searchfunc <- function(last, zj, i, as, stdv, za, zb, tol, bs, delta){
  maxnn <- 50; upper <- zb[i - 1]*stdv[i,2]
  if(bs == TRUE) upper <- za[i - 1]*stdv[i,2]
  del <- 10
  qout <- prob(xq = upper, last = last, zj, k = i, stdv = stdv,
                           bs = bs, delta = delta)
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
        qout <- prob(xq = upper, last = last, zj, k = i, stdv = stdv,
                                 bs = bs, delta = delta)
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
        qout <- prob(xq = upper, last = last, zj, k = i, stdv = stdv,
                                 bs = bs, delta = delta)
        if( qout >= as - tol){
          break
        }
      }
    }
  }
  return(upper/stdv[i,2])
}

# beta_boundary ----
#' @importFrom stats qnorm
beta_boundary <- function(inf_frac, beta, side, alpha_boundaries,
                          zninf = -20, tol = 1e-15,
                          rm_bs = 0,
                          es_beta = "esOF", delta = NULL, r = 18,
                          design_R = NULL, warp_root = NULL){
  nn <- length(inf_frac); lnn <- 5000; h <- 0.05

  if(is.null(delta)) delta = abs(qnorm(alpha_boundaries$alpha/side)+qnorm(beta)) 
  alpha_bound <- alpha_boundaries$alpha_ubound

  nn <- length(inf_frac)
  za <- numeric(length = nn); zb <- numeric(length = nn)
  stdv <- numeric(length = nn)
  ya <- numeric(length = nn); yb <- numeric(length = nn)
  last <- numeric(length = lnn)

  org_inf_frac <- inf_frac
  if(!is.null(warp_root)){
    org_inf_frac <- inf_frac*warp_root
  }

  beta_timing <- inf_frac
  if(!is.null(design_R)){
    beta_timing <- inf_frac/design_R
    org_inf_frac <- inf_frac
    if(max(org_inf_frac) < design_R){
      org_inf_frac <- c(org_inf_frac, design_R)
    }
    beta_timing <- c(beta_timing[beta_timing < 1],1)
    nn <- length(beta_timing)
    }


  if(sum(beta_timing > 1) > 0){
    beta_timing <- c(beta_timing[beta_timing < 1],1)
    nn <- length(beta_timing)
  }

  # calculate the beta spending
  if(rm_bs != 0){
    beta_timing <- c(rep(0, rm_bs),beta_timing[-c(1:rm_bs)])
  }
  
  beta_spend <- switch(which(es_beta == c("esOF", "esPoc", "HSDC", "rho")),
                        esOF(beta/side, beta_timing), esPoc(beta/side,beta_timing),
                        HSDC(beta/side, beta_timing, gamma),
                        rho(beta/side, beta_timing, rho))

  info <- sd_inf(timing = org_inf_frac) 

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
    za[1] <- qnorm(beta_spend$as_incr[1], mean = info$sd_proc[1]*delta, sd = 1)
    ya[1] <- za[1] * info$sd_incr[1]
  }

  zb <- alpha_bound
  yb <- zb*info$sd_proc

  zj_wj <- z_n_w(r = r, info = info, za = za, zb = zb, i = 1,
                 delta = delta)
  zj <- zj_wj$zj
  wj <- zj_wj$wj

  for(i in 2:nn){

    if(i == 2){
      last <- init_int(wj = wj, zj = zj, delta = delta, stdv = info$sd_incr)
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
      za[i] <- searchfunc(last = last,
                          i = i, as = beta_spend$as_incr[i],
                          stdv = matrix(unlist(info), ncol = 2),
                          za = za, zb = zb, tol = tol, bs =TRUE,
                          delta = delta, zj = zj)
      ya[i] <- za[i]*info$sd_proc[i]
    }

    if(i != nn){
      zj_wj_up <- z_n_w(r = r, info = info, za = za, zb = zb, i = i,
                        delta = delta)

      last <- recur_int(k = i, stdv = matrix(unlist(info),ncol = 2),
                        zj = zj, last = last, zj_up = zj_wj_up$zj,
                        wj_up = zj_wj_up$wj, delta = delta, bs = FALSE)

      zj <- zj_wj_up$zj
      wj <- zj_wj_up$wj
    }
  }
  testDrift <- 0 + abs(ya[length(inf_frac)])
  ret1 <- org_inf_frac
  ret2 <- info$sd_proc*delta
  return(list(ret1 = ret1, ret2 = ret2, as_incr = beta_spend$as_incr,
              as_cum = beta_spend$as_cum, za = za, zb = zb,
              ya = ya, yb = yb, drift = testDrift, delta = delta, info = info))
}

# define - z_n_w ----
# compute and update weights and z values
z_n_w <- function(r, info, za, zb, i, delta){
  j <- 1:(6*r-1)
  xi <- delta*info$sd_incr[i]+(j < r)*(-3-4*log(r/j)) +
    (r <= j & j <= 5*r)*(-3+3*(j-r)/(2*r))+ (5*r < j)*(3+4*log(r/(6*r-j)))

  if(sum(xi < za[i]) > 0){
    indi <- max(which(xi < za[i]))
    xi <- xi[indi:length(xi)]
    xi[1] <- za[i]
  }
  if(sum(xi > zb[i]) > 0){
    indi <- min(which(xi > zb[i]))
    xi <- xi[1:indi]
    xi[indi] <- zb[i]
  }

  m <- length(xi)*2-1
  zj <- numeric(m)

  zj[seq(1,m,2)] <- xi
  zj[seq(2,m-1,2)] <- (xi[seq(1,length(xi)-1,1)]+xi[seq(2,length(xi),1)])/2

  ij <- 1:m

  wj <- numeric(m)
  for(k in ij){
    if(k == 1){
      wj[k] <- (1/6)*(zj[3]-zj[1])
    } else if(k %in% seq(3,m-2,2)){
      wj[k] <- (1/6)*(zj[k+2]-zj[k-2])
    } else if(k %in% seq(2,m-1,2)){
      wj[k] <- (4/6)*(zj[k+1]-zj[k-1])
    } else{
      wj[k] <- (1/6)*(zj[m]-zj[m-2])
    }
  }
  return(list(zj = zj, wj = wj))
}

## old TSA
fcab_old <- function(last, nint, yam1, ybm1, h, x,
                 stdv, delta){ # fcad stands for?
  nlim <- 5000
  f <- numeric(length = nlim)
  
  for(fcon in 1:(nint+1)){
    grid2 <- yam1 + h*(fcon-1)
    f[fcon] <- last[fcon] * gfunc((x-grid2)/stdv, delta = delta)/stdv
  }
  return(trap_old(f = f, n = nint, h = h))
}

other_old <- function(ya, yb, i, stdv, h, last, nints){
  nlim <- 5000
  fn <- numeric(length = nlim)
  hh <- (yb[i] - ya[i])/nints[i]
  hlast <- (yb[i-1] - ya[i-1])/nints[i-1]
  
  for(j in 1:(nints[i]+1)){
    grid1 <- ya[i] + hh*(j-1) # renamed to grid1 from grid
    fn[j] <- fcab_old(last = last, nint = nints[i-1], yam1 = ya[i-1], ybm1 = yb[i-1],
                  h = hlast, x = grid1, stdv = stdv, delta = 0)
  }
  last <- fn
  return(last)
}

trap_old <- function(f, n, h){
  sum1 = f[1] # rename from sum to sum1
  
  for(j in 2:n){
    sum1 <- sum1 + 2*f[j]
  }
  sum1 <- sum1 + f[n + 1]
  
  return(h/2*sum1)
}

qpos_old <- function(xq, last, nint, yam1, ybm1, stdv){
  nwork <- 5000
  fun1 <- numeric(length = nwork)
  hlast <- (ybm1 - yam1)/nint
  
  for(j in 1:(nint+1)){
    grid3 <- yam1 + hlast * (j-1)
    fun1[j] <- last[j]*(1 - pnorm((xq - grid3)/stdv, mean= 0, sd = 1))
  }
  return(trap_old(f = fun1, n = nint, h = hlast))
}

gfunc <- function(x, delta){
  exp(-0.5*(x - delta)^2)/sqrt(2*pi)
}

first_old <- function(ya, yb, h, stdv, nints, delta, lnn){
  hh <- (yb - ya)/nints[1]
  
  last <- numeric(length = lnn)
  
  for(j in 1:(nints[1])+1){
    last[j] <- gfunc(x = (ya+hh*(j-1))/stdv, delta = delta)/stdv
  }
  return(last)
}


# alphas <- function(alpha, side, ti, tol = 1e-13){
#   alphaSpendCum <- numeric(length(ti))
#   alphaSpendDelta <- numeric(length(ti))
#   
#   for(i in 1:length(ti)){
#     alphaSpendCum[i] <- 2*(1 - pnorm(qnorm(1-alpha/side/2)/sqrt(ti[i]), mean = 0, sd = 1))
#     alphaSpendCum[i] <- side*alphaSpendCum[i]
#     if(i == 1){
#       alphaSpendDelta[i] <- alphaSpendCum[i]
#     } else{
#       alphaSpendDelta[i] <- alphaSpendCum[i] - alphaSpendCum[i-1]
#     }
#     if(alphaSpendDelta[i] < tol){
#       alphaSpendDelta[i] <- 0 }
#   }
#   return(list(alphaValuesCum = alphaSpendCum, alphaValuesDelta = alphaSpendDelta))
# }


betas_Obf <- function(use, nn, beta, informationFractions, tol = 1e-13){
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

searchfunc_old <- function(last, nints, i, valSF, stdv, ya, yb){
  maxnn <- max(c(length(nints), 50)); upper <- yb[i - 1]
  del <- 10
  eps = 1e-7
  qout <- qpos_old(xq = upper, last = last, nint = nints[i-1], yam1 = ya[i-1],
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
        qout <- qpos_old(xq = upper, last = last, nint = nints[i-1], yam1 = ya[i-1],
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
        qout <- qpos_old(xq = upper, last = last, nint = nints[i-1], yam1 = ya[i-1],
                     ybm1 = yb[i-1], stdv = stdv)
        if( qout >= valSF - eps){
          break
        }
      }
    }
  }
  return(yb)
}

## old TSA functions (translated from java)
getInnerWedge <- function(informationFractions, beta, bsInf, delta, side, fakeIFY,
                          zninf = -20, tol = 1e-13){
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
  
  outbeta = betas_Obf(use, nn, beta, informationFractions = informationFractions)
  outsd = sdfunc(ti = informationFractions)
  
  if(outbeta$betaValuesDelta[1] <= 0 || outbeta$betaValuesDelta[1] >= beta){
    outbeta$betaValuesDelta[1] <- pmin(beta, outbeta$betaValuesDelta[1])
    outbeta$betaValuesDelta[1] <- pmax(0, outbeta$betaValuesDelta[1])
  }
  
  if(outbeta$betaValuesDelta[1] < tol){
    za[1] <- zninf
    ya[1] <- za[1]*outsd$sdincr[1]
  } else if(outbeta$betaValuesDelta[1] == beta){
    za[1] <- 0
    ya[1] <- za[1]*outsd$sdincr[1]
  } else {
    za[1] <- qnorm(outbeta$betaValuesDelta[1]) # the first is easy
    ya[1] <- za[1] * outsd$sdincr[1]
  }
  
  if(side == 1){
    zb[1] <- -zninf
    yb[1] <- zb[1]* outsd$sdincr[1]
  } else {
    zb[1] <- -za[1]
    yb[1] <- -ya[1]
  }
  
  nints[1] <- round((abs(yb[1] - ya[1])/h*outsd$sdincr[1])) + 1
  for(i in 2:nn){
    if(i == 2){
      last <- first_old(ya = ya[1], yb = yb[1], h = h, stdv = outsd$sdincr[1],
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
      yb <- searchfunc_old(last = last, nints = nints,
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
      last <- other_old(ya, yb, i, outsd$sdincr[i], h, last, nints)
    }
  }
  testDrift <- fakeIFY + abs(ya[length(informationFractions)])
  ret1 <- informationFractions
  ret2 <- za + outsd$sdproc*testDrift
  return(list(ret1 = ret1, ret2 = ret2, betaValuesDelta = outbeta$betaValuesDelta, za = za, zb = zb,
              ya = ya, yb = yb, drift = testDrift))
}