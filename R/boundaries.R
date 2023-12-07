#' @title Boundaries for group sequential designs
#' @description
#' Calculates alpha- and potentially beta-spending boundaries for group sequential designs for meta-analysis. Should be used for exploring how the different arguments affect the sequential design. The function is not intended to be used individually for Trial Sequential Analysis. For this purpose, we recommend RTSA().
#' 
#'
#' @param timing Expected timings of interim analyses and final analysis as a vector consisting of values from 0 to 1. 
#' @param alpha The level of type I error as a percentage, the default is 0.05 corresponding to 5\%.
#' @param beta The level of type II error as a percentage, the default is 0.1 corresponding to 10\%.
#' @param side Whether a 1- or 2-sided hypothesis test is used. Defaults to 2. Options are 1 or 2.
#' @param futility Futility boundaries added to design. Options are: none, non-binding and binding. Default is "none".
#' @param es_alpha The error spending function for alpha-spending. Options are: "esOF" (Lan & DeMets version of O'Brien-Fleming boundaries), "esPoc" (Lan & DeMets version of Pocock boundaries), "HSDC" (Hwang Sihi and DeCani) and "rho" (rho family). Defaults to "esOF".
#' @param es_beta The error spending function for beta-spending. For options see es_alpha. Defaults to NULL.
#' @param type Whether the boundaries are used for design or analysis. We recommend only to use the boundaries() function with type equal to design. Defaults to design.
#' @param design_R If type is analysis, a scalar for achieving the right amount of power is required. It is recommended not to use the boundaires() function with the setting type equal to analysis. Defaults to NULL. 
#' @param tol Tolerance level for numerical integration. Defaults to 1e-09.
#'
#' @returns A \code{boundaries} object which includes:
#' \item{inf_frac}{Timing of interim analyses and final analysis. Potentially modified if \code{type = "analysis"}.}
#' \item{org_inf_frac}{Original timing. If \code{type = "design"}.}
#' \item{alpha_ubound}{Upper alpha-spending boundaries}
#' \item{alpha_lbound}{Lower alpha-spending boundaries}
#' \item{alpha}{As input}
#' \item{alpha_spend}{List of cumulative and incremental spending}
#' \item{delta}{Drift parameter}
#' \item{design_R}{If \code{type = "analysis"} it is the scalar for correct power in the design. Else NULL.}
#' \item{info}{List of the information as the squareroot of the information increments and the squareroot of the cumulative information}
#' \item{beta_ubound}{Upper beta-spending boundaries}
#' \item{beta_lbound}{Lower beta-spending boundaries}
#' \item{root}{Scalar for achieving correct power}
#' \item{beta_spend}{List of cumulative and incremental spending}
#' \item{pwr}{List of probabilities for rejecting the null under the sample size settings being true at each analysis and the sum.}
#' \item{tIe}{List of probabilities for type-I-error at each analysis and the sum}
#' \item{side}{As input}
#' \item{beta}{As input}
#' \item{es_alpha}{As input}
#' \item{es_beta}{As input}
#' \item{type}{As input}
#' \item{futility}{As input}
#'
#' @examples
#' boundaries(timing = c(0.25, 0.5, 0.75, 1), alpha = 0.05, beta = 0.1,
#'  side = 2, futility = "non-binding", es_alpha = "esOF", es_beta = "esOF")
#'
#' @export
boundaries <- function(timing, alpha = 0.05, zninf, beta = 0.1, side = 2,
                       futility = "none", es_alpha = "esOF", es_beta = NULL,
                       type = "design", design_R = NULL, tol = 1e-09){

  n_max <- 15 # finding root for power
  nn_max <- 50 # finding root for type-I-error

  if(futility == "none"){
    es_beta <- NULL
  }

  if(max(timing) < 1 & type == "design"){
    timing <- c(timing, 1)
  }
  
  # calculate the initial alpha boundaries
  boundout <- alpha_boundary(inf_frac = timing,
                            side = side, alpha = alpha,
                             zninf = zninf, beta = beta,
                                    es_alpha = es_alpha,
                                    type = type, design_R = design_R)

  # calculate the boundaries depending on the design:
  # 1-sided without futility or with either non-binding or binding futility
  # 2-sided without futility or with either non-binding or binding futility
  if(side == 1){
    if(futility == "none"){
      if(type == "design"){
        lb <- rep(-20,length(timing)) # dummy lower boundaries

        root <- uniroot(right_power, lower = 0.9, upper = 1.2,
                        za = lb, zb = boundout$alpha_ubound, delta = boundout$delta,
                        timing = timing, beta = beta)$root

        info <- sd_inf(timing*root)
        pwr <- ma_power(zb = boundout$alpha_ubound, za = lb,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound, za = lb,
                               delta = 0,
                               info = info)

        boundout$beta_ubound <- NA
        boundout$beta_lbound <- NA
        boundout$alpha_lbound <- NA
        boundout$root <- root
      } else {
        # type is analysis
        if(is.null(design_R)){
          stop("Root is needed to calculate analysis")
        }

        lb <- rep(-20,length(timing)) # dummy lower boundaries

        info <- sd_inf(timing)
        pwr <- ma_power(zb = boundout$alpha_ubound, za = lb,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound, za = lb,
                               delta = 0,
                               info = info)

        boundout$beta_ubound <- NA
        boundout$beta_lbound <- NA
        boundout$alpha_lbound <- NA

        boundout$root <- design_R
      }
    }
    else if(futility == "non-binding"){
      if(type == "design"){
        lb <- beta_boundary(inf_frac = timing, beta = beta, side = 1,
                                   alpha_boundaries = boundout,
                                   es_beta = es_beta)

        root <- uniroot(inf_warp, lower = 0.9, upper = 1.5,
                        alpha_boundaries = boundout,
                        tol = 1e-09, side = 1, beta = beta, timing = timing,
                        es_beta = es_beta)$root

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   es_beta = es_beta,
                                   warp_root = root)

        info <- sd_inf(timing*root)

        pwr <- ma_power(zb = boundout$alpha_ubound, za = lb$za,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound, za = lb$za,
                               delta = 0,
                               info = info)

        boundout$root <- root
        boundout$beta_lbound <- lb$za
        boundout$beta_ubound <- NA
        boundout$alpha_lbound <- NA

        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)
      } else {
        if(is.null(design_R)){
          stop("Root is needed to calculate analysis")
        }

        lb <- beta_boundary(inf_frac = timing,
                                   beta = beta,
                                   side = 1,
                                   alpha_boundaries = boundout,
                                   design_R = design_R)

        boundout$beta_lbound <- lb$za
        boundout$beta_ubound <- NA
        boundout$alpha_lbound <- NA

        if(boundout$beta_lbound[length(boundout$alpha_ubound)] >
           boundout$alpha_ubound[length(boundout$alpha_ubound)]){
          boundout$beta_lbound[length(boundout$alpha_ubound)] =
            boundout$alpha_ubound[length(boundout$alpha_ubound)]
        }

        info <- sd_inf(timing)
        pwr <- ma_power(zb = boundout$alpha_ubound, za = lb$za,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound, za = lb$za,
                               delta = 0,
                               info = info)

        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)
        boundout$root <- design_R

      }
    } else if(futility == "binding"){
      if(type == "design"){
        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   es_beta = es_beta)

        # get then the upper boundaries again
        lb$alpha_ubound <- -lb$za
        lb$alpha <- alpha

        ub <- beta_boundary(inf_frac = timing, beta = alpha,
                                   side = 1, alpha_boundaries = lb,
                                   delta = 0, es_beta = es_alpha)

        ub$alpha_ubound <- -ub$za
        ub$alpha <- alpha

        upperRoot <- 0.95
        n_itr <- 1
        while(n_itr <= nn_max){
        root <- try(uniroot(inf_warp, lower = upperRoot-0.1, upper = upperRoot, alpha_boundaries = ub,
                        tol = 1e-09, beta = beta, timing = timing,
                        es_beta= es_beta)$root, TRUE)
        if(inherits(root,"try-error")){
            upperRoot <- upperRoot + 0.1
            n_itr <- n_itr + 1
        } else {
          break
          }
        }

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = ub,
                                   es_beta = es_beta,
                                   warp_root = root)

        info <- sd_inf(timing*root)
        tIe <- ma_power(zb = -ub$za, za = lb$za, delta = 0,
                               info = info)

        while(abs(tIe[[2]] - alpha) > 1e-6){

          lb$alpha_ubound <- -lb$za
          lb$alpha <- alpha

          ub <- beta_boundary(inf_frac = timing, beta = alpha,
                                     side = 1, alpha_boundaries = lb,
                                     delta = 0, es_beta = es_alpha)

          ub$alpha_ubound <- -ub$za
          ub$alpha <- alpha

          upperRoot <- 1.05
          n_itr <- 1
          while(n_itr <= n_max){
          root <- try(uniroot(inf_warp, lower = upperRoot - 0.1, upper = upperRoot, alpha_boundaries = ub,
                          tol = 1e-09, beta = beta, timing = timing,
                          es_beta = es_beta)$root, TRUE)
          if(inherits(root,"try-error")){
            upperRoot <- upperRoot + 0.1
            n_itr <- n_itr + 1
          } else {
            break
          }
          }

          lb <- beta_boundary(inf_frac = timing, beta = beta,
                                     side = 1, alpha_boundaries = ub,
                                     delta = NULL, es_beta = es_beta,
                                     warp_root = root)

          info <- sd_inf(timing*root)
          tIe <- ma_power(zb = -ub$za, za = lb$za, delta = 0,
                                 info = info)
          pwr <- ma_power(zb = -ub$za, za = lb$za, delta = boundout$delta,
                                 info = info)
          tIe
        }

        boundout$root <- root
        boundout$alpha_ubound <- -ub$za
        boundout$beta_lbound <- lb$za
        boundout$beta_ubound <- NA
        boundout$alpha_lbound <- NA

        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)
      } else {
        # type is analysis

        if(is.null(design_R)){
          stop("Root is needed to calculate analysis")
        }

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   es_beta = es_beta,
                                   design_R = design_R)

        # get then the upper boundaries again
        lb$alpha_ubound <- -lb$za
        lb$alpha <- alpha

        ub <- beta_boundary(inf_frac = timing, beta = alpha,
                                   side = 1, alpha_boundaries = lb,
                                   delta = 0, es_beta = es_alpha,
                                   design_R = design_R)

        ub$alpha_ubound <- -ub$za
        ub$alpha <- alpha

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = ub,
                                   es_beta = es_beta,
                                   design_R = design_R)

        if(lb$za[length(ub$za)] >
           -ub$za[length(ub$za)]){
          lb$za[length(ub$za)] =
            -ub$za[length(ub$za)]
        }

        info <- sd_inf(timing)
        tIe <- ma_power(zb = -ub$za, za = lb$za, delta = 0,
                               info = info)
        pwr <- ma_power(zb = -ub$za, za = lb$za, delta = boundout$delta,
                               info = info)

        boundout$alpha_ubound <- -ub$za
        boundout$beta_lbound <- lb$za
        boundout$beta_ubound <- NA
        boundout$alpha_lbound <- NA

        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)

        boundout$root <- design_R
      }
    }
  } else {
    if(futility == "none"){
      if(type == "design"){
        
        upperRoot <- 1.05
        n_itr <- 1
        while(n_itr <= n_max){
          root <- try(uniroot(right_power, lower = upperRoot-0.1, upper = upperRoot,
                              za = boundout$alpha_lbound, zb = boundout$alpha_ubound,
                              delta = boundout$delta, timing = timing,
                              beta = beta, tol = 1e-13)$root,TRUE)
          if(inherits(root,"try-error")){
            upperRoot <- upperRoot + 0.1
            n_itr <- n_itr + 1
          } else {
            break
          }
        }

        info <- sd_inf(timing*root)

        pwr <- ma_power(zb = boundout$alpha_ubound, za = boundout$alpha_lbound,
                        delta = boundout$delta,
                        info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound, za = boundout$alpha_lbound,
                        delta = 0,
                        info = info)

        boundout$beta_ubound <- NA
        boundout$beta_lbound <- NA

        boundout$root <- root
      } else {
        # type is analysis
        if(is.null(design_R)){
          stop("Root is needed to calculate analysis")
        }

        info <- sd_inf(timing)
        pwr <- ma_power(zb = boundout$alpha_ubound,
                               za = boundout$alpha_lbound,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound,
                               za = boundout$alpha_lbound,
                               delta = 0,
                               info = info)

        boundout$beta_ubound <- NA
        boundout$beta_lbound <- NA

        boundout$root <- design_R

      }
    }
    else if(futility == "non-binding"){
      if(type == "design"){
        delta <- abs(qnorm(alpha/side)+qnorm(beta))

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   delta = delta, es_beta = es_beta)

        upperRoot <- .95
        n_itr <- 1
        while(n_itr <= nn_max){
        root <- try(uniroot(inf_warp, lower = upperRoot-0.02, upper = upperRoot, alpha_boundaries = boundout,
                        tol = 1e-09, side = 1, delta = delta, timing = timing,
                        es_beta = es_beta, beta = beta)$root,TRUE)
        if(inherits(root,"try-error")){
          upperRoot <- upperRoot + 0.02
          n_itr <- n_itr + 1
        } else {
          break
        }
        }
        
        if(inherits(root,"try-error")){
          stop("Non-binding futility boundaries could not be computed. Consider setting futility to 'none' or change some of the design parameters.")
        }

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   delta = delta, warp_root = root,
                                   es_beta = es_beta)
      
        n_itr <- 1
        upperRoot <- 0.95
        while(n_itr <= nn_max){
        root <- try(uniroot(inf_warp, lower = upperRoot-0.05, upper = upperRoot, alpha_boundaries = boundout,
                        tol = 1e-09, side = 1, rm_bs = sum(lb$za < 0),
                        delta = delta, timing = timing, es_beta = es_beta,
                        beta = beta)$root,TRUE)
        if(inherits(root,"try-error")){
          upperRoot <- upperRoot + 0.05
          n_itr <- n_itr + 1
        } else {
          break
        }
      }

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   delta = delta, rm_bs = sum(lb$za < 0),
                                   warp_root = root, es_beta = es_beta)

        info <- sd_inf(timing*root)

        boundout$beta_ubound <- c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20])
        boundout$beta_lbound <- c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20])

        boundout$root <- root

        pwr <- ma_power(zb = boundout$alpha_ubound, za = lb$za,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound, za = lb$za,
                               delta = 0,
                               info = info)

        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)

      } else {
        #if type is analysis

        lb <- beta_boundary(inf_frac = timing,
                                   beta = beta,
                                   side = 1,
                                   alpha_boundaries = boundout,
                                   delta = boundout$delta,
                                   design_R = design_R)

        lb <- beta_boundary(inf_frac = timing,
                                   beta = beta,
                                   side = 1,
                                   alpha_boundaries = boundout,
                                   delta = boundout$delta,
                                   design_R = design_R,
                                   rm_bs = sum(lb$za < 0))
        
        

        lb <- beta_boundary(inf_frac = timing,
                                   beta = beta,
                                   side = 1,
                                   alpha_boundaries = boundout,
                                   delta = boundout$delta,
                                   design_R = design_R,
                                   rm_bs = sum(lb$za < 0))

        boundout$beta_ubound <- c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20])
        boundout$beta_lbound <- c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20])


        if(boundout$beta_ubound[length(boundout$alpha_ubound)] >
           boundout$alpha_ubound[length(boundout$alpha_ubound)]){
          boundout$beta_ubound[length(boundout$alpha_ubound)] =
            boundout$alpha_ubound[length(boundout$alpha_ubound)]
          boundout$beta_lbound[length(boundout$alpha_ubound)] =
            boundout$alpha_lbound[length(boundout$alpha_ubound)]
        }

        info <- sd_inf(timing)
        pwr <- ma_power(zb = boundout$alpha_ubound,
                               za = boundout$alpha_lbound,
                               zd = boundout$beta_ubound,
                               zc = boundout$beta_lbound,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound,
                               za = boundout$alpha_lbound,
                               zd = boundout$beta_ubound,
                               zc = boundout$beta_lbound,
                               delta = 0,
                               info = info)

        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)

        boundout$root <- design_R

      }
    }
    else if(futility == "binding"){
      if(type == "design"){
        delta <- abs(qnorm(alpha/side)+qnorm(beta))
        
        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   delta = delta, es_beta = es_beta)

        # get then the upper boundaries again
        lb$alpha_ubound <- -lb$za
        lb$alpha <- alpha

        ub <- beta_boundary(inf_frac = timing, beta = alpha,
                                   side = 2, alpha_boundaries = lb,
                                   delta = 0, es_beta = es_alpha)

        ub$alpha_ubound <- -ub$za
        ub$alpha <- alpha

        upperRoot <- 1.1
        n_itr <- 1
        while(n_itr <= n_max){
        root <- try(uniroot(inf_warp, lower = upperRoot - 0.1, upper = upperRoot, alpha_boundaries = ub,
                        tol = 1e-09, side = 1, delta = delta, timing = timing,
                        es_beta = es_beta, beta = beta)$root,TRUE)
        if(inherits(root,"try-error")){
          upperRoot <- upperRoot + 0.1
          n_itr <- n_itr + 1
        } else {
          break
        }
        }

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = ub,
                                   delta = delta,
                                   warp_root = root,es_beta = es_beta)
        if(length(timing)-1 == sum(lb$za < 0)){
          stop("It is not possible to compute the design with futility.")
        }
        
        upperRoot <- 1.05
        n_itr <- 1
        while(n_itr <= n_max){
        root <- try(uniroot(inf_warp, lower = upperRoot - 0.1, upper = upperRoot, alpha_boundaries = ub,
                        tol = 1e-09, rm_bs = sum(lb$za < 0), timing = timing,
                        es_beta = es_beta, beta = beta, side = 1,
                        delta = delta)$root,TRUE)
        if(inherits(root,"try-error")){
          upperRoot <- upperRoot + 0.1
          n_itr <- n_itr + 1
        } else {
          break
        }
        }

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = ub,
                                   delta = delta, rm_bs = sum(lb$za < 0),
                                   warp_root = root, es_beta = es_beta)

        info = sd_inf(timing*root)
        tIe <- ma_power(zb = -ub$za, delta = 0,
                               info = info, za = ub$za,
                               zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                               zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))
        pre_tIe <- tIe

        nN_max <- 50
        i <- 0
        while(abs(tIe[[2]]*2 - alpha) > 1e-6 & i < nn_max){
          lb$alpha_ubound <- -lb$za
          lb$alpha <- alpha

          ub <- beta_boundary(inf_frac = timing, beta = alpha, # root
                                     side = 2,alpha_boundaries = lb,
                                     delta = 0, es_beta = es_alpha)

          ub$alpha_ubound <- -ub$za
          ub$alpha <- alpha

          upperRoot <- 1.1
          root_max <- 10
          n_itr <- 1
          while(n_itr < root_max){
          root <- try(uniroot(inf_warp, lower = upperRoot - 0.1, upper = upperRoot, alpha_boundaries = ub,
                          side = 1,
                          tol = 1e-09, rm_bs = sum(lb$za < 0),
                          timing = timing, es_beta = es_beta, beta = beta,
                          delta = delta)$root,TRUE)
          if(inherits(root,"try-error")){
            upperRoot <- upperRoot + 0.1
            n_itr <- n_itr + 1
          } else {
            break
          }
          }

          lb <- beta_boundary(inf_frac = timing, beta = beta,
                                     side = 1, alpha_boundaries = ub,
                                     delta = delta, rm_bs = sum(lb$za < 0),
                                     warp_root = root, es_beta = es_beta)

          info = sd_inf(timing*root)

          tIe <- ma_power(zb = -ub$za, delta = 0,
                                 info = info, za = ub$za,
                                 zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                                 zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))
          pwr <- ma_power(zb = -ub$za, delta = delta,
                                 info = info, za = ub$za,
                                 zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                                 zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))

          if(abs(tIe[[2]]-pre_tIe[[2]]) < 1e-8){
            i <- nn_max
          }
          i <- i + 1
          pre_tIe <- tIe
        }

        if(i == nn_max){
          no_con <- function(x){
            nnn_max <- 20
            i <- 0
            while(abs(tIe[[2]]*2 - alpha) > 1e-7 & i < nnn_max){
              lb$alpha_ubound <- -lb$za*x
              lb$alpha <- alpha

              ub <- beta_boundary(inf_frac = timing*root, beta = alpha,
                                         side = 2, alpha_boundaries = lb,
                                         delta = 0, es_beta = es_alpha)

              ub$alpha_ubound <- -ub$za
              ub$alpha <- alpha

              upperRoot <- 1.1
              root_max <- 5
              n_itr <- 1
              while(n_itr < root_max){
              root <- try(uniroot(inf_warp, lower = upperRoot - 0.1, upper = upperRoot, alpha_boundaries = ub,
                              tol = 1e-09, rm_bs = sum(lb$za < 0),
                              timing = timing,es_beta = es_beta, beta = beta,
                              delta = delta)$root,TRUE)
              if(inherits(root,"try-error")){
                upperRoot <- upperRoot + 0.1
                n_itr <- n_itr + 1
              } else {
                break
              }
              }

              lb <- beta_boundary(inf_frac = timing*root, beta = beta,
                                         side = 1, alpha_boundaries = ub,
                                         delta = delta, rm_bs = sum(lb$za < 0),
                                         warp_root = root)

              info = sd_inf(timing*root)
              tIe <- ma_power(zb = -ub$za, delta = 0,
                                     info = info, za = ub$za,
                                     zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                                     zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))
              i <- i + 1
            }
            return(tIe[[2]] - alpha/side)
          }

          x <- seq(-1,1.6,by = 0.2)
          out <- 0
          for(i in x){
            trytest <- try(no_con(x=i), TRUE)
            if(!(inherits(trytest,"try-error"))){
              out <- c(out, no_con(x=i))
              if(sum(diff(sign(out)) == 2) > 0){
                break
              }
            } else {
              out <- c(out, 0)
            }
            print(out)
          }
          lower <- c(0,x)[max(which(diff(sign(out)) == 2))]

          if(is.na(lower)){
            stop("No futility boundaries were found. Select non-binding futility or
             no futility boundaries")
          }

          upper <- lower + 0.2
          root_x <- uniroot(no_con, upper = upper, lower = lower,tol = 1e-08)$root

          lb$alpha_ubound <- -lb$za*root_x
          lb$alpha <- alpha

          ub <- beta_boundary(inf_frac = timing*root, beta = alpha,
                                     side = 2,
                                     alpha_boundaries = lb,
                                     delta = 0, es_beta = es_alpha)

          ub$alpha_ubound <- -ub$za
          ub$alpha <- alpha

          upperRoot <- 1.1
          root_max <- 5
          n_itr <- 1
          while(n_itr < root_max){
          root <- try(uniroot(inf_warp, lower = upperRoot - 0.1, upper = upperRoot, alpha_boundaries = ub,
                          tol = 1e-08, rm_bs = sum(lb$za < 0), timing = timing,
                          es_beta = es_beta, beta = beta, delta = delta)$root,TRUE)
          if(inherits(root,"try-error")){
            upperRoot <- upperRoot + 0.1
            n_itr <- n_itr + 1
          } else {
            break
          }
          }

          lb <- beta_boundary(inf_frac = timing, beta = beta,
                                     side = 1, alpha_boundaries = ub,
                                     delta = delta, rm_bs = sum(lb$za < 0),
                                     warp_root = root)

          info <- sd_inf(timing*root)
          tIe <- ma_power(zb = -ub$za, delta = 0,
                                 info = info, za = ub$za,
                                 zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                                 zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))
          pwr <- ma_power(zb = -ub$za, delta = boundout$delta,
                                 info = info, za = ub$za,
                                 zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                                 zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))
        } else {

          upperRoot <- 1.1
          n_itr <- 1
          while(n_itr <= n_max){
          root <- try(uniroot(inf_warp, lower = upperRoot - 0.1, upper = upperRoot, alpha_boundaries = ub,
                          tol = 1e-08, rm_bs = sum(lb$za < 0), timing = timing,
                          es_beta = es_beta, beta = beta, side = 1,
                          delta = delta)$root,TRUE)
          if(inherits(root,"try-error")){
            upperRoot <- upperRoot + 0.1
            n_itr <- n_itr + 1
          } else {
            break
          }
        }


          lb <- beta_boundary(inf_frac = timing, beta = beta,
                                     side = 1, delta = delta,
                                     alpha_boundaries = ub,
                                     rm_bs = sum(lb$za < 0),
                                     warp_root = root, es_beta = es_beta)

          if(sum(lb$za[-length(timing)] == -20) == (length(timing)-1)){
            lb$za[length(timing)] <- -ub$za[length(timing)]
          }

          root <- uniroot(right_power, lower = 0.9, upper = 1.4, tol = 1e-09,
                          za = ub$za, zb = -ub$za,
                          zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                          zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]),
                          delta = boundout$delta, timing = timing, beta = beta)$root

          info = sd_inf(timing*root)
          tIe <- ma_power(zb = -ub$za, delta = 0,
                                 info = info, za = ub$za,
                                 zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                                 zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))

          pwr <- ma_power(zb = -ub$za, delta = boundout$delta,
                                 info = info, za = ub$za,
                                 zc = c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20]),
                                 zd = c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20]))
        }
        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)
        boundout$alpha_ubound <- -ub$za
        boundout$alpha_lbound <- ub$za
        boundout$beta_ubound <- c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20])
        boundout$beta_lbound <- c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20])
        boundout$root <- root

      }  else {
        # if type is analysis

        delta <- abs(qnorm(alpha/side)+qnorm(beta))

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = boundout,
                                   delta = delta,
                                   design_R = design_R)

        # get then the upper boundaries again
        lb$alpha_ubound <- -lb$za
        lb$alpha <- alpha

        ub <- beta_boundary(inf_frac = timing, beta = alpha,
                                   side = 2, alpha_boundaries = lb,
                                   delta = 0,
                                   design_R = design_R,
                                   es_beta = es_alpha)

        ub$alpha_ubound <- -ub$za
        ub$alpha <- alpha

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = ub,
                                   delta = delta,
                                   design_R = design_R)

        lb <- beta_boundary(inf_frac = timing, beta = beta,
                                   side = 1, alpha_boundaries = ub,
                                   delta = delta, rm_bs = sum(lb$za < 0),
                                   design_R = design_R)

        boundout$beta_spend <- list(bs_incr = lb$as_incr, bs_cum = lb$as_cum)
        boundout$alpha_ubound <- -ub$za
        boundout$alpha_lbound <- ub$za
        boundout$beta_ubound <- c(rep(NA,sum(abs(lb$za) == 20)), lb$za[abs(lb$za) < 20])
        boundout$beta_lbound <- c(rep(NA,sum(abs(lb$za) == 20)), -lb$za[abs(lb$za) < 20])

        if(!(sum(is.na(boundout$beta_ubound)) == length(boundout$beta_ubound))){
          if(boundout$beta_ubound[length(boundout$alpha_ubound)] >
             boundout$alpha_ubound[length(boundout$alpha_ubound)]){
            boundout$beta_ubound[length(boundout$alpha_ubound)] =
              boundout$alpha_ubound[length(boundout$alpha_ubound)]
            boundout$beta_lbound[length(boundout$alpha_ubound)] =
              boundout$alpha_lbound[length(boundout$alpha_ubound)]
          }
        }

        info <- sd_inf(timing)
        pwr <- ma_power(zb = boundout$alpha_ubound,
                               za = boundout$alpha_lbound,
                               zd = boundout$beta_ubound,
                               zc = boundout$beta_lbound,
                               delta = boundout$delta,
                               info = info)
        tIe <- ma_power(zb = boundout$alpha_ubound,
                               za = boundout$alpha_lbound,
                               zd = boundout$beta_ubound,
                               zc = boundout$beta_lbound,
                               delta = 0,
                               info = info)

        boundout$root <- design_R

      }
    }
  }
  
  out <- c(boundout,list(pwr = pwr, tIe = tIe, side = side,
              beta = beta, es_alpha = es_alpha,
              es_beta = es_beta, type = type, futility = futility, design_R = design_R))
  class(out) <- "boundaries"
  return(out)
}


# FUNCTION | print boundaries ----
#' @method print boundaries
#' @export
print.boundaries <- function(x, ...) {
  cat(paste0("Boundaries for a ", x$side, "-sided design with a type-I-error of ",
             x$alpha, " and type-II-error of ", x$beta,".\n"))
  cat(paste0("Futility is set to: ", x$futility, ". Alpha-spending function: ",
             x$es_alpha, ".\n", "Beta-spending function: ", x$es_beta, ".\n\n"))
  cat("Timing and boundaries:\n")
  if(x$side == 1){
    if(!is.null(x$design_R)){
      sma_timing <- x$inf_frac
    } else {
      sma_timing <- x$inf_frac*x$root
    }
    y <- data.frame("SMA_Timing" = sma_timing,
                    "Upper" = x$alpha_ubound)
    if(x$futility != "none"){
      y <- cbind(y, data.frame("FutLower" = x$beta_lbound))
    }
  } else {
    if(!is.null(x$design_R)){
      sma_timing <- x$inf_frac
    } else {
      sma_timing <- x$inf_frac*x$root
    }
    y <- data.frame("SMA_Timing" = sma_timing,
                    "Upper" = x$alpha_ubound,
                    "Lower" = x$alpha_lbound)
    if(x$futility != "none"){
      y <- cbind(y, data.frame("FutUpper" = x$beta_ubound,
                               "FutLower" = x$beta_lbound))
    }
  }
  print(round(y,3), row.names = FALSE)
}

#' @title Plot of boundaries for group sequential designs
#'
#' @param x boundaries object
#' @param theme Whether the theme is "classic" or "aussie"
#' @param ... Other arguments to plot.boundaries
#'
#' @return Plot. Either a plot for two- or one-sided testing.
#' @export
#'
#' @importFrom ggplot2 ggplot coord_cartesian geom_hline theme_bw geom_vline geom_line geom_point aes theme element_blank geom_ribbon xlab ylab scale_x_continuous expansion scale_y_continuous scale_fill_identity scale_colour_manual ggtitle geom_segment geom_label scale_y_reverse sec_axis theme_classic element_text margin scale_fill_manual scale_fill_discrete guides guide_legend
#'
#' @examples
#' bounds <- boundaries(timing = c(0.5,0.75, 1), alpha  = 0.025, beta = 0.2,
#' side = 1, futility = "none", es_alpha = "esOF")
#' plot(x = bounds)
#'
plot.boundaries = function(x, theme = "classic", ...){
  plot.RTSA(x, type = "classic", theme = theme)
}

