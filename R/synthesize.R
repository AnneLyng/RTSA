#' synthesize
#'
#' Synthesize trials in a meta-analysis
#'
#' @param y list-object. Preferably a "synthPrepped" object from the metaPrepare function.
#' @param sign TBA
#' @param fixedStudy TRUE or FALSE. When using GLM methods for the meta-analysis, should the study effect be fixed-effect. Defaults to TRUE.
#' @param hksj TRUE or FALSE. When conducting a random-effects model should the Harting-Knapp-Sidik-Jonkman adjustment be used. Else is DerSimonian-Laird used. Defaults to FALSE.
#'
#' @return A list is returned with the a set of the following items:
#' \item{fw}{Weights for the fixed-effect model}
#' \item{peF}{Vector containing fixed-effect model results - pooled effect size, lower 95\% confidence limit, upper 95\% confidence limit, z-value, p-value, log pooled effect size, variance of estimate}
#' \item{rwR}{Weights for the random-effects model}
#' \item{peR}{Vector containing random-effects model results - pooled effect size, lower 95\% confidence limit, upper 95\% confidence limit, z-value (DL) or t-value (HKSJ), p-value, variance of estimate}
#' \item{Q}{vector containing Q-measure, degrees of freedom (df) and p-value for Q}
#' \item{U}{vector containing tau2 (random effect variance estimate), H-measure, I2 (inconsistency)
#' and D2 (diversity)}
#'
#' @importFrom stats binomial glm coef vcov qnorm pnorm pchisq pt qt
#' @export
#'
#' @examples
#' data(perioOxy)
#' m.prepare = with(perioOxy, metaPrepare(outcome = "RR", eI = eI,
#'  nI = nI, eC = eC, nC = nC, method = "MH"))
#' RTSA::synthesize(m.prepare)
#'
synthesize <- function(y,
           sign = NULL,
           fixedStudy = TRUE,
           hksj = FALSE) {
    if (class(y) != "synthPrepped") {
      warning('sig object is not of synthPrepped-class, results may be inaccurate')
    }

    w <- y$w   # collect objects
    sig <- y$sig
    te <- y$te
    pe <- y$pe
    eI <- y$eI
    eC <- y$eC
    nI <- y$nI
    nC <- y$nC
    df <- length(w) - 1
    if (length(w) == 1)
      df <- 1 #check up on this

    if (y$method == "GLM") {
      bi <- nI - eI
      di <- nC - eC

      grpOut <- cbind(xi = c(rbind(eI, eC)), mi = c(rbind(bi, di)))
      k <- length(eI)
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
          lme4::glmer(grpOut ~ -1 + eff + trial + (sv - 1 |
                                                     trial),
                      family = binomial,
                      nAGQ = 7)
      }

      esR <- glmRandomFit@beta[1]
      seR <- sqrt(vcov(glmRandomFit)[1, 1])
      if (fixedStudy == TRUE) {
        tau2 <- lme4::VarCorr(glmRandomFit)[[1]][1]
      } else {
        tau2 <- lme4::VarCorr(glmRandomFit)[[2]][1]
        sigma2 <- lme4::VarCorr(glmRandomFit)[[1]][1]
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
      # NOT GLM

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
        # method the same for OR and RR
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

        pval <- (1 - pnorm(abs(zval))) * 2

      } else {
        lpeF <- sum(log(te) * rw) # fixed effect log pooled estimate
        uci <- exp(lpeF + 1.96 * sqrt(vw))
        lci <- exp(lpeF - 1.96 * sqrt(vw))
        peF <- exp(lpeF)
        if (is.null(sign)) {
          zval <- lpeF / sqrt(vw)
        } else {
          zval <- sign * lpeF / sqrt(vw)
        }
        pval <- (1 - pnorm(abs(zval))) * 2

      }

      if (y$method != "cont") {
        w <- 1 / (sig ^ 2)
        Q <- sum(w * log(te) ^ 2) - (sum(w * log(te))) ^ 2 / sum(w)
        U <- sum(w) - sum(w ^ 2) / sum(w)
        tau2 <-
          ifelse(Q > df, (Q - df) / U, 0) # DerSimonian-Laird estimate
        pQ <- pchisq(Q, df, lower.tail = FALSE)
        if (!is.na(Q) & Q / df <= 0) {
          H <- 0
        } else if (is.na(Q)) {
          H <- NA
        } else {
          H <- sqrt(Q / df)
        }
        I2 <- ifelse((Q - df) / Q >= 0, (Q - df) / Q, 0)

        # random effect
        if (tau2 != 0) {
          # cal. random effect weights and belonging pooled estimate
          wR <- 1 / (sig ^ 2 + tau2)
          vwR <- 1 / sum(wR) # variance of pooled effect (random)
          rwR <- wR * vwR
          teR <- sum(te * wR) / sum(wR)
          leR <- sum(log(te) * rwR)
          peR <- exp(leR)
          if (hksj == TRUE) {
            vwR <- 1 / df * sum(wR * (log(te) - leR) ^ 2 / sum(wR))
            if (is.null(sign)) {
              zvalR <- leR / sqrt(vwR)
            } else {
              zvalR <- sign * leR / sqrt(vwR)
            }
            pvalR <- 2 * pt(abs(zvalR), df = df, lower.tail = FALSE)
            lciR <- exp(leR - qt(1 - 0.05 / 2, df = df) * sqrt(vwR))
            uciR <- exp(leR + qt(1 - 0.05 / 2, df = df) * sqrt(vwR))
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

          vw <- 1 / sum(w)
          D2 <- 1 - vw / vwR
          synth <-
            list(
              fw = round(rw * 100, 1),
              peF = c(peF, lci, uci, zval, pval, lpeF, vw),
              rwR = rwR * 100,
              peR = c(peR, lciR, uciR, zvalR, pvalR, vwR),
              Q = c(Q, df, pQ),
              U = c(tau2, H, I2, D2)
            )
          class(synth) <- "synthesized"
          return(synth)
        } else {
          synth <-
            list(
              peF = c(peF, lci, uci, zval, pval, lpeF, vw),
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
