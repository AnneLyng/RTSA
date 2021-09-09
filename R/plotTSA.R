#' plot.RTSA
#'
#' @param x RTSA object
#' @param ... Other arguments to plot.RTSA
#'
#' @return Plot. Either a plot for two sided testing or one-sided
#' @export
#'
#' @examples
#' data(perioOxy)
#' count <- cumsum(perioOxy$nI+perioOxy$nC)
#' timing <- c(count/5540,1)
#' mp = metaPrepare(outcome = "RR", eI = perioOxy$eI, nI = perioOxy$nI,
#' eC = perioOxy$eC, nC = perioOxy$nC, method = "IV")
#' RTSAout = RTSA(timing = timing, synth = mp, anaTimes = c(4,5,7,8),
#' side = 2, alpha = 0.05)
#' plot(x = RTSAout)
#'
plot.RTSA = function(x, ...){
  if(x$side == 2){
    plot(x = c(0,x$boundout$informationFractions), y =
           c(8, x$boundout$alpha.boundaries.upper), xlim = c(0,1),
         ylim = c(-8,8), type = "l", col = "red",
         xlab = "Information fraction",
         ylab = "Z-value")
    points(x = c(0,x$boundout$informationFractions), y =
             c(8, x$boundout$alpha.boundaries.upper), pch = 20)
    lines(x = c(0,x$boundout$informationFractions), y =
            c(-8, x$boundout$alpha.boundaries.lower), col = "red")
    points(x = c(0,x$boundout$informationFractions), y =
             c(-8, x$boundout$alpha.boundaries.lower), pch = 20)
    abline(h = 0)
    abline(h = 1.96, col = "red")
    abline(h = -1.96, col = "red")
    lines(x = c(0,x$boundout$informationFractions[-4]), y =
            c(0, x$zvalues[2,]), col = "blue")
    points(x = c(0,x$boundout$informationFractions[-4]), y =
             c(0, x$zvalues[2,]), pch = 20)
    lines(x = c(0,x$boundout$informationFractions[-4]), y =
            c(0, x$zvalues[1,]), col = "blue", lty = 2)
    points(x = c(0,x$boundout$informationFractions[-4]), y =
             c(0, x$zvalues[1,]), pch = 20)
  } else {
    plot(x = c(0,x$boundout$informationFractions), y =
           c(8, x$boundout$alpha.boundaries.upper), xlim = c(0,1),
         ylim = c(-8,8), type = "l", col = "red",
         xlab = "Information fraction",
         ylab = "Z-value")
    points(x = c(0,x$boundout$informationFractions), y =
             c(8, x$boundout$alpha.boundaries.upper), pch = 20)
    abline(h = 0)
    abline(h = 1.96, col = "red")
    lines(x = c(0,x$boundout$informationFractions[-4]), y =
            c(0, x$zvalues[2,]), col = "blue")
    points(x = c(0,x$boundout$informationFractions[-4]), y =
             c(0, x$zvalues[2,]), pch = 20)
    lines(x = c(0,x$boundout$informationFractions[-4]), y =
            c(0, x$zvalues[1,]), col = "blue", lty = 2)
    points(x = c(0,x$boundout$informationFractions[-4]), y =
             c(0, x$zvalues[1,]), pch = 20)
  }
}
