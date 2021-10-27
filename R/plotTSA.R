#' plot.RTSA
#'
#' @param x RTSA object
#' @param ylim y axis limits
#' @param ... Other arguments to plot.RTSA
#'
#' @return Plot. Either a plot for two sided testing or one-sided
#' @export
#'
#' @importFrom graphics abline lines points
#'
#' @examples
#' data(perioOxy)
#' RTSAout = RTSA(data = perioOxy, outcome = "RR", mc = 0.9)
#' plot(x = RTSAout)
#'
plot.RTSA = function(x, ylim = c(-8,8),...){
  if(x$side == 2){
    plot(x = c(0,x$boundout$informationFractions), y =
           c(20, x$boundout$alpha.boundaries.upper), xlim = c(0,1),
         ylim = ylim, type = "l", col = "red",
         xlab = "Information fraction",
         ylab = "Z-value")
    points(x = c(0,x$boundout$informationFractions), y =
             c(20, x$boundout$alpha.boundaries.upper), pch = 20)
    lines(x = c(0,x$boundout$informationFractions), y =
            c(-20, x$boundout$alpha.boundaries.lower), col = "red")
    points(x = c(0,x$boundout$informationFractions), y =
             c(-20, x$boundout$alpha.boundaries.lower), pch = 20)
    abline(h = 0)
    abline(h = 1.96, col = "red")
    abline(h = -1.96, col = "red")
    if(x$boundout$informationFractions[length(x$boundout$informationFractions)] ==1 ){
      index = -length(x$boundout$informationFractions)
    } else {
      index = 1:length(x$boundout$informationFractions)
    }
    lines(x = c(0,x$boundout$informationFractions[index]), y =
            c(0, x$zvalues[2,]), col = "blue")
    points(x = c(0,x$boundout$informationFractions[index]), y =
             c(0, x$zvalues[2,]), pch = 20)
    lines(x = c(0,x$boundout$informationFractions[index]), y =
            c(0, x$zvalues[1,]), col = "blue", lty = 2)
    points(x = c(0,x$boundout$informationFractions[index]), y =
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
