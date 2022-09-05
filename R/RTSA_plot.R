#' plot.RTSA
#'
#' @param x RTSA object
#' @param ... Other arguments to plot.RTSA
#'
#' @return Plot. Either a plot for two sided testing or one-sided
#' @export
#'
#' @importFrom ggplot2 ggplot coord_cartesian geom_hline theme_bw geom_vline geom_line geom_point aes theme element_blank geom_ribbon xlab ylab scale_x_continuous expansion scale_y_continuous scale_fill_identity scale_colour_manual ggtitle
#'
#' @examples
#' data(perioOxy)
#' outRTSA <- RTSA(data = perioOxy, outcome = "RR", mc = 0.7, side = 2, alpha = 0.05, fixed = FALSE)
#' plot(x = outRTSA)
#'
plot.RTSA = function(x, ...){

  if(is.null(dim(x$zvalues))){
    indi <- 1
    zval <- x$zvalues
  } else {
    indi <- 2
    zval <- x$zvalues[2,]
  }

  andat <- data.frame(zval = c(0,zval),
                      infFrac = c(0,x$orgTiming))
  bounddat <- rbind(c(0,20,-20),data.frame(do.call(cbind, x$boundout)))
  if(max(andat$infFrac) > max(bounddat$informationFractions)){
    bounddat <- rbind(bounddat, c(max(andat$infFrac), round(qnorm(1-x$alpha/x$side),2),
                                  -round(qnorm(1-x$alpha/x$side),2)))
  }
  gplot <- ggplot(data = bounddat) +
    coord_cartesian(ylim=c(-10, 10), xlim = c(0,max(x$orgTiming,1.1))) + theme_bw() +
    geom_hline(yintercept = c(round(qnorm(1-x$alpha/x$side),2),
                              -round(qnorm(1-x$alpha/x$side),2), 0),
               cex = 0.25, col = c("green", "green", "gray")) +
    geom_vline(xintercept = 1, cex = 0.25, col = "black") +
    geom_line(aes(x = informationFractions, y = alpha.boundaries.upper, col = "red"),
              cex = 0.25) +
    geom_point(aes(x = informationFractions, y = alpha.boundaries.upper), col = "red", cex = 1.25) +
    {if(x$side == 2)geom_line(aes(x = informationFractions, y = alpha.boundaries.lower), col = "red",
                                    cex = 0.25)} +
    {if(x$side == 2)geom_point(aes(x = informationFractions,
                                         y = alpha.boundaries.lower), col = "red", cex = 1.25)} +
    geom_line(data = andat, aes(x = infFrac,
                                y = zval, col = "blue"), cex = 0.25) +
    geom_point(data = andat, aes(x = infFrac,
                                 y = zval), cex = 1.25) +
    theme(panel.border = element_blank()) +
    {if(x$side == 2)geom_ribbon(aes(ymax = alpha.boundaries.lower,
                                          ymin = -20,
                                          x = informationFractions, fill = "red"), alpha = 0.25)} +
    geom_ribbon(aes(ymax = alpha.boundaries.upper,
                    ymin = 20,
                    x = informationFractions, fill = "green"), alpha = 0.25) +
    xlab("Information fraction") + ylab("Cumulative Z-score")  +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    scale_fill_identity(name = 'Stopping areas', guide = 'legend',labels = c('Benefit', 'Harm')) +
    theme(legend.position = 'bottom',  legend.box="vertical") +
    scale_colour_manual(name = 'Lines', values =c('blue'='blue', 'green'='green','red' = 'red', 'black' = 'black'),
                        labels = c('Cumulative Z-score', 'Naive stopping boundaries', 'Sequential boundaries', 'RIS reached')) +
    ggtitle(paste0("TSA plot - ", ifelse(x$side == 1, "one-sided ",
                                         "two-sided "), ifelse(indi == 1, "fixed-effect ",
                                                               "random-effects "), "model"))
  gplot
}
