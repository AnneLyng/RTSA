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

  if(x$side == 2 | (x$side == 1 & !is.null(x$futility))){
    bounddat <- rbind(c(0,-20,20),data.frame(inf_frac = x$boundout$inf_frac, alpha_ubound = -x$boundout$alpha_ubound,
                                             alpha_lbound = -x$boundout$alpha_lbound))
  } else {
    bounddat <- rbind(c(0,-20),data.frame(inf_frac = x$boundout$inf_frac, alpha_ubound = -x$boundout$alpha_ubound))
  }

  if(max(andat$infFrac) > max(bounddat$inf_frac)){
    bounddat <- rbind(bounddat, c(max(andat$infFrac), round(qnorm(1-x$alpha/x$side),2),
                                  -round(qnorm(1-x$alpha/x$side),2)))
  }

  bounddat$RIS <- round(c(0,x$boundout$inf_frac)*x$adjRIS)
  bounddat$RIS[1:(nrow(bounddat)-2)] <- NA
  bounddat$RIS[nrow(bounddat)-1] <- x$AIS
  bounddat$labels <- c(rep("",nrow(bounddat)-2),"AIS","DARIS")


ggplot(data = bounddat) +
    coord_cartesian(xlim = c(0,max(inf_frac,1.1)),
                    ylim =c(ifelse("alpha_lbound" %in% colnames(bounddat), bounddat$alpha_lbound[bounddat$alpha_lbound < 20], 2.5), -10)) +
    geom_segment(x=0,xend=1, y=round(qnorm(1-x$alpha/x$side),2),
      yend = round(qnorm(1-x$alpha/x$side),2), cex = 0.25,
      col = "#006400", linetype="dashed") +
    geom_segment(x=0,xend=1, y=-round(qnorm(1-x$alpha/x$side),2),
      yend = -round(qnorm(1-x$alpha/x$side),2), cex = 0.25,
      col = "#006400", linetype="dashed") +
    geom_segment(x=0,xend=1, y=0, yend = 0,
      cex = 0.25, col = "gray", linetype="solid") +
    geom_vline(xintercept = 1, cex = 0.25, col = "black") +
    geom_line(aes(x = inf_frac, y = alpha_ubound, col = "red"), cex = 0.25) +
    geom_point(aes(x = inf_frac, y = alpha_ubound), col = "black", cex = 1) +
  {if(x$side == 2)
    geom_line(aes(x = inf_frac, y = alpha_lbound), col = "red", cex = 0.25)
    } +
  {if(x$side == 1 & !is.null(x$futility))
    geom_line(aes(x = inf_frac, y = alpha_lbound), col = "red", cex = 0.25,
              linetype = "dashed")
  } +
  {if(x$side == 2 | (x$side == 1 & !is.null(x$futility)))
      geom_point(aes(x = inf_frac,y = alpha_lbound), col = "black", cex = 1)
  } +
    geom_line(data = andat, aes(x = infFrac,y = zval, col = "blue"), cex = 0.25) +
    geom_point(data = andat, aes(x = infFrac,y = zval), cex = 1.25) +
    geom_label(aes(x=inf_frac, y=-10, label=RIS),vjust=1,label.size=NA) +
    scale_x_continuous(expand = expansion(0), name="Information fraction",
                       sec.axis = sec_axis(~., breaks=bounddat$inf_frac,labels=bounddat$labels)) +
    scale_y_reverse(expand = expansion(0), "Cumulative Z-score") +
    theme(legend.position = 'bottom',  legend.box="vertical") +
    scale_colour_manual(name = ' ', values =c('blue'='blue', '#006400'='#006400','red' = 'red', 'black' = 'black'),
                        labels = c('Z-score', 'Naive boundaries', 'Sequential boundaries', 'DARIS')) +
  scale_linetype_manual(values = c("solid","dashed", "solid", "solid")) +
    # labs(title=title, subtitle="Proportion in control group, ") +
    theme_classic() +
    theme(legend.position="bottom",
          legend.title = element_blank(),
          legend.margin = margin(),
          legend.text = element_text(size=9),
          plot.title = element_text(hjust = 0.5),
          plot.title.position = "plot",
          plot.subtitle = element_text(hjust=0.5),
          axis.title.y = element_text(color="black",size=10),
          axis.title.x = element_text(color="black",size=10),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          axis.title.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.line.x.top = element_blank())
}


