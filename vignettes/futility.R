## -----------------------------------------------------------------------------
library(RTSA)
bound_none <- boundaries(timing = c(0.5,0.75, 1), alpha  = 0.025, beta = 0.2, 
                    side = 1, futility = "none", es_alpha = "esOF")
bound_bind <- boundaries(timing = c(0.5,0.75, 1), alpha  = 0.025, beta = 0.2, 
                    side = 1, futility = "binding", es_alpha = "esOF",
                    es_beta = "esOF")
bound_none2 <- boundaries(timing = c(0.5,0.75, 1), alpha  = 0.05, beta = 0.2, 
                    side = 2, futility = "none", es_alpha = "esOF")
bound_bind2 <- boundaries(timing = c(0.5,0.75, 1), alpha  = 0.05, beta = 0.2, 
                    side = 2, futility = "binding", es_alpha = "esOF",
                    es_beta = "esOF")

## ----fig1, warning=FALSE, fig.height=8, fig.cap="Boundaries for benefit and harm (red lines) and boundaries for futility (blue lines) are shown on the 4 plots.", fig.width=10, echo=FALSE----
library(ggplot2)
library(gridExtra)
grid.arrange(plot(bound_none)+ggtitle("a) One-sided - no futility"),
             plot(bound_bind)+ggtitle("b) One-sided - bind. futility"),
             plot(bound_none2)+ggtitle("c) Two-sided - no futility"),
             plot(bound_bind2)+ggtitle("d) Two-sided - bind. futility"),
             ncol = 2)

## -----------------------------------------------------------------------------
bound_none

## -----------------------------------------------------------------------------
bound_bind

## ----fig3, warning=FALSE, fig.height=4, fig.cap="Boundaries for benefit and harm (red lines) and boundaries for futility (blue lines) are shown on the 2 plots.", fig.width=10, echo=FALSE----
grid.arrange(plot(bound_bind)+ggtitle("One-sided - binding futility") + 
               geom_line(aes(x = c(0,0.5,0.75,1)*1.075, y = c(0,1.1, 1.0, 2.1)), line_width = 0.25) + 
               geom_point(aes(x = c(0,0.5,0.75,1)*1.075, y = c(0,1.1, 1.0, 2.1))),
             plot(bound_bind)+ggtitle("One-sided - binding futility") + 
               geom_line(aes(x = c(0,0.5,0.75,1)*1.075, y = c(0,1.7, 2, 2.7)), line_width = 0.25) + 
               geom_point(aes(x = c(0,0.5,0.75,1)*1.075, y = c(0,1.7, 2, 2.7))),
             ncol = 2)

## -----------------------------------------------------------------------------
bound_nbind <- boundaries(timing = c(0.5,0.75, 1), alpha  = 0.025, beta = 0.2, 
                    side = 1, futility = "non-binding", es_alpha = "esOF",
                    es_beta = "esOF")
bound_nbind2 <- boundaries(timing = c(0.5,0.75, 1), alpha  = 0.05, beta = 0.2, 
                    side = 2, futility = "non-binding", es_alpha = "esOF",
                    es_beta = "esOF")

## ----fig2, warning=FALSE, fig.height=8, fig.cap="Boundaries for benefit and harm (red lines) and boundaries for futility (blue lines) are shown on the 4 plots.", fig.width=10, echo=FALSE----
grid.arrange(plot(bound_none)+ggtitle("a) One-sided - no futility"),
             plot(bound_nbind)+ggtitle("b) One-sided - non-bind. futility"),
             plot(bound_none2)+ggtitle("c) Two-sided - no futility"),
             plot(bound_nbind2)+ggtitle("d) Two-sided - non-bind. futility"),
             ncol = 2)

## -----------------------------------------------------------------------------
bound_nbind

