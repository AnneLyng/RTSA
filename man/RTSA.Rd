% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RTSA.R
\name{RTSA}
\alias{RTSA}
\alias{print.RTSA}
\title{R version of Trial Sequential Analysis. Used for designing and analysing sequential meta-analyses.}
\usage{
RTSA(
  type = "design",
  outcome = NULL,
  side = 2,
  alpha = 0.05,
  beta = 0.1,
  futility = "none",
  es_alpha = "esOF",
  es_beta = NULL,
  timing = NULL,
  data = NULL,
  design = NULL,
  ana_times = NULL,
  fixed = FALSE,
  mc = NULL,
  RRR = NULL,
  sd_mc = NULL,
  pC = NULL,
  weights = "MH",
  re_method = "DL_HKSJ",
  tau_ci_method = "BJ",
  gamma = NULL,
  rho = NULL,
  study = NULL,
  cont_vartype = "equal",
  zero_adj = 0.5,
  tau2 = NULL,
  I2 = NULL,
  D2 = NULL,
  trials = NULL,
  final_analysis = NULL,
  inf_type = "sw",
  conf_level = 0.95,
  random_adj = "tau2",
  power_adj = TRUE,
  ...
)
}
\arguments{
\item{type}{Type of RTSA. Options are "design" or "analysis".}

\item{outcome}{Outcome metric. Options are: RR (risk ratio/relative risk), OR (odds ratio), RD (risk difference) and MD (mean difference).}

\item{side}{Whether a 1- or 2-sided hypothesis test is used. Options are 1 or 2. Default is 2.}

\item{alpha}{The level of type I error as a percentage, the default is 0.05 corresponding to 5\%.}

\item{beta}{The level of type II error as a percentage, the default is 0.1 corresponding to 10\%.}

\item{futility}{Futility boundaries added to design. Options are: none, non-binding and binding. Default is "none".}

\item{es_alpha}{The spending function for alpha-spending. Options are: esOF (Lan & DeMets version of O'Brien-Fleming), esPoc (Lan & DeMets version of Pocock), HSDC (Hwang Sihi and DeCani) and rho (rho family).}

\item{es_beta}{The spending function for beta-spending. For options see es_alpha.}

\item{timing}{Expected timings of interim analyses when type = "design". Defaults to NULL.}

\item{data}{A data.frame containing the study results. The data set must containing a specific set of columns. These are respectively `eI` (events in intervention group), `eC` (events in control group), `nC` (participants intervention group) or `nI` (participants control group) for discrete data, or, `mI` (mean intervention group), `mC` (mean control group), `sdI` (standard error intervention group), `sdC` (standard error control group),`nC` (participants intervention group) and `nI` (participants control group)  for continuous outcomes. Preferable also a `study` column as an indicator of study.}

\item{design}{RTSA object where type is design.}

\item{ana_times}{An optional vector of analysis times. Used if the sequential analysis is not done for all studies included in the meta-analysis.}

\item{fixed}{Should only a fixed-effect meta-analysis be computed. Default is FALSE.}

\item{mc}{Minimal clinical relevant outcome value}

\item{RRR}{Relative risk reduction. Used for binary outcomes with outcome metric RR. Argument mc can be used instead. Must be a value between 0 and 1.}

\item{sd_mc}{The expected standard deviation. Used for sample size calculation for mean differences.}

\item{pC}{The expected probability of event in the control group. Used for sample size calculation for binary outcomes.}

\item{weights}{Weighting method options include IV (inverse-variance) and MH (Mantel-Haenszel). Defaults to MH.}

\item{re_method}{Method for calculating the estimate of heterogeneity, tau^2, and the random-effects meta-analysis variance. Options are "DL" for DerSimonian-Laird and "DL_HKSJ" for the Hartung-Knapp-Sidik-Jonkman adjustment of the DerSimonian-Laird estimator.}

\item{tau_ci_method}{Method for calculating confidence intervals for the estimated heterogeneity tau^2. Options are "QP" for Q-profiling and "BJ" for Biggelstaff ....}

\item{gamma}{Parameter for the HSDC error spending function.}

\item{rho}{Parameter for the rho family error spending function.}

\item{study}{An optional vector of study names and perhaps year of study. Defaults to NULL.}

\item{cont_vartype}{For mean difference outcomes, do we expect the variance in the different groups to be "equal" or "non-equal".}

\item{zero_adj}{Zero adjustment. Options for now is 0.5.}

\item{tau2}{Heterogeneity estimate. Used for sample and trial size calculation. Defaults to NULL.}

\item{I2}{Inconsistency estimate. Used for sample and trial size calculation. Defaults to NULL.}

\item{D2}{Diversity estimate. Used for sample and trial size calculation. Defaults to NULL.}

\item{trials}{Number of anticipated extra trials. Used for heterogeneity adjustment by tau2.}

\item{final_analysis}{Whether or not the current analysis is the final analysis.}

\item{inf_type}{Stopping time confidence interval. Options for now is sw (stage-wise).}

\item{conf_level}{Confidence level on stopping time confidence interval.}

\item{random_adj}{The sample size adjustment based on presence of heterogeneity. Options are "D2" (Diversity), "I2" (Inconsistency) and "tau2" (the heterogeneity estimate). Default is "tau2".}

\item{power_adj}{Whether the sample size should be adjusted by the sequential design. Defaults to TRUE.}

\item{...}{other arguments}
}
\value{
A RTSA object, a list of five elements:
\item{settings}{A list containing all of the settings used in the \code{RTSA} call. See Arguments.}
\item{ris}{List containing sample and trial size calculations for a non-sequential meta-analysis. See documentation for \code{ris} function.}
\item{bounds}{List of stopping boundaries, timing of trials and more. See documentation for \code{boundaries} function.}
\item{results}{List of 3 to 7 elements. \code{AIS} Achieved information size. \code{RIS} Fixed-effect required information size for a non-sequential meta-analysis. \code{SMA_RIS} RIS adjusted for sequential analysis. \code{HARIS} Heterogeneity adjusted required information size for a non-sequential meta-analysis. \code{SMA_HARIS} HARIS adjusted for sequential analysis. \code{results_df} a data.frame of inference, see documentation for \code{inference} function.  \code{seq_inf} a list of conditional inference, see documentation for \code{inference} function. \code{metaanalysis} A metaanalysis object, see documentation for \code{metaanalysis} function. \code{design_df} a data.frame containing the stopping boundaries and timings from the design.}
\item{warnings}{List of warnings}
}
\description{
R version of Trial Sequential Analysis. Used for designing and analysing sequential meta-analyses.
}
\examples{
\dontrun{
### Retrospective sequential meta-analysis:
# A RRR of 20\% is expected which gives mc = 1 - RRR = 0.8. 
# No futility boundaries
data(perioOxy)
RTSA(type = "analysis", data = perioOxy, outcome = "RR", mc = 0.8, side = 2,
 alpha = 0.05, beta = 0.2, es_alpha = "esOF")

# Set binding futility boundaries
# And use Lan and DeMets' version of Pocock stopping boundaries 
RTSA(type = "analysis", data = perioOxy, outcome = "RR", mc = 0.8, side = 2,
 alpha = 0.05, beta = 0.2, es_alpha = "esOF", futility = "binding",
 es_beta = "esPoc", random_adj = "D2")

# Set non-binding futility boundaries
RTSA(type = "analysis", data = perioOxy, outcome = "RR", mc = 0.8, side = 2,
 alpha = 0.05, beta = 0.2, es_alpha = "esOF", futility = "non-binding",
 es_beta = "esOF")
 
### Design a prospective sequential meta-analysis
# For continuous data without expected heterogeneity
RTSA(type = "design", outcome = "MD", mc = 5, sd_mc = 10, side = 1, 
timing = c(0.33, 0.66, 1), fixed = TRUE,
alpha = 0.025, beta = 0.1, es_alpha = "esOF", futility = "non-binding", 
es_beta = "esPoc")

# For binary outcome
RTSA(type = "design", outcome = "RR", mc = 0.75, side = 1, 
timing = c(0.33, 0.66, 1), pC = 0.1, D2 = 0.1, 
alpha = 0.025, beta = 0.2, es_alpha = "esOF", futility = "non-binding", 
es_beta = "esOF")

# extract sample size calculation
out_rtsa <-  RTSA(type = "design", outcome = "RR", mc = 0.75, side = 1, 
timing = c(0.33, 0.66, 1), pC = 0.1, D2 = 0.1, 
alpha = 0.025, beta = 0.2, es_alpha = "esOF", futility = "non-binding", 
es_beta = "esOF")
out_rtsa$ris

# plot the design
plot(out_rtsa)

# update the design with data as it accumulates (here toy-data)
fake_data <- data.frame(eI = c(10,10), eC = c(13, 11), nI = c(750, 750),
nC = c(750,750))
RTSA(type = "analysis", design = out_rtsa, data = fake_data)

# plot the analysis
an_rtsa <- RTSA(type = "analysis", design = out_rtsa, data = fake_data)
plot(an_rtsa)
}
}
