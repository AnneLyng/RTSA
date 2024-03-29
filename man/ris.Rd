% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ris.R
\name{ris}
\alias{ris}
\alias{print.ris}
\title{Calculate required sample and trials size.}
\usage{
ris(
  outcome,
  mc,
  side = 2,
  alpha = 0.05,
  beta = 0.1,
  fixed = TRUE,
  sd_mc = NULL,
  pC = NULL,
  p1 = NULL,
  ma = NULL,
  tau2 = NULL,
  I2 = NULL,
  D2 = NULL,
  type = "prospective",
  trials = NULL,
  RTSA = FALSE,
  ...
)
}
\arguments{
\item{outcome}{Choose between: "MD" (mean difference), "RR" (relative risk), "OR" (odds ratio) or "RD" (risk difference).}

\item{mc}{Minimum clinical relevant effect. For "OR" or "RR" set to natural scale, not log scale.}

\item{side}{Test type. Set to 1 or 2 depending on the test being 1- or 2-sided.}

\item{alpha}{The level of type I error as a percentage, the default is 0.05 corresponding to 5\%.}

\item{beta}{The level of type II error as a percentage, the default is 0.1 corresponding to 10\%.}

\item{fixed}{Should sample size be based on a fixed-effect (TRUE) or random-effects (FALSE) model. Defaults to TRUE.}

\item{sd_mc}{Standard deviation of estimated effect. Only needed when outcome type is "MD".}

\item{pC}{Probability of event in control group. Only needed when outcome type is "OR", "RR" or "RD".}

\item{p1}{Probability of event in treatment group. Only needed when outcome type is "RD".}

\item{ma}{An optional \code{metaanalysis} object. Required for retrospective sample size calculations.}

\item{tau2}{The value of the heterogeneity. Use when estimating the sample size under a random effects model. If data is provided, the estimated heterogeneity is used instead.}

\item{I2}{Optional argument. Inconsistency.}

\item{D2}{Optional argument. Diversity.}

\item{type}{Whehter the type of calculaiton is for "prospective" meta-analysis or "retrospective" meta-analysis. If the type is retrospective, one should add a meta-analysis object to the function. See argument ma.}

\item{trials}{Optional numeric argument. If one is interested in a specific number of trials.}

\item{RTSA}{Whether the ris function was called via the RTSA function. Purely operational argument.}

\item{...}{additional arguments}
}
\value{
A list of up to 6 elements:
\item{settings}{A list containing the arguments provided to the \code{ris} function.}
\item{NF}{The total number of required participants in a fixed-effect meta-analysis if type is prospective. Contains a list if the type is retrospective, where \code{NF} is the additional required number of participants and \code{NF_full} is the total required number of participants.}
\item{NR_tau}{A list containing: \code{minTrial} the minimum number of trials. \code{nPax} a matrix containing four possible number of trials with the number of participants per trial and total number of participants. \code{tau2} the estimate used for the calculation. Might contain \code{NR_tau_ll} and \code{NR_tau_ul} which contain the same three elements. \code{NR_tau_ll} is based on the lower value in the confidence interval of tau2. \code{NR_tau_ul} is based on the upper value in the confidence interval for tau2. If the type is prospective the numbers are the total required. If the type is retrospective the numbers are the additional required.}
\item{NR_D2}{The total number of required participants in a random-effects meta-analysis adjusted by diversity (\code{D2}) if type is prospective. Contains a list if the type is retrospective, where \code{NR_D2} is the additional required number of participants and \code{NR_D2_full} is the total required number of participants.}
\item{NR_I2}{The total number of required participants in a random-effects meta-analysis adjusted by inconsistency (\code{I2}) if type is prospective. Contains a list if the type is retrospective, where \code{NR_I2} is the additional required number of participants and \code{NR_I2_full} is the total required number of participants.}
}
\description{
Calculate required sample and trials size.
}
\examples{
# Sample and trial size calculation for prospective meta-analysis
ris(outcome = "RR", mc = 0.8, pC = 0.12, fixed = TRUE, alpha = 0.05,
beta = 0.1, side = 2)

# Additional sample and trial size calculation for retrospective meta-analysis 
# It is calculated directly from the metaanalysis() function
data("perioOxy")
ma <- metaanalysis(outcome = "RR", data = perioOxy, mc = 0.8, beta = 0.2)
ma$ris
# Or by using the two functions in sequence
ma <- metaanalysis(outcome = "RR", data = perioOxy)
ris(outcome = "RR", mc = 0.8, ma = ma, type = "retrospective", fixed = FALSE,
 beta = 0.2, alpha = 0.05, side = 2)
}
