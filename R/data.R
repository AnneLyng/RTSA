#' Dataset of RCTs investigating the effect of 80\% perioperative oxygen vs. 30
#' -35\% perioperative oxygen on surgical site infection.
#'
#' A dataset containing data on seven trials which includes their number of events per treatment group, where intervention is 80\% oxygen and control is 30-35\% oxygen, number of participants in each treatment group and the year of the trial.
#'
#' @format A data frame with 7 rows and 6 variables:
#' \describe{
#'   \item{study}{Name of first author of the trial}
#'   \item{eI}{Number of events in the intervention group (80\% oxygen)}
#'   \item{nI}{Number of pax in the intervention group (80\% oxygen)}
#'   \item{eC}{Number of events in the control group (30-35\% oxygen)}
#'   \item{nC}{Number of pax in the control group (30-35\% oxygen)}
#' }
"perioOxy"

#' A dataset containing data on ...
#'
#' @format A data frame with 6 rows and 5 variables:
#' \describe{
#'   \item{study}{Name of first author of the trial}
#'   \item{eI}{Number of events in the intervention group}
#'   \item{nI}{Number of pax in the intervention group}
#'   \item{eC}{Number of events in the control group}
#'   \item{nC}{Number of pax in the control group}
#' }
"coronary"

#' Early supported discharge services
#'
#' A dataset containing data on the length of hospital stay when receiving early
#' supported discharge (ESD) service versus conventional care. The outcome is
#' length of initial hospital stay counted in days.
#'
#' @format A data frame with 9 studies and 8 variables:
#'
#' @details
#' \itemize{
#'   \item study. Name of the city of the study
#'   \item year. Year of the trial
#'   \item mI. Mean duration at hospital in intervention (ESD) group
#'   \item mC. Mean duration at hospital in control group
#'   \item sdI. Standard deviation of intervention (ESD) estimate
#'   \item sdC. Standard deviation of control estimate
#'   \item nI. Number of participants in the intervention (ESD) group
#'   \item nC. Number of participants in the control group
#' }
#'
#' @references Fearon P, Langhorne P. Services for reducing duration of hospital care for acute stroke patients. Cochrane Database of Systematic Reviews 2012, Issue 9. Art. No.: CD000443. DOI: 10.1002/14651858.CD000443.pub3. Accessed 17 October 2022.
"eds"

