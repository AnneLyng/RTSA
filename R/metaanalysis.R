#' metaanalysis
#'
#' @param data Data.frame. If data is provided as a data set. Dataset must then containing arguments for
#'  meta-analysis. Either `study`, `eI`, `eC`, `nC` or `nI` for discrete data, or, `study`, `mI`, `mC`, `sdI` and `sdC` for cont.
#'  See details. If a data frame is not provided, the columns must be set in the function.
#' @param outcome Outcome metric for the studies. Choose between: cont, RR, RD or OR
#' @param vartype Variance type for continuous outcomes. Choices are "equal" or "non-equal". Defaults to "equal".
#' @param method Method for calculating weights. Options include: MH (Mantel-Haenzel), Inverse variance weighting (IV)
#'  or GLM
#' @param fixedStudy TRUE or FALSE. For using GLM methods for the meta-analysis, should the study effect be
#' fixed-effect. Defaults to TRUE.
#' @param hksj TRUE or FALSE. Should the Hartung-Knapp-Sidik-Jonkman adjustment be used to the random-effects
#'  meta-analysis. Defaults to FALSE.
#' @param sign TBA.
#' @param ... Additional variables. See Details.
#'
#' @return A list object containing two data frames.
#' 1. studyResults contains information about the individual studies
#' 2. metaResults contains information about the meta-analysis (traditional, non-sequential)
#' @export
#' @aliases print.metaanalysis
#'
#' @examples
#' data(perioOxy)
#' metaanalysis(outcome = "RR", data = perioOxy, study = perioOxy$trial)
metaanalysis <- function(data = NULL,
                         outcome = "RR",
                         vartype = "equal",
                         method = "MH",
                         fixedStudy = TRUE,
                         hksj = FALSE,
                         sign = NULL,
                         ...) {

  #skal vi lige have ordnet når vi er enige om trial eller study.
  study <- NULL
  # check if the correct outcome metric is choosen
  if (!(outcome %in% c("cont", "RD", "RR", "OR"))) {
    stop("Outcome must be either: cont, RD, RR or OR.")
  }

  # check if the input is correct for outcome == "cont"
  if (outcome == "cont") {
    if (!is.null(data)) {
      mI = data$mI
      mC = data$mC
      sdI = data$sdI
      sdC = data$sdC
    }
  }

  # figure out if data is in data set or needs to be specified
  if (outcome != "cont") {
    if (!is.null(data)) {
      eI = data$eI
      nI = data$nI
      eC = data$eC
      nC = data$nC
    }
  }

  if (outcome == "cont") {
    mp = metaPrepare(
      outcome = "cont",
      mI = mI,
      mC = mC,
      sdI = sdI,
      sdC = sdC,
      vartype = vartype,
      method = method
    )
  } else {
    mp = metaPrepare(
      outcome = outcome,
      eI = eI,
      eC = eC,
      nI = nI,
      nC = nC,
      vartype = vartype,
      method = method
    )
  }

  sy = synthesize(y = mp, sign = sign, fixedStudy = fixedStudy, hksj = hksj)

  # create an output object
  # Her skal vi nok blive enige om at være konsekvente - Trial eller study.
  if(!is.null(data$trial)) colnames(data)[colnames(data) == "trial"] <- "study"
  if (is.null(study) & is.null(data$study)) {
    trial = 1:length(mp$te)
    message("Provide a study vector to name studies.")
  } else {
    if (is.null(study))
      trial = data$study
    if (is.null(data$study))
      trial = study
  }

  studyResults = data.frame(
    study = trial,
    "ES" = mp$te,
    "stdError" = sqrt(mp$sig),
    "lowerCI" = mp$lower,
    "upperCI" = mp$upper,
    weightFixed = mp$w,
    weightRandom = ifelse(is.null(sy$rwR), NA, sy$rwR)
  )
  colnames(studyResults)[2] = c(outcome)

  if(!is.null(sy$peR)){
    metaResults = data.frame(
      type = c("Fixed", "Random"),
      "ES" = c(sy$peF[1], sy$peR[1]),
      "stdError" = c(sqrt(sy$peF[7]), sqrt(sy$peR[6])),
      "lowerCI" = c(sy$peF[2], sy$peR[2]),
      "upperCI" = c(sy$peF[3], sy$peR[3]),
      "pValue" = c(sy$peF[5], sy$peR[5])
    )} else {
      metaResults = data.frame(
        type = c("Fixed"),
        "ES" = c(sy$peF[1]),
        "stdError" = c(sqrt(sy$peF[7])),
        "lowerCI" = c(sy$peF[2]),
        "upperCI" = c(sy$peF[3]),
        "pValue" = c(sy$peF[5]))
    }
  colnames(metaResults)[2] = c(outcome)

  out = list(studyResults = studyResults, metaResults = metaResults, metaPrepare = mp, synthesize = sy)
  class(out) <- "metaanalysis"
  return(out)
}


#' @method print metaanalysis
#' @export
print.metaanalysis <- function(x,...){
  cat("Individual trial results: \n \n")
  y <- x$studyResults
  print(y)
  cat("\n Metanalysis results: \n \n")
  y <- x$metaResults
  print(y)
  invisible(x)
  if(testo$metaPrepare$method == "GLM"){
    cat("\n NB. Only fixed-effect is analysed since method is GLM")
  }

}

