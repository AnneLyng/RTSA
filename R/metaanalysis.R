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
#' @param study XXX.
#' @param mI,mC,sdI,sdC See details.
#' @param eI,nI,eC,nC See details.
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
<<<<<<< HEAD
#' data(eds)
#' metaanalysis(outcome = "cont", data = eds, study = eds$study)
=======

>>>>>>> 16e9e491d1421b57742269701d323482c9ca689e
metaanalysis <- function(data = NULL,
                         outcome = "RR",
                         vartype = "unequal",
                         method = "MH",
                         fixedStudy = TRUE,
                         hksj = FALSE,
                         sign = NULL,
                         study = NULL,
                         mI = NULL, mC = NULL,
                         sdI = NULL, sdC = NULL,
                         eI = NULL, nI = NULL,
                         eC = NULL, nC = NULL,
                         ...) {

  # check if the correct outcome metric is choosen
  if (!(outcome %in% c("cont", "RD", "RR", "OR"))) {
    stop("Outcome must be either: cont, RD, RR or OR.")
  }

<<<<<<< HEAD
  # check if the input is correct for outcome == "cont"
  if (outcome == "cont") {
    if (!is.null(data)) {
      mI = data$mI
      mC = data$mC
      sdI = data$sdI
      sdC = data$sdC
      nI = data$nI
      nC = data$nC
    }
  }

  # figure out if data is in data set or needs to be specified
  if (outcome != "cont") {
    if (!is.null(data)) {
      eI = data$eI
      nI = data$nI
      eC = data$eC
      nC = data$nC
=======
  #dataframe to variables.
  if(!is.null(data)){
    if(any(colnames(data) == "study")) study = data$study
    if(outcome == "cont"){
      if(any(colnames(data) == "mI")) mI = data$mI
      if(any(colnames(data) == "mC")) mC = data$mC
      if(any(colnames(data) == "sdI")) sdI = data$sdI
      if(any(colnames(data) == "sdC")) sdC = data$sdC
    }else{
      if(any(colnames(data) == "eI")) eI = data$eI
      if(any(colnames(data) == "nI")) nI = data$nI
      if(any(colnames(data) == "eC")) eC = data$eC
      if(any(colnames(data) == "nC")) nC = data$nC
>>>>>>> 16e9e491d1421b57742269701d323482c9ca689e
    }
  }

  #metaprepare
  if(outcome == "cont") {
    mp = metaPrepare(
      outcome = "cont",
<<<<<<< HEAD
      mI = mI,
      mC = mC,
      sdI = sdI,
      sdC = sdC,
      nI = nI,
      nC = nC,
      vartype = vartype,
      method = method
=======
      mI = mI,mC = mC,
      sdI = sdI,sdC = sdC,
      vartype = vartype, method = method
>>>>>>> 16e9e491d1421b57742269701d323482c9ca689e
    )
  }else{
    mp = metaPrepare(
      outcome = outcome,
      eI = eI,eC = eC,
      nI = nI,nC = nC,
      vartype = vartype,
      method = method
    )
  }

  #synthesize
  sy = synthesize(y = mp, sign = sign, fixedStudy = fixedStudy, hksj = hksj)


  #create an output object

  nonevent <- NULL
  missing_vec <- NULL
  if(is.null(study)){
    missing_vec <- 1
    study = 1:length(mp$te)
    if(length(study) != length(eC)){
      nonevent <- "A study"
    }
  }else{
    if(!is.null(mp$nonevent)){
      nonevent <- study[mp$nonevent]
      study <- study[-mp$nonevent]
    }
  }

  studyResults = data.frame(
    study = study,
    "ES" = mp$te,
    "stdError" = sqrt(mp$sig),
    "lowerCI" = mp$lower,
    "upperCI" = mp$upper,
    weightFixed = mp$w,
    weightRandom = ifelse(is.null(sy$rwR), NA, sy$rwR)
  )
  colnames(studyResults)[2] = outcome

  if(!is.null(sy$peR)){
    metaResults = data.frame(
      type = c("Fixed", "Random"),
      "ES" = c(sy$peF[1], sy$peR[1]),
      "stdError" = c(sqrt(sy$peF[7]), sqrt(sy$peR[6])),
      "lowerCI" = c(sy$peF[2], sy$peR[2]),
      "upperCI" = c(sy$peF[3], sy$peR[3]),
      "pValue" = c(sy$peF[5], sy$peR[5])
    )
  }else{
    metaResults = data.frame(
      type = c("Fixed"),
      "ES" = c(sy$peF[1]),
      "stdError" = c(sqrt(sy$peF[7])),
      "lowerCI" = c(sy$peF[2]),
      "upperCI" = c(sy$peF[3]),
      "pValue" = c(sy$peF[5]))
  }
  colnames(metaResults)[2] = outcome

  out = list(studyResults = studyResults, metaResults = metaResults,
             metaPrepare = mp, synthesize = sy, nonevent = nonevent,
             missing_vec = missing_vec)
  class(out) <- "metaanalysis"
  return(out)
}


#' @method print metaanalysis
#' @export
print.metaanalysis <- function(x,...){
  cat("Individual trial results: \n \n")
  y <- x$studyResults
  print(y)
  cat("\n Non-sequential metaanalysis results: \n \n")
  y <- x$metaResults
  print(y)
  invisible(x)
  if(x$metaPrepare$method == "GLM"){
    message("\n NB. Only fixed-effect is analysed since method is GLM")
  }
  #ZERO TRIAL
  if(!is.null(x$missing_vec)){
    message("\n NB. Please provide a study vector to name studies.")
  }
  if(!is.null(x$nonevent)){
    message(paste("\n NB.",x$nonevent,"was excluded from the analysis due to zero events, consider changing outcome to RD"))
  }

}

