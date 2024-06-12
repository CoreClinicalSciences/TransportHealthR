#' @title Transportability analysis using g-computation
#' 
#' @description
#' Estimates the coefficients of a marginal structural model (MSM) using g-computation in a generalizability or transportability analysis. In particular, the estimators should be unbiased for the coefficients in the superpopulation or the target population, respectively.
#' 
#' @param msmFormula A formula for the MSM to be fitted, which usually includes the outcome, the treatment and any effect modifiers.
#' @param preparedModel A \code{transportGCPreparedModel} object. This is obtained by using the \code{transportGCPreparedModel} function to fit an outcome model using the study data.
#' @param targetData A target dataset.
#' @param bootstrapNum Number of bootstrap datasets to simulate to obtain robust variance estimators
#' 
#' @details
#' The expected workflow is as follows:
#' 
#' \enumerate{
#'  \item{A researcher who wants to perform a generalizability/transportability analysis collects data from the target population.}
#'  \item{They then request the owner of the study data from which they wish to generalize/transport to provide an outcome model fitted using the study data.}
#'  \item{The owner of the study data runs the \code{transportGCPreparedModel} function on the study data to obtain a \code{transportGCPreparedModel} object which contains the fitted outcome model.}
#'  \item{The owner of the study data provides the \code{transportGCPreparedModel} object to the researcher, perhaps in a \code{.rds} file.}
#'  \item{The researcher uses this function and the provided \code{transportGCPreparedModel} object to perform the analysis using g-computation.}
#' }
#' 
#' Since model-fitting objects in \code{R} often contain the data used to fit the model, the \code{transportGCPreparedModel} function wipes this data in the model-fitting object and keeps additional information about the name of the response variable, the name of the treatment variable and the levels of treatment. This is to comply with government regulations regarding access and integration of data sources from different countries.
#' 
#' The MSM-fitting functions do not provide correct standard errors as-is. Bootstrap is used to calculate robust variance estimates of the MSM coefficient estimators. Note that these standard errors are only valid conditional on the observed study data because it is not possible to resample the study data when access to it is restricted.
#'
#' @return
#' A \code{transportGC} object with the following components:
#' * \code{msm}: Raw model fit object of the MSM. Its class will be the same as that of \code{outcomeModel} in the provided \code{transportGCPreparedModel} object.
#' 
#' @export
#' 
#' @md
transportGC <- function (msmFormula,
                         preparedModel,
                         targetData, bootstrapNum = 500) {
  transportGCResult <- transportGCFit(msmFormula,
                                      preparedModel,
                                      targetData)
  
  bootstrapEstimates <- t(sapply(1:bootstrapNum,
                                 function (x) {
                                   targetBoot <- targetData[sample.int(n = nrow(targetData), replace = T), ]
                                   
                                   resultBoot <- transportGCFit(msmFormula,
                                                                preparedModel,
                                                                targetData = targetBoot)
                                   
                                   return(c(resultBoot$msm$coefficients, resultBoot$msm$zeta))
                                 }))
  
  varMatrix <- var(bootstrapEstimates)
  colnames(varMatrix) <- rownames(varMatrix) <- c(names(transportGCResult$msm$coefficients), names(transportGCResult$msm$zeta))
  transportGCResult$msm$var <- varMatrix
  
  return(transportGCResult)
}

transportGCFit <- function (msmFormula,
                         preparedModel,
                         targetData) {
  if (!is.transportGCPreparedModel(preparedModel)) stop("Please provide a transportGCPreparedModel object.")
  if (preparedModel$response != all.vars(msmFormula)[1]) stop("Response variable does not match transportGCPreparedModel object.")
  
  response <- preparedModel$response
  treatment <- preparedModel$treatment
  outcomeModel <- preparedModel$outcomeModel
  
  targetDataReferenceList <- list()
  for (level in preparedModel$treatmentLevels) {
    targetDataReferenceList[[level]] <- targetData
    targetDataReferenceList[[level]][[treatment]] <- level
  }
  
  targetDataCounterfactualFrame <- Reduce(function(d1, d2) rbind(d1, d2), targetDataReferenceList)
  targetDataCounterfactualFrame[[treatment]] <- as.factor(targetDataCounterfactualFrame[[treatment]])
  
  # TODO: adapt to coxph since it doesn't give survival time predictions - see Lee (2024) (transporting survival of hiv...) and Chen & Tsiatis (2001)
  if (!inherits(preparedModel$outcomeModel, "coxph")) {
    targetDataCounterfactualFrame[[response]] <- stats::predict(outcomeModel,
                                                       newdata = targetDataCounterfactualFrame,
                                                       type = "response")
    
    if (inherits(outcomeModel, "glm"))
      msm <- glm(msmFormula, family = outcomeModel$family, data = targetDataCounterfactualFrame)
    else if (inherits(outcomeModel, "survreg"))
      msm <- survival::survreg(msmFormula, data = targetDataCounterfactualFrame)
    else if (inherits(outcomeModel, "polr"))
      msm <- update(outcomeModel, .formula = msmFormula, data = targetDataCounterfactualFrame, evaluate = T)
  }
  else {
    # Still very WIP
    counterfactualSurvCurves <- survival::survfit(preparedModel$outcomeModel, newdata = targetDataCounterfactualFrame)
    treatmentSurvCurves <- aggregate(counterfactualSurvCurve, targetDataCounterfactualFrame[[preparedModel$treatment]])
  }
  
  transportGCResult <- list(msm = msm)
  
  class(transportGCResult) <- "transportGC"
  
  return(transportGCResult)
  
}