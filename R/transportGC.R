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
#' * \code{bootstrapNum}: Integer indicating number of bootstrap datasets simulated to calculate robust variance estimators.
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
  
  if (!preparedModel$wipe) {
    studyData <- preparedModel$outcomeModel$data
    treatmentLevels <- preparedModel$treatmentLevels
    treatmentGroupData <- list()
    treatmentVector <- as.character(studyData[[preparedModel$treatment]])
    for (level in treatmentLevels) {
      treatmentGroupData[[level]] <- studyData[treatmentVector == level, ]
    }
    nLevels <- sapply(treatmentGroupData, nrow)
    names(nLevels) <- treatmentLevels
  }
  
  nTarget <- nrow(targetData)
  
  bootstrapEstimates <- t(sapply(1:bootstrapNum,
                                 function (x) {
                                   if (preparedModel$wipe) {
                                    targetBoot <- targetData[sample.int(n = nTarget, replace = T), ]
                                   
                                    resultBoot <- transportGCFit(msmFormula,
                                                                preparedModel,
                                                                targetData = targetBoot)
                                   } else {
                                      treatmentGroupBoot <- list()
                                      for (level in treatmentLevels) {
                                        nSample <- nLevels[level]
                                        treatmentGroupBoot[[level]] <- treatmentGroupData[[level]][sample.int(nSample, replace = T), ]
                                      }
                                      for (level in treatmentLevels) {
                                        if (!exists("studyBoot")) studyBoot <- treatmentGroupBoot[[level]]
                                        else studyBoot <- rbind(studyBoot, treatmentGroupBoot[[level]])
                                      }
                                     
                                      targetBoot <- targetData[sample.int(n = nTarget, replace = T), ]
                                     
                                      preparedBoot <- transportGCPreparedModel(preparedModel$outcomeModel$formula,
                                                                               response = preparedModel$response,
                                                                               treatment = preparedModel$treatment,
                                                                               treatmentLevels = preparedModel$treatmentLevels,
                                                                               family = preparedModel$family,
                                                                               method = preparedModel$method,
                                                                               studyData = studyBoot,
                                                                               wipe = F)
                                      resultBoot <- transportGCFit(msmFormula,
                                                                  preparedBoot,
                                                                  targetData = targetBoot)
                                   }
                                   
                                   return(c(resultBoot$msm$coefficients, resultBoot$msm$zeta))
                                 }))
  
  varMatrix <- stats::var(bootstrapEstimates)
  colnames(varMatrix) <- rownames(varMatrix) <- c(names(transportGCResult$msm$coefficients), names(transportGCResult$msm$zeta))
  transportGCResult$msm$var <- varMatrix
  
  transportGCResult$bootstrapNum <- bootstrapNum
  
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
  if (!inherits(outcomeModel, "coxph")) {
    targetDataCounterfactualFrame[[response]] <- stats::predict(outcomeModel,
                                                       newdata = targetDataCounterfactualFrame,
                                                       type = "response")
    
    if (inherits(outcomeModel, "glm"))
      msm <- stats::glm(msmFormula, data = targetDataCounterfactualFrame, family = outcomeModel$family)
    else if (inherits(outcomeModel, "survreg"))
      msm <- survival::survreg(msmFormula, data = targetDataCounterfactualFrame)
    else if (inherits(outcomeModel, "polr"))
      msm <- MASS::polr(msmFormula, data = targetDataCounterfactualFrame, Hess = T, method = outcomeModel$method)
  }
  else {
    # Still very WIP
    counterfactualSurvCurves <- survival::survfit(preparedModel$outcomeModel, newdata = targetDataCounterfactualFrame)
    treatmentSurvCurves <- stats::aggregate(counterfactualSurvCurves, targetDataCounterfactualFrame[[preparedModel$treatment]])
  }
  
  transportGCResult <- list(msm = msm,
                            preparedModel = preparedModel,
                            data = targetData)
  
  class(transportGCResult) <- "transportGC"
  
  return(transportGCResult)
  
}

#' @title Summarize results of a fitted MSM object using g-computation
#' 
#' @description
#' Returns summary object which contains summary objects for the MSM and the outcome model, as well as information about response and treatment variables. In the MSM summary object, the correct variance estimators are calculated.
#' 
#' @rdname summary.transportGC
#'
#' @param object Result from \code{transportGC} function
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' The \code{summary.transportGC} function returns a \code{summary.transportGC} object containing the following components:
#' * \code{msmSummary}: Summary object of MSM with correct variance estimates. 
#' * \code{preparedModelSummary}: Summary object of outcome model, provided only for information. No conclusions should be drawn from the outcome model.
#' * \code{response}: String indicating response variable name.
#' * \code{treatment}: String indicating treatment variable name.
#' * \code{treatmentLevels}: Vector of strings indicating levels of treatment variable
#' 
#' @export
summary.transportGC <- function (object, ...) {
  transportGCResult <- object
  
  preparedModel <- transportGCResult$preparedModel
  
  preparedModelSummary <- summary(preparedModel$outcomeModel)
  response <- preparedModel$response
  treatment <- preparedModel$treatment
  treatmentLevels <- preparedModel$treatmentLevels
  
  # If model is glm, calculate and replace correct SEs
  
  msm <- transportGCResult$msm
  
  msmSummary <- summary(msm)
  
  if (inherits(msmSummary, "summary.glm")) {
    if (!is.null(msm$var)) msmSummary$cov.scaled <- msm$var
    msmSummary$cov.unscaled <- msmSummary$cov.scaled / msmSummary$dispersion
    msmSummary$coefficients[, 2] <- sqrt(diag(msmSummary$cov.scaled))
    msmSummary$coefficients[, 3] <- msmSummary$coefficients[, 1] / msmSummary$coefficients[, 2]
    if (msmSummary$family$family == "gaussian") msmSummary$coefficients[, 4] <- 2 * stats::pt(abs(msmSummary$coefficients[, 3]), msmSummary$df[2], lower.tail = F)
    else msmSummary$coefficients[, 4] <- 2 * stats::pnorm(abs(msmSummary$coefficients[, 3]), lower.tail = F)
  }
  
  # Same for polr
  
  if (inherits(msmSummary, "summary.polr")) {
    if (!is.null(msm$var)) msmSummary$coefficients[, 2] <- sqrt(diag(msm$var))
    msmSummary$coefficients[, 3] <- msmSummary$coefficients[, 1] / msmSummary$coefficients[, 2]
  }
  
  summaryTransportGC <- list(msmSummary = msmSummary,
                    preparedModelSummary = preparedModelSummary,
                    response = response,
                    treatment = treatment,
                    treatmentLevels = treatmentLevels)
  
  class(summaryTransportGC) <- "summary.transportGC"
  
  return(summaryTransportGC)
}

#' @rdname summary.transportGC
#'
#' @param x \code{summary.transportGC} object.
#' @param out Output stream.
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @export
#'
#'
print.summary.transportGC <- function (x, out = stdout(), ...) {
  summaryTransportGC <- x
  
  write(paste0("Response: ", summaryTransportGC$response), out)
  write(paste0("Treatment: ", summaryTransportGC$treatment), out)
  
  write("Fitted outcome model:", out)
  print(summaryTransportGC$preparedModelSummary, out)
  
  write("Fitted MSM:", out)
  print(summaryTransportGC$msmSummary, out)
}

#' @title Visually represent results of transportability analysis using g-computation
#' 
#' @description
#' This function is a wrapper for \code{modelsummary::modelplot} to plot the coefficient estimates in a transportability analysis using g-computation. Note that the correct variance estimates are used in this function.
#' 
#' @param x Result from \code{transportGC} function.
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' A \code{ggplot} object showing the estimates and confidence intervals of the MSM coefficients.
#' 
#' @export
plot.transportGC <- function (x, ...) {
  transportGCResult <- x
  resultPlot <- modelsummary::modelplot(transportGCResult$msm, vcov = list(transportGCResult$msm$var))
  
  return(resultPlot)
}

#' @title Check validity of g-computation result object
#' 
#' @description
#' A simple helper function that validates whether the components of the given \code{transportGC} object are of the correct types.
#' 
#' @param transportGCResult Result object from \code{transportGC} function
#'
#' @return A boolean indicating whether all components of \code{transportGC} object have the correct types.
#' 
#' @export
#'
is.transportGC <- function (transportGCResult) {
  return((inherits(transportGCResult$msm, "glm") | inherits(transportGCResult$msm, "survreg") | inherits(transportGCResult$msm, "coxph") | inherits(transportGCResult$msm, "polr")) &
           is.transportGCPreparedModel(transportGCResult$preparedModel) &
           is.data.frame(transportGCResult$data))
}