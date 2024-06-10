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
                                   
                                   return(resultBoot$msm$coefficients)
                                 }))
  
  varMatrix <- var(bootstrapEstimates)
  colnames(varMatrix) <- rownames(varMatrix) <- names(transportGCResult$msm$coefficients)
  transportGCResult$msm$var <- varMatrix
  
  return(transportGCResult)
}

transportGCFit <- function (msmFormula,
                         preparedModel,
                         targetData) {
  if (!is.transportGCPreparedModel(preparedModel)) stop("Please provide a transportGCPreparedModel object.")
  if (preparedModel$response != all.vars(msmFormula)[1]) stop("Response variable does not match transportGCPreparedModel object.")
  
  targetDataReferenceList <- list()
  for (level in preparedModel$treatmentLevels) {
    targetDataReferenceList[[level]] <- targetData
    targetDataReferenceList[[level]][[preparedModel$treatment]] <- level
  }
  
  targetDataCounterfactualFrame <- Reduce(function(d1, d2) rbind(d1, d2), targetDataReferenceList)
  targetDataCounterfactualFrame[[treatment]] <- as.factor(targetDataCounterfactualFrame[[treatment]])
  
  # TODO: adapt to coxph since it doesn't give survival time predictions - see Lee (2024) (transporting survival of hiv...) and Chen & Tsiatis (2001)
  if (!inherits(preparedModel$outcomeModel, "coxph")) {
    targetDataCounterfactualFrame[[response]] <- stats::predict(preparedModel$outcomeModel,
                                                       newdata = targetDataCounterfactualFrame,
                                                       type = "response")
    
    if (inherits(preparedModel$outcomeModel, "glm"))
      msm <- glm(msmFormula, family = preparedModel$outcomeModel$family, data = targetDataCounterfactualFrame)
    else msm <- survival::survreg(msmFormula, data = targetDataCounterfactualFrame)
  }
  else {
    counterfactualSurvCurves <- list()
    for (level in preparedModel$treatmentLevels) {
      # Very WIP
      treatmentLevelSurvCurves <- survival::survfit(preparedModel$outcomeModel,
                                                    newdata = targetDataReferenceList[[level]])
      
    }
  }
  
  transportGCResult <- list(msm = msm)
  
  class(transportGCResult) <- "transportGC"
  
  return(transportGCResult)
  
}