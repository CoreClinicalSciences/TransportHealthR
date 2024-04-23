#' Obtain transported causal effect estimates using IOPW
#'
#' @param msmFormula 
#' @param propScoreMod 
#' @param participationMod Either a formula or a glm object. If a formula and the "data" argument is a list, then add participation as the left-hand side of the formula.
#' @param data Either a single data frame containing merged study and target datasets, or a list containing the study dataset and the target dataset. Note that if participationModel is a glm object, the datasets would have been merged, so provide the merged dataset containing response, treatment, covariates controlled for in the original study, study participation and effect modifiers if this is the case. Make sure to code treatment and participation as 0-1 or TRUE-FALSE, with 1 and TRUE representing treatment group and study data, respectively
#'
#'
#'
#' @return
#' @export
#'
#' @examples
transportIP <- function(msmModel,
                        propensityScoreModel = NULL,
                        participationModel = NULL,
                        propensityWeights = NULL,
                        participationWeights = NULL,
                        treatment = NULL,
                        participation = NULL,
                        response = NULL,
                        family, data, transport) {
  
  if (is.null(treatment)) treatment <- all.terms(ifelse(is.formula(propensityScoreModel), propensityScoreModel, propensityScoreModel$formula))[1]
  if (is.null(participation)) participation <- all.terms(ifelse(is.formula(participationModel), participationModel, participationModel$formula))[1]
  if (is.null(response)) response <- all.terms(msmFormula)[1]
  
  studyData <- NULL
  targetData <- NULL
  treatmentIndex <- NULL
  participationIndex <- NULL
  
  if (!is.data.frame(data)) {
    if (response %in% names(data[[1]])) {
      studyData <- data[[1]]
      targetData <- data[[2]]
    }
    else{
      studyData <- data[[2]]
      targetData <- data[[1]]
    }
  } else {
    treatmentIndex <- which(treatment %in% names(data))
    participationIndex <- which(participation %in% names(data))
  }
  

  if (class(propensityScoreModel) == "formula") {
    # Fit model ourselves then reassign it to propensityScoreModel
    propensityScoreModel <- glm(propensityScoreModel, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                                            data[,participationIndex] == T,]), family = binomial())
  }
  
  if (class(participationModel) == "formula") {
    # Fit model ourselves then reassign it to participationModel
    if (!is.data.frame(data)) {
      effectModifiers <- all.vars(participationModel)[-1]
      studyParticipationData <- studyData[, names(studyData) %in% c(effectModifiers, "participation")]
      targetParticipationData <- targetData[, names(studyData) %in% c(effectModifiers, "participation")]
      if (!(participation %in% names(studyParticipationData))) studyParticipationData$participation <- 1
      if (!(participation %in% names(targetParticipationData))) targetParticipationData$participation <- 0
      participationData <- rbind(studyParticipationData, targetParticipationData)
    }
    participationModel <- glm(participationModel, data = ifelse(!is.data.frame(data), participationData, data), family = binomial())
  }
  
  stopifnot(class(propensityScoreModel) == "glm", class(participationModel) == "glm")
  
  warning("Custom propensity weights are being used. Please ensure that these weights are meaningful.",
          !is.null(propensityWeights))
  
  warning("Custom participation weights are being used. Please ensure that these weights are meaningful.",
          !is.null(participationWeights))
  
  warning("Both propensity model and custom weights are provided, using custom weights.",
          !is.null(propensityScoreModel) & !is.null(propensityWeights))
  
  warning("Both participation model and custom weights are provided, using custom weights.",
          !is.null(participationModel) & !is.null(participationWeights))
  
  propensityWeights <- ifelse(is.null(propensityWeights), obtainWeights(propensityScoreModel, type = "probability"), propensityWeights)
  participationWeights <- ifelse(is.null(participationWeights), obtainWeights(participationModel, type = ifelse(transport, "odds", "probability")), participationWeights)
  
  finalWeights <-  propensityWeights * participationWeights
              
  if (!(treatment %in% all.terms(msmFormula)[-1])) stop("Treatment is not included in MSM.")
  
  if (family == "coxph") {
    model <- survival::coxph(msmFormula, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                     data[,participationIndex] == T,]), weight = finalWeights)
  } else if (family == "survreg") {
    model <- survival::survreg(msmFormula, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                       data[,participationIndex] == T,]), weight = finalWeights)
  } else {
    model <- glm(msmFormula, family = family, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                                    data[,participationIndex] == T,]), weight = finalWeights)
  }
  
  transportIPResult <- list(msm = model,
                            propensityScoreModel = propensityScoreModel,
                            participationModel = participationModel,
                            propensityWeights = propensityWeights,
                            participationWeights = participationWeights,
                            finalWeights = finalWeights,
                            treatment = treatment,
                            participation = participation,
                            data = data)
  
  class(transportIPResult) <- "transportIP"
  
  return(transportIPResult)
}

obtainWeights <- function(model, type = c("probability", "odds")) {
  type <- match.arg(type, c("probability", "odds"))
  
  if (type == "probability") {
    return(ifelse(model$y == T | model$y == 1, 1 / model$fitted.values, 1 / (1 - model$fitted.values)))
  } else {
    participationProb <- predict(model, newdata = model[model$y == 1 | model$y == T])
    return((1-participationProb) / participationProb)
  }
}

#' Summarize results of a fitted MSM using the IOPW approach
#'
#' @param transportIPResult Result from transportIP function
#'
#' @return
#' @export
#'
#' @examples
summary.transportIP <- function(transportIPResult, covariates = NULL, effectModifiers = NULL, out = stdout()) {
  treatment <- transportIPResult$treatment
  participation <- transportIPResult$participation
  
  # Todo: print covariate balance for propensity model
  if (!is.null(transportIPResult$propensityScoreModel)) {
    propensityScoreModel <- transportIPResult$propensityScoreModel
    propensityFormula <- propensityScoreModel$formula
    studyData <- propensityScoreModel$data
    if (is.null(covariates)) covariates <- all.vars(propensityFormula)[-1]
  } else {
    data <- transportIPResult$data
    if (is.data.frame(data)) {
      participationIndex <- which(names(data) == participation)
      studyData <- data[data[, participationIndex] == 1 |
                        data[, participationIndex] == T, ]
    } else {
      if (transportIPResult$response %in% names(data[[1]])) {
        studyData <- data[[1]]
        targetData <- data[[2]]
      }
      else{
        studyData <- data[[2]]
        targetData <- data[[1]]
      }
    }
  }
  propensityWeights <- transportIPResult$propensityWeights
  
  propensityBalanceTables <- lapply(covariates, function (covariate) tidysmd::tidy_smd(studyData,
                                                                                       as.name(covariate),
                                                                                       as.name(treatment),
                                                                                       propensityWeights,
                                                                                       include_observed = T))
  propensityBalance <- data.table::rbindlist(propensityBalanceTables)
  
  write("SMDs of covariates between treatments before and after weighting:", out)
  print(propensityBalance, out)
  
  # Todo: print covariate balance for participation model
  
  if (!is.null(transportIPResult$participationModel)) {
    participationModel <- transportIPResult$participationModel
    participationFormula <- participationModel$formula
    participationData <- participationModel$data
    if (is.null(effectModifiers)) effectModifiers <- all.vars(participationFormula)[-1]
  }
  participation <- transportIPResult$participation
  participationData <- transportIPResult$data
  participationWeights <- transportIPResult$participationWeights
  
  participationBalanceTables <- lapply(effectModifiers, function (effectModifier) tidysmd::tidy_smd(participationData,
                                                                                                    as.name(effectModifier),
                                                                                                    as.name(participation),
                                                                                                    participationWeights,
                                                                                                    include_observed = T))
  participationBalance <- data.table::rbindlist(participationBalanceTables)
  
  write("SMDs of effect modifiers between study and target populations before and after weighting", out)
  print(participationBalance, out)
  
  # Todo: print results of MSM
  write("MSM results:", out)
  summary(transportIPResult$msm)
}