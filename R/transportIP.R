#' @title Transportability analysis using IOPW
#' 
#' @description
#' Estimates the coefficients of a marginal structural model (MSM) using IP weighting in a generalizability or transportability analysis. In particular, the estimators should be unbiased for the coefficients in the superpopulation or the target population, respectively.
#' 
#'
#' @param msmFormula A formula for the MSM to be fitted, which usually includes the outcome, the treatment and any effect modifiers.
#' @param propensityScoreModel Either a formula or a \code{glm} object representing the model for treatment assignment given covariates.
#' @param participationModel Either a formula or a \code{glm} object representing the model for study participation given effect modifiers. If a formula and the "data" argument is a list, then add participation as the left-hand side of the formula.
#' @param propensityWeights Vector of custom weights balancing covariates between treatments. Providing them will override the formula or model provided by \code{propensityScoreModel}. This vector should have as any entries as the sample size of the study data.
#' @param participationWeights Vector of custom weights balancing effect modifiers between study and target populations. Providing them will override the formula or model provided by \code{participationModel}. This vector should have as any entries as the sample size of the study data.
#' @param treatment String indicating name of treatment variable. If \code{NULL}, it will be auto-detected from \code{propensityScoreModel} if provided; otherwise it will remain \code{NULL}. Note that when using custom weights, \code{treatment} should be provided so that \code{summary.transportIP} and \code{plot.transportIP} works.
#' @param participation String indicating name of participation variable. If \code{NULL}, it will be auto-detected from \code{participationModel} if provided; otherwise it will remain \code{NULL}. Note that when using custom weights, \code{participation} should be provided so that \code{summary.transportIP} and \code{plot.transportIP} works.
#' @param response String indicating name of response variable. If \code{NULL}, it will be auto-detected form \code{msmFormula}.
#' @param family Either a \code{family} function as used for \code{glm}, or one of \code{c("coxph", "survreg")}.
#' @param data Either a single data frame containing merged study and target datasets, or a list containing the study dataset and the target dataset. Note that if participationModel is a glm object, the datasets would have been merged, so provide the merged dataset containing response, treatment, covariates controlled for in the original study, study participation and effect modifiers if this is the case. Make sure to code treatment and participation as 0-1 or TRUE-FALSE, with 1 and TRUE representing treatment group and study data, respectively.
#' @param transport A boolean indicating whether a generalizability analysis (false) or a transportability analysis (true) is done.
#'
#' @details
#' The function fits models of treatment assignment and study participation in order to calculate the weights used to fit the MSM. For each of these models, if a formula is provided, logistic regression is used by default. If a \code{glm} object is provided, the function extracts the necessary weights from the object. The function does not support other weighting methods, so if they are required, provide custom weights.
#' 
#' The MSM-fitting functions do not provide correct standard errors as-is. \code{sandwich::vcovBS} is used to calculate robust bootstrap variance estimators of the parameter estimators. The function replaces the variance component in \code{summary.glm}, \code{coxph} and \code{survreg} with the robust variance estimators directly. This does not seem to behave well with \code{predict.glm} yet, but prediction is not of primary interest in a generalizability or transportability analysis.
#'
#' @return
#' A \code{transportIP} object containing the following components:
#' * \code{msm}: Raw model fit object for MSM of class \code{glm}, \code{survreg} and \code{coxph}, with the correct variance estimators appropriately replaced
#' * \code{propensityScoreModel}: Model of treatment assignment, \code{NULL} if not provided and custom propensity weights are used
#' * \code{participationModel}: Model of study participation, \code{NULL} if not provided and custom propensity weights are used
#' * \code{propensityWeights}: Propensity weights used
#' * \code{participationWeights}: Participation weights used
#' * \code{finalWeights}: Weights used to fit MSM
#' * \code{customPropensity}: Boolean indicating whether custom propensity weights are used
#' * \code{customParticipation}: Boolean indicating whether custom participation weights are used
#' * \code{treatment}: String indicating variable name of treatment
#' * \code{participation}: String indicating variable name of participation
#' * \code{data}: Data provided in \code{data} argument. Either a list containing study data and target data or a data frame containing both.
#' 
#' @export
#'
#' @md
transportIP <- function(msmFormula,
                        propensityScoreModel = NULL,
                        participationModel = NULL,
                        propensityWeights = NULL,
                        participationWeights = NULL,
                        treatment = NULL,
                        participation = NULL,
                        response = NULL,
                        family, data, transport = T) {
  
  # Auto-detect response, treatment and participation if provided
  if (is.null(response)) response <- all.vars(msmFormula)[1]

  if (is.null(treatment) & !is.null(propensityScoreModel)) {
    if (is.glm(propensityScoreModel)) treatment <- all.vars(propensityScoreModel$formula)[1]
    else treatment <- all.vars(propensityScoreModel)[1]
  }

  if (is.null(participation) & !is.null(participationModel)) {
    if (is.glm(participationModel)) participation <- all.vars(participationModel$formula)[1]
    else participation <- all.vars(participationModel)[1]
  }
  
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
    if (!is.null(treatment)) treatmentIndex <- which(names(data) == treatment)
    if (!is.null(participation)) participationIndex <- which(names(data) == participation)
  }
  

  if (inherits(propensityScoreModel, "formula")) {
    # Fit model ourselves then reassign it to propensityScoreModel
    if (!is.data.frame(data)) propensityScoreModel <- glm(propensityScoreModel, data = studyData, family = binomial())
    else propensityScoreModel <- glm(propensityScoreModel, data = data[data[,participationIndex] == 1 | data[,participationIndex] == T,], family = binomial())
  }
  
  if (inherits(participationModel, "formula")) {
    # Fit model ourselves then reassign it to participationModel
    if (!is.data.frame(data)) {
      # Handle case when study and target data are not yet merged
      effectModifiers <- all.vars(participationModel)[-1]
      studyParticipationData <- studyData[, names(studyData) %in% c(effectModifiers, "participation")]
      targetParticipationData <- targetData[, names(targetData) %in% c(effectModifiers, "participation")]
      if (!(participation %in% names(studyParticipationData))) {
        studyParticipationData$participation <- 1
        participation <- "participation"
      }
      if (!(participation %in% names(targetParticipationData))) {
        targetParticipationData$participation <- 0
        participation <- "participation"
      }
      participationData <- rbind(studyParticipationData, targetParticipationData)
      participationIndex <- which(names(participationData) == participation)
      participationModel <- glm(participationModel, data = participationData, family = binomial())
    }
    else participationModel <- glm(participationModel, data = data, family = binomial())
  }
  
  stopifnot(is.glm(propensityScoreModel) | is.null(propensityScoreModel), is.glm(participationModel) | is.null(participationModel))
  
  # Track custom weights
  
  customPropensity <- F
  if (!is.null(propensityWeights)) {
    warning("Custom propensity weights are being used. Please ensure that these weights are meaningful.")
    customPropensity <- T
  }
  
  customParticipation <- F
  if (!is.null(participationWeights)){
    warning("Custom participation weights are being used. Please ensure that these weights are meaningful.")
    customParticipation <- T
  }
  
  if (!is.null(propensityScoreModel) & !is.null(propensityWeights)) warning("Both propensity model and custom weights are provided, using custom weights.")
  
  if (!is.null(participationModel) & !is.null(participationWeights)) warning("Both participation model and custom weights are provided, using custom weights.")
  
  # Calculate actual weights to be used
  
  if (is.null(propensityWeights)) propensityWeights <- obtainWeights(propensityScoreModel, type = "probability")
  if (is.null(participationWeights)) participationWeights <- obtainWeights(participationModel, type = ifelse(transport, "odds", "probability"))
  
  # Makeshift solution to account for generalizability analysis
  
  if (!transport & length(participationWeights) > length(propensityWeights))
    participationWeights <- participationWeights[participationModel$data[, participationIndex] == 1 |
                                                   participationModel$data[, participationIndex] == 1]
  
  finalWeights <-  propensityWeights * participationWeights
              
  if (!customPropensity) if (!(treatment %in% all.vars(msmFormula)[-1])) stop("Treatment is not included in MSM.")
  
  # Fit MSM
  # Variance of coeffs are corrected in this method for coxph and survreg
  # Variance of coeffs are corrected in summary.transportIP for glm
  
  # Extract study data because data frames don't behave well with ifelse
  if (!is.data.frame(data)) toAnalyze <- studyData
  else toAnalyze <- data[data[,participationIndex] == 1 | data[,participationIndex] == T,]
  
  # The model fitting functions require weights to be part of the data frame
  toAnalyze$finalWeights <- finalWeights
  
  if (is.character(family)) {
    if (family == "coxph") {
    model <- survival::coxph(msmFormula, data = toAnalyze, weight = finalWeights)
    model$var <- sandwich::vcovBS(model)
    } else if (family == "survreg") {
    model <- survival::survreg(msmFormula, data = toAnalyze, weight = finalWeights)
    model$var <- sandwich::vcovBS(model)
  }} else {
    model <- glm(msmFormula, family = family, data = toAnalyze, weight = finalWeights)
  }
  
  # NOTE: this model object by itself has wrong SEs. The right ones are calculated by sandwich::vcovBS. Should we handle this here or in summary.transportIP?
  # Answer: for glm objects, it should be handled in summary.transportIP. For survreg and coxph, it should be handled in this function
  
  transportIPResult <- list(msm = model,
                            propensityScoreModel = propensityScoreModel,
                            participationModel = participationModel,
                            propensityWeights = propensityWeights,
                            participationWeights = participationWeights,
                            finalWeights = finalWeights,
                            customPropensity = customPropensity,
                            customParticipation = customParticipation,
                            treatment = treatment,
                            participation = participation,
                            data = data,
                            transport = transport)
  
  class(transportIPResult) <- "transportIP"
  
  return(transportIPResult)
}

# Helper function that extracts weights from models
obtainWeights <- function(model, type = c("probability", "odds")) {
  type <- match.arg(type, c("probability", "odds"))
  
  if (type == "probability") {
    return(ifelse(model$y == T | model$y == 1, 1 / model$fitted.values, 1 / (1 - model$fitted.values)))
  } else {
    participationProb <- predict(model, newdata = model$data[model$y == 1 | model$y == T,], type = "response")
    return((1-participationProb) / participationProb)
  }
}

# Helper function that detects glms
is.glm <- function(x) inherits(x, "glm")

#' @title Summarize results of a fitted MSM using the IOPW approach
#' 
#' @description
#' Returns summary object which contains a summary of the fitted MSM as well as pre- and post-weighting standardized mean differences (SMDs).
#' 
#' @rdname summary.transportIP
#'
#' @param transportIPResult Result from transportIP function
#' @param covariates Vector of strings indicating names of covariates in propensity model
#' @param effectModifiers Vector of strings indicating names of effect modifiers in participation model
#'
#' @return
#' The \code{summary.transportIP} function returns a \code{summary.transportIP} object containing the following components:
#' * \code{propensitySMD}: Table of unweighted and weighted SMDs of covariates between treatment groups. Only propensity weights are used.
#' * \code{participationSMD}: Table of unweighted and weighted SMDs of effect modifiers between study data and target data. Only participation weights are used.
#' * \code{msmSummary}: Summary object of model object for MSM. The correct variance estimators are set here for \code{glm}, whereas they are set in \code{transportIP} for \code{survreg} and \code{coxph}.
#'
#' @export
summary.transportIP <- function(transportIPResult, covariates = NULL, effectModifiers = NULL) {
  
  
  treatment <- transportIPResult$treatment
  participation <- transportIPResult$participation
  
  # Calculate SMDs for covariates
  
  # Extract study data and weights
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
  treatmentIndex <- which(names(studyData) == treatment)
  
  # Calculate SMDs for covariates (should optimize)
  prePropensitySMD <- sapply(covariates, function (covariate) as.double(smd::smd(x = studyData[, which(names(studyData) == covariate)],
                                                                       g = studyData[, treatmentIndex])$estimate))
  prePropensityBalance <- data.frame(variable = covariates, smd = prePropensitySMD, method = rep("Observed", length(covariates)))
  postPropensitySMD <- sapply(covariates, function (covariate) as.double(smd::smd(x = studyData[, which(names(studyData) == covariate)],
                                                                    g = studyData[, treatmentIndex],
                                                                    w = propensityWeights)$estimate))
  postPropensityBalance <- data.frame(variable = covariates, smd = postPropensitySMD, method = rep("Weighted", length(covariates)))
  
  propensityBalance <- rbind(prePropensityBalance, postPropensityBalance)
  
  # Calculate SMDs for effect modifiers
  
  if (!is.null(transportIPResult$participationModel)) {
    participationModel <- transportIPResult$participationModel
    participationFormula <- participationModel$formula
    participationData <- participationModel$data
    participationIndex <- which(names(participationData) == participation)
    if (is.null(effectModifiers)) effectModifiers <- all.vars(participationFormula)[-1]
  } else {
    if (is.data.frame(data)) {
      participationData <- data
    } else {
      # This only happens when user provides custom participation weights and separate study and target data
      studyDataEM <- studyData[, names(studyData) %in% effectModifiers]
      studyDataEM$participation <- 1
      targetDataEM <- targetData[, names(targetData) %in% effectModifiers]
      targetDataEM$participation <- 0
      participationData <- rbind(studyDataEM, targetDataEM)
      participation <- "participation"
      participationIndex <- which(names(participationData) == participation)
    }
  }
  # Only addresses transportability - need to implement generalizability later
  if (transportIPResult$transport) {
    participationWeights <- rep(1, nrow(participationData))
    participationWeights[participationData[, participationIndex] == 1 | participationData[, participationIndex] == T] <- transportIPResult$participationWeights
  } else if (!is.null(participationModel)) {
    participationProb <- predict(participationModel, newdata = participationData, type = "response")
    ifelse(participationData[, participationIndex] == 1 | participationData[, participationIndex] == T, 1 / participationProb, 1 / (1 - participationProb))
  } else {
    warning("Custom participation weights used for generalizability analysis. Defaulting to transportability analysis because I don't know what weights you use for target data.")
    participationWeights <- rep(1, nrow(participationData))
    participationWeights[participationData[, participationIndex] == 1 | participationData[, participationIndex] == T] <- transportIPResult$participationWeights
  }
  
  # Should also optimize
  preParticipationSMD <- sapply(effectModifiers, function (effectModifier) as.double(smd::smd(x = participationData[, which(names(participationData) == effectModifier)],
                                                                                 g = participationData[, participationIndex])$estimate))
  preParticipationBalance <- data.frame(variable = effectModifiers, smd = preParticipationSMD, method = rep("Observed", length(effectModifiers)))
  postParticipationSMD <- sapply(effectModifiers, function (effectModifier) as.double(smd::smd(x = participationData[, which(names(participationData) == effectModifier)],
                                                                                              g = participationData[, participationIndex],
                                                                                              w = participationWeights)$estimate))
  postParticipationBalance <- data.frame(variable = effectModifiers, smd = postParticipationSMD, method = rep("Weighted", length(effectModifiers)))
  participationBalance <- rbind(preParticipationBalance, postParticipationBalance)
  
  # If model is glm, calculate and replace correct SEs
  
  msm <- transportIPResult$msm
  
  msmSummary <- summary(msm)
  
  if (inherits(msmSummary, "summary.glm")) {
    msmSummary$cov.unscaled <- sandwich::vcovBS(msm)
    msmSummary$cov.scaled <- msmSummary$cov.unscaled / msmSummary$dispersion
  }
  
  summaryTransportIP <- list(propensitySMD = propensityBalance,
                             participationSMD = participationBalance,
                             msmSummary = msmSummary)
  
  class(summaryTransportIP) <- "summary.transportIP"
  
  return(summaryTransportIP)
}

#' @rdname summary.transportIP
#'
#' @param summaryTransportIP Result from transportIP function.
#' @param out Output stream.
#'
#' @export
#'
#' 
print.summary.transportIP <- function(summaryTransportIP, out = stdout()) {
  write("SMDs of covariates between treatments before and after weighting:", out)
  print(propensitySMD, out)
  
  write("SMDs of effect modifiers between study and target populations before and after weighting", out)
  print(participationSMD, out)
  
  write("MSM results:", out)
  summary(transportIPResult$msm)
}

#' @title Plot graphs relevant to transportability analysis using IOPW
#'
#' @description
#' Plot graphs for assessment of covariate balance and results in a IOPW analysis. This function currently supports mirrored histograms, SMD plots and model coefficient plots.
#'
#' @param transportIPResult Result from transportIP function
#' @param type One of \code{"propensityHist", "propensitySMD", "participationHist", "participationSMD", "msm"}. \code{Hist} produces mirrored histograms of estimated probability of treatment between treatment groups (for \code{propensity}), or of estimated probability of participation between study and target data (for \code{participation}). \code{SMD} produces SMD plots of covariates between treatment groups (for \code{propensity}) or effect modifiers between study and target data (for \code{participation}). \code{msm} produces plots showing confidence intervals for the model coefficients, which should have the correct standard errors.
#'
#' @return
#' A \code{ggplot} object which contains the desired plot.
#'
#' @export
#'
#'
plot.transportIP <- function(transportIPResult, type = "propensityHist") {
  summaryTransportIP <- summary(transportIPResult)
  resultPlot <- NULL
  
  # Match argument
  type <- match.arg(type, c("propensityHist", "propensitySMD", "participationHist", "participationSMD"))
  
  if (type == "propensityHist") {
    # Mirrored histogram of propensity scores
    if (!customPropensity) {
      propensityModel <- transportIPResult$propensityModel
      studyData <- propensityModel$data
      studyData$propensityScore <- propensityModel$fitted.values
      treatment <- transportIPResult$treatment
      resultPlot <- ggplot2::ggplot(studyData, aes(propensityScore)) + halfmoon::geom_mirror_histogram(aes(group = as.name(treatment)))
    } else {
      stop("Custom propensity weights were used. Please plot your previously estimated propensity scores using the halfmoon package, if desired.")
    }
  } else if (type == "propensitySMD") {
    # SMD plot of covariates
      propensitySMD <- summaryTransportIP$propensitySMD
      resultPlot <- ggplot2::ggplot(propensitySMD, aes(x = variable, y = smd, group = method, color = method)) + 
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0.1) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.key = element_blank())
    } else if (type == "participationHist") {
    # Mirrored histogram of probability of participation
      if (!customParticipation) {
        participationModel <- transportIPResult$participationModel
        allData <- participationModel$data
        allData$participationScore <- participationModel$fitted.values
        participation <- transportIPResult$participation
        resultPlot <- ggplot2::ggplot(allData, aes(participationScore)) + halfmoon::geom_mirror_histogram(aes(group = as.name(participation)))
      } else {
        stop("Custom participation weights were used. Please plot your previously estimated participation scores using the halfmoon package, if desired.")
      }
    } else if (type == "participationSMD") {
    # SMD plots of effect modifiers
      participationSMD <- summaryTransportIP$participationSMD
      resultPlot <- ggplot2::ggplot(participationSMD, aes(x = variable, y = smd, group = method, color = method)) + 
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0.1) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.key = element_blank())
    } else if (type == "msm") {
      resultPlot <- modelsummary::modelplot(transportIPResult$msm, vcov = sandwich::vcovBS)
    }
  
  return(resultPlot)
}

#' @title Check validity of IOPW result object
#' 
#' @description
#' A simple helper function that validates whether the components of the given \code{transportIP} object are of the correct types.
#' 
#' @param transportIPResult Result object from \code{transportIP} function
#'
#' @return A boolean indicating whether all components of \code{transportIP} object have the correct types.
#' 
#' @export
#'
is.transportIP <- function (transportIPResult) {
  return((inherits(transportIPResult$msm, "glm") | inherits(transportIPResult$msm, "coxph") | inherits(transportIPResult$msm, "survreg")) &
           (inherits(transportIPResult$propensityScoreModel, "glm") | is.null(transportIPResult$propensityScoreModel)) & 
           (inherits(transportIPResult$participationModel, "glm") | is.null(transportIPResult$participationModel)) &
           inherits(transportIPResult$propensityWeights, "numeric") & 
           inherits(transportIPResult$participationWeights, "numeric") &
           inherits(transportIPResult$finalWeights, "numeric") &
           inherits(transportIPResult$customPropensity, "logical") &
           inherits(transportIPResult$customParticipation, "logical") &
           (inherits(transportIPResult$treatment, "character") | is.null(transportIPResult$treatment)) &
           (inherits(transportIPResult$participation, "character") | is.null(transportIPResult$participation)) &
           (inherits(transportIPResult$data, "data.frame") | length(transportIPResult$data) == 2) &
           inherits(transportIPResult, "transportIP"))
}