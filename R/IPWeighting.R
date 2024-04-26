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
  if (is.null(response)) response <- all.terms(msmFormula)[1]
  
  if (is.null(treatment) & !is.null(propensityScoreModel)) treatment <- all.terms(ifelse(is.formula(propensityScoreModel), propensityScoreModel, propensityScoreModel$formula))[1]
  if (is.null(participation) & !is.null(participationModel)) participation <- all.terms(ifelse(is.formula(participationModel), participationModel, participationModel$formula))[1]
  
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
    if (!is.null(treatment)) treatmentIndex <- which(treatment %in% names(data))
    if (!is.null(participation)) participationIndex <- which(participation %in% names(data))
  }
  

  if (class(propensityScoreModel) == "formula") {
    # Fit model ourselves then reassign it to propensityScoreModel
    propensityScoreModel <- glm(propensityScoreModel, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                                            data[,participationIndex] == T,]), family = binomial())
  }
  
  if (class(participationModel) == "formula") {
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
    }
    participationModel <- glm(participationModel, data = ifelse(!is.data.frame(data), participationData, data), family = binomial())
  }
  
  stopifnot(class(propensityScoreModel) == "glm" | is.null(propensityScoreModel), class(participationModel) == "glm" | is.null(participationModel))
  
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
  
  warning("Both propensity model and custom weights are provided, using custom weights.",
          !is.null(propensityScoreModel) & !is.null(propensityWeights))
  
  warning("Both participation model and custom weights are provided, using custom weights.",
          !is.null(participationModel) & !is.null(participationWeights))
  
  # Calculate actual weights to be used
  
  propensityWeights <- ifelse(is.null(propensityWeights), obtainWeights(propensityScoreModel, type = "probability"), propensityWeights)
  participationWeights <- ifelse(is.null(participationWeights), obtainWeights(participationModel, type = ifelse(transport, "odds", "probability")), participationWeights)
  
  # Makeshift solution to account for generalizability analysis
  
  if (!transport & length(participationWeights) > length(propensityWeights))
    participationWeights <- participationWeights[participationModel$data[, participationIndex] == 1 |
                                                   participationModel$data[, participationIndex] == 1]
  
  finalWeights <-  propensityWeights * participationWeights
              
  if (!(treatment %in% all.terms(msmFormula)[-1])) stop("Treatment is not included in MSM.")
  
  # Fit MSM
  # Variance of coeffs are corrected in this method for coxph and survreg
  # Variance of coeffs are corrected in summary.transportIP for glm
  
  if (family == "coxph") {
    model <- survival::coxph(msmFormula, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                     data[,participationIndex] == T,]), weight = finalWeights)
    model$var <- sandwich::vcovBS(model)
  } else if (family == "survreg") {
    model <- survival::survreg(msmFormula, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                       data[,participationIndex] == T,]), weight = finalWeights)
    model$var <- sandwich::vcovBS(model)
  } else {
    model <- glm(msmFormula, family = family, data = ifelse(!is.data.frame(data), studyData, data[data[,participationIndex] == 1 | 
                                                                                                    data[,participationIndex] == T,]), weight = finalWeights)
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
                            data = data)
  
  class(transportIPResult) <- "transportIP"
  
  return(transportIPResult)
}

# Helper function that extracts weights from models
obtainWeights <- function(model, type = c("probability", "odds")) {
  type <- match.arg(type, c("probability", "odds"))
  
  if (type == "probability") {
    return(ifelse(model$y == T | model$y == 1, 1 / model$fitted.values, 1 / (1 - model$fitted.values)))
  } else {
    participationProb <- predict(model, newdata = model$data[model$y == 1 | model$y == T,])
    return((1-participationProb) / participationProb)
  }
}

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
  
  # Calculate SMDs for covariates (should optimize)
  propensityBalanceTables <- lapply(covariates, function (covariate) tidysmd::tidy_smd(studyData,
                                                                                       as.name(covariate),
                                                                                       as.name(treatment),
                                                                                       propensityWeights,
                                                                                       include_observed = T))
  propensityBalance <- data.table::rbindlist(propensityBalanceTables)
  
  # Calculate SMDs for effect modifiers
  
  if (!is.null(transportIPResult$participationModel)) {
    participationModel <- transportIPResult$participationModel
    participationFormula <- participationModel$formula
    participationData <- participationModel$data
    if (is.null(effectModifiers)) effectModifiers <- all.vars(participationFormula)[-1]
  }
  participation <- transportIPResult$participation
  participationData <- transportIPResult$data
  participationWeights <- transportIPResult$participationWeights
  
  # Should also optimize
  participationBalanceTables <- lapply(effectModifiers, function (effectModifier) tidysmd::tidy_smd(participationData,
                                                                                                    as.name(effectModifier),
                                                                                                    as.name(participation),
                                                                                                    participationWeights,
                                                                                                    include_observed = T))
  participationBalance <- data.table::rbindlist(participationBalanceTables)
  
  # If model is glm, calculate and replace correct SEs
  
  msm <- transportIPResult$msm
  
  msmSummary <- summary(msm)
  
  if (class(msmSummary) == "summary.glm") {
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