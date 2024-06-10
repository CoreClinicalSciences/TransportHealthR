#' @title Prepare an outcome model object for \code{transportGC}
#'
#' @description
#' An outcome model needs to be fitted using the study data to perform transportability analysis using g-computation. However, \code{glm} and functions in \code{survival} typically contain the data used to fit the model in their respective result objects. This function removes most components containing the data from these result objects to comply with data sharing regulations. The party with sole access to the study data may use only this function and provide the results object (possibly in a .rds file) to others who request it.
#' 
#' @param outcomeModel Either a formula or a \code{glm}, \code{survreg} or \code{coxph} object representing the outcome model.
#' @param response String indicating name of response variable. If \code{NULL}, it will be auto-detected from \code{outcomeModel}.
#' @param treatment String indicating name of treatment variable. This argument is required.
#' @param treatmentLevels Vector of strings indicating levels of treatment variable in the study data. If \code{NULL}, it will be auto-detected using \code{treatment} and \code{studyData}.
#' @param family Either a family function as used for \code{glm}, or one of \code{c("coxph", "survreg")}. Only required if \code{outcomeModel} is a formula.
#' @param studyData Data frame of the study data.
#'
#' @return
#' 
#' A \code{transportGCPreparedModel} object containing the following components:
#' * \code{outcomeModel}: The fitted outcome model with all components containing the data used to fit the model removed
#' * \code{response}: String indicating name of response variable
#' * \code{treatment}: String indicating name of treatment variable
#' * \code{treatmentLevels}: Vector of strings indicating levels of treatment variable
#' 
#' @export
#' 
#' @md
transportGCPreparedModel <- function(outcomeModel,
                             response,
                             treatment,
                             treatmentLevels = NULL,
                             family = stats::gaussian,
                             studyData = NULL) {
  
  # Extract response variable
  if (is.null(response)) {
    if (inherits(outcomeModel, "formula")) response <- all.vars(outcomeModel)[1]
    else response <- all.vars(outcomeModel$formula)[1]
  }
  
  # Extract treatment variable information. There is no way of detecting treatment from the outcome model formula because it often has covariates in it as well.
  if (is.null(treatment)) stop("Treatment variable name is not specified.")

  treatmentLevels <- levels(studyData[[treatment]])
  
  # If not already a model object, fit the outcome model ourselves. Recall that only the party that has access to the study data should run this function and provide its output to the other party.
  if (inherits(outcomeModel, "formula")) {
    if (is.character(family)) {
      if (family == "coxph") {
        outcomeModel <- survival::coxph(outcomeModel, data = studyData)
      } else if (family == "survreg") {
        outcomeModel <- survival::survreg(outcomeModel, data = studyData)
      }
    } else {
        outcomeModel <- stats::glm(outcomeModel, family = family, data = studyData)
    }
  }
  
  # Erase study data from outcome model object
  if (inherits(outcomeModel, "glm")) {
    outcomeModel$residuals <- outcomeModel$fitted.values <- outcomeModel$linear.predictors <-
      outcomeModel$weights <- outcomeModel$y <- outcomeModel$x <- outcomeModel$model <-
      outcomeModel$data <- outcomeModel$offset <- outcomeModel$xlevels <- NULL
  } else if (inherits(outcomeModel, "survreg")) {
    outcomeModel$linear.predictors <- outcomeModel$means <- outcomeModel$y <- outcomeModel$x <-
      outcomeModel$model <- NULL
  } else if (inherits(outcomeModel, "coxph")) {
    outcomeModel$linear.predictors <- outcomeModel$residuals <- outcomeModel$n <- outcomeModel$nevent <-
      outcomeModel$concordance <- outcomeModel$means <- outcomeModel$y <- outcomeModel$x <-
      outcomeModel$model <- NULL
  }
  
  preparedModel <- list(outcomeModel= outcomeModel,
                        response = response,
                        treatment = treatment,
                        treatmentLevels = treatmentLevels)
  
  class(preparedModel) <- "transportGCPreparedModel"
  
  return(preparedModel)
}

#' @title Check validity of prepared model object for g-computation
#' 
#' @description
#' A simple helper function that validates whether the components of the given \code{transportGCPreparedModel} object are of the correct types.
#' 
#' @param preparedModel Output object from \code{transportGCPreparedModel} function
#'
#' @return
#' A boolean indicating whether all components of \code{transportGCPreparedModel} object have the correct types.
#' 
#' @export
is.transportGCPreparedModel <- function (preparedModel) {
  return((inherits(preparedModel$outcomeModel, "glm") | inherits(preparedModel$outcomeModel, "survreg") | inherits(preparedModel$outcomeModel, "coxph")) &
           is.character(response) & is.character(treatment) & is.character(treatmentLevels))
}