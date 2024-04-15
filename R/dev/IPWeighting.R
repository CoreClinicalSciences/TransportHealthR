#' Obtain transported causal effect estimates using IOPW
#'
#' @param msmFormula 
#' @param propScoreMod 
#' @param participationMod 
#'
#' @return
#' @export
#'
#' @examples
transportIP <- function(msmModel,
                        propensityScoreModel = NULL,
                        participationModel = NULL,
                        propensityWeights = NULL,
                        family, data, transport) {
  
  if (class(propensityScoreModel) == "formula") {
    # Fit model ourselves then reassign it to propensityScoreModel
    propensityScoreModel <- glm(propensityScoreModel, data = data, family = binomial())
  }
  
  if (class(participationModel) == "formula") {
    # Fit model ourselves then reassign it to participationModel
    participationModel <- glm(participationModel, data = data, family = binomial())
  }
  
  weights <- ifelse(is.null(propensityWeights), obtainWeights(propensityScoreModel, type = "probability"), propensityWeights) *
              obtainWeights(participationModel, type = ifelse(transport, "odds", "probability"))
  
  # Adapt to survival data as well
  model <- glm(msmFormula, family = family, data = data)
  
  transportIPResult <- list(msm = model,
                            propScoreModel = propScoreModel,
                            participationModel = participationModel,
                            weights = weights)
  
  return(transportIPResult)
}

obtainWeights <- function(model, type = c("probability", "odds")) {
  type <- match.arg(type, c("probability", "odds"))
  
  if (type == "probability") {
    return(1 / model$fitted.values)
  } else {
    return((1 - model$fitted.values) / model$fitted.values)
  }
}