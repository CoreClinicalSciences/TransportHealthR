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
transportIP <- function(msmFormula,
                        propensityScoreModel,
                        participationModel,
                        family, data, transport) {
  
  weights <- obtainWeights(propScoreModel, type = "probability") *
              obtainWeights(participationModel, type = ifelse(transport, "odds", "probability"))
  
  model <- glm(msmFormula, family = family, data = data)
  
  transportIPResult <- list(msm = model,
                            propScoreModel = propScoreModel,
                            participationModel = participationModel)
  
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