test_that("Scenario 1: separate study and target data, formula provided for propensityScoreModel and participationModel", {
  set.seed(20240429)
  data <- generateTestData()
  expect_no_warning(testResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                                               participationModel = participation ~ stress + med2,
                                               family = gaussian,
                                               data = data,
                                               transport = T))
  
  expect_true(is.transportIP(testResult))
  expect_true(!testResult$customPropensity & !testResult$customParticipation)
  expect_false(is.data.frame(testResult$data))
  
  expect_no_warning(testSummary <- summary(testResult))
  
  expect_true(inherits(testSummary, "summary.transportIP"))
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$participationSMD, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensityHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensitySMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationSMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "msm"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  # Test if function works properly when study data and target data are provided in opposite order
  
  swapData <- list()
  
  swapData[[1]] <- data$targetData
  swapData[[2]] <- data$studyData
  
  expect_no_warning(testSwapResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                              propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                                              participationModel = participation ~ stress + med2,
                                              family = gaussian,
                                              data = swapData,
                                              transport = T))
  
  expect_true(is.transportIP(testSwapResult))
  expect_true(!testSwapResult$customPropensity & !testSwapResult$customParticipation)
  expect_false(is.data.frame(testSwapResult$data))
  
  expect_no_warning(testSwapSummary <- summary(testSwapResult))
  
  expect_equal(testResult$msm$coefficients, testSwapResult$msm$coefficients)
})

test_that("Scenario 2: merged study and target data, formula provided for propensityScoreModel and participationModel", {
  set.seed(20240429)
  data <- generateTestData()
  
  studyData <- data$studyData
  targetData <- data$targetData
  
  studyData$participation <- 1
  
  targetData$sysBloodPressure <- targetData$med1 <- NA
  
  targetData$participation <- 0
  
  mergedData <- rbind(studyData, targetData)
  
  expect_no_warning(testResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                                               participationModel = participation ~ stress + med2,
                                               family = gaussian,
                                               data = mergedData,
                                               transport = T))
  
  expect_no_warning(testSummary <- summary(testResult))
  
  expect_true(is.transportIP(testResult))
  expect_true(!testResult$customPropensity & !testResult$customParticipation)
  expect_true(is.data.frame(testResult$data))
  
  expect_true(inherits(testSummary, "summary.transportIP"))
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$participationSMD, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  compareResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                               participationModel = participation ~ stress + med2,
                               family = gaussian,
                               data = data,
                               transport = T)
  
  compareSummary <- summary(compareResult)
  
  expect_equal(testResult$msm$coefficients, compareResult$msm$coefficients)
  
  expect_equal(testSummary$propensitySMD$smd, compareSummary$propensitySMD$smd)
  expect_equal(testSummary$participationSMD$smd, compareSummary$participationSMD$smd)
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensityHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensitySMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationSMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "msm"))
  expect_true(ggplot2::is.ggplot(testPlot))
})

test_that("Scenario 3: glm provided for propensityScoreModel and participationModel (so merged data is provided)", {
  set.seed(20240429)
  data <- generateTestData()
  
  studyData <- data$studyData
  targetData <- data$targetData
  
  studyData$participation <- 1
  
  targetData$sysBloodPressure <- targetData$med1 <- NA
  
  targetData$participation <- 0
  
  mergedData <- rbind(studyData, targetData)
  
  propensityScoreModel <- glm(med1 ~ sex + percentBodyFat + stress, data = studyData, family = binomial())
  
  participationModel <- glm(participation ~ stress + med2, data = mergedData, family = binomial())
  
  expect_no_warning(testResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                              propensityScoreModel = propensityScoreModel,
                                              participationModel = participationModel,
                                              family = gaussian,
                                              data = mergedData,
                                              transport = T))
  
  expect_no_warning(testSummary <- summary(testResult))
  
  expect_true(is.transportIP(testResult))
  expect_true(!testResult$customPropensity & !testResult$customParticipation)
  expect_true(is.data.frame(testResult$data))
  
  expect_true(inherits(testSummary, "summary.transportIP"))
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$participationSMD, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  compareResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                               participationModel = participation ~ stress + med2,
                               family = gaussian,
                               data = data,
                               transport = T)
  
  compareSummary <- summary(compareResult)
  
  expect_equal(testResult$msm$coefficients, compareResult$msm$coefficients)
  
  expect_equal(testSummary$propensitySMD$smd, compareSummary$propensitySMD$smd)
  expect_equal(testSummary$participationSMD$smd, compareSummary$participationSMD$smd)
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensityHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensitySMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationSMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "msm"))
  expect_true(ggplot2::is.ggplot(testPlot))
})

test_that("Scenario 4: separate study and target data, glm provided for propensityScoreModel, but not participationModel", {
  set.seed(20240429)
  data <- generateTestData()
  
  studyData <- data$studyData
  targetData <- data$targetData
  
  studyData$participation <- 1
  
  targetData$sysBloodPressure <- targetData$med1 <- NA
  
  targetData$participation <- 0
  
  mergedData <- rbind(studyData, targetData)
  
  propensityScoreModel <- glm(med1 ~ sex + percentBodyFat + stress, data = studyData, family = binomial())
  
  expect_no_warning(testResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                              propensityScoreModel = propensityScoreModel,
                                              participationModel = participation ~ stress + med2,
                                              family = gaussian,
                                              data = data,
                                              transport = T))
  
  expect_no_warning(testSummary <- summary(testResult))
  
  expect_true(is.transportIP(testResult))
  expect_true(!testResult$customPropensity & !testResult$customParticipation)
  expect_false(is.data.frame(testResult$data))
  
  expect_true(inherits(testSummary, "summary.transportIP"))
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$participationSMD, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  compareResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                               participationModel = participation ~ stress + med2,
                               family = gaussian,
                               data = data,
                               transport = T)
  
  compareSummary <- summary(compareResult)
  
  expect_equal(testResult$msm$coefficients, compareResult$msm$coefficients)
  
  expect_equal(testSummary$propensitySMD$smd, compareSummary$propensitySMD$smd)
  expect_equal(testSummary$participationSMD$smd, compareSummary$participationSMD$smd)
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensityHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensitySMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationSMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "msm"))
  expect_true(ggplot2::is.ggplot(testPlot))
})

test_that("Scenario 5: custom weights", {
  set.seed(20240429)
  data <- generateTestData()
  
  expect_warning(testResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                              propensityWeights = rep(1, nrow(data$studyData)),
                                              participationWeights = rep(1, nrow(data$studyData)),
                                              treatment = "med1",
                                              family = gaussian,
                                              data = data,
                                              transport = T)) |> expect_warning() |> expect_warning()
  #browser()
  expect_no_warning(testSummary <- summary(testResult, covariates = c("sex", "percentBodyFat", "stress"), effectModifiers = c("stress", "med2")))
  
  expect_true(is.transportIP(testResult))
  expect_true(testResult$customPropensity & testResult$customParticipation)
  expect_false(is.data.frame(testResult$data))
  expect_true(all(testResult$propensityWeights == 1))
  expect_true(all(testResult$participationWeights == 1))
  
  expect_true(inherits(testSummary, "summary.transportIP"))
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$participationSMD, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  expect_error(testPlot <- plot(testResult, type = "propensityHist", covariates = c("sex", "percentBodyFat", "stress"), effectModifiers = c("stress", "med2")))
  
  expect_no_warning(testPlot <- plot(testResult, type = "propensitySMD", covariates = c("sex", "percentBodyFat", "stress"), effectModifiers = c("stress", "med2")))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_error(testPlot <- plot(testResult, type = "participationHist", covariates = c("sex", "percentBodyFat", "stress"), effectModifiers = c("stress", "med2")))
  
  expect_no_warning(testPlot <- plot(testResult, type = "participationSMD", covariates = c("sex", "percentBodyFat", "stress"), effectModifiers = c("stress", "med2")))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_warning(testPlot <- plot(testResult, type = "msm", covariates = c("sex", "percentBodyFat", "stress"), effectModifiers = c("stress", "med2")))
  expect_true(ggplot2::is.ggplot(testPlot))
})
