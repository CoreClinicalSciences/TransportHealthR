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
  
  expect_true(is.transportIP(testResult))
  expect_true(!testResult$customPropensity & !testResult$customParticipation)
  expect_true(is.data.frame(testResult$data))
  
  compareResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                               participationModel = participation ~ stress + med2,
                               family = gaussian,
                               data = data,
                               transport = T)
  
  expect_equal(testResult$msm$coefficients, compareResult$msm$coefficients)
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
  
  expect_true(is.transportIP(testResult))
  expect_true(!testResult$customPropensity & !testResult$customParticipation)
  expect_true(is.data.frame(testResult$data))
  
  compareResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                               participationModel = participation ~ stress + med2,
                               family = gaussian,
                               data = data,
                               transport = T)
  
  expect_equal(testResult$msm$coefficients, compareResult$msm$coefficients)
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
                                              participationModel = participation ~ med2 + stress,
                                              family = gaussian,
                                              data = data,
                                              transport = T))
  
  expect_true(is.transportIP(testResult))
  expect_true(!testResult$customPropensity & !testResult$customParticipation)
  expect_false(is.data.frame(testResult$data))
  
  compareResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                               participationModel = participation ~ stress + med2,
                               family = gaussian,
                               data = data,
                               transport = T)
  
  expect_equal(testResult$msm$coefficients, compareResult$msm$coefficients)
})

test_that("Scenario 5: custom weights", {
  set.seed(20240429)
  data <- generateTestData()
  
  expect_warning(testResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                              propensityWeights = rep(1, nrow(data$studyData)),
                                              participationWeights = rep(1, nrow(data$studyData)),
                                              family = gaussian,
                                              data = data,
                                              transport = T))
  
  expect_true(is.transportIP(testResult))
  expect_true(testResult$customPropensity & testResult$customParticipation)
  expect_false(is.data.frame(testResult$data))
  expect_true(all(testResult$propensityWeights == 1))
  expect_true(all(testResult$participationWeights == 1))
})

test_that("Testing summary.transportIP", {
  set.seed(20240429)
  data <- generateTestData()
  
  testResult <- transportIP(msmFormula = sysBloodPressure ~ med1,
                                              propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                                              participationModel = participation ~ stress + med2,
                                              family = gaussian,
                                              data = data,
                                              transport = T)
  browser()
  testSummary <- summary(testResult)
  
  expect_true(inherits(testSummary, "summary.transportIP"))
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$participationSMD, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
})
