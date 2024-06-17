test_that("Outcome model provided as formula", {
  testData <- generateTestData()
  expect_no_error(preparedModel <- transportGCPreparedModel(outcomeModel = sysBloodPressure ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                                            treatment = "med1",
                                            family = stats::gaussian,
                                            studyData = testData$studyData))
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC(sysBloodPressure ~ med1,
                                   preparedModel,
                                   testData$targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
  
  expect_true(inherits(transportGCResult$msm, "glm"))
  
  expect_no_error(preparedModel <- transportGCPreparedModel(outcomeModel = sysBloodPressure ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                                                            treatment = "med1",
                                                            family = stats::gaussian,
                                                            studyData = testData$studyData,
                                                            wipe = F))
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC(sysBloodPressure ~ med1,
                                                   preparedModel,
                                                   testData$targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
  
  expect_true(inherits(transportGCResult$msm, "glm"))
})
