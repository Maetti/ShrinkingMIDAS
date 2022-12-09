library(midasr)
sDir <- "inst/simulation"
i <- 1
nTrain <- 200

## get simulated data
nSimRunSave <- helper_create_number_name(i)

sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

lData <- readRDS(sDirInputLoad)

y_train <- lData$model_data$y_train
y_test <- lData$model_data$y_test

x_train <- lData$x_raw[1:(6 + 3*nTrain), ]
x_test <- lData$x_raw[(6 + 3*nTrain + 1):nrow(lData$x_raw), ]

## prepare data
y_train2 <- as.matrix(c(NA, y_train, NA))
sFormula <- as.formula(paste0("y_train2 ~ ", paste0("mls(x_", 1:50, ", 0:5, 3, nealmon)", collapse = " + "), " - 1"))
lStart <- lapply(1:50, function(x) {c(1, -0.5)})
names(lStart) <- paste0("x_", 1:50)

for (i in 1:50) {
      assign(x = paste0("x_", i), x_train[, i])
      assign(x = paste0("x_test_", i), x_test[, i])
}

## estimate model
md_midas <- midasr::midas_r(formula = sFormula, start = lStart)


## start with tests
test_that("test nealmon function", {

      mX1 <- midasr::mls(x_train[, 1], 0:5, 3, midasr::nealmon)[2:201, ]
      mX2 <- lData$x_align[[1]][1:200, ]

      testthat::expect_equivalent(mX1, mX2)

})


test_that("midas design matrix and y train", {

      ## check y response
      checkY <- md_midas$model[, 1]
      checkYTrue <- lData$model_data$y_train

      all(checkY == checkYTrue)

      ## check design matrix
      checkDesign <- unname(md_midas$model[, -1])
      checkXTrue <- lData$model_data$x_train
      dim(checkDesign)
      dim(checkXTrue)
      checkDesign[1:5, 1:6]
      checkXTrue[1:5, 1:6]

      testthat::expect_equivalent(checkDesign, checkXTrue, na.rm = TRUE)


})


test_that("midas coefficient check", {

      ## check if coefficient are the same
      # - yes they are -> can use midas_coefficients on the whole test matrix
      nAlmon1 <- md_midas$coefficients[1:2]

      lCoefTest <- vector("list", length = length(md_midas$coefficients) / 2)
      for (i in 1:(length(md_midas$coefficients) / 2)) {
            lCoefTest[[i]] <- nealmon(p = md_midas$coefficients[(2*i - 1):(2 * i)], d = 6)
      }
      nCoefTest <- unlist(lCoefTest)

      testthat::expect_equivalent(nCoefTest, md_midas$midas_coefficients)

})


test_that("prediction check", {

      ## prediction method #1
      mPred1 <- lData$model_data$x_test %*% md_midas$midas_coefficients

      ## prediction method #2
      lTestData <- lapply(1:50, function(x) {lData$x_raw[601:753, x]})
      names(lTestData) <- paste0("x_", 1:50)
      mPred2 <- forecast(md_midas, newdata = lTestData)


      round(cbind(mPred1, mPred2$mean[-1], lData$model_data$y_test), 4)
      testthat::expect_equivalent(round(mPred1, 5), round(mPred2$mean[-1], 5))

})
