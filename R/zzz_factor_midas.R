zzz_factor_midas <- function() {


      #   ____________________________________________________________________________
      #   Set Inputs                                                              ####

      library(dplyr)
      library(ggplot2)
      library(midasr)

      ## set logging file
      logger::log_appender(logger::appender_tee(here::here("inst", "simulation", "logging", "model_check",
                                                           glue::glue("{gsub('-', '_', Sys.Date())}_model_check_log"))))
      logger::log_info("Setting input variables")

      ## simulations
      nSimulation <- 15
      vSeed <- 1001:(1001 + nSimulation)

      bForceNewData <- FALSE


      ## specifiy general amount
      nY <- 250
      nTrain <- 200

      ## predictor variable setup
      vBeta <- c(1.2, -0.5, 3, 1.1, 0.5, -2.5, -1.5, 0.8, 3.5, rep(0, 41))
      nK = length(vBeta) ## amount of predictors
      nFreq = 3
      nLag = 6
      nP = 3

      ## data generating process variables
      nMu = 0.1
      nRho = 0.5
      nVar = 1
      nWithin = 0.5
      nBetween = 0

      ## data generation for response
      .sPoly = "beta"
      nT1 = 1
      nT2 = 3


      bolTrue <- vBeta != 0
      vTrueValue <- vBeta[bolTrue]

      nG <- length(vBeta)
      nGroupSize <- rep(nLag, nG) # group Index

      # dgp_exp_almon_lag(nLag = vLag, 0.005, -0.05)

      ## set seed
      nSeed <- 1001:(1001 + nSimulation)

      ## direction
      sDir <- "inst/simulation"
      # sDir <- here::here("inst/data/simulation")

      nSim <- here::here(sDir, "input")
      nSim <- nSim[1:10]


      #   ____________________________________________________________________________
      #   Create Models (if necessary)                                            ####

      logger::log_info("Creating Stan models")

      lCMDmodels <- vector("list", length(stanmodels))
      dirStan <- here::here("inst/stan")
      sStanFiles <- dir(dirStan)[grep("stan", dir(dirStan))]

      for (i in seq_along(sStanFiles)) {
            sStanModel <- glue::glue("{dirStan}/{sStanFiles[i]}")
            lCMDmodels[[i]] <- cmdstanr::cmdstan_model(stan_file = sStanModel)
            #stanc_options = list("01"))
      }

      sModelName <- stringr::str_to_title(gsub("_", " ", stringr::str_extract(sStanFiles, "[^\\.]+")))

      names(lCMDmodels) <- sModelName

      logger::log_info("Stan models creation done: {length(sModelName)} models")

      lCMDmodels <- lCMDmodels[-which(names(lCMDmodels) == "Horseshoe Plus Wo Group")]
      sModelName <- sModelName[!(sModelName == "Horseshoe Plus Wo Group")]

      #   ____________________________________________________________________________
      #   Prepare Posterior Beta                                                  ####


      #   ____________________________________________________________________________
      #   Create Models (if necessary)                                            ####

      logger::log_info("Creating Stan models")

      lCMDmodels <- vector("list", length(stanmodels))
      dirStan <- here::here("inst/stan")
      sStanFiles <- dir(dirStan)[grep("stan", dir(dirStan))]

      for (i in seq_along(sStanFiles)) {
            sStanModel <- glue::glue("{dirStan}/{sStanFiles[i]}")
            lCMDmodels[[i]] <- cmdstanr::cmdstan_model(stan_file = sStanModel)
            #stanc_options = list("01"))
      }

      sModelName <- stringr::str_to_title(gsub("_", " ", stringr::str_extract(sStanFiles, "[^\\.]+")))

      names(lCMDmodels) <- sModelName

      logger::log_info("Stan models creation done: {length(sModelName)} models")

      lCMDmodels <- lCMDmodels[-which(names(lCMDmodels) == "Horseshoe Plus Wo Group")]
      sModelName <- sModelName[!(sModelName == "Horseshoe Plus Wo Group")]

      #   ____________________________________________________________________________
      #   Prepare Posterior Beta                                                  ####

      arrRawOut <- arrow::open_dataset(here::here("inst", "simulation",
                                                  "output", "02_raw_extracted"),
                                       partitioning = c("model", "simulation"))

      ## models to compare
      vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
      lSeriesModel <- vector("list", length = length(vMDcompare))






      i <- 1
      nSimRunSave <- helper_create_number_name(i)

      sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

      lData <- readRDS(sDirInputLoad)

      y_train <- lData$model_data$y_train
      x_train <- lData$x_raw[1:(6 + 3*nTrain), ]

      x_test <- lData$x_raw[601:753, ]
      y_test <- lData$model_data$y_test

      dfX <- lData$x_raw

      dfXTrain <- lData$x_raw[1:(6 + 3*nTrain), ]


      ## playing with PCA
      dfXAdj <- apply(dfX, 2, function(x) x - mean(x))
      mdPCA <- prcomp(dfXAdj, center = FALSE, scale. = FALSE)

      mEigenVector <- mdPCA$rotation

      ## use first 5 PCAs
      dfFeature <- dfXAdj %*% mEigenVector[, 1:5]

      dfXBack <- dfFeature %*% t(mEigenVector[, 1:5])


      kmax <- 50

      kInd <- 5

      dfFeature <- dfXAdj %*% mEigenVector[, 1:kInd]

      lCP <- vector("numeric", length = kmax)


      for (kInd in 1:kmax) {

         dfFeature <- dfXAdj %*% mEigenVector[, 1:kInd, drop = FALSE]
         mdReg <- lm(dfXAdj ~ dfFeature - 1)

         nN <- nrow(dfFeature)
         nT <- ncol(dfXAdj)

         nVKF <- sum(apply(mdReg$residuals, 2, function(x) sum(x^2))) / (nN * nT)
         lCP[kInd] <- log(nVKF) + kInd*( (nN + nT) / (nN * nT)) * log(min(nN, nT))
      }

      unlist(lCP)
      which.max(unlist(lCP))


      #   ____________________________________________________________________________
      #   PCA with training data                                                  ####

      ## playing with PCA

      ## prep training data
      dfXTrain <- lData$x_raw[1:(6 + 3 * nTrain), ]
      dfXAdj <- apply(dfXTrain, 2, function(x) x - mean(x))
      mdPCA <- prcomp(dfXAdj, center = FALSE, scale. = FALSE)

      mEigenVector <- mdPCA$rotation
      kInd <- 15

      ## use first 5 PCAs
      maFeature <- dfXAdj %*% mEigenVector[, 1:kInd]

      lFeature <- dgp_pred_align(nTrain, maFeature, nFreq = 3, nLag = 6)
      maTest <- dgp_pred_align(50, dfXTest, nFreq = 3, nLag = 6)
      maTest <- do.call("cbind", maTest)


      ## model
      maFeatureAligned <- do.call("cbind", lFeature)
      dfFeature <- as.data.frame(maFeatureAligned)
      colnames(dfFeature) <- paste0("x_", 1:ncol(dfFeature))
      mdFactor <- lm(y_train ~ . - 1, data = dfFeature)

      ## prep test data
      dfXAll <- lData$x_raw
      dfXAdjAll <- apply(dfXAll, 2, function(x) x - mean(x))

      ## PCA
      mdPCA_all <- prcomp(dfXAdjAll, center = FALSE, scale. = FALSE)
      mEigenVectorAll <- mdPCA_all$rotation
      maFeatureAll <- dfXAdjAll %*% mEigenVectorAll[, 1:kInd]

      ## align test data
      lFeatureAll <- dgp_pred_align(nY, maFeatureAll, nFreq = 3, nLag = 6)
      maTransAll <- do.call("cbind", lFeatureAll)
      maTest <- maTransAll[201:250, ]
      dfTest <- as.data.frame(maTest)
      colnames(dfTest) <- paste0("x_", 1:ncol(dfTest))

      predict.lm(mdFactor, dfTest)

}

