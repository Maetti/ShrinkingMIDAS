masterPosteriorChecking <- function() {

      #   ____________________________________________________________________________
      #   Set Inputs                                                              ####

      library(dplyr)
      library(ggplot2)

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



      #   ____________________________________________________________________________
      #   Prediction                                                              ####

      for (i in 1:nSimulation) {


            ## get simulated data
            nSimRunSave <- helper_create_number_name(i)

            sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

            lData <- readRDS(sDirInputLoad)

            y_train <- lData$model_data$y_train
            y_test <- lData$model_data$y_test

            x_train <- lData$x_raw[1:(6 + 3*nTrain), ]
            x_test <- lData$x_raw[(nTrain + 1):nY, ]

            ## MIDAS WITH ALMON LAG

            # It is possible to re-estimate the NLS problem with the different algorithm using as starting
            # values the final solution of the previous algorithm. For example it is known, that the default
            # algorithm in nls is sensitive to starting values. So first we can use the standard Nelder-Mead
            # algorithm to find “more feasible” starting values and then use nls to get the final result:

            # In order to improve the convergence it is possible to use user defined gradient functions. To
            # use them it is necessary to define the gradient function of the restriction. For example for the
            # nealmon restriction the gradient function is defined in the following way:

            midasr::mls(x_train[, 1], 0:5, 3, midasr::nealmon)

            midasr::fmls(x_train[, 1], k = 5, m = 3, midasr::nealmon)

            head(midasr::mls(x_train[, 1], 0:5, 3, midasr::nealmon), n = 10)

            mX1 <- midasr::mls(x_train[, 1], 0:5, 3, midasr::nealmon)[2:201, ]
            mX2 <- lData$x_align[[1]][1:200, ]

            all(mX1 == mX2)

            eq_r <- midasr::midas_r(y_train ~ midasr::mls(x_train, 0:5, 3, midasr::nealmon), start = list(x_train = c(1, -0.5)))

            ## U-MIDAS



            ## FACTOR MIDAS



            ## SPARSE GROUP LASSO WITH U-MIDAS


            ## BAYESIAN SHRINKAGE MODELS

            tblBayesPred <-
                  arrRawOut %>%
                  dplyr::filter(simulation  == i, parameter == "y_pred") %>%
                  dplyr::select(simulation, model, key, value) %>%
                  dplyr::collect()

            tblBayesPred %>%
                  dplyr::group_by(simulation, model, key) %>%
                  dplyr::summarise(median = median(value, na.rm = TRUE)) %>%
                  dplyr::ungroup()





}

