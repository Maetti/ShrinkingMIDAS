masterPosteriorChecking <- function() {

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



      #   ____________________________________________________________________________
      #   Prediction                                                              ####

      lResults <- vector("list", length = nSimulation)

      for (i in 1:nSimulation) {

         logger::log_info("Simulation: {i}")

            ## get simulated data
            logger::log_info("Loading simulation data...")

            nSimRunSave <- helper_create_number_name(i)

            sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

            lData <- readRDS(sDirInputLoad)

            y_train <- lData$model_data$y_train
            x_train <- lData$x_raw[1:(6 + 3*nTrain), ]

            x_test <- lData$x_raw[(nTrain + 1):nY, ]
            y_test <- lData$model_data$y_test




            ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
            ### MIDAS WITH ALMON POLYNOMIAL                                             ####


            logger::log_info("MIDAS with Almon Polynomial")

            # It is possible to re-estimate the NLS problem with the different algorithm using as starting
            # values the final solution of the previous algorithm. For example it is known, that the default
            # algorithm in nls is sensitive to starting values. So first we can use the standard Nelder-Mead
            # algorithm to find “more feasible” starting values and then use nls to get the final result:

            # In order to improve the convergence it is possible to use user defined gradient functions. To
            # use them it is necessary to define the gradient function of the restriction. For example for the
            # nealmon restriction the gradient function is defined in the following way:

            ## prepare data
            y_train2 <- as.matrix(c(NA, y_train, NA))
            sFormula <- as.formula(paste0("y_train2 ~ ",
                                          paste0("mls(x_", 1:50, ", 0:5, 3, nealmon)", collapse = " + "), " - 1"))
            lStart <- lapply(1:50, function(x) {c(1, -0.5)})
            names(lStart) <- paste0("x_", 1:50)

            for (j in 1:50) {
               assign(x = paste0("x_", j), x_train[, j])
               assign(x = paste0("x_test_", j), x_test[, j])
            }

            ## estimate model
            md_midas <- midasr::midas_r(formula = sFormula, start = lStart)

            ## check design matrix
            checkY <- md_midas$model[, 1]
            checkYTrue <- lData$model_data$y_train

            if (all(checkY == checkYTrue)) {
               logger::log_success("y predict as expected")
            } else {
               logger::log_error("y predict not as expected")
            }

            ## check design matrix
            checkDesign <- unname(md_midas$model[, -1])
            checkXTrue <- lData$model_data$x_train

            if (all(checkDesign == checkXTrue, na.rm = TRUE)) {
               logger::log_success("design matrix as expected")
            } else {
               logger::log_error("design matrix not as expected")
            }

            ## prediction
            pred_midas <- lData$model_data$x_test %*% md_midas$midas_coefficients
            colnames(pred_midas) <- "midas"


            ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
            ### UMIDAS                                                                  ####

            logger::log_info("U-MIDAS")

            dfTrainUmidas <- as.data.frame(cbind(lData$model_data$y_train, lData$model_data$x_train))
            colnames(dfTrainUmidas) <- c("y", paste0("x", 1:300))

            md_umidas <- lm(y ~ . - 1, data = dfTrainUmidas)
            md_umidas$coefficients[is.na(md_umidas$coefficients)] <- 0

            pred_umidas <- lData$model_data$x_test %*% md_umidas$coefficients

            dfTestUmidas <- as.data.frame(lData$model_data$x_test)
            colnames(dfTestUmidas) <- paste0("x", 1:300)
            pred_umidas2 <- predict(md_umidas, dfTestUmidas)

            if (all(pred_umidas2 == pred_umidas)) {
               logger::log_success("predicting u-midas")
            } else {
               logger::log_error("predicting u-midas strange")
            }

            colnames(pred_umidas) <- "u_midas"


            ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
            ### FACTOR MIDAS                                                            ####

            logger::log_info("Factor MIDAS")

            pred_factor_midas <- factor_midas(lData)
            pred_factor_midas <- as.matrix(pred_factor_midas)
            colnames(pred_factor_midas) <- "factor_midas"


            ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
            ### SPARSE GROUP LASSO WITH U-MIDAS                                         ####

            logger::log_info("Sparse Group Lasso with U-MIDAS")

            nGroup <- sort(rep(1:50, 6))
            md_sgl = SGL::SGL(data = list(y = lData$model_data$y_train, x = lData$model_data$x_train),
                           index = nGroup, standardize = TRUE,
                           type = "linear")

            pred_sgl <- SGL::predictSGL(x = md_sgl, newX = lData$model_data$x_test, lam = 8)

            colnames(pred_sgl) <- "sgl"

            ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
            ### BAYESIAN SHRINKAGE MODELS                                               ####

            logger::log_info("Bayesian Models")

            tblBayesPred <-
                  arrRawOut %>%
                  dplyr::filter(simulation  == i, parameter == "y_pred") %>%
                  dplyr::select(simulation, model, key, value) %>%
                  dplyr::collect()

            tblBayesRes <-
               tblBayesPred %>%
                  dplyr::group_by(simulation, model, key) %>%
                  dplyr::summarise(value = median(value, na.rm = TRUE)) %>%
                  dplyr::ungroup()


         #   ____________________________________________________________________________
         #   Combining Results                                                       ####

            logger::log_info("Combining Results")

            dfResults <- as.data.frame(cbind(pred_midas, pred_umidas, pred_factor_midas, pred_sgl))
            dfResults$key <- 1:nrow(dfResults)
            dfResults$simulation <- i

            dfResults <-
               dfResults %>%
               dplyr::select(simulation, key, dplyr::everything()) %>%
               tidyr::pivot_longer(names_to = "model", values_to = "value", -c(1, 2))


            tblResOut <-
               dplyr::bind_rows(tblBayesRes, dfResults) %>%
               dplyr::arrange(model, key)

            lResults[[i]] <- tblResOut


      } ## end of for loop


      tblResults <- dplyr::bind_rows(lResults)


      ## true y
      lTrueY <- vector("list", length = nSimulation)
      for (i in 1:nSimulation) {

         logger::log_info("Simulation: {i}")

         ## get simulated data
         logger::log_info("Loading simulation data...")

         nSimRunSave <- helper_create_number_name(i)

         sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

         lData <- readRDS(sDirInputLoad)

         y_test <- lData$model_data$y_test
         dfTrue <- as.data.frame(y_test)
         colnames(dfTrue) <- "value"
         dfTrue$simulation <- i
         dfTrue$key <- 1:nrow(dfTrue)
         dfTrue$model <- "y_test"
         lTrueY[[i]] <- dfTrue
      }

      tblTrueY <- dplyr::bind_rows(lTrueY)


      tblResults <- dplyr::bind_rows(tblResults, tblTrueY) %>% dplyr::arrange(simulation, model, key)
      # ## PLOT TESTING

      lPlot <- vector("list", length = nSimulation)
      for (i in seq_along(lPlot)) {
            nSimPlot <- i
            tblPlot <- tblResults %>% dplyr::filter(simulation == nSimPlot, model != "u_midas")
            dfTrue <- tblTrueY %>% dplyr::filter(simulation == nSimPlot)

            # tblPlot %>%
            #    ggplot2::ggplot(aes(x = key, y = value, color = model)) +
            #    ggplot2::geom_line() +
            #    ggplot2::geom_point(data = dfTrue, aes(x = key, y = value))

            lPlot[[i]] <-
            tblPlot %>%
               ggplot2::ggplot(aes(x = key, y = value, color = model)) +
               ggplot2::geom_line() +
               ggplot2::facet_wrap(~model) +
               ggplot2::annotate(geom = 'point', x = dfTrue$key, y = dfTrue$value) +
               ggplot2::annotate(geom = 'line', x = dfTrue$key, y = dfTrue$value) +
               ggplot2::ggtitle(label = paste0("Simulation: ", i))
      }



   #   ____________________________________________________________________________
   #   Calculate RMSE                                                          ####

      tblTrueY <- tblTrueY %>% dplyr::rename(y_true = value)
      tblRMSEPrep <- dplyr::left_join(tblResults, tblTrueY[, c("y_true", "simulation", "key")], by = c("simulation", "key"))
      nTest <- lData$model_data$nY_test

      tblRMSE <-
         tblRMSEPrep %>%
            dplyr::group_by(simulation, model) %>%
            dplyr::summarise(
               mspe = sqrt(sum((value - y_true)^2) / nTest)
            )


      tblRMSE %>%
         dplyr::filter(model != "u_midas") %>%
         ggplot2::ggplot(aes(x = simulation, y = mspe, color = model)) +
         ggplot2::geom_line()

      tblRMSE %>%
         dplyr::filter(model != "u_midas") %>%
         ggplot2::ggplot(aes(x = model, y = mspe)) +
         ggplot2::geom_boxplot()


   #   ____________________________________________________________________________
   #   Calculate RMSFE                                                         ####

      tblMidas <-
         tblResults %>%
         dplyr::filter(model == "midas") %>%
         dplyr::rename(midas = value) %>%
         dplyr::select(-model)

      tblRMSFEPrep <-
         dplyr::left_join(tblResults, tblTrueY[, c("y_true", "simulation", "key")], by = c("simulation", "key")) %>%
         dplyr::left_join(., tblMidas, by = c("simulation", "key"))


      tblRMSFE <-
         tblRMSFEPrep %>%s
            dplyr::group_by(simulation, model) %>%
            dplyr::summarise(
               sum_model = sum((value - y_true)^2),
               sum_midas = sum((midas - y_true)^2),
               rmsfe = sqrt(sum_model / sum_midas)
            ) %>%
         dplyr::ungroup()

      tblRMSFE %>%
         dplyr::filter(model != "u_midas") %>%
         ggplot2::ggplot(aes(x = simulation, y = rmsfe, color = model)) +
         ggplot2::geom_line()

      tblRMSFE %>%
         dplyr::filter(model != "u_midas") %>%
         ggplot2::ggplot(aes(x = model, y = rmsfe)) +
         ggplot2::geom_boxplot()


      tblRMSFE %>%
         dplyr::group_by(model) %>%
         dplyr::summarise(
            min_rmsfe = min(rmsfe),
            max_rmsfe = max(rmsfe),
            mean_rmsfe = mean(rmsfe),
            median = median(rmsfe)
         )



}

