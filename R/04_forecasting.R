masterForecasting <- function() {

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
      nSimulation <- 100
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

      ## combining results and saveing
      tblForecasting <- dplyr::bind_rows(lResults)

      saveRDS(tblForecasting, file = here::here("inst", "simulation", "output", "05_model_checking", "tblForecasting.rds"))



      #   ____________________________________________________________________________
      #   Forecasting Analyzing                                                   ####

      ## read data
      tblForecasting <- readRDS(file = here::here("inst", "simulation", "output", "05_model_checking", "tblForecasting.rds"))
      tblResults <- tblForecasting


      #   ____________________________________________________________________________
      #   Prepare Data                                                            ####

      ## true y
      lTrueY <- vector("list", length = nSimulation)
      for (i in 1:nSimulation) {

         logger::log_info("Simulation: {i}")

         ## get simulated data
         # logger::log_info("Loading simulation data...")

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


      ## prep data for RMSFE and RRMSFE Calculation
      tblTrueY <- tblTrueY %>% dplyr::rename(y_true = value)

      tblMidas <-
         tblResults %>%
            dplyr::filter(model == "midas") %>%
            dplyr::rename(midas = value) %>%
            dplyr::select(-model)

      tblRMSFEPrep <-
         dplyr::left_join(tblResults, tblTrueY[, c("y_true", "simulation", "key")], by = c("simulation", "key")) %>%
         dplyr::left_join(., tblMidas, by = c("simulation", "key"))



      #   ____________________________________________________________________________
      #   Calculate Forecasting Metrics                                           ####

      nTest <- lData$model_data$nY_test

      tblRMSFE <-
         tblRMSFEPrep %>%
            dplyr::group_by(simulation, model) %>%
            dplyr::summarise(
               mspe = sqrt(sum((value - y_true)^2) / nTest),
               sum_model = sum((value - y_true)^2),
               sum_midas = sum((midas - y_true)^2),
               rmsfe = sqrt(sum_model / sum_midas)
            ) %>%
         dplyr::ungroup()


      dfNameMatching <-
         data.frame(
            "full_model" = c("Horseshoe Plus","Group Lasso",  "Horseshoe", "MIDAS", "Sparse Group Lasso", "Factor MIDAS", "U MIDAS", "Y Test"),
            "model" = c("horseshoe_group_plus", "group_lasso_hierarchical", "horseshoe_group", "midas", "sgl", "factor_midas", "u_midas", "y_test")
         )


      tblRMSFE <- dplyr::left_join(tblRMSFE, dfNameMatching, by = "model")

      tblRMSFE$model <- factor(tblRMSFE$model, levels = dfNameMatching$model)
      tblRMSFE$full_model <- factor(tblRMSFE$full_model, levels = dfNameMatching$full_model)

      #   ____________________________________________________________________________
      #   Tables                                                                  ####

      ## Overview Tables

      fc_overview_table(tblRMSFE, .sMetric = "mspe", sCaption = "RMSFE")

      fc_overview_table(tblRMSFE, .sMetric = "rmsfe", sCaption = "Relative RMSFE")


      ## % of beats model comparison
      maRMSFE <-
            tblRMSFE %>%
            dplyr::filter(model != "y_test") %>%
               dplyr::select(simulation, model, mspe) %>%
               tidyr::pivot_wider(names_from = "model", values_from = "mspe") %>%
               dplyr::select(-simulation) %>%
               as.matrix() %>%
               round(., 5)

      maRMSFE <- maRMSFE[, dfNameMatching$model[-nrow(dfNameMatching)]]


      lCompMatrix <- vector("list", length = ncol(maRMSFE))
      for (i in 1:ncol(maRMSFE)) {

         vComp <- maRMSFE[, i]
         lCompMatrix[[i]] <- apply(apply(maRMSFE, 2, function(x) {vComp <= x}), 2, sum)

      }

      tblComp <- dplyr::bind_rows(lCompMatrix)
      tblComp <- as.matrix(tblComp)
      diag(tblComp) <- NA_real_
      tblComp <- dplyr::as_tibble(tblComp)

      colnames(tblComp) <- dfNameMatching$full_model
      tblComp$model <- colnames(tblComp)

      tblComp <-
         tblComp %>%
            dplyr::select(model, dplyr::everything()) %>%
            # dplyr::mutate_if(is.numeric, ~as.character(.)) %>%
            force()

      tblPlot <-
         tblComp %>%
            tidyr::pivot_longer(names_to = "model_comp", values_to = "value", -1)

      tblPlot$model <- factor(tblPlot$model, levels = rev(dfNameMatching$full_model))
      tblPlot$model_comp <- factor(tblPlot$model_comp, levels = dfNameMatching$full_model)


      tblPlotDiag <-
         tblPlot %>%
            dplyr::filter(model == model_comp)

      gpHeatForecast <-
         tblPlot %>%
            ggplot2::ggplot(ggplot2::aes(model_comp, model, fill = value)) +
               ggplot2::geom_tile(color = "lightgrey",
                                  lwd = .5,
                                  linetype = 1) +
               ggplot2::geom_tile(data = tblPlotDiag, fill = "white") +
               ggplot2::geom_text(ggplot2::aes(label = value)) +
               ggplot2::scale_fill_gradient2(
                  name = "", # changes legend title
                  low = "#e5f5f9",
                  mid = "#99d8c9",
                  high = "#2ca25f",
                  limit = c(0, 100),
                  space = "Lab",
                  guide = "colourbar"
               ) +
               # ggplot2::scale_fill_manual(values = c(0 = "#67a9cf", 100 = "white"),
               #                            name = "") +
               # ggplot2::facet_wrap(~model, nrow = 3) +
               ggplot2::scale_x_discrete(position = "top") +
               ggplot2::xlab(label = "") +
               ggplot2::ylab(label = "") +
               theme_custom_thesis() +
               ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 0))
               # ggplot2::theme(legend.position = "none")


      sForecastHeat <- here::here("inst/simulation/output/04_plots/03_forecasting/forecasting_heatmap.pdf")

      ggplot2::ggsave(sForecastHeat, plot = gpHeatForecast, width = 9, height = 6)
      file.copy(from = sForecastHeat,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/08_forecasting/forecasting_heatmap.pdf",
                overwrite = TRUE)

      #   ____________________________________________________________________________
      #   Plots                                                                   ####

      ## Boxplot per Model

      gpBoxRMSFE <-
         tblRMSFE %>%
            dplyr::filter(!(model %in% c("u_midas", "y_test"))) %>%
               ggplot2::ggplot(aes(x = full_model, y = mspe, fill = full_model)) +
               ggplot2::geom_boxplot(alpha = 0.8,          # Fill transparency
                                     outlier.colour = 1,
                                     show.legend = F) + # Outlier color
               # ggplot2::geom_jitter(color = "lightgrey", alpha = 0.5) +
               ggplot2::xlab(label = "") +
               ggplot2::ylab(label = "") +
               ggplot2::theme(legend.position = "none") +
               ggplot2::scale_fill_manual(values = rev(RColorBrewer::brewer.pal(6, "Blues"))) +  # Fill color
               ggplot2::theme(legend.position = "none") +
               theme_custom_thesis()

      sgpBoxRMSFE <- here::here("inst/simulation/output/04_plots/03_forecasting/boxplot_rmsfe.pdf")
      ggplot2::ggsave(sgpBoxRMSFE, plot = gpBoxRMSFE, width = 9, height = 6)
      file.copy(from = sgpBoxRMSFE,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/08_forecasting/boxplot_rmsfe.pdf",
                overwrite = TRUE)


      gpBoxRRMSFE <-
         tblRMSFE %>%
            dplyr::filter(!(model %in% c("u_midas", "y_test"))) %>%
               ggplot2::ggplot(aes(x = full_model, y = rmsfe, fill = full_model)) +
               ggplot2::geom_boxplot(alpha = 0.8,          # Fill transparency
                                     outlier.colour = 1,
                                     show.legend = F) + # Outlier color
               # ggplot2::geom_jitter(color = "lightgrey", alpha = 0.5) +
               ggplot2::geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
               ggplot2::xlab(label = "") +
               ggplot2::ylab(label = "") +
               ggplot2::theme(legend.position = "none") +
               ggplot2::scale_fill_manual(values = rev(RColorBrewer::brewer.pal(6, "Blues"))) +  # Fill color
               ggplot2::theme(legend.position = "none") +
               theme_custom_thesis()

      sgpBoxRRMSFE <- here::here("inst/simulation/output/04_plots/03_forecasting/boxplot_rrmsfe.pdf")
      ggplot2::ggsave(sgpBoxRRMSFE, plot = gpBoxRRMSFE, width = 9, height = 6)
      file.copy(from = sgpBoxRRMSFE,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/08_forecasting/boxplot_rrmsfe.pdf",
                overwrite = TRUE)



      #   ____________________________________________________________________________
      #   Line Plot                                                               ####


      gpLineRMSFE <-
         tblRMSFE %>%
            dplyr::filter(!(model %in% c("u_midas", "y_test"))) %>%
               ggplot2::ggplot(ggplot2::aes(x = simulation, y = mspe, color = full_model)) +
               ggplot2::geom_line() +
               ggplot2::scale_color_manual(values = rev(RColorBrewer::brewer.pal(6, "Dark2"))) +  # Fill color
               ggplot2::theme(legend.position = "none") +
               ggplot2::labs(x = "Simulation", y = "", color = "Model") +
               theme_custom_thesis()

      sgpLineRM <- here::here("inst/simulation/output/04_plots/03_forecasting/line_rmsfe.pdf")
      ggplot2::ggsave(sgpLineRM, plot = gpLineRMSFE, width = 9, height = 6)
      file.copy(from = sgpLineRM,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/08_forecasting/line_rmsfe.pdf",
                overwrite = TRUE)



      gpLineRRMSFE <-
         tblRMSFE %>%
            dplyr::filter(!(model %in% c("u_midas", "y_test"))) %>%
            ggplot2::ggplot(ggplot2::aes(x = simulation, y = rmsfe, color = full_model)) +
            ggplot2::geom_line() +
            ggplot2::scale_color_manual(values = rev(RColorBrewer::brewer.pal(6, "Dark2"))) +  # Fill color
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = "Simulation", y = "", color = "Model") +
            theme_custom_thesis()

      sgpLineRRM <- here::here("inst/simulation/output/04_plots/03_forecasting/line_rrmsfe.pdf")
      ggplot2::ggsave(sgpLineRRM, plot = gpLineRRMSFE, width = 9, height = 6)
      file.copy(from = sgpLineRRM,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/08_forecasting/line_rrmsfe.pdf",
                overwrite = TRUE)

}

