masterModelChecking <- function() {

      #   ____________________________________________________________________________
      #   Set Inputs                                                              ####

      library(dplyr)
      library(ggplot2)

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
      #   Diagnostic Preparation                                                  ####

      ## all true Y values
      lYTrue_raw <- vector("list", nSimulation)
      for (i in 1:nSimulation) {


            ## get simulated data
            nSimRunSave <- helper_create_number_name(i)

            sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

            lData <- readRDS(sDirInputLoad)
            lDataStanInput <- lData$model_data

            lYTrue_raw[[i]] <- c(lDataStanInput$y_train, lDataStanInput$y_test)
      }

      tblYTrueRaw <-
         lYTrue_raw %>%
            purrr::map2(., 1:nSimulation, function(x, y) {
               dplyr::tibble(simulation = y, key = 1:length(x), value = x)
            }) %>%
            dplyr::bind_rows()

      rm(list = c("lData", "lDataStanInput"))
      gc()

      #   ____________________________________________________________________________
      #   Model Checking                                                          ####

      ## TODO: 30.07.2022 stopped here

      arrRawOut <- arrow::open_dataset(here::here("inst", "simulation",
                                                  "output", "02_raw_extracted"),
                                       partitioning = c("model", "simulation"))



      ## TODO: 02.08.2022

      ## overall statistic
      # lOverallStat <- vector("list", length = nSimulation)
      #
      # lOverallStat[[i]] <- data.frame("simulation" = i,
      #                               "mean" = mean(tblYRep$value),
      #                               "median" = median(tblYRep$value),
      #                               "sd" = sd(tblYRep$value),
      #                               "min" = min(tblYRep$value),
      #                               "max" = max(tblYRep$value))

      ## TODO: 9.8.2022 plots
      # tblMean <- tblYperSeries %>% dplyr::filter(stat == "mean")
      # plot(density(tblMean$value))



      ##  ............................................................................
      ##  Overall Series per each Simulation                                      ####


      ## 9.8.2022 create and plot for all simulation and models
      ## - Boxplots and true value for mean, median, max, min, sd

      ## calculate statistics for true Y
      tblYTrue <-
         tblYTrueRaw %>%
            dplyr::group_by(simulation) %>%
            dplyr::summarise(mean = mean(value),
                             median = median(value),
                             sd = sd(value),
                             max = max(value),
                             min = min(value)) %>%
            dplyr::ungroup()

      ## models to compare
      vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
      lSeriesModel <- vector("list", length = length(vMDcompare))

      for (md in seq_along(vMDcompare)) {

         sModel <- vMDcompare[[md]]
         logger::log_info("----- {sModel} ------")

         lOverallStat <- vector("list", length = nSimulation)
         for (i in 1:nSimulation) {

            logger::log_info("Simulation: {i}/{nSimulation}")
            ## load data per simulation
            ## true y
            tblTrueStat <-
               tblYTrue %>%
               dplyr::filter(simulation == i)

            ## posterior y_rep
            tblYRep <-
               arrRawOut %>%
               dplyr::filter(model == sModel,
                             parameter == "y_rep",
                             simulation == i) %>%
               dplyr::collect()



            tblTrueStat <- dplyr::tibble("key" = colnames(tblTrueStat), "true" = as.matrix(tblTrueStat)[1, ])

            ## Compute replicated statistics and combine with true model stat
            lOverallStat[[i]] <-
               tblYRep %>%
                  dplyr::group_by(run) %>%
                  dplyr::summarise(mean = mean(value),
                                   median = median(value),
                                   sd = sd(value),
                                   max = max(value),
                                   min = min(value)) %>%
                  dplyr::ungroup() %>%
                  tidyr::pivot_longer(names_to = "stat", values_to = "value", -1) %>%
                  dplyr::left_join(., tblTrueStat, by = c("stat" = "key")) %>%
                  dplyr::group_by(stat) %>%
                  dplyr::mutate(
                     p_value = sum(value <= true) / n()
                  ) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(simulation = i) %>%
                  dplyr::select(simulation, dplyr::everything())

            rm(list = c("tblYRep"))
            gc()

         }

         lSeriesModel[[md]] <-
            dplyr::bind_rows(lOverallStat) %>%
            dplyr::mutate(model = sModel)

      }


      tblSeriesModel <-
         lSeriesModel %>%
            dplyr::bind_rows() %>%
            dplyr::mutate(model = dplyr::case_when(
               model == "group_lasso_hierarchical" ~ "Group Lasso",
               model == "horseshoe_group_plus" ~ "Horseshoe Plus",
               model == "horseshoe_group" ~ "Horseshoe"
            ))

      saveRDS(tblSeriesModel, file = here::here("inst", "simulation", "output", "05_model_checking", "tblSeriesModel.rds"))



      #   ____________________________________________________________________________
      #   Analyzing Model Fitting to Y                                            ####

      tblSeriesModel <- readRDS(file = here::here("inst", "simulation", "output", "05_model_checking", "tblSeriesModel.rds"))

      ## plotting single
      # lSeriesModel[[3]] %>%
      #    moc_yrep_stat_violin(., "median",
      #                         sTitle = "Mean of Replicated Y") +
      #    theme_custom_thesis()

      tblSeriesModel %>%
         dplyr::filter(model == "Group Lasso") %>%
            moc_yrep_stat_violin(., "median", sTitle = "Mean of Replicated Y") +
            theme_custom_thesis()


      ## plot and save
      dfSeriesStatSave <-
         data.frame(
            "stat" = c("mean", "median", "min", "max", "sd")
         )

      dfSeriesStatSave$title <- paste0(stringr::str_to_title(dfSeriesStatSave$stat), " of Replicated Y per Simulation")
      dfSeriesStatSave$file <- here::here("inst/simulation/output/04_plots/01_yrep_check/01_statistics", paste0(dfSeriesStatSave$stat, ".pdf"))

      for (i in 1:nrow(dfSeriesStatSave)) {

         logger::log_info("Plot for: {dfSeriesStatSave$stat[[i]]}")

         gpSave <-
            tblSeriesModel %>%
            moc_yrep_stat_violin_facet(., as.character(dfSeriesStatSave[[i, 1]]),
                                       sTitle = dfSeriesStatSave[[i, 2]], sYlab = "") +
            theme_custom_thesis(base_size = 12)

         ggplot2::ggsave(dfSeriesStatSave[[i, 3]], plot = gpSave, width = 9, height = 6)
      }


      rm(list = c("tblSeriesModel", "gpSave", "tblPlot", "tblStat"))
      gc()

      ## TODO: finish ~facet_wrap per model for all stats and custom theme
      ## should be align such that one plot per row
      ## * next p-value summary
      ## * MSE to y_rep
      ## * log_likelihood
      ## ENDED here on 12.8.2022

      ##  ............................................................................
      ##  F(x) distribution for each y_i aggregated per Simulation             ####


      ## models to compare
      vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
      lECDF <- vector("list", length = length(vMDcompare))

      for (md in seq_along(vMDcompare)) {

         sModel <- vMDcompare[[md]]
         logger::log_info("----- {sModel} ------")

         lOverallStat <- vector("list", length = nSimulation)
         for (i in 1:nSimulation) {

            logger::log_info("Simulation: {i}/{nSimulation}")
            ## load data per simulation
            tblYTrue <-
               tblYTrueRaw %>%
                  dplyr::filter(simulation == i) %>%
                  dplyr::arrange(key) %>%
                  dplyr::rename(ytrue = value)

            ## posterior y_rep
            tblYRep <-
               arrRawOut %>%
               dplyr::filter(model == sModel,
                             parameter == "y_rep",
                             simulation == i) %>%
               dplyr::collect()


            ## compute F(x) for each y_i
            lOverallStat[[i]] <-
               dplyr::left_join(tblYRep, tblYTrue, by = c("simulation", "key")) %>%
                  dplyr::group_by(key) %>%
                  dplyr::summarise(
                     ecdf = sum(ytrue <= value) / n()
                  ) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(simulation = i)

         }

         lECDF[[md]] <-
            dplyr::bind_rows(lOverallStat) %>%
            dplyr::mutate(model = sModel) %>%
            dplyr::select(model, simulation, key, ecdf)

      }

      tblECDF <-
         lECDF %>%
         dplyr::bind_rows() %>%
         dplyr::mutate(model = dplyr::case_when(
            model == "group_lasso_hierarchical" ~ "Group Lasso",
            model == "horseshoe_group_plus" ~ "Horseshoe Plus",
            model == "horseshoe_group" ~ "Horseshoe"
         )) %>%
         dplyr::mutate(simulation = as.factor(simulation))


      logger::log_info("Saving: {tblECDF}")
      saveRDS(tblECDF,
              file = here::here("inst", "simulation", "output", "05_model_checking", "tblECDF.rds"))


      ##  ............................................................................
      ##  Percentage of y_i outside xx% replication                               ####

      tblECDF <- readRDS(here::here("inst", "simulation", "output", "05_model_checking", "tblECDF.rds"))

      ## use calculation from F(x)/ECDF
      dfIntervalSave <-
         data.frame(
            "interval" = c(50, 75, 90, 95)
         )

      dfIntervalSave$title <- paste0("Intervall check for ", dfIntervalSave$interval, "% interval")

      dfIntervalSave$file <- here::here("inst/simulation/output/04_plots/01_yrep_check/02_interval",
                                          paste0(dfIntervalSave$interval, ".pdf"))

      for (i in 1:nrow(dfIntervalSave)) {

         gpSave <-
            tblECDF %>%
            moc_yrep_intervall_cover(., nIntervall = dfIntervalSave$interval[[i]] / 100,
                                     sTitle = dfIntervalSave$title[[i]], sYlab = "Percentage") +
            theme_custom_thesis(base_size = 12)

         ggplot2::ggsave(dfIntervalSave$file[[i]], plot = gpSave, width = 9, height = 6)
      }

      ##  ............................................................................
      ##  Boxplot of empirical CDF                                                ####

      gpBoxCDF <- moc_yrep_ecdf_per_y(tblStat = tblECDF, sTitle = "", sYlab = "", sXlab = "Simulation")

      ggplot2::ggsave(here::here("inst/simulation/output/04_plots/01_yrep_check/04_cdf_boxplot/cdf_boxplot.pdf"),
             plot = gpBoxCDF, width = 9, height = 6)

      ##  ............................................................................
      ##  Best/Worst Simulations per Model                                        ####


      ## models to compare
      vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
      lBestWorst <- vector("list", length = length(vMDcompare))

      for (md in seq_along(vMDcompare)) {

         sModel <- vMDcompare[[md]]
         logger::log_info("----- {sModel} ------")

         lOverallStat <- vector("list", length = nSimulation)
         for (i in 1:nSimulation) {

            logger::log_info("Simulation: {i}/{nSimulation}")
            ## load data per simulation
            tblYTrue <-
               tblYTrueRaw %>%
               dplyr::filter(simulation == i) %>%
               dplyr::arrange(key) %>%
               dplyr::rename(ytrue = value)

            ## posterior y_rep
            tblYRep <-
               arrRawOut %>%
               dplyr::filter(model == sModel,
                             parameter == "y_rep",
                             simulation == i) %>%
               dplyr::collect()


            ## compute mean/median + intervall and mse
            lOverallStat[[i]] <-
               tblYRep %>%
                  dplyr::group_by(key) %>%
                  dplyr::summarise(mean = mean(value),
                                   median = median(value),
                                   q05 = quantile(value, 0.05),
                                   q95 = quantile(value, 0.95)) %>%
                  dplyr::left_join(., tblYTrue, by = c("key")) %>%
                  dplyr::mutate(mse = sqrt(sum((mean - ytrue)^2)))

         }

         lBestWorst[[md]] <-
            dplyr::bind_rows(lOverallStat) %>%
            dplyr::mutate(model = sModel)

      }

      tblBestWorst <-
         lBestWorst %>%
         dplyr::bind_rows() %>%
         dplyr::mutate(model = dplyr::case_when(
            model == "group_lasso_hierarchical" ~ "Group Lasso",
            model == "horseshoe_group_plus" ~ "Horseshoe Plus",
            model == "horseshoe_group" ~ "Horseshoe"
         ))

      logger::log_info("Saving: {tblECDF}")
      saveRDS(tblBestWorst,
              file = here::here("inst", "simulation", "output", "05_model_checking", "tblBestWorst.rds"))


      ## plotting

      tblBestWorst <- readRDS(here::here("inst", "simulation", "output", "05_model_checking", "tblBestWorst.rds"))

      # checkmate::matchArg(sModel, choices = c("Group Lasso", "Horseshoe Plus", "Horseshoe"))
      lHorsePlus <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Horseshoe Plus")
      lHorse <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Horseshoe")
      lLasso <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Group Lasso")

      lBestWorstPlot <- list(
         "horseshoe_plus" = lHorsePlus,
         "horseshoe" = lHorse,
         "lasso" = lLasso
      )

      for (i in 1:length(lBestWorstPlot)) {

         gpBest <- lBestWorstPlot[[i]]$best
         gpWorst <- lBestWorstPlot[[i]]$worst

         sFileSave <- paste0(here::here("inst/simulation/output/04_plots/01_yrep_check/03_best_worst/"),
                             names(lBestWorstPlot)[[i]], c("_best.pdf", "_worst.pdf"))

         ggplot2::ggsave(sFileSave[[1]], plot = gpBest, width = 9, height = 6)
         ggplot2::ggsave(sFileSave[[2]], plot = gpWorst, width = 9, height = 6)
      }


}
