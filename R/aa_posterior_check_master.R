masterPosteriorChecking <- function() {

      #   ____________________________________________________________________________
      #   Set Inputs                                                              ####

      library(dplyr)

      ## set logging file
      logger::log_appender(logger::appender_tee(here::here("inst", "simulation", "logging", "model_check",
                                                           glue::glue("{gsub('-', '_', Sys.Date())}_model_check_log"))))
      logger::log_info("Setting input variables")

      ## simulations
      nSimulation <- 10
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



      #   ____________________________________________________________________________
      #   Prepare Posterior Beta                                                  ####

      arrRawOut <- arrow::open_dataset(here::here("inst", "simulation",
                                                  "output", "02_raw_extracted"),
                                       partitioning = c("model", "simulation"))

      ## true beta prep
      tblBetaTrue <- data.frame(
            parameter = "beta",
            key = seq_along(vBeta),
            beta_true = vBeta,
            bolTrue = bolTrue
      )

      ## models to compare
      vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
      lSeriesModel <- vector("list", length = length(vMDcompare))


      tblBeta <-
         arrRawOut %>%
               dplyr::filter(parameter == "beta") %>%  #model == "horseshoe_group_vector") %>%
               dplyr::collect()


      tblConf <-
            tblBeta %>%
                  dplyr::group_by(model, simulation) %>%
                  tidyr::nest() %>%
                  dplyr::mutate(data = purrr::map(data, function(x, bolTrue) {
                        x %>%
                              tidyr::pivot_wider(names_from = "key", values_from = "value") %>%
                              dplyr::select(-run, -chain, -parameter) %>%
                              as.matrix() %>%
                              md_check_beta_create_confusion_matrix_foo(., bolTrue)
                  }, bolTrue)) %>%
                  tidyr::unnest(cols = c(data)) %>%
                  dplyr::ungroup()

      tblConf %>%
         tidyr::pivot_longer(names_to = "")


      ## old function test
      md_check_beta_mse(tblBeta, tblBetaTrue)

      md_check_beta_tpfp_table(tblBeta, tblBetaTrue)
      md_check_beta_true_false_mcc(tblConf)


      md_check_beta_tp_fp_plot(tblConf, .sType = "TP")

      lBetaPlot1 <- md_check_beta_distribution(tblBeta, tblBetaTrue, bolTrueOnly = TRUE, sWrap = "model")
      lBetaPlot1[[1]]
      lBetaPlot1[[2]]

      lBetaModelPlot <- md_check_beta_distribution_per_coef(tblBeta, tblBetaTrue, bolTrueOnly = TRUE)
      lBetaModelPlot[[1]]
      lBetaModelPlot[[2]]



   #   ____________________________________________________________________________
   #   Plots                                                                   ####

      ## MSE
      tblMSE <- md_check_beta_mse(tblBeta, tblBetaTrue)

      tblMSE %>%
         tidyr::pivot_longer(names_to = "key", values_to = "value", -c(1, 2)) %>%
         foo_simulation_as_factor() %>%
         dplyr::filter(key %in% c("mse_var", "mse_bias")) %>%
            ggplot2::ggplot(., ggplot2::aes(x = simulation, y = value, fill = key)) +
               ggplot2::geom_bar(position="stack", stat="identity") +
               ggplot2::facet_wrap(~model, nrow = 3)

      ## MCC
      tblMCC <- md_check_beta_true_false_mcc(tblConf)

      tblMCC %>%
         foo_simulation_as_factor() %>%
            ggplot2::ggplot(., ggplot2::aes(x = simulation, y = mcc)) +
               ggplot2::geom_bar(position="stack", stat="identity") +
               ggplot2::facet_wrap(~model, nrow = 3)


      ##
      md_check_boxplot_tp_fp_mcc(tblMCC, sStat = "mcc")
      md_check_boxplot_tp_fp_mcc(tblMCC, sStat = "tpr")
      md_check_boxplot_tp_fp_mcc(tblMCC, sStat = "fnr")


      ## TPR (sensitivity) and TNR (specificity)


      ## beta coeffsys
      moc_tp_fp_plot(tblConf, .sType = "TP")
      moc_tp_fp_table(lOut$beta, bolTrue)
      moc_beta_distr_plot(lOut$beta, bolTrue, vTrueValue)
      moc_mse_table(lOut$beta, vBeta)
      moc_matthes_corr(tblConf)



      #   ____________________________________________________________________________
      #   Table of Posterior Statistics                                           ####






#   ____________________________________________________________________________
#   Old Stuff                                                               ####


      ### old function test




      maTest <-
      dfTes %>%
            dplyr::arrange(run, chain, key) %>%
            dplyr::mutate(parameter = paste0(parameter, "_", key)) %>%
            dplyr::select(parameter, value, run)

      sSortBeta <- unique(maTest$parameter)

      maTest <- maTest %>%
            dplyr::group_by(parameter, run) %>%
            dplyr::summarise(mean = mean(value)) %>%
            tidyr::pivot_wider(names_from = "parameter", values_from = "mean")

      maTest[, sSortBeta]

      maTest <- as.matrix(maTest[, sSortBeta])
      bayesplot::mcmc_areas(maTest)





}


#   ____________________________________________________________________________
#   From other file....not needed                                           ####

#       for (md in seq_along(vMDcompare)) {
#
#             sModel <- vMDcompare[[md]]
#             logger::log_info("----- {sModel} ------")
#
#             lOverallStat <- vector("list", length = nSimulation)
#             for (i in 1:nSimulation) {
#
#                   logger::log_info("Simulation: {i}/{nSimulation}")
#                   ## load data per simulation
#                   ## true y
#                   tblTrueStat <-
#                         tblYTrue %>%
#                         dplyr::filter(simulation == i)
#
#                   ## posterior y_rep
#                   tblYRep <-
#                         arrRawOut %>%
#                         dplyr::filter(model == sModel,
#                                       parameter == "y_rep",
#                                       simulation == i) %>%
#                         dplyr::collect()
#
#
#
#                   tblTrueStat <- dplyr::tibble("key" = colnames(tblTrueStat), "true" = as.matrix(tblTrueStat)[1, ])
#
#                   ## Compute replicated statistics and combine with true model stat
#                   lOverallStat[[i]] <-
#                         tblYRep %>%
#                         dplyr::group_by(run) %>%
#                         dplyr::summarise(mean = mean(value),
#                                          median = median(value),
#                                          sd = sd(value),
#                                          max = max(value),
#                                          min = min(value)) %>%
#                         dplyr::ungroup() %>%
#                         tidyr::pivot_longer(names_to = "stat", values_to = "value", -1) %>%
#                         dplyr::left_join(., tblTrueStat, by = c("stat" = "key")) %>%
#                         dplyr::group_by(stat) %>%
#                         dplyr::mutate(
#                               p_value = sum(value <= true) / n()
#                         ) %>%
#                         dplyr::ungroup() %>%
#                         dplyr::mutate(simulation = i) %>%
#                         dplyr::select(simulation, dplyr::everything())
#
#             }
#
#             lSeriesModel[[md]] <-
#                   dplyr::bind_rows(lOverallStat) %>%
#                   dplyr::mutate(model = sModel)
#
#       }
#
#       tblSeriesModel <-
#             lSeriesModel %>%
#             dplyr::bind_rows() %>%
#             dplyr::mutate(model = dplyr::case_when(
#                   model == "group_lasso_hierarchical" ~ "Group Lasso",
#                   model == "horseshoe_group_plus" ~ "Horseshoe Plus",
#                   model == "horseshoe_group_vector" ~ "Horseshoe"
#             ))
#
#       ## plotting single
#       lSeriesModel[[3]] %>%
#             moc_yrep_stat_violin(., "median",
#                                  sTitle = "Mean of Replicated Y") +
#             theme_custom_thesis()
#
#
#       ## plot and save
#       dfSeriesStatSave <-
#             data.frame(
#                   "stat" = c("mean", "median", "min", "max", "sd")
#             )
#
#       dfSeriesStatSave$title <- paste0(stringr::str_to_title(dfSeriesStatSave$stat), " of Replicated Y per Simulation")
#       dfSeriesStatSave$file <- here::here("inst/simulation/output/04_plots/01_yrep_check/01_statistics", paste0(dfSeriesStatSave$stat, ".pdf"))
#
#       for (i in 1:nrow(dfSeriesStatSave)) {
#
#             gpSave <-
#                   tblSeriesModel %>%
#                   moc_yrep_stat_violin_facet(., as.character(dfSeriesStatSave[[i, 1]]),
#                                              sTitle = dfSeriesStatSave[[i, 2]], sYlab = "") +
#                   theme_custom_thesis(base_size = 12)
#
#             ggsave(dfSeriesStatSave[[i, 3]], plot = gpSave, width = 9, height = 6)
#       }
#
#       file.copy("/home/matthias/Schreibtisch/SoSe 21/ShrinkingMidas/inst/simulation/output/04_plots/01_yrep_check/01_statistics",
#                 "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/", recursive=TRUE)
#
#
#
#       tblSeriesModel %>%
#             moc_yrep_stat_violin_facet(., "mean",
#                                        sTitle = "Mean of Replicated Y per Simulation", sYlab = "") +
#             theme_custom_thesis()
#
#
#
#       ### test for 100 simulation
#       # lTest <- list()
#       # for (i in 1:10) {
#       #    lTest[[i]] <- tblOverall
#       # }
#       #
#       # tblTest <- bind_rows(lTest)
#       # tblTest$simulation <- sort(rep(1:100, 20000))
#       #
#       # tblTest %>%
#       #    moc_yrep_stat_violin(., "median",
#       #                         sTitle = "Mean of Replicated Y") +
#       #    theme_custom_thesis()
#
#       ## TODO: finish ~facet_wrap per model for all stats and custom theme
#       ## should be align such that one plot per row
#       ## * next p-value summary
#       ## * MSE to y_rep
#       ## * log_likelihood
#       ## ENDED here on 12.8.2022
#
#
#
#       ##  ............................................................................
#       ##  F(x) distribution for each y_i aggregated per Simulation             ####
#
#
#       ## models to compare
#       vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
#       lECDF <- vector("list", length = length(vMDcompare))
#
#       for (md in seq_along(vMDcompare)) {
#
#             sModel <- vMDcompare[[md]]
#             logger::log_info("----- {sModel} ------")
#
#             lOverallStat <- vector("list", length = nSimulation)
#             for (i in 1:nSimulation) {
#
#                   logger::log_info("Simulation: {i}/{nSimulation}")
#                   ## load data per simulation
#                   tblYTrue <-
#                         tblYTrueRaw %>%
#                         dplyr::filter(simulation == i) %>%
#                         dplyr::arrange(key) %>%
#                         dplyr::rename(ytrue = value)
#
#                   ## posterior y_rep
#                   tblYRep <-
#                         arrRawOut %>%
#                         dplyr::filter(model == sModel,
#                                       parameter == "y_rep",
#                                       simulation == i) %>%
#                         dplyr::collect()
#
#
#                   ## compute F(x) for each y_i
#                   lOverallStat[[i]] <-
#                         dplyr::left_join(tblYRep, tblYTrue, by = c("simulation", "key")) %>%
#                         dplyr::group_by(key) %>%
#                         dplyr::summarise(
#                               ecdf = sum(ytrue <= value) / n()
#                         ) %>%
#                         dplyr::ungroup() %>%
#                         dplyr::mutate(simulation = i)
#
#             }
#
#             lECDF[[md]] <-
#                   dplyr::bind_rows(lOverallStat) %>%
#                   dplyr::mutate(model = sModel) %>%
#                   dplyr::select(model, simulation, key, ecdf)
#
#       }
#
#       tblECDF <-
#             lECDF %>%
#             dplyr::bind_rows() %>%
#             dplyr::mutate(model = dplyr::case_when(
#                   model == "group_lasso_hierarchical" ~ "Group Lasso",
#                   model == "horseshoe_group_plus" ~ "Horseshoe Plus",
#                   model == "horseshoe_group_vector" ~ "Horseshoe"
#             )) %>%
#             dplyr::mutate(simulation = as.factor(simulation))
#
#
#       ## plotting
#       tblECDF %>%
#             moc_yrep_ecdf_per_y() +
#             theme_custom_thesis()
#
#
#
#
#       ##  ............................................................................
#       ##  Percentage of y_i outside xx% replication                               ####
#
#       ## use calculation from F(x)/ECDF
#       tblECDF %>%
#             moc_yrep_intervall_cover(., nIntervall = 0.90,
#                                      sTitle = "Intervall check for 90% intervall", sYlab = "Percentage") +
#             theme_custom_thesis()
#
#
#       dfIntervalSave <-
#             data.frame(
#                   "interval" = c(50, 75, 90, 95)
#             )
#
#       dfIntervalSave$title <- paste0("Intervall check for ", dfIntervalSave$interval, "% interval")
#
#       dfIntervalSave$file <- here::here("inst/simulation/output/04_plots/01_yrep_check/02_interval",
#                                         paste0(dfIntervalSave$interval, ".pdf"))
#
#       for (i in 1:nrow(dfIntervalSave)) {
#
#             gpSave <-
#                   tblECDF %>%
#                   moc_yrep_intervall_cover(., nIntervall = dfIntervalSave$interval[[i]] / 100,
#                                            sTitle = dfIntervalSave$title[[i]], sYlab = "Percentage") +
#                   theme_custom_thesis(base_size = 12)
#
#             ggsave(dfIntervalSave$file[[i]], plot = gpSave, width = 9, height = 6)
#       }
#
#       file.copy("/home/matthias/Schreibtisch/SoSe 21/ShrinkingMidas/inst/simulation/output/04_plots/01_yrep_check/02_interval",
#                 "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/", recursive=TRUE)
#
#
#
#       ##  ............................................................................
#       ##  Best/Worst Simulations per Model                                        ####
#
#
#       ## models to compare
#       vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
#       lBestWorst <- vector("list", length = length(vMDcompare))
#
#       for (md in seq_along(vMDcompare)) {
#
#             sModel <- vMDcompare[[md]]
#             logger::log_info("----- {sModel} ------")
#
#             lOverallStat <- vector("list", length = nSimulation)
#             for (i in 1:nSimulation) {
#
#                   logger::log_info("Simulation: {i}/{nSimulation}")
#                   ## load data per simulation
#                   tblYTrue <-
#                         tblYTrueRaw %>%
#                         dplyr::filter(simulation == i) %>%
#                         dplyr::arrange(key) %>%
#                         dplyr::rename(ytrue = value)
#
#                   ## posterior y_rep
#                   tblYRep <-
#                         arrRawOut %>%
#                         dplyr::filter(model == sModel,
#                                       parameter == "y_rep",
#                                       simulation == i) %>%
#                         dplyr::collect()
#
#
#                   ## compute mean/median + intervall and mse
#                   lOverallStat[[i]] <-
#                         tblYRep %>%
#                         dplyr::group_by(key) %>%
#                         dplyr::summarise(mean = mean(value),
#                                          median = median(value),
#                                          q05 = quantile(value, 0.05),
#                                          q95 = quantile(value, 0.95)) %>%
#                         dplyr::left_join(., tblYTrue, by = c("key")) %>%
#                         dplyr::mutate(mse = sqrt(sum((mean - ytrue)^2)))
#
#             }
#
#             lBestWorst[[md]] <-
#                   dplyr::bind_rows(lOverallStat) %>%
#                   dplyr::mutate(model = sModel)
#
#       }
#
#       tblBestWorst <-
#             lBestWorst %>%
#             dplyr::bind_rows() %>%
#             dplyr::mutate(model = dplyr::case_when(
#                   model == "group_lasso_hierarchical" ~ "Group Lasso",
#                   model == "horseshoe_group_plus" ~ "Horseshoe Plus",
#                   model == "horseshoe_group_vector" ~ "Horseshoe"
#             ))
#
#
#       ## plotting
#
#       # checkmate::matchArg(sModel, choices = c("Group Lasso", "Horseshoe Plus", "Horseshoe"))
#       lHorsePlus <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Horseshoe Plus")
#       lHorse <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Horseshoe")
#       lLasso <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Group Lasso")
#
#       lBestWorstPlot <- list(
#             "horseshoe_plus" = lHorsePlus,
#             "horseshoe" = lHorse,
#             "lasso" = lLasso
#       )
#
#       for (i in 1:length(lBestWorstPlot)) {
#
#             gpBest <- lBestWorstPlot[[i]]$best
#             gpWorst <- lBestWorstPlot[[i]]$worst
#
#             sFileSave <- paste0(here::here("inst/simulation/output/04_plots/01_yrep_check/03_best_worst/"),
#                                 names(lBestWorstPlot)[[i]], c("_best.pdf", "_worst.pdf"))
#
#             ggsave(sFileSave[[1]], plot = gpBest, width = 9, height = 6)
#             ggsave(sFileSave[[2]], plot = gpWorst, width = 9, height = 6)
#       }
#
#       file.copy("/home/matthias/Schreibtisch/SoSe 21/ShrinkingMidas/inst/simulation/output/04_plots/01_yrep_check/03_best_worst",
#                 "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/", recursive=TRUE)
#
#
#
# }
#
# #
# #
# #       tblYTrue <- tblYTrue %>% dplyr::rename(true = value)
# #
# #       tblYperYRep <-
# #          tblYRep %>%
# #             dplyr::group_by(key) %>%
# #             dplyr::summarise(mean = mean(value),
# #                              median = median(value),
# #                              sd = sd(value),
# #                              max = max(value),
# #                              min = min(value)) %>%
# #             dplyr::ungroup() %>%
# #             tidyr::pivot_longer(names_to = "stat", values_to = "value", -1) %>%
# #             dplyr::left_join(., tblYTrue, by = c("key"))
# #
# #
# #
# #
# #
# #       # horseshoe_group_plus, group_lasso_hierarchical
# #
# #       tblBeta <-
# #             arrRawOut %>%
# #             dplyr::filter(model == "group_lasso_hierarchical",
# #                           parameter == "beta",
# #                           key %in% c(1:10)) %>%
# #             dplyr::collect()
# #
# #
# #       df1 <- tblBeta %>% dplyr::filter(key == 1, simulation == 8)
# #       plot(density(df1$value))
# #
# #       tblBeta %>%
# #             dplyr::group_by(key, simulation) %>%
# #             dplyr::summarise(
# #                   mean = mean(value),
# #                   sd = sd(value),
# #                   median = median(value),
# #                   q005 = quantile(value, 0.05),
# #                   q095 = quantile(value, 0.95)
# #             ) %>%
# #             View()
# #
# #       tblTheta <-
# #             arrRawOut %>%
# #             dplyr::filter(model == "group_lasso_hierarchical", parameter == "theta",
# #                           key %in% c(1:6)) %>%
# #             dplyr::collect()
# #
# #
# #       df1 <- tblTheta %>% dplyr::filter(key == 1, simulation == 8)
# #       plot(density(df1$value))
# #
# #       tblTheta %>%
# #             dplyr::group_by(key, simulation) %>%
# #             dplyr::summarise(
# #                   mean = mean(value),
# #                   sd = sd(value),
# #                   median = median(value),
# #                   q005 = quantile(value, 0.05),
# #                   q095 = quantile(value, 0.95)
# #             ) %>%
# #             dplyr::arrange(simulation) %>%
# #             View()
# #
# #       dfT <-
# #             tblTheta %>%
# #             dplyr::filter(simulation == 1) %>%
# #             dplyr::group_by(key) %>%
# #             dplyr::summarise(
# #                   mean = mean(value),
# #                   sd = sd(value),
# #                   median = median(value),
# #                   q005 = quantile(value, 0.05),
# #                   q095 = quantile(value, 0.95)
# #             )
# #
# #       nTheta <- dfT$mean
# #
# #       #   ____________________________________________________________________________
# #       #   Old Stuff                                                               ####
# #
# #
# #       lCombined <- purrr::compact(lCombined)
# #
# #       ## getting Beta Coef from Theta
# #       dfTheta <-
# #             lCombined[[2]]$theta %>%
# #             dplyr::filter(chain == "chain:1") %>%
# #             dplyr::select(run, parameter, value) %>%
# #             tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
# #             dplyr::select(-run)
# #
# #
# #       ## Group Index
# #       vIndGroup <- c(0, cumsum(nGroupSize))
# #
# #       ## Lag Matrix
# #       # mQ <- prep_model_lag_matrix(vP, vLag, .bEndpoint = FALSE)
# #
# #       maTheta <- as.matrix(dfTheta)
# #       lBeta <- vector("list", length(nGroupSize))
# #       lWeigthing <- vector("list", length(nGroupSize))
# #       for (i in 1:(length(vIndGroup) - 1)) {
# #             sWPrep <- maTheta[, (vIndGroup[[i]] + 1):vIndGroup[[i + 1]]] # %*% mQ
# #             lWeigthing[[i]] <- sWPrep # / rowSums(sWPrep)
# #             lBeta[[i]] <- rowSums(sWPrep)
# #       }
# #
# #       maBeta <- do.call("cbind", lBeta)
# #       maWeight <- do.call("cbind", lWeigthing)
# #
# #       apply(maBeta, 2, mean)[1:9]
# #       plot(apply(maBeta, 2, mean))
# #       abline(h = 0)
# #       # 1.2, -0.5, 3, 1.1, 0.5, -2.5, -1.5, 0.8, 3.5
# #
# #       ## Weighting Check
# #       lWOutPlot <- vector("list", length(lWeigthing))
# #       for (i in seq_along(lWeigthing)) {
# #
# #             sM <- apply(lWeigthing[[i]], 2, mean)
# #             sD <- apply(lWeigthing[[i]], 2, sd)
# #
# #             s1 <- sM
# #             s2 <- sM + sD
# #             s3 <- sM - sD
# #
# #             lWOutPlot[[i]] <- rbind(s1, s2, s3)
# #       }
# #
# #
# #       #   ____________________________________________________________________________
# #       #   Weighting Plot                                                          ####
# #
# #       ## can be missleading since the theta parameters are divided by beta coeff (which may be close to zero)
# #
# #       tblThetaPlot <-
# #             purrr::map2(lWeigthing, 1:length(lWeigthing),
# #                         function(x, y) {
# #                               tibble::tibble(
# #                                     "theta" = y,
# #                                     "para" = 1:length(apply(x, 2, mean)),
# #                                     "mean" = apply(x, 2, mean),
# #                                     "sdUp" = apply(x, 2, mean) + apply(x, 2, sd),
# #                                     "sdDown" = apply(x, 2, mean) - apply(x, 2, sd)
# #                               )
# #                         }) %>%
# #             dplyr::bind_rows() %>%
# #             tidyr::pivot_longer(-c(1, 2), names_to = "key", values_to = "value")
# #
# #
# #       lGPlot <- vector("list", 50)
# #       for (i in 1:50) {
# #             lGPlot[[i]] <-
# #                   tblThetaPlot %>%
# #                   dplyr::filter(theta == i) %>%
# #                   ggplot2::ggplot(ggplot2::aes(x = para, y = value, group = key)) +
# #                   ggplot2::geom_line()
# #       }
# #
# #       # > apply(maBeta[, 1:9], 2, mean)
# #       # 0.93082954 -0.17315322  2.93918057  0.79205462  0.05328455 -1.90212491 -1.01271238  0.56329044  3.31076842
# #       # 1.04221105 -0.17941081  3.01592710  0.84181671  0.05302292 -1.99368471 -1.09861901  0.85483801  3.44890155
# #
# #       tblThetaPlot %>%
# #             dplyr::filter(theta %in% c(1:9)) %>%
# #             ggplot2::ggplot(ggplot2::aes(x = para, y = value, group = key)) +
# #             ggplot2::geom_line() +
# #             ggplot2::facet_wrap(~theta)
# #
# #       ## should look like that
# #       yExpAlmon <- dgp_exp_almon_lag(nLag = vLag, nT1, nT2)
# #       plot(1:vLag, yExpAlmon, type = "l")
# #
# #       ## End of getting Beta Coef
# #
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   UMIDAS                                                                  ####
# #
# #       xfull <- do.call("cbind", x1$x_align[[1]])
# #
# #       regMod1 <- lm(x1$y ~ xfull - 1)
# #       modCoef1 <- regMod1$coefficients
# #       modCoef1 <- matrix(modCoef1, nrow = 6, byrow = F)
# #
# #       modBeta1 <- rowSums(t(modCoef1))
# #       modWeight1 <- t(modCoef1) / rowSums(t(modCoef1))
# #
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   MIDAS Package                                                           ####
# #
# #       library(midasr)
# #       xMidas <- x1$x_raw[[1]][1:600, 1]
# #       for (i in 1:vK) {
# #             assign(x = paste0("xmidas_", i), value = x1$x_raw[[1]][1:600, i])
# #       }
# #
# #       # eq_r <- midas_r(x1$y ~ mls(xMidas, 0:7, 3, nealmon) - 1, start = list(xMidas = c(1, -0.5)))
# #
# #       sFullMidas <- paste0("eq_r <- midas_r(x1$y ~ ", paste0("mls(", paste0("xmidas_", 1:vK), ", 0:7, 3, nealmon)", collapse = " + "), " - 1, start = list(", paste0(paste0("xmidas_", 1:vK), " = c(1, -0.5)", collapse = ", "), "))")
# #
# #       eval(parse(text = sFullMidas))
# #
# #       summary(eq_r)
# #       maNealmon <- matrix(nrow = 50, ncol = 6)
# #       for (i in 1:50) {
# #             maNealmon[i, ] <- nealmon(p = c(eq_r$coefficients[2 * i - 1], eq_r$coefficients[2 * i]), 6)
# #       }
# #       maNealmon <- maNealmon
# #
# #       maNealBeta <- rowSums(maNealmon)
# #       maNealWeight <- maNealmon / maNealBeta
# #
# #       rowSums(maNealWeight)
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   UMIDAS Lasso                                                            ####
# #       library(lars)
# #
# #       modLasso <- lars::lars(x = xfull, y = x1$y, normalize = FALSE, intercept = FALSE, type = "lasso")
# #
# #       maLasso <- matrix(coef(modLasso)[200, ], ncol = 6, byrow = TRUE)
# #
# #       modBeta1 <- rowSums(maLasso)
# #       modWeight1 <- maLasso / rowSums(maLasso)
# #
# #       #   ____________________________________________________________________________
# #       #   Distribution of Theta                                                   ####
# #
# #
# #
# #       rm(comb_fit_model)
# #       # names(lCombined) <- sModelName
# #       names(lCombined) <- sModelName[1:2]
# #
# #       lOut <- moc_data_combine_all(lCombined)
# #
# #
# #       ## connection to files
# #       dirDGP <- here::here("inst/data/01_dgp/")
# #
# #       dirExtract <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_extract/")
# #       l1 <- lapply(dirExtract, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
# #       arrExtract <- arrow::open_dataset(l1)
# #
# #
# #       dirSummary <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_summary/")
# #       l1 <- lapply(dirSummary, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
# #       arrSummary <- arrow::open_dataset(l1)
# #
# #       ## load data
# #       tblTheta <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "theta") %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, -parameter_2)
# #
# #       tblBeta <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "beta") %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, -parameter_2)
# #
# #       moc_create_confusion_matrix(tblBeta, bolTrue)
# #
# #       ## all true Y values
# #       fSims <- dir(dirDGP)
# #       lYTrue_raw <- vector("list", length(fSims))
# #       for (i in seq_along(nSimulation)) {
# #             sFile <- dir(paste0(dirDGP, "/", fSims[[i]], "/raw_input"))
# #             s1 <- readRDS(paste0(dirDGP, "/", fSims[[i]], "/raw_input/", sFile))
# #             lYTrue_raw[[i]] <- data.frame(simulation = i,
# #                                           key = seq(1, s1$nY_train +  s1$nY_test),
# #                                           value = c(s1$y_train, s1$y_test))
# #       }
# #       tblYTrue <- dplyr::bind_rows(lYTrue_raw) %>% dplyr::as_tibble()
# #
# #
# #       tblYPred <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "y_pred") %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, parameter_2)
# #
# #
# #       tblYRep <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "y_rep", simulation == 1) %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, parameter_2)
# #
# #       tbl1 <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "y_rep", simulation == 1) %>%
# #             dplyr::collect() %>%
# #             dplyr::group_by(model, parameter_1) %>%
# #             dplyr::filter(value == min(value)) %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, parameter_2)
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   Model Checking                                                          ####
# #
# #       ## Note: maybe only read in only every 4 iteration
# #
# #       ## diagnose
# #       ## TODO:
# #       ## * min/max of rhat and neff for each parameter and model simulation
# #       ## * line plot with % divergence per simulation and group by model
# #       ## * test statistic plot from bayesplot package
# #       ## * test stat plot for with simulation on x and p-value on y, group by model
# #       moc_diagnose_plot(lOut$rhat)
# #
# #
# #       ## beta coeffsys
# #       moc_tp_fp_plot(lOut$confusion, .sType = "TP")
# #       moc_tp_fp_table(lOut$beta, bolTrue)
# #       moc_beta_distr_plot(lOut$beta, bolTrue, vTrueValue)
# #       moc_mse_table(lOut$beta, vBeta)
# #       moc_matthes_corr(lOut$confusion)
# #
# #       ## y rep
# #       moc_yrep_plot(lOut$yrep, "mean", nY)
# #       moc_yrep_plot_per_simulation(lOut$yrep, 2)
# #       moc_yrep_95_covered(lOut$yrep)
# #
# #
# #       ## Table of overall fit with p values for yrep
# #       tbl1 %>%
# #             dplyr::group_by(model, iteration) %>%
# #             dplyr::summarise(min = min(value), max = max(value), mean = mean(value))
# #
# #       ## Visualize Fit for each simulation
# #       tbl1 <-
# #             tblYRep %>%
# #             dplyr::filter(simulation == 1) %>%
# #             dplyr::select(iteration, value, model, parameter) %>%
# #             split(., tblYRep$model) %>%
# #             purrr::map(., function(x) {
# #                   x %>%
# #                         dplyr::select(-model) %>%
# #                         tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
# #                         dplyr::select(-iteration) %>%
# #                         as.matrix() %>%
# #                         unname()
# #             })
# #
# #       vGroup <- rep(names(tbl1), purrr::map_dbl(tbl1, ncol))
# #
# #       tbl1 <- tbl1 %>% purrr::reduce(., dplyr::bind_cols) %>% as.matrix() %>% unname()
# #
# #       s1 <- tblYTrain %>% dplyr::filter(simulation == 1)
# #       s2 <- rep(s1$value, 5)
# #
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "min")
# #       q25 <- function(y) quantile(y, 0.25)
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "min")
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "max")
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "mean")
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "sd")
# #
# #       q25 <- function(y) quantile(y, 0.25)
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "q25")
# #
# #
# #       ppc_ribbon_grouped(s2, tbl1, group = vGroup)
# #
# #       ppc_intervals_grouped(s2, tbl1, group = vGroup)
# #       ppc_scatter_avg_grouped(s2, tbl1, group = vGroup)
# #       ## log_lik
# #
# # }
# #
# #
# #
# #
# # foo_fit_model <- function(sName, lStanModel, lData) {
# #       tmp_envir <- new.env(parent = baseenv())
# #       tmp_envir$model <- lStanModel[[sName]]
# #       print(tmp_envir$model)
# #       tmp_envir$fit <- rstan::sampling(tmp_envir$model, data = lData)
# #
# #       tmp_envir$fit
# # }
# #
# #
# # foo_allModel_combi <- function(x, y) {
# #       lOut <- list(x)
# #       names(lOut) <- y
# #       lOut
# # }
# #
# #
# # create_save_model_input <- function(lInput) {
# #
# #       for (i in seq_along(lInput)) {
# #             saveRDS(lInput[[i]], paste0(here::here("inst/data/input/"), i, ".rds"))
# #       }
# #
# # }





masterPosteriorChecking <- function() {

   #   ____________________________________________________________________________
   #   Set Inputs                                                              ####

   library(dplyr)

   ## set logging file
   logger::log_appender(logger::appender_tee(here::here("inst", "simulation", "logging", "model_check",
                                                        glue::glue("{gsub('-', '_', Sys.Date())}_model_check_log"))))
   logger::log_info("Setting input variables")

   ## simulations
   nSimulation <- 10
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



   #   ____________________________________________________________________________
   #   Prepare Posterior Beta                                                  ####

   arrRawOut <- arrow::open_dataset(here::here("inst", "simulation",
                                               "output", "02_raw_extracted"),
                                    partitioning = c("model", "simulation"))

   ## true beta prep
   tblBetaTrue <- data.frame(
      parameter = "beta",
      key = seq_along(vBeta),
      beta_true = vBeta,
      bolTrue = bolTrue
   )

   ## models to compare
   vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
   lSeriesModel <- vector("list", length = length(vMDcompare))


   tblBeta <-
      arrRawOut %>%
      dplyr::filter(parameter == "beta") %>%  #model == "horseshoe_group_vector") %>%
      dplyr::collect()


   tblConf <-
      tblBeta %>%
      dplyr::group_by(model, simulation) %>%
      tidyr::nest() %>%
      dplyr::mutate(data = purrr::map(data, function(x, bolTrue) {
         x %>%
            tidyr::pivot_wider(names_from = "key", values_from = "value") %>%
            dplyr::select(-run, -chain, -parameter) %>%
            as.matrix() %>%
            md_check_beta_create_confusion_matrix_foo(., bolTrue)
      }, bolTrue)) %>%
      tidyr::unnest(cols = c(data)) %>%
      dplyr::ungroup()

   tblConf %>%
      tidyr::pivot_longer(names_to = "")


   ## old function test
   md_check_beta_mse(tblBeta, tblBetaTrue)

   md_check_beta_tpfp_table(tblBeta, tblBetaTrue)
   md_check_beta_true_false_mcc(tblConf)


   md_check_beta_tp_fp_plot(tblConf, .sType = "TP")

   lBetaPlot1 <- md_check_beta_distribution(tblBeta, tblBetaTrue, bolTrueOnly = TRUE, sWrap = "model")
   lBetaPlot1[[1]]
   lBetaPlot1[[2]]

   lBetaModelPlot <- md_check_beta_distribution_per_coef(tblBeta, tblBetaTrue, bolTrueOnly = TRUE)
   lBetaModelPlot[[1]]
   lBetaModelPlot[[2]]



   #   ____________________________________________________________________________
   #   Table of Posterior Statistics                                           ####





   #   ____________________________________________________________________________
   #   Plots                                                                   ####

   ## MSE
   tblMSE <- md_check_beta_mse(tblBeta, tblBetaTrue)

   tblMSE %>%
      tidyr::pivot_longer(names_to = "key", values_to = "value", -c(1, 2)) %>%
      foo_simulation_as_factor() %>%
      dplyr::filter(key %in% c("mse_var", "mse_bias")) %>%
      ggplot2::ggplot(., ggplot2::aes(x = simulation, y = value, fill = key)) +
      ggplot2::geom_bar(position="stack", stat="identity") +
      ggplot2::facet_wrap(~model, nrow = 3)

   ## MCC
   tblMCC <- md_check_beta_true_false_mcc(tblConf)

   tblMCC %>%
      foo_simulation_as_factor() %>%
      ggplot2::ggplot(., ggplot2::aes(x = simulation, y = mcc)) +
      ggplot2::geom_bar(position="stack", stat="identity") +
      ggplot2::facet_wrap(~model, nrow = 3)


   ## TPR (sensitivity) and TNR (specificity)


   ## beta coeffsys
   moc_tp_fp_plot(tblConf, .sType = "TP")
   moc_tp_fp_table(lOut$beta, bolTrue)
   moc_beta_distr_plot(lOut$beta, bolTrue, vTrueValue)
   moc_mse_table(lOut$beta, vBeta)
   moc_matthes_corr(tblConf)


   ### old function test




   maTest <-
      dfTes %>%
      dplyr::arrange(run, chain, key) %>%
      dplyr::mutate(parameter = paste0(parameter, "_", key)) %>%
      dplyr::select(parameter, value, run)

   sSortBeta <- unique(maTest$parameter)

   maTest <- maTest %>%
      dplyr::group_by(parameter, run) %>%
      dplyr::summarise(mean = mean(value)) %>%
      tidyr::pivot_wider(names_from = "parameter", values_from = "mean")

   maTest[, sSortBeta]

   maTest <- as.matrix(maTest[, sSortBeta])
   bayesplot::mcmc_areas(maTest)





}


#   ____________________________________________________________________________
#   From other file....not needed                                           ####

#       for (md in seq_along(vMDcompare)) {
#
#             sModel <- vMDcompare[[md]]
#             logger::log_info("----- {sModel} ------")
#
#             lOverallStat <- vector("list", length = nSimulation)
#             for (i in 1:nSimulation) {
#
#                   logger::log_info("Simulation: {i}/{nSimulation}")
#                   ## load data per simulation
#                   ## true y
#                   tblTrueStat <-
#                         tblYTrue %>%
#                         dplyr::filter(simulation == i)
#
#                   ## posterior y_rep
#                   tblYRep <-
#                         arrRawOut %>%
#                         dplyr::filter(model == sModel,
#                                       parameter == "y_rep",
#                                       simulation == i) %>%
#                         dplyr::collect()
#
#
#
#                   tblTrueStat <- dplyr::tibble("key" = colnames(tblTrueStat), "true" = as.matrix(tblTrueStat)[1, ])
#
#                   ## Compute replicated statistics and combine with true model stat
#                   lOverallStat[[i]] <-
#                         tblYRep %>%
#                         dplyr::group_by(run) %>%
#                         dplyr::summarise(mean = mean(value),
#                                          median = median(value),
#                                          sd = sd(value),
#                                          max = max(value),
#                                          min = min(value)) %>%
#                         dplyr::ungroup() %>%
#                         tidyr::pivot_longer(names_to = "stat", values_to = "value", -1) %>%
#                         dplyr::left_join(., tblTrueStat, by = c("stat" = "key")) %>%
#                         dplyr::group_by(stat) %>%
#                         dplyr::mutate(
#                               p_value = sum(value <= true) / n()
#                         ) %>%
#                         dplyr::ungroup() %>%
#                         dplyr::mutate(simulation = i) %>%
#                         dplyr::select(simulation, dplyr::everything())
#
#             }
#
#             lSeriesModel[[md]] <-
#                   dplyr::bind_rows(lOverallStat) %>%
#                   dplyr::mutate(model = sModel)
#
#       }
#
#       tblSeriesModel <-
#             lSeriesModel %>%
#             dplyr::bind_rows() %>%
#             dplyr::mutate(model = dplyr::case_when(
#                   model == "group_lasso_hierarchical" ~ "Group Lasso",
#                   model == "horseshoe_group_plus" ~ "Horseshoe Plus",
#                   model == "horseshoe_group_vector" ~ "Horseshoe"
#             ))
#
#       ## plotting single
#       lSeriesModel[[3]] %>%
#             moc_yrep_stat_violin(., "median",
#                                  sTitle = "Mean of Replicated Y") +
#             theme_custom_thesis()
#
#
#       ## plot and save
#       dfSeriesStatSave <-
#             data.frame(
#                   "stat" = c("mean", "median", "min", "max", "sd")
#             )
#
#       dfSeriesStatSave$title <- paste0(stringr::str_to_title(dfSeriesStatSave$stat), " of Replicated Y per Simulation")
#       dfSeriesStatSave$file <- here::here("inst/simulation/output/04_plots/01_yrep_check/01_statistics", paste0(dfSeriesStatSave$stat, ".pdf"))
#
#       for (i in 1:nrow(dfSeriesStatSave)) {
#
#             gpSave <-
#                   tblSeriesModel %>%
#                   moc_yrep_stat_violin_facet(., as.character(dfSeriesStatSave[[i, 1]]),
#                                              sTitle = dfSeriesStatSave[[i, 2]], sYlab = "") +
#                   theme_custom_thesis(base_size = 12)
#
#             ggsave(dfSeriesStatSave[[i, 3]], plot = gpSave, width = 9, height = 6)
#       }
#
#       file.copy("/home/matthias/Schreibtisch/SoSe 21/ShrinkingMidas/inst/simulation/output/04_plots/01_yrep_check/01_statistics",
#                 "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/", recursive=TRUE)
#
#
#
#       tblSeriesModel %>%
#             moc_yrep_stat_violin_facet(., "mean",
#                                        sTitle = "Mean of Replicated Y per Simulation", sYlab = "") +
#             theme_custom_thesis()
#
#
#
#       ### test for 100 simulation
#       # lTest <- list()
#       # for (i in 1:10) {
#       #    lTest[[i]] <- tblOverall
#       # }
#       #
#       # tblTest <- bind_rows(lTest)
#       # tblTest$simulation <- sort(rep(1:100, 20000))
#       #
#       # tblTest %>%
#       #    moc_yrep_stat_violin(., "median",
#       #                         sTitle = "Mean of Replicated Y") +
#       #    theme_custom_thesis()
#
#       ## TODO: finish ~facet_wrap per model for all stats and custom theme
#       ## should be align such that one plot per row
#       ## * next p-value summary
#       ## * MSE to y_rep
#       ## * log_likelihood
#       ## ENDED here on 12.8.2022
#
#
#
#       ##  ............................................................................
#       ##  F(x) distribution for each y_i aggregated per Simulation             ####
#
#
#       ## models to compare
#       vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
#       lECDF <- vector("list", length = length(vMDcompare))
#
#       for (md in seq_along(vMDcompare)) {
#
#             sModel <- vMDcompare[[md]]
#             logger::log_info("----- {sModel} ------")
#
#             lOverallStat <- vector("list", length = nSimulation)
#             for (i in 1:nSimulation) {
#
#                   logger::log_info("Simulation: {i}/{nSimulation}")
#                   ## load data per simulation
#                   tblYTrue <-
#                         tblYTrueRaw %>%
#                         dplyr::filter(simulation == i) %>%
#                         dplyr::arrange(key) %>%
#                         dplyr::rename(ytrue = value)
#
#                   ## posterior y_rep
#                   tblYRep <-
#                         arrRawOut %>%
#                         dplyr::filter(model == sModel,
#                                       parameter == "y_rep",
#                                       simulation == i) %>%
#                         dplyr::collect()
#
#
#                   ## compute F(x) for each y_i
#                   lOverallStat[[i]] <-
#                         dplyr::left_join(tblYRep, tblYTrue, by = c("simulation", "key")) %>%
#                         dplyr::group_by(key) %>%
#                         dplyr::summarise(
#                               ecdf = sum(ytrue <= value) / n()
#                         ) %>%
#                         dplyr::ungroup() %>%
#                         dplyr::mutate(simulation = i)
#
#             }
#
#             lECDF[[md]] <-
#                   dplyr::bind_rows(lOverallStat) %>%
#                   dplyr::mutate(model = sModel) %>%
#                   dplyr::select(model, simulation, key, ecdf)
#
#       }
#
#       tblECDF <-
#             lECDF %>%
#             dplyr::bind_rows() %>%
#             dplyr::mutate(model = dplyr::case_when(
#                   model == "group_lasso_hierarchical" ~ "Group Lasso",
#                   model == "horseshoe_group_plus" ~ "Horseshoe Plus",
#                   model == "horseshoe_group_vector" ~ "Horseshoe"
#             )) %>%
#             dplyr::mutate(simulation = as.factor(simulation))
#
#
#       ## plotting
#       tblECDF %>%
#             moc_yrep_ecdf_per_y() +
#             theme_custom_thesis()
#
#
#
#
#       ##  ............................................................................
#       ##  Percentage of y_i outside xx% replication                               ####
#
#       ## use calculation from F(x)/ECDF
#       tblECDF %>%
#             moc_yrep_intervall_cover(., nIntervall = 0.90,
#                                      sTitle = "Intervall check for 90% intervall", sYlab = "Percentage") +
#             theme_custom_thesis()
#
#
#       dfIntervalSave <-
#             data.frame(
#                   "interval" = c(50, 75, 90, 95)
#             )
#
#       dfIntervalSave$title <- paste0("Intervall check for ", dfIntervalSave$interval, "% interval")
#
#       dfIntervalSave$file <- here::here("inst/simulation/output/04_plots/01_yrep_check/02_interval",
#                                         paste0(dfIntervalSave$interval, ".pdf"))
#
#       for (i in 1:nrow(dfIntervalSave)) {
#
#             gpSave <-
#                   tblECDF %>%
#                   moc_yrep_intervall_cover(., nIntervall = dfIntervalSave$interval[[i]] / 100,
#                                            sTitle = dfIntervalSave$title[[i]], sYlab = "Percentage") +
#                   theme_custom_thesis(base_size = 12)
#
#             ggsave(dfIntervalSave$file[[i]], plot = gpSave, width = 9, height = 6)
#       }
#
#       file.copy("/home/matthias/Schreibtisch/SoSe 21/ShrinkingMidas/inst/simulation/output/04_plots/01_yrep_check/02_interval",
#                 "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/", recursive=TRUE)
#
#
#
#       ##  ............................................................................
#       ##  Best/Worst Simulations per Model                                        ####
#
#
#       ## models to compare
#       vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
#       lBestWorst <- vector("list", length = length(vMDcompare))
#
#       for (md in seq_along(vMDcompare)) {
#
#             sModel <- vMDcompare[[md]]
#             logger::log_info("----- {sModel} ------")
#
#             lOverallStat <- vector("list", length = nSimulation)
#             for (i in 1:nSimulation) {
#
#                   logger::log_info("Simulation: {i}/{nSimulation}")
#                   ## load data per simulation
#                   tblYTrue <-
#                         tblYTrueRaw %>%
#                         dplyr::filter(simulation == i) %>%
#                         dplyr::arrange(key) %>%
#                         dplyr::rename(ytrue = value)
#
#                   ## posterior y_rep
#                   tblYRep <-
#                         arrRawOut %>%
#                         dplyr::filter(model == sModel,
#                                       parameter == "y_rep",
#                                       simulation == i) %>%
#                         dplyr::collect()
#
#
#                   ## compute mean/median + intervall and mse
#                   lOverallStat[[i]] <-
#                         tblYRep %>%
#                         dplyr::group_by(key) %>%
#                         dplyr::summarise(mean = mean(value),
#                                          median = median(value),
#                                          q05 = quantile(value, 0.05),
#                                          q95 = quantile(value, 0.95)) %>%
#                         dplyr::left_join(., tblYTrue, by = c("key")) %>%
#                         dplyr::mutate(mse = sqrt(sum((mean - ytrue)^2)))
#
#             }
#
#             lBestWorst[[md]] <-
#                   dplyr::bind_rows(lOverallStat) %>%
#                   dplyr::mutate(model = sModel)
#
#       }
#
#       tblBestWorst <-
#             lBestWorst %>%
#             dplyr::bind_rows() %>%
#             dplyr::mutate(model = dplyr::case_when(
#                   model == "group_lasso_hierarchical" ~ "Group Lasso",
#                   model == "horseshoe_group_plus" ~ "Horseshoe Plus",
#                   model == "horseshoe_group_vector" ~ "Horseshoe"
#             ))
#
#
#       ## plotting
#
#       # checkmate::matchArg(sModel, choices = c("Group Lasso", "Horseshoe Plus", "Horseshoe"))
#       lHorsePlus <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Horseshoe Plus")
#       lHorse <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Horseshoe")
#       lLasso <- moc_yrep_best_worst_lines(tblBestWorst, sModel = "Group Lasso")
#
#       lBestWorstPlot <- list(
#             "horseshoe_plus" = lHorsePlus,
#             "horseshoe" = lHorse,
#             "lasso" = lLasso
#       )
#
#       for (i in 1:length(lBestWorstPlot)) {
#
#             gpBest <- lBestWorstPlot[[i]]$best
#             gpWorst <- lBestWorstPlot[[i]]$worst
#
#             sFileSave <- paste0(here::here("inst/simulation/output/04_plots/01_yrep_check/03_best_worst/"),
#                                 names(lBestWorstPlot)[[i]], c("_best.pdf", "_worst.pdf"))
#
#             ggsave(sFileSave[[1]], plot = gpBest, width = 9, height = 6)
#             ggsave(sFileSave[[2]], plot = gpWorst, width = 9, height = 6)
#       }
#
#       file.copy("/home/matthias/Schreibtisch/SoSe 21/ShrinkingMidas/inst/simulation/output/04_plots/01_yrep_check/03_best_worst",
#                 "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/", recursive=TRUE)
#
#
#
# }
#
# #
# #
# #       tblYTrue <- tblYTrue %>% dplyr::rename(true = value)
# #
# #       tblYperYRep <-
# #          tblYRep %>%
# #             dplyr::group_by(key) %>%
# #             dplyr::summarise(mean = mean(value),
# #                              median = median(value),
# #                              sd = sd(value),
# #                              max = max(value),
# #                              min = min(value)) %>%
# #             dplyr::ungroup() %>%
# #             tidyr::pivot_longer(names_to = "stat", values_to = "value", -1) %>%
# #             dplyr::left_join(., tblYTrue, by = c("key"))
# #
# #
# #
# #
# #
# #       # horseshoe_group_plus, group_lasso_hierarchical
# #
# #       tblBeta <-
# #             arrRawOut %>%
# #             dplyr::filter(model == "group_lasso_hierarchical",
# #                           parameter == "beta",
# #                           key %in% c(1:10)) %>%
# #             dplyr::collect()
# #
# #
# #       df1 <- tblBeta %>% dplyr::filter(key == 1, simulation == 8)
# #       plot(density(df1$value))
# #
# #       tblBeta %>%
# #             dplyr::group_by(key, simulation) %>%
# #             dplyr::summarise(
# #                   mean = mean(value),
# #                   sd = sd(value),
# #                   median = median(value),
# #                   q005 = quantile(value, 0.05),
# #                   q095 = quantile(value, 0.95)
# #             ) %>%
# #             View()
# #
# #       tblTheta <-
# #             arrRawOut %>%
# #             dplyr::filter(model == "group_lasso_hierarchical", parameter == "theta",
# #                           key %in% c(1:6)) %>%
# #             dplyr::collect()
# #
# #
# #       df1 <- tblTheta %>% dplyr::filter(key == 1, simulation == 8)
# #       plot(density(df1$value))
# #
# #       tblTheta %>%
# #             dplyr::group_by(key, simulation) %>%
# #             dplyr::summarise(
# #                   mean = mean(value),
# #                   sd = sd(value),
# #                   median = median(value),
# #                   q005 = quantile(value, 0.05),
# #                   q095 = quantile(value, 0.95)
# #             ) %>%
# #             dplyr::arrange(simulation) %>%
# #             View()
# #
# #       dfT <-
# #             tblTheta %>%
# #             dplyr::filter(simulation == 1) %>%
# #             dplyr::group_by(key) %>%
# #             dplyr::summarise(
# #                   mean = mean(value),
# #                   sd = sd(value),
# #                   median = median(value),
# #                   q005 = quantile(value, 0.05),
# #                   q095 = quantile(value, 0.95)
# #             )
# #
# #       nTheta <- dfT$mean
# #
# #       #   ____________________________________________________________________________
# #       #   Old Stuff                                                               ####
# #
# #
# #       lCombined <- purrr::compact(lCombined)
# #
# #       ## getting Beta Coef from Theta
# #       dfTheta <-
# #             lCombined[[2]]$theta %>%
# #             dplyr::filter(chain == "chain:1") %>%
# #             dplyr::select(run, parameter, value) %>%
# #             tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
# #             dplyr::select(-run)
# #
# #
# #       ## Group Index
# #       vIndGroup <- c(0, cumsum(nGroupSize))
# #
# #       ## Lag Matrix
# #       # mQ <- prep_model_lag_matrix(vP, vLag, .bEndpoint = FALSE)
# #
# #       maTheta <- as.matrix(dfTheta)
# #       lBeta <- vector("list", length(nGroupSize))
# #       lWeigthing <- vector("list", length(nGroupSize))
# #       for (i in 1:(length(vIndGroup) - 1)) {
# #             sWPrep <- maTheta[, (vIndGroup[[i]] + 1):vIndGroup[[i + 1]]] # %*% mQ
# #             lWeigthing[[i]] <- sWPrep # / rowSums(sWPrep)
# #             lBeta[[i]] <- rowSums(sWPrep)
# #       }
# #
# #       maBeta <- do.call("cbind", lBeta)
# #       maWeight <- do.call("cbind", lWeigthing)
# #
# #       apply(maBeta, 2, mean)[1:9]
# #       plot(apply(maBeta, 2, mean))
# #       abline(h = 0)
# #       # 1.2, -0.5, 3, 1.1, 0.5, -2.5, -1.5, 0.8, 3.5
# #
# #       ## Weighting Check
# #       lWOutPlot <- vector("list", length(lWeigthing))
# #       for (i in seq_along(lWeigthing)) {
# #
# #             sM <- apply(lWeigthing[[i]], 2, mean)
# #             sD <- apply(lWeigthing[[i]], 2, sd)
# #
# #             s1 <- sM
# #             s2 <- sM + sD
# #             s3 <- sM - sD
# #
# #             lWOutPlot[[i]] <- rbind(s1, s2, s3)
# #       }
# #
# #
# #       #   ____________________________________________________________________________
# #       #   Weighting Plot                                                          ####
# #
# #       ## can be missleading since the theta parameters are divided by beta coeff (which may be close to zero)
# #
# #       tblThetaPlot <-
# #             purrr::map2(lWeigthing, 1:length(lWeigthing),
# #                         function(x, y) {
# #                               tibble::tibble(
# #                                     "theta" = y,
# #                                     "para" = 1:length(apply(x, 2, mean)),
# #                                     "mean" = apply(x, 2, mean),
# #                                     "sdUp" = apply(x, 2, mean) + apply(x, 2, sd),
# #                                     "sdDown" = apply(x, 2, mean) - apply(x, 2, sd)
# #                               )
# #                         }) %>%
# #             dplyr::bind_rows() %>%
# #             tidyr::pivot_longer(-c(1, 2), names_to = "key", values_to = "value")
# #
# #
# #       lGPlot <- vector("list", 50)
# #       for (i in 1:50) {
# #             lGPlot[[i]] <-
# #                   tblThetaPlot %>%
# #                   dplyr::filter(theta == i) %>%
# #                   ggplot2::ggplot(ggplot2::aes(x = para, y = value, group = key)) +
# #                   ggplot2::geom_line()
# #       }
# #
# #       # > apply(maBeta[, 1:9], 2, mean)
# #       # 0.93082954 -0.17315322  2.93918057  0.79205462  0.05328455 -1.90212491 -1.01271238  0.56329044  3.31076842
# #       # 1.04221105 -0.17941081  3.01592710  0.84181671  0.05302292 -1.99368471 -1.09861901  0.85483801  3.44890155
# #
# #       tblThetaPlot %>%
# #             dplyr::filter(theta %in% c(1:9)) %>%
# #             ggplot2::ggplot(ggplot2::aes(x = para, y = value, group = key)) +
# #             ggplot2::geom_line() +
# #             ggplot2::facet_wrap(~theta)
# #
# #       ## should look like that
# #       yExpAlmon <- dgp_exp_almon_lag(nLag = vLag, nT1, nT2)
# #       plot(1:vLag, yExpAlmon, type = "l")
# #
# #       ## End of getting Beta Coef
# #
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   UMIDAS                                                                  ####
# #
# #       xfull <- do.call("cbind", x1$x_align[[1]])
# #
# #       regMod1 <- lm(x1$y ~ xfull - 1)
# #       modCoef1 <- regMod1$coefficients
# #       modCoef1 <- matrix(modCoef1, nrow = 6, byrow = F)
# #
# #       modBeta1 <- rowSums(t(modCoef1))
# #       modWeight1 <- t(modCoef1) / rowSums(t(modCoef1))
# #
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   MIDAS Package                                                           ####
# #
# #       library(midasr)
# #       xMidas <- x1$x_raw[[1]][1:600, 1]
# #       for (i in 1:vK) {
# #             assign(x = paste0("xmidas_", i), value = x1$x_raw[[1]][1:600, i])
# #       }
# #
# #       # eq_r <- midas_r(x1$y ~ mls(xMidas, 0:7, 3, nealmon) - 1, start = list(xMidas = c(1, -0.5)))
# #
# #       sFullMidas <- paste0("eq_r <- midas_r(x1$y ~ ", paste0("mls(", paste0("xmidas_", 1:vK), ", 0:7, 3, nealmon)", collapse = " + "), " - 1, start = list(", paste0(paste0("xmidas_", 1:vK), " = c(1, -0.5)", collapse = ", "), "))")
# #
# #       eval(parse(text = sFullMidas))
# #
# #       summary(eq_r)
# #       maNealmon <- matrix(nrow = 50, ncol = 6)
# #       for (i in 1:50) {
# #             maNealmon[i, ] <- nealmon(p = c(eq_r$coefficients[2 * i - 1], eq_r$coefficients[2 * i]), 6)
# #       }
# #       maNealmon <- maNealmon
# #
# #       maNealBeta <- rowSums(maNealmon)
# #       maNealWeight <- maNealmon / maNealBeta
# #
# #       rowSums(maNealWeight)
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   UMIDAS Lasso                                                            ####
# #       library(lars)
# #
# #       modLasso <- lars::lars(x = xfull, y = x1$y, normalize = FALSE, intercept = FALSE, type = "lasso")
# #
# #       maLasso <- matrix(coef(modLasso)[200, ], ncol = 6, byrow = TRUE)
# #
# #       modBeta1 <- rowSums(maLasso)
# #       modWeight1 <- maLasso / rowSums(maLasso)
# #
# #       #   ____________________________________________________________________________
# #       #   Distribution of Theta                                                   ####
# #
# #
# #
# #       rm(comb_fit_model)
# #       # names(lCombined) <- sModelName
# #       names(lCombined) <- sModelName[1:2]
# #
# #       lOut <- moc_data_combine_all(lCombined)
# #
# #
# #       ## connection to files
# #       dirDGP <- here::here("inst/data/01_dgp/")
# #
# #       dirExtract <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_extract/")
# #       l1 <- lapply(dirExtract, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
# #       arrExtract <- arrow::open_dataset(l1)
# #
# #
# #       dirSummary <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_summary/")
# #       l1 <- lapply(dirSummary, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
# #       arrSummary <- arrow::open_dataset(l1)
# #
# #       ## load data
# #       tblTheta <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "theta") %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, -parameter_2)
# #
# #       tblBeta <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "beta") %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, -parameter_2)
# #
# #       moc_create_confusion_matrix(tblBeta, bolTrue)
# #
# #       ## all true Y values
# #       fSims <- dir(dirDGP)
# #       lYTrue_raw <- vector("list", length(fSims))
# #       for (i in seq_along(nSimulation)) {
# #             sFile <- dir(paste0(dirDGP, "/", fSims[[i]], "/raw_input"))
# #             s1 <- readRDS(paste0(dirDGP, "/", fSims[[i]], "/raw_input/", sFile))
# #             lYTrue_raw[[i]] <- data.frame(simulation = i,
# #                                           key = seq(1, s1$nY_train +  s1$nY_test),
# #                                           value = c(s1$y_train, s1$y_test))
# #       }
# #       tblYTrue <- dplyr::bind_rows(lYTrue_raw) %>% dplyr::as_tibble()
# #
# #
# #       tblYPred <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "y_pred") %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, parameter_2)
# #
# #
# #       tblYRep <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "y_rep", simulation == 1) %>%
# #             dplyr::collect() %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, parameter_2)
# #
# #       tbl1 <-
# #             arrExtract %>%
# #             dplyr::filter(parameter_1 == "y_rep", simulation == 1) %>%
# #             dplyr::collect() %>%
# #             dplyr::group_by(model, parameter_1) %>%
# #             dplyr::filter(value == min(value)) %>%
# #             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
# #             dplyr::select(-parameter_1, parameter_2)
# #
# #
# #
# #       #   ____________________________________________________________________________
# #       #   Model Checking                                                          ####
# #
# #       ## Note: maybe only read in only every 4 iteration
# #
# #       ## diagnose
# #       ## TODO:
# #       ## * min/max of rhat and neff for each parameter and model simulation
# #       ## * line plot with % divergence per simulation and group by model
# #       ## * test statistic plot from bayesplot package
# #       ## * test stat plot for with simulation on x and p-value on y, group by model
# #       moc_diagnose_plot(lOut$rhat)
# #
# #
# #       ## beta coeffsys
# #       moc_tp_fp_plot(lOut$confusion, .sType = "TP")
# #       moc_tp_fp_table(lOut$beta, bolTrue)
# #       moc_beta_distr_plot(lOut$beta, bolTrue, vTrueValue)
# #       moc_mse_table(lOut$beta, vBeta)
# #       moc_matthes_corr(lOut$confusion)
# #
# #       ## y rep
# #       moc_yrep_plot(lOut$yrep, "mean", nY)
# #       moc_yrep_plot_per_simulation(lOut$yrep, 2)
# #       moc_yrep_95_covered(lOut$yrep)
# #
# #
# #       ## Table of overall fit with p values for yrep
# #       tbl1 %>%
# #             dplyr::group_by(model, iteration) %>%
# #             dplyr::summarise(min = min(value), max = max(value), mean = mean(value))
# #
# #       ## Visualize Fit for each simulation
# #       tbl1 <-
# #             tblYRep %>%
# #             dplyr::filter(simulation == 1) %>%
# #             dplyr::select(iteration, value, model, parameter) %>%
# #             split(., tblYRep$model) %>%
# #             purrr::map(., function(x) {
# #                   x %>%
# #                         dplyr::select(-model) %>%
# #                         tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
# #                         dplyr::select(-iteration) %>%
# #                         as.matrix() %>%
# #                         unname()
# #             })
# #
# #       vGroup <- rep(names(tbl1), purrr::map_dbl(tbl1, ncol))
# #
# #       tbl1 <- tbl1 %>% purrr::reduce(., dplyr::bind_cols) %>% as.matrix() %>% unname()
# #
# #       s1 <- tblYTrain %>% dplyr::filter(simulation == 1)
# #       s2 <- rep(s1$value, 5)
# #
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "min")
# #       q25 <- function(y) quantile(y, 0.25)
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "min")
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "max")
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "mean")
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "sd")
# #
# #       q25 <- function(y) quantile(y, 0.25)
# #       ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "q25")
# #
# #
# #       ppc_ribbon_grouped(s2, tbl1, group = vGroup)
# #
# #       ppc_intervals_grouped(s2, tbl1, group = vGroup)
# #       ppc_scatter_avg_grouped(s2, tbl1, group = vGroup)
# #       ## log_lik
# #
# # }
# #
# #
# #
# #
# # foo_fit_model <- function(sName, lStanModel, lData) {
# #       tmp_envir <- new.env(parent = baseenv())
# #       tmp_envir$model <- lStanModel[[sName]]
# #       print(tmp_envir$model)
# #       tmp_envir$fit <- rstan::sampling(tmp_envir$model, data = lData)
# #
# #       tmp_envir$fit
# # }
# #
# #
# # foo_allModel_combi <- function(x, y) {
# #       lOut <- list(x)
# #       names(lOut) <- y
# #       lOut
# # }
# #
# #
# # create_save_model_input <- function(lInput) {
# #
# #       for (i in seq_along(lInput)) {
# #             saveRDS(lInput[[i]], paste0(here::here("inst/data/input/"), i, ".rds"))
# #       }
# #
# # }













