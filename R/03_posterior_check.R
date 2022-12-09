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

      ## iterate over models in order to fit data into memory

      ## 11.10.2022
      lMSEAll <- lConfAll <- lIncluAll <- vector("list", length = length(vMDcompare))
      for (md in seq_along(vMDcompare)) {

         sModel <- vMDcompare[[md]]
         logger::log_info("----- {sModel} ------")

         lMSE <- lConf <- lInclu <- vector("list", length = nSimulation)
         for (i in 1:nSimulation) {

            logger::log_info("Simulation: {i}/{nSimulation}")

            tblBeta <-
               arrRawOut %>%
                  dplyr::filter(parameter == "beta" &
                                   model == sModel &
                                   simulation == i) %>%
                  dplyr::collect()


            ## MSE
            lMSE[[i]] <- md_check_beta_mse(tblBeta, tblBetaTrue)

            ## TP/FP/TN/FN/MCC
            tblConfInclusion <-
               tblBeta %>%
                  dplyr::group_by(model, simulation) %>%
                  tidyr::nest() %>%
                  ## data prep
                  dplyr::mutate(data_prep = purrr::map(data, function(x) {
                     x %>%
                        tidyr::pivot_wider(names_from = "key", values_from = "value") %>%
                        dplyr::select(-run, -chain, -parameter) %>%
                        as.matrix()
                  })) %>%
                  ## confusion matrix
                  dplyr::mutate(conf = purrr::map(data_prep, function(x, bolTrue) {
                     x %>% md_check_beta_create_confusion_matrix_foo(., bolTrue)
                  }, bolTrue)) %>%
                  ## true/false per each beta coefficient
                  dplyr::mutate(inclusion = purrr::map(data_prep, function(x) {
                     k <- as.data.frame(apply(x, 2, function(x) {findInterval(0, quantile(x, probs = c(0.05, 0.95))) == 0L }))
                     colnames(k) <- "included"
                     k$beta <- 1:nrow(k)
                     k
                  })) %>%
                  dplyr::ungroup()

            lConf[[i]] <-
               tblConfInclusion %>%
                  dplyr::select(model, simulation, conf) %>%
                  tidyr::unnest(cols = c(conf))

            lInclu[[i]] <-
               tblConfInclusion %>%
                  dplyr::select(model, simulation, inclusion) %>%
                  tidyr::unnest(cols = c(inclusion))

            rm(tblBeta)

         }

         lMSEAll[[md]] <- dplyr::bind_rows(lMSE)
         lConfAll[[md]] <- dplyr::bind_rows(lConf)
         lIncluAll[[md]] <- dplyr::bind_rows(lInclu)


      }


      tblMSE <- dplyr::bind_rows(lMSEAll)
      tblConf <- dplyr::bind_rows(lConfAll)
      tblInclu <- dplyr::bind_rows(lIncluAll)

      saveRDS(tblMSE, file = here::here("inst", "simulation", "output", "05_model_checking", "tblMSE.rds"))
      saveRDS(tblConf, file = here::here("inst", "simulation", "output", "05_model_checking", "tblConf.rds"))
      saveRDS(tblInclu, file = here::here("inst", "simulation", "output", "05_model_checking", "tblInclu.rds"))



      tblMSE <- readRDS(file = here::here("inst", "simulation", "output", "05_model_checking", "tblMSE.rds"))
      tblConf <- readRDS(file = here::here("inst", "simulation", "output", "05_model_checking", "tblConf.rds"))
      tblInclu <- readRDS(file = here::here("inst", "simulation", "output", "05_model_checking", "tblInclu.rds"))


      tblMCC <- md_check_beta_true_false_mcc(tblConf)



      # tblBeta <-
      #    arrRawOut %>%
      #          dplyr::filter(parameter == "beta") %>%  #model == "horseshoe_group_vector") %>%
      #          dplyr::collect()
      #
      #
      # tblConf <-
      #       tblBeta %>%
      #             dplyr::group_by(model, simulation) %>%
      #             tidyr::nest() %>%
      #             dplyr::mutate(data = purrr::map(data, function(x, bolTrue) {
      #                   x %>%
      #                         tidyr::pivot_wider(names_from = "key", values_from = "value") %>%
      #                         dplyr::select(-run, -chain, -parameter) %>%
      #                         as.matrix() %>%
      #                         md_check_beta_create_confusion_matrix_foo(., bolTrue)
      #             }, bolTrue)) %>%
      #             tidyr::unnest(cols = c(data)) %>%
      #             dplyr::ungroup()
      #
      # tblConf %>%
      #    tidyr::pivot_longer(names_to = "")


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
      # tblMSE <- md_check_beta_mse(tblBeta, tblBetaTrue)

      x <- tblMSE
      for (i in 1:9) {
         x <- rbind(x, tblMSE)
      }

      x <-
      x %>%
         dplyr::arrange(model)
      x$simulation <- rep(1:100, 3)


      ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
      ### MSE                                                                     ####

      gpMSE <- tblMSE %>% md_check_posterior_mse_beta_plot()

      sMSEPath <- here::here("inst/simulation/output/04_plots/02_posterior_beta/01_mse/beta_mse.pdf")
      ggplot2::ggsave(sMSEPath, plot = gpMSE, width = 9, height = 6)
      file.copy(from = sMSEPath,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/04_mse", overwrite = TRUE)


      gpMSEBoxplot <- tblMSE %>% md_check_posterior_mse_boxplot()

      sMSEBoxPath <- here::here("inst/simulation/output/04_plots/02_posterior_beta/01_mse/beta_mse_boxplot.pdf")
      ggplot2::ggsave(sMSEBoxPath, plot = gpMSEBoxplot, width = 9, height = 6)
      file.copy(from = sMSEBoxPath,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/04_mse", overwrite = TRUE)


      ## MCC
      tblMCC <- md_check_beta_true_false_mcc(tblConf)

      tblMCC %>%
         foo_simulation_as_factor() %>%
            ggplot2::ggplot(., ggplot2::aes(x = simulation, y = mcc)) +
               ggplot2::geom_bar(position = "stack", stat = "identity") +
               ggplot2::facet_wrap(~model, nrow = 3)


      md_check_posterior_mcc_plot(tblMCC, sStat = "tpr")

      ## Violine Plot of MCC/TPR/FPR


      ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
      ### Confusion Violine                                                       ####

      vVioline <- c("mcc", "tpr", "fpr")
      for (i in seq_along(vVioline)) {

         sViolinePath <- here::here("inst/simulation/output/04_plots/02_posterior_beta/02_confusion_violine",
                                    paste0(vVioline[[i]], "_violine.pdf"))

         gpVioline <- md_check_boxplot_tp_fp_mcc(tblMCC, sStat = vVioline[[i]])

         ggplot2::ggsave(sViolinePath, plot = gpVioline, width = 9, height = 6)
         file.copy(from = sViolinePath,
                   to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/05_confusion_violine", overwrite = TRUE)

      }

      md_check_boxplot_tp_fp_mcc(tblMCC, sStat = "mcc")
      md_check_boxplot_tp_fp_mcc(tblMCC, sStat = "tpr")
      md_check_boxplot_tp_fp_mcc(tblMCC, sStat = "fpr")

      ## Inclusion Plot
      tbl1 <-
         tblInclu %>%
            # dplyr::filter(model == "group_lasso_hierarchical") %>%
            dplyr::mutate(beta = factor(beta, levels = sort(unique(beta), decreasing = TRUE)))


      tbl1 <-
         tblInclu %>%
         dplyr::filter(model == "group_lasso_hierarchical") %>%
         dplyr::mutate(beta = factor(beta, levels = sort(unique(beta), decreasing = TRUE))) %>%
         foo_simulation_as_factor()

      # x <- tbl1
      # for (i in 1:9) {
      #    x <- rbind(x, tbl1)
      # }
      #
      # x <-
      #    x %>%
      #    dplyr::arrange(simulation)
      # x$simulation <- sort(rep(1:100, 50))

      tbl1 %>%
         ggplot2::ggplot(ggplot2::aes(simulation, beta, fill = included)) +
            ggplot2::geom_tile(color = "lightgrey",
                      lwd = .5,
                      linetype = 1) +
            ggplot2::scale_fill_manual(values = c("TRUE" = "#67a9cf", "FALSE" = "white"),
                                       name = "") +
            ggplot2::coord_fixed() +
            # ggplot2::facet_wrap(~model, nrow = 3) +
            theme_custom_thesis() +
            ggplot2::theme(legend.position = "none")



      lHeatmap <- vector("list", length = length(vMDcompare))
      for (i in seq_along(vMDcompare)) {

         tblHeat <-
            tblInclu %>%
               dplyr::filter(model == vMDcompare[[i]])

         gpHeat <- md_check_posterior_heatmap_plot(tblHeat)

         sHeatPath <- here::here("inst/simulation/output/04_plots/02_posterior_beta/03_heatmap",
                                                          paste0(vMDcompare[[i]], "_heatmap.pdf"))

         ggplot2::ggsave(sHeatPath, plot = gpHeat, width = 9, height = 6)
         file.copy(from = sHeatPath,
                   to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/06_beta_heatmap", overwrite = TRUE)
      }

      #   ____________________________________________________________________________
      #   Table of Posterior Statistics                                           ####

      tblMean <-
         dplyr::left_join(tblMSE, tblMCC, by = c("model", "simulation")) %>%
            dplyr::group_by(model) %>%
            dplyr::summarise(
               MSE = mean(mse),
               VAR = mean(mse_var),
               BIAS = mean(mse_bias),
               MCC = mean(mcc),
               TPR = mean(tpr),
               FNR = mean(fpr)
            ) %>%
            dplyr::mutate_if(is.numeric, . %>% round(., 3)) %>%
            dplyr::mutate(model = stringr::str_to_title(gsub("_", " ", model)))

      Hmisc::latex(tblMean, cdec = c(0, 3, 3, 3, 3, 3, 3), na.blank = TRUE,
                   booktabs = TRUE, table.env=FALSE, center = "center", file="", title="", rowname=NULL,
                   caption = "Simulation: MSE Estimation and Selection Accuracy")

}
