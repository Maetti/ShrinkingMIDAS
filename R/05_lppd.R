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
      #   Load lppd                                                               ####


      arrRawOut <- arrow::open_dataset(here::here("inst", "simulation",
                                                  "output", "02_raw_extracted"),
                                       partitioning = c("model", "simulation"))


      dfObs <-
         data.frame(
            "run" = sort(rep(1:4000, 4)),
            "chain" = rep(1:4, 4000),
            "obs" = 1:16000
         )

      lOut <- vector("list", length = nSimulation)

      for (i in 1:nSimulation) {

         logger::log_info("----- simulation: {i}/{nSimulation} -----")

         ## load simulated data
         logger::log_info("load simulated data")
         nSimRunSave <- helper_create_number_name(i)

         sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))
         lData <- readRDS(sDirInputLoad)

         ## get mean and sd for lppd
         logger::log_info("load theta...")
         tblTheta <-
            arrRawOut %>%
               dplyr::filter(simulation == i, parameter == "theta") %>%
               dplyr::collect() %>%
               dplyr::mutate(run = as.numeric(run),
                             chain = as.numeric(chain))


         tblTheta <- dplyr::left_join(tblTheta, dfObs, by = c("run", "chain"))


         logger::log_info("load sigma...")
         tblSigma <-
            arrRawOut %>%
               dplyr::filter(simulation == 1, parameter == "sigma") %>%
               dplyr::collect() %>%
               dplyr::mutate(run = as.numeric(run),
                             chain = as.numeric(chain)) %>%
               dplyr::arrange(key, model, run, chain)


         tblSigma <- dplyr::left_join(tblSigma, dfObs, by = c("run", "chain"))

         nYTest <- lData$model_data$y_test
         maXTest <- lData$model_data$x_test
         sModelNamesRaw <- tblTheta$model %>% unique()

         dfOut <-
            data.frame(
               "model" = sModelNamesRaw,
               "value" = rep(0, length(sModelNamesRaw))
            )

         for (j in seq_along(sModelNamesRaw)) {

            logger::log_info("model: {sModelNamesRaw[[j]]}")
            tblModel <-
               tblTheta %>%
                  dplyr::filter(model == sModelNamesRaw[[j]]) %>%
                  dplyr::select(key, value, obs) %>%
                  tidyr::pivot_wider(values_from = value, names_from = obs)

            maTheta <- as.matrix(tblModel[, -1])

            maMean <- maXTest %*% as.matrix(maTheta[, -1])


            vSigma <-
               tblSigma %>%
                  dplyr::filter(model == sModelNamesRaw[[j]]) %>%
                  dplyr::select(value) %>%
                  unlist()


            vLogCalc <- vector("numeric", length = length(nYTest))

            for (k in seq_along(nYTest)) {
               y_i <- nYTest[[k]]
               vLogCalc[[k]] <- log(mean(dnorm(y_i, maMean[k, ], sd = vSigma)))
            }

            dfOut[dfOut$model == sModelNamesRaw[[j]], "value"] <- sum(vLogCalc)

         }

         logger::log_info("append data")
         dfOut$simulation <- i
         lOut[[i]] <- dfOut


      }


      tblLPPD <- dplyr::bind_rows(lOut)

      saveRDS(tblLPPD, file = here::here("inst", "simulation", "output", "05_model_checking", "tblLPPD.rds"))



   #   ____________________________________________________________________________
   #   LPPD                                                                    ####

      tblLPPD <- readRDS(here::here("inst", "simulation", "output", "05_model_checking", "tblLPPD.rds"))
      tblLPPD <-
         dplyr::as_tibble(tblLPPD) %>%
         dplyr::mutate(model = dplyr::case_when(
            model == "group_lasso_hierarchical" ~ "Group Lasso Hierarchical",
            model == "horseshoe_group" ~ "Horseshoe Group",
            model == "horseshoe_group_plus" ~ "Horseshoe Group Plus"
         ))



      ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
      ### Table                                                                   ####

      tblTable <-
         tblLPPD %>%
         dplyr::group_by(model) %>%
         dplyr::summarise(
            Mean = mean(value),
            SD = sd(value),
            Min = min(value),
            Max = max(value),
            IQR = IQR(value)
         ) %>%
         dplyr::mutate_if(is.numeric, ~round(., 1))


      Hmisc::latex(tblTable, cdec = c(0, 1, 1, 1, 1, 1), na.blank = TRUE,
                   booktabs = TRUE, table.env=FALSE, center = "center", file="", title="", rowname=NULL,
                   caption = "Summary log pointwise predictive densitiy")

   ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
   ### Heatmap                                                                 ####

      maLPPD <-
         tblLPPD %>%
            tidyr::pivot_wider(names_from = "model", values_from = "value")
      maLPPD <- maLPPD[, -1]

      lCompMatrix <- vector("list", length = ncol(maLPPD))
      for (i in 1:ncol(maLPPD)) {

         vComp <- maLPPD[, i]
         lCompMatrix[[i]] <- apply(apply(maLPPD, 2, function(x) {vComp >= x}), 2, sum)

      }

      tblComp <- dplyr::bind_rows(lCompMatrix)
      tblComp <- as.matrix(tblComp)
      diag(tblComp) <- NA_real_
      tblComp <- dplyr::as_tibble(tblComp)

      tblComp$model <- colnames(tblComp)

      tblComp <-
         tblComp %>%
         dplyr::select(model, dplyr::everything()) %>%
         # dplyr::mutate_if(is.numeric, ~as.character(.)) %>%
         force()

      tblPlot <-
         tblComp %>%
         tidyr::pivot_longer(names_to = "model_comp", values_to = "value", -1)

      tblPlot$model <- factor(tblPlot$model, levels = rev(colnames(maLPPD)))
      tblPlot$model_comp <- factor(tblPlot$model_comp, levels = colnames(maLPPD))


      tblPlotDiag <-
         tblPlot %>%
         dplyr::filter(model == model_comp)

      gpHeatLPPD <-
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


      sForecastHeat <- here::here("inst/simulation/output/04_plots/04_lppd/lppd_heatmap.pdf")

      ggplot2::ggsave(sForecastHeat, plot = gpHeatLPPD, width = 9, height = 6)
      file.copy(from = sForecastHeat,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/09_lppd/lppd_heatmap.pdf",
                overwrite = TRUE)



   ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
   ### Boxplot                                                                 ####

      gpBoxLPPD <-
         tblLPPD %>%
            ggplot2::ggplot(aes(x = model, y = value, fill = model)) +
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

      sLPPDBoxplot <- here::here("inst/simulation/output/04_plots/04_lppd/lppd_boxplot.pdf")

      ggplot2::ggsave(sLPPDBoxplot, plot = gpBoxLPPD, width = 9, height = 6)
      file.copy(from = sLPPDBoxplot,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/09_lppd/lppd_boxplot.pdf",
                overwrite = TRUE)

   ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
   ### Linechart                                                               ####

      gpLineLPPD <-
         tblLPPD %>%
            ggplot2::ggplot(ggplot2::aes(x = simulation, y = value, color = model)) +
            ggplot2::geom_line() +
            ggplot2::scale_color_manual(values = rev(RColorBrewer::brewer.pal(6, "Dark2"))) +  # Fill color
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = "Simulation", y = "", color = "Model") +
            theme_custom_thesis()


      sLPPDLine <- here::here("inst/simulation/output/04_plots/04_lppd/lppd_line.pdf")

      ggplot2::ggsave(sLPPDLine, plot = gpLineLPPD, width = 9, height = 6)
      file.copy(from = sLPPDLine,
                to = "/home/matthias/Schreibtisch/SoSe 21/magister arbeit/paper/plot/09_lppd/lppd_line.pdf",
                overwrite = TRUE)


}

