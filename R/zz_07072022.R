#' Title
#'
#' This function is just a wrapper so the code does not get triggered when package is load
#'
#' @return
#' @export
#'
#' @examples
fooSimulation <- function() {

      #   ____________________________________________________________________________
      #   Set Inputs                                                              ####

      library(dplyr)

      ## simulations
      nSimulation <- 10
      vSeed <- 1001:(1001 + nSimulation)


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
      nGroupSize <- rep(nP + 1, lBlock[[1]]) # group Index

      # dgp_exp_almon_lag(nLag = vLag, 0.005, -0.05)

      ## set seed
      nSeed <- 1001:(1001 + nSimulation)

      ## direction
      sDir <- here::here("inst/data/newData_01022022/")

      nSim <- dir(paste0(sDir, "input"))
      nSim <- nSim[1:10]

      #   ____________________________________________________________________________
      #   Creating Data                                                           ####

      for (i in 1:nSimulation) {

            logger::log_info("Run: {i}")
            x1 <- dgp_create_full(nY = nY, vBeta = vBeta,
                                  nK = nK, nFreq = nFreq, nLag = nLag,
                                  nMu = nMu, lRho = lRho,
                                  nVar = nVar, nWithin = nWithin, nBetween = nBetween,
                                  nT1 = nT1, nT2 = nT2,
                                  .sSeed = nSeed[[i]])

            # x2 <- create_predictor_lag_poly(x1[["x_align"]],
            #                                 vLag = vLag, vP = vP,
            #                                 .sPolyMatrix = "almon", .bEndpoint = FALSE)

            lModelInput <- create_model_input(maY = x1[["y"]], dfDataX = x2[[1]],
                                              nTrain = nTrain, nG = lBlock[[1]],
                                              nGroupSize = nGroupSize)

            sDataSave <- list(
                  "x_raw" = x1$x_raw[[1]],
                  "x_align" = x1$x_align[[1]],
                  "y" = x1$y,
                  "model_data" = lModelInput)

            # saveRDS(sDataSave,
            #         paste0(sDir, "/input/simulatedData_", i, "_seed_", nSeed[[i]], ".rds"))
      }




      #   ____________________________________________________________________________
      #   Create Models (if necessary)                                            ####

      lCMDmodels <- vector("list", length(stanmodels))
      dirStan <- here::here("inst/stan")
      sStanFiles <- dir(dirStan)[grep("stan", dir(dirStan))]

      for (i in seq_along(sStanFiles)) {
            sStanModel <- glue::glue("{dirStan}/{sStanFiles[i]}")
            lCMDmodels[[i]] <- cmdstanr::cmdstan_model(stan_file = sStanModel)
      }

      sModelName <- stringr::str_to_title(gsub("_", " ", stringr::str_extract(sStanFiles, "[^\\.]+")))

      names(lCMDmodels) <- sModelName



      #   ____________________________________________________________________________
      #   Fit Models                                                              ####

      # TODO: 30.10.2020 stopped here...
      # * for each simulation create full output
      #  - directory structure
      #  - stan object
      #  - extract parquet
      #  - summary parquet

      start_time <- Sys.time()
      lTimes <- vector("list", length(sModelName))
      for (i in 1:nSimulation) {

            logger::log_info("----------------- Run: {i}/{length(nSim)} -----------------")

            dirSim_input <- paste0(sDir, "input")
            lData <- readRDS(paste0(dirSim_input, "/", dir(dirSim_input)[[i]]))
            lData <- lData$model_data

            for (j in seq_along(lCMDmodels)) {
                  sModel <- names(lCMDmodels)[[j]]
                  logger::log_info("{sModel} (run {i})")

                  sDirStanModel <- paste0(sDir, "output/",
                                          tolower(gsub(" ", "_", sModel)), "/",
                                          tolower(gsub(" ", "_", sModel)), "_", i, ".rds")

                  #if (file.exists(sDirStanModel)) {

                  #logger::log_info("File already exists - skipped to next...")
                  #next()

                  #} else {

                  start_time_sampling <- Sys.time()


                  md_stan <- lCMDmodels[[j]]

                  lStanObj <- md_stan$sample(
                        data = lData,
                        seed = 123,
                        chains = 4,
                        parallel_chains = 4,
                        iter_sampling = 4000, iter_warmup = 4000,
                        save_warmup = FALSE,
                        adapt_delta = 0.99, max_treedepth = 10
                  )

                  stanfit <- rstan::read_stan_csv(lStanObj$output_files())
                  saveRDS(stanfit, sDirStanModel)

                  rm(lStanObj)
                  rm(stanfit)

                  end_time_sampling <- Sys.time()

                  lTimes[[j]] <- end_time_sampling - start_time_sampling

                  #}

                  # lStanObj$save_object(file = paste0(here::here("inst/data/simulation_output/"),
                  #                                    tolower(gsub(" ", "_", sModel)), "/",
                  #                                    tolower(gsub(" ", "_", sModel)), "_", i, ".rds"))


            }
            rm(lData)
            Sys.sleep(5)
      }





      #   ____________________________________________________________________________
      #   Diagnostic Preparation                                                  ####

      ## all true Y values
      lYTrue_raw <- vector("list", length(nSim))
      for (i in seq_along(nSim)) {
            s1 <- readRDS(paste0(sDir, "input/", nSim[[i]]))
            lYTrue_raw[[i]] <- c(s1$y_train, s1$y_test)
      }




      ## combine data from simulation for each model
      lCombined <- vector("list", length(sModelName))


      for (i in seq_along(sModelName)) {
            # for (i in seq_along(sModelName[1:2])) {
            # i <- 1
            print(sModelName[[i]])

            sDirOutput <- paste0(sDir, "output/")
            sFileName <- dir(paste0(sDirOutput, tolower(gsub(" ", "_", sModelName[[i]]))))

            comb_fit_model <- lapply(paste0(sDirOutput, tolower(gsub(" ", "_", sModelName[[i]])), "/", sFileName), readRDS)

            nSimRegex <- as.numeric(unlist(stringr::str_extract_all(sFileName, "\\(?[0-9]+\\)?")))
            lYTrue_sort <- lYTrue_raw[nSimRegex]
            lYTrue_sort <- lapply(lYTrue_sort, function(x, nSplit) {list("train" = x[1:nSplit], "test" = x[(nSplit+1): length(x)])},
                                  nSplit = nTrain)
            lYTrue_sort <- purrr::transpose(lYTrue_sort)

            lCombined[[i]] <- moc_data_prep_model(lStanObj = comb_fit_model, lYTrue = lYTrue_sort,
                                                  nGroupSize = nGroupSize, vLag = vLag, vP = vP,
                                                  .sPolyMatrix = "almon", .bEndpoint = FALSE, bolTrue = bolTrue)
      }

      lCombined <- purrr::compact(lCombined)

      ## getting Beta Coef from Theta
      dfTheta <-
            lCombined[[4]]$theta %>%
            dplyr::filter(chain == "chain:1") %>%
            dplyr::select(run, parameter, value) %>%
            tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
            dplyr::select(-run)


      ## Group Index
      vIndGroup <- c(0, cumsum(nGroupSize))

      ## Lag Matrix
      mQ <- prep_model_lag_matrix(vP, vLag, .bEndpoint = FALSE)

      maTheta <- as.matrix(dfTheta)
      lBeta <- vector("list", length(nGroupSize))
      lWeigthing <- vector("list", length(nGroupSize))
      for (i in 1:(length(vIndGroup) - 1)) {
            sWPrep <- maTheta[, (vIndGroup[[i]] + 1):vIndGroup[[i + 1]]] %*% mQ
            lWeigthing[[i]] <- sWPrep / rowSums(sWPrep)
            lBeta[[i]] <- rowSums(sWPrep)
      }

      maBeta <- do.call("cbind", lBeta)
      maWeight <- do.call("cbind", lWeigthing)

      apply(maBeta, 2, mean)[1:9]
      plot(apply(maBeta, 2, mean))
      abline(h = 0)
      # 1.2, -0.5, 3, 1.1, 0.5, -2.5, -1.5, 0.8, 3.5

      ## Weighting Check
      lWOutPlot <- vector("list", length(lWeigthing))
      for (i in seq_along(lWeigthing)) {

            sM <- apply(lWeigthing[[i]], 2, mean)
            sD <- apply(lWeigthing[[i]], 2, sd)

            s1 <- sM
            s2 <- sM + sD
            s3 <- sM - sD

            lWOutPlot[[i]] <- rbind(s1, s2, s3)
      }


      #   ____________________________________________________________________________
      #   Weighting Plot                                                          ####

      tblThetaPlot <-
            purrr::map2(lWeigthing, 1:length(lWeigthing),
                        function(x, y) {
                              tibble::tibble(
                                    "theta" = y,
                                    "para" = 1:length(apply(x, 2, mean)),
                                    "mean" = apply(x, 2, mean),
                                    "sdUp" = apply(x, 2, mean) + apply(x, 2, sd),
                                    "sdDown" = apply(x, 2, mean) - apply(x, 2, sd)
                              )
                        }) %>%
            dplyr::bind_rows() %>%
            tidyr::pivot_longer(-c(1, 2), names_to = "key", values_to = "value")


      lGPlot <- vector("list", 9)
      for (i in 1:9) {
            lGPlot[[i]] <-
                  tblThetaPlot %>%
                  dplyr::filter(theta == i) %>%
                  ggplot2::ggplot(ggplot2::aes(x = para, y = value, group = key)) +
                  ggplot2::geom_line()
      }

      # > apply(maBeta[, 1:9], 2, mean)
      # 0.93082954 -0.17315322  2.93918057  0.79205462  0.05328455 -1.90212491 -1.01271238  0.56329044  3.31076842
      # 1.04221105 -0.17941081  3.01592710  0.84181671  0.05302292 -1.99368471 -1.09861901  0.85483801  3.44890155

      tblThetaPlot %>%
            dplyr::filter(theta %in% c(1:9)) %>%
            ggplot2::ggplot(ggplot2::aes(x = para, y = value, group = key)) +
            ggplot2::geom_line() +
            ggplot2::facet_wrap(~theta)

      ## should look like that
      yExpAlmon <- dgp_exp_almon_lag(nLag = vLag, nT1, nT2)
      plot(1:vLag, yExpAlmon, type = "l")

      ## End of getting Beta Coef




      #   ____________________________________________________________________________
      #   UMIDAS                                                                  ####

      xfull <- do.call("cbind", x1$x_align[[1]])

      regMod1 <- lm(x1$y ~ xfull - 1)
      modCoef1 <- regMod1$coefficients
      modCoef1 <- matrix(modCoef1, nrow = 6, byrow = F)

      modBeta1 <- rowSums(t(modCoef1))
      modWeight1 <- t(modCoef1) / rowSums(t(modCoef1))




      #   ____________________________________________________________________________
      #   MIDAS Package                                                           ####

      library(midasr)
      xMidas <- x1$x_raw[[1]][1:600, 1]
      for (i in 1:vK) {
            assign(x = paste0("xmidas_", i), value = x1$x_raw[[1]][1:600, i])
      }

      # eq_r <- midas_r(x1$y ~ mls(xMidas, 0:7, 3, nealmon) - 1, start = list(xMidas = c(1, -0.5)))

      sFullMidas <- paste0("eq_r <- midas_r(x1$y ~ ", paste0("mls(", paste0("xmidas_", 1:vK), ", 0:7, 3, nealmon)", collapse = " + "), " - 1, start = list(", paste0(paste0("xmidas_", 1:vK), " = c(1, -0.5)", collapse = ", "), "))")

      eval(parse(text = sFullMidas))

      summary(eq_r)
      maNealmon <- matrix(nrow = 50, ncol = 6)
      for (i in 1:50) {
            maNealmon[i, ] <- nealmon(p = c(eq_r$coefficients[2 * i - 1], eq_r$coefficients[2 * i]), 6)
      }
      maNealmon <- maNealmon

      maNealBeta <- rowSums(maNealmon)
      maNealWeight <- maNealmon / maNealBeta

      rowSums(maNealWeight)



      #   ____________________________________________________________________________
      #   UMIDAS Lasso                                                            ####
      library(lars)

      modLasso <- lars::lars(x = xfull, y = x1$y, normalize = FALSE, intercept = FALSE, type = "lasso")

      maLasso <- matrix(coef(modLasso)[200, ], ncol = 6, byrow = TRUE)

      modBeta1 <- rowSums(maLasso)
      modWeight1 <- maLasso / rowSums(maLasso)

      #   ____________________________________________________________________________
      #   Distribution of Theta                                                   ####



      rm(comb_fit_model)
      # names(lCombined) <- sModelName
      names(lCombined) <- sModelName[1:2]

      lOut <- moc_data_combine_all(lCombined)


      ## connection to files
      dirDGP <- here::here("inst/data/01_dgp/")

      dirExtract <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_extract/")
      l1 <- lapply(dirExtract, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
      arrExtract <- arrow::open_dataset(l1)


      dirSummary <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_summary/")
      l1 <- lapply(dirSummary, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
      arrSummary <- arrow::open_dataset(l1)

      ## load data
      tblTheta <-
            arrExtract %>%
            dplyr::filter(parameter_1 == "theta") %>%
            dplyr::collect() %>%
            dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
            dplyr::select(-parameter_1, -parameter_2)

      tblBeta <-
            arrExtract %>%
            dplyr::filter(parameter_1 == "beta") %>%
            dplyr::collect() %>%
            dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
            dplyr::select(-parameter_1, -parameter_2)

      moc_create_confusion_matrix(tblBeta, bolTrue)

      ## all true Y values
      fSims <- dir(dirDGP)
      lYTrue_raw <- vector("list", length(fSims))
      for (i in seq_along(nSim)) {
            sFile <- dir(paste0(dirDGP, "/", fSims[[i]], "/raw_input"))
            s1 <- readRDS(paste0(dirDGP, "/", fSims[[i]], "/raw_input/", sFile))
            lYTrue_raw[[i]] <- data.frame(simulation = i,
                                          key = seq(1, s1$nY_train +  s1$nY_test),
                                          value = c(s1$y_train, s1$y_test))
      }
      tblYTrue <- dplyr::bind_rows(lYTrue_raw) %>% dplyr::as_tibble()


      tblYPred <-
            arrExtract %>%
            dplyr::filter(parameter_1 == "y_pred") %>%
            dplyr::collect() %>%
            dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
            dplyr::select(-parameter_1, parameter_2)


      tblYRep <-
            arrExtract %>%
            dplyr::filter(parameter_1 == "y_rep", simulation == 1) %>%
            dplyr::collect() %>%
            dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
            dplyr::select(-parameter_1, parameter_2)

      tbl1 <-
            arrExtract %>%
            dplyr::filter(parameter_1 == "y_rep", simulation == 1) %>%
            dplyr::collect() %>%
            dplyr::group_by(model, parameter_1) %>%
            dplyr::filter(value == min(value)) %>%
            dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
            dplyr::select(-parameter_1, parameter_2)



      #   ____________________________________________________________________________
      #   Model Checking                                                          ####

      ## Note: maybe only read in only every 4 iteration

      ## diagnose
      ## TODO:
      ## * min/max of rhat and neff for each parameter and model simulation
      ## * line plot with % divergence per simulation and group by model
      ## * test statistic plot from bayesplot package
      ## * test stat plot for with simulation on x and p-value on y, group by model
      moc_diagnose_plot(lOut$rhat)


      ## beta coeffsys
      moc_tp_fp_plot(lOut$confusion, .sType = "TP")
      moc_tp_fp_table(lOut$beta, bolTrue)
      moc_beta_distr_plot(lOut$beta, bolTrue, vTrueValue)
      moc_mse_table(lOut$beta, vBeta)
      moc_matthes_corr(lOut$confusion)

      ## y rep
      moc_yrep_plot(lOut$yrep, "mean", nY)
      moc_yrep_plot_per_simulation(lOut$yrep, 2)
      moc_yrep_95_covered(lOut$yrep)


      ## Table of overall fit with p values for yrep
      tbl1 %>%
            dplyr::group_by(model, iteration) %>%
            dplyr::summarise(min = min(value), max = max(value), mean = mean(value))

      ## Visualize Fit for each simulation
      tbl1 <-
            tblYRep %>%
            dplyr::filter(simulation == 1) %>%
            dplyr::select(iteration, value, model, parameter) %>%
            split(., tblYRep$model) %>%
            purrr::map(., function(x) {
                  x %>%
                        dplyr::select(-model) %>%
                        tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
                        dplyr::select(-iteration) %>%
                        as.matrix() %>%
                        unname()
            })

      vGroup <- rep(names(tbl1), purrr::map_dbl(tbl1, ncol))

      tbl1 <- tbl1 %>% purrr::reduce(., dplyr::bind_cols) %>% as.matrix() %>% unname()

      s1 <- tblYTrain %>% dplyr::filter(simulation == 1)
      s2 <- rep(s1$value, 5)

      ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "min")
      q25 <- function(y) quantile(y, 0.25)
      ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "min")
      ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "max")
      ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "mean")
      ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "sd")

      q25 <- function(y) quantile(y, 0.25)
      ppc_stat_grouped(s2, tbl1, group = vGroup, stat = "q25")


      ppc_ribbon_grouped(s2, tbl1, group = vGroup)

      ppc_intervals_grouped(s2, tbl1, group = vGroup)
      ppc_scatter_avg_grouped(s2, tbl1, group = vGroup)
      ## log_lik

}




foo_fit_model <- function(sName, lStanModel, lData) {
      tmp_envir <- new.env(parent = baseenv())
      tmp_envir$model <- lStanModel[[sName]]
      print(tmp_envir$model)
      tmp_envir$fit <- rstan::sampling(tmp_envir$model, data = lData)

      tmp_envir$fit
}


foo_allModel_combi <- function(x, y) {
      lOut <- list(x)
      names(lOut) <- y
      lOut
}


create_save_model_input <- function(lInput) {

      for (i in seq_along(lInput)) {
            saveRDS(lInput[[i]], paste0(here::here("inst/data/input/"), i, ".rds"))
      }


}
