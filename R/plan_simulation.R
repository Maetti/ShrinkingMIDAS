#' Simulation Plan
#' This function contains the drake plan for the Monte-Carlo Simulation
#'
#' @return
#' @export
#'
#' @examples
get_simulation_plan <- function() {


      library(dplyr)
      # library(rstan)
      library(cmdstanr)
      library(bayesplot)
      parallel::detectCores()
      options(mc.cores = 4)

      ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
      ### Input specification                                                     ####
      # vBeta <- c(0, 0.3, 0.5, 0, 0.3, 0.5, 0, 0, 0.8, rep(0, 21))
      vBeta <- c(0, 0.4, 1.2, 0, 0, 0.8, 1, 0, 0, 0)
      bolTrue <- vBeta != 0
      vTrueValue <- vBeta[bolTrue]
      nY <- 200
      nTrain <- 150
      vK = c(10)
      lBlock = list(c(10))
      vFreq = 3
      vLag = 6
      nBlockNames = "m"
      nMu = 0.1
      lRho = list(c(0.5))
      nVar = 1
      nWithin = 0.5
      nBetween = 0
      nT1 = 0.0005
      nT2 = -0.00007


      ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
      ### Model Input                                                             ####

      nP <- 3
      nGroupSize <- rep(nP + 1, lBlock[[1]]) # group Index


      # lStanModel <- list(
      #    stanmodels$adaptive_group_lasso,
      #    stanmodels$adaptive_group_lasso_ss,
      #    stanmodels$group_lasso_hierarchical
      # )

      lStanModel <- vector("list", length(stanmodels))
      for (i in seq_along(stanmodels)) {
      # for (i in 1:2) {
         lStanModel[[i]] <- stanmodels[[i]]
      }

      ## fitting rstan models
      nSim <- dir(here::here("inst/data/input/"))
      nSim <- nSim[1:10]
      # lStanModel[[1]] <- stanmodels$adaptive_group_lasso
      # lStanModel[[2]] <- stanmodels$horseshoe_group
      # lStanModel[[3]] <- stanmodels$adaptive_group_lasso_ss

      # lStanModel <- purrr::compact(lStanModel)
      # sModelName <- c("Group Lasso", "Horseshoe")
      # names(lStanModel) <- sModelName
      names(lStanModel) <- stringr::str_to_title(gsub("_", " ", names(stanmodels)))
      # names(lStanModel)[2] <- "Adaptive Group Lasso SS"
      sModelName <- names(lStanModel)
      ##
      # drake::readd(lSimData)
      # drake::readd(lPredictor)
      # drake::readd(fit_model)
      # drake::readd(comb_summary)

      # lStanObj <- drake::readd(comb_summary_group)
      # l2 <- drake::readd(comb_summary_group)
      # lOut <- drake::readd(allModelSplit)

      kSeed <- 1000:1050

      p1 <- drake::drake_plan(
            ## Monte Carlo Simulation
            # vSeed = drake::target(seq(1000, 1001, 1)),
            lSimData = drake::target(dgp_create_full(nY = nY, vBeta = vBeta,
                                                     vK = vK, lBlock = lBlock, vFreq = vFreq, vLag = vLag,
                                                     nBlockNames = nBlockNames,
                                                     nMu = nMu, lRho = lRho,
                                                     nVar = nVar, nWithin = nWithin, nBetween = nBetween,
                                                     nT1 = nT1, nT2 = nT2,
                                                     .sSeed = nSeed),
                                     transform = map(nSeed = !!kSeed)
                              ),
            ## prepare input data for model estimation
            lPredictor = drake::target(create_predictor_lag_poly(lSimData[["x_align"]],
                                                                 vLag = vLag, vP = nP,
                                                                 .sPolyMatrix = "almon", .bEndpoint = FALSE),
                                       transform = map(lSimData)
                              ),
            lModelInput = drake::target(create_model_input(maY = lSimData[["y"]], dfDataX = lPredictor[[1]],
                                                           nTrain = nTrain,
                                                           nG = lBlock[[1]], nGroupSize = nGroupSize),
                                        transform = map(lSimData, lPredictor)
                              ),
            saveData = drake::target(saveRDS(lModelInput,
                                             file = !!paste0(here::here("inst/data/input/"), .id_chr, ".rds")),
                                     transform = map(lModelInput))
            # fit_model = drake::target(rstan::sampling(lStanModel, data = lModelInput),
            #                           # dynamic = map(lModelInput)
            #                           transform = cross(lStanModel, lModelInput),
            #                           # dynamic = cross(lStanModel, lModelInput)
            #                   ),
            # fit_model = drake::target(foo_fit_model(sName, lStanModel = lStanModel, lData = lModelInput),
            #                           transform = cross(sName = !!sModelName, lModelInput)
            #                      ),
            # ## combine results per model
            # comb_fit_model = drake::target(list(fit_model),
            #                              transform = combine(fit_model, .by = sName)),
            # check_model_data = drake::target(moc_data_prep_model(comb_fit_model,
            #                                                      nGroupSize = nGroupSize, vLag = vLag, vP = nP,
            #                                                      .sPolyMatrix = "almon", .bEndpoint = FALSE,
            #                                                      bolTrue = bolTrue),
            #                               transform = map(comb_fit_model)),
            # rename_model = drake::target(foo_allModel_combi(check_model_data, y),
            #                              transform = map(check_model_data, y = !!sModelName)),
            # allModel = drake::target(rename_model, transform = combine(rename_model)),
            # allModelSplit = drake::target(moc_data_combine_all(allModel))
            # combine_model = drake::target(rstan::extract(fit_model), dynamic = combine(fit_model, .by = stanmodels))
            # model_checking_plot = drake::target()
      )

      plot(p1)
      drake::vis_drake_graph(p1)

      drake::make(p1, garbage_collection = TRUE, memory_strategy = "autoclean")

      # for (i in seq_along(sModelName)) {
      #    dir.create(paste0(here::here("inst/data/simulation_output/"), tolower(gsub(" ", "_", sModelName[[i]]))))
      # }

      sLogFile <- here::here("inst/log", gsub("-", "_", Sys.Date()))
      logger::log_appender(logger::appender_tee(sLogFile))

      logger::log_info("Create Models")

      lCMDmodels <- vector("list", length(sModelName))

      dirStan <- here::here("inst/stan")
      sStanFiles <- dir(dirStan)[grep("stan", dir(dirStan))]

      for (i in seq_along(sStanFiles)) {
         sStanModel <- glue::glue("{dirStan}/{sStanFiles[i]}")
         lCMDmodels[[i]] <- cmdstanr::cmdstan_model(stan_file = sStanModel)
      }

      names(lCMDmodels) <- stringr::str_to_title(gsub("_", " ", stringr::str_extract(sStanFiles, "[^\\.]+")))

      sModelName <- names(lCMDmodels)

      nSim <- nSim[1:10]

      start_time <- Sys.time()
      for (i in seq_along(nSim)) {
         logger::log_info("----------------- Run: {i}/{length(nSim)} -----------------")
         lData <- readRDS(here::here("inst/data/input", nSim[[i]]))

         for (j in seq_along(lCMDmodels)) {
            sModel <- names(lCMDmodels)[[j]]
            logger::log_info("{sModel} (run {i})")

            sDirStanModel <- paste0(here::here("inst/data/simulation_output/"),
                                    tolower(gsub(" ", "_", sModel)), "/",
                                    tolower(gsub(" ", "_", sModel)), "_", i, ".rds")

            if (file.exists(sDirStanModel)) {

               logger::log_info("File already exists - skipped to next...")
               next()

            } else {

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

            }

            # lStanObj$save_object(file = paste0(here::here("inst/data/simulation_output/"),
            #                                    tolower(gsub(" ", "_", sModel)), "/",
            #                                    tolower(gsub(" ", "_", sModel)), "_", i, ".rds"))


         }
         rm(lData)
         Sys.sleep(5)
      }

      end_time <- Sys.time()
      end_time - start_time


      # for (i in seq_along(nSim)) {
      #    logger::log_info("----------------- Run: {i} -----------------")
      #    lData <- readRDS(here::here("inst/data/input", nSim[[i]]))
      #
      #    for (j in seq_along(sModelName)) {
      #       sModel <- sModelName[[j]]
      #       logger::log_info("{sModel} (run {i})")
      #
      #       lStanObj <- rstan::sampling(lStanModel[[sModel]], data = lData, iter = 4000)
      #       saveRDS(lStanObj, paste0(here::here("inst/data/simulation_output/"), tolower(gsub(" ", "_", sModel)), "/", tolower(gsub(" ", "_", sModel)), "_", i, ".rds"))
      #
      #       rm(lStanObj)
      #    }
      #    rm(lData)
      #    Sys.sleep(5)
      # }



   #   ____________________________________________________________________________
   #   Diagnostic Preparation                                                  ####

      ## all true Y values
      lYTrue_raw <- vector("list", length(nSim))
      for (i in seq_along(nSim)) {
         s1 <- readRDS(here::here("inst/data/input", nSim[[i]]))
         lYTrue_raw[[i]] <- c(s1$y_train, s1$y_test)
      }




      ## combine data from simulation for each model
      lCombined <- vector("list", length(sModelName))

      # for (i in seq_along(sModelName)) {
      for (i in seq_along(sModelName[1:2])) {
         print(sModelName[[i]])

         sFileName <- dir(paste0(here::here("inst/data/simulation_output/"),
                                 tolower(gsub(" ", "_", sModelName[[i]]))))

         comb_fit_model <- lapply(here::here("inst/data/simulation_output",
                                             tolower(gsub(" ", "_", sModelName[[i]])), sFileName),
                                  readRDS)

         nSimRegex <- as.numeric(unlist(stringr::str_extract_all(sFileName, "\\(?[0-9]+\\)?")))
         lYTrue_sort <- lYTrue_raw[nSimRegex]
         lYTrue_sort <- lapply(lYTrue_sort, function(x, nSplit) {list("train" = x[1:nSplit], "test" = x[(nSplit+1): length(x)])}, nSplit = nTrain)
         lYTrue_sort <- purrr::transpose(lYTrue_sort)

         lCombined[[i]] <- moc_data_prep_model(lStanObj = comb_fit_model, lYTrue = lYTrue_sort,
                                               nGroupSize = nGroupSize, vLag = vLag, vP = nP,
                                               .sPolyMatrix = "almon", .bEndpoint = FALSE, bolTrue = bolTrue)

      }

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
