#' Title
#'
#' This function is just a wrapper so the code does not get triggered when package is load
#'
#' @return
#' @export
#'
#' @examples
masterSimulation <- function() {

      #   ____________________________________________________________________________
      #   Set Inputs                                                              ####

      library(dplyr)

      ## set logging file
      logger::log_appender(logger::appender_tee(here::here("inst", "simulation", "logging", "simulation",
                                                           glue::glue("{gsub('-', '_', Sys.Date())}_simulation_log"))))
      logger::log_info("Setting input variables")

      ## simulations
      nSimulation <- 25
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
      # nSim <- nSim[1:10]

      #   ____________________________________________________________________________
      #   Creating Data                                                           ####

      logger::log_info("Simulating input data for: {nSimulation} simulations")
      for (i in 1:nSimulation) {

            ## check if simulated file already exists
            ## save input data
            if (i < 10) {
               nSimRunSave <- paste0("00", i)
            } else if (i >= 10 & i < 100) {
               nSimRunSave <- paste0("0", i)
            } else {
               nSimRunSave <- as.character(i)
            }

            sDirSave <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

            if (!file.exists(sDirSave) || bForceNewData) {

               logger::log_info("Create data for simulation: {nSimRunSave}")

               lRawData <- dgp_create_full(nY = nY, vBeta = vBeta,
                                     nK = nK, nFreq = nFreq, nLag = nLag,
                                     nMu = nMu, nRho = nRho,
                                     nVar = nVar, nWithin = nWithin, nBetween = nBetween,
                                     nT1 = nT1, nT2 = nT2,
                                     .sSeed = nSeed[[i]])

               # x2 <- create_predictor_lag_poly(x1[["x_align"]],
               #                                 vLag = vLag, vP = vP,
               #                                 .sPolyMatrix = "almon", .bEndpoint = FALSE)

               dfX <- do.call("cbind", lRawData$x_align)

               lModelInput <- create_model_input(maY = lRawData[["y"]], dfDataX = dfX,
                                                 nTrain = nTrain, nG = nG,
                                                 nGroupSize = nGroupSize)

               lDataSave <- list(
                     "x_raw" = lRawData$x_raw,
                     "x_align" = lRawData$x_align,
                     "y" = lRawData$y,
                     "model_data" = lModelInput)





               # helper_createFolder(sDirInput)

               ## !!OUTCOMMENT!!
               saveRDS(lDataSave, file = sDirSave)
            }

            # saveRDS(sDataSave,
            #         paste0(sDir, "/input/simulatedData_", i, "_seed_", nSeed[[i]], ".rds"))
      }


      logger::log_info("Done with input data generation")

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
      #   Fit Models                                                              ####

      # TODO: 19.07.2022 stopped here...such a shame
      # * for each simulation create full output
      #  - directory structure
      #  - stan object
      #  - extract parquet
      #  - summary parquet

      logger::log_info("Starting model fitting for each generated data")

      lTimes <- vector("list", length(sModelName))

      for (i in 1:nSimulation) {

            logger::log_info("----------------- Run: {i}/{nSimulation} -----------------")

            ## get simulated data
            if (i < 10) {
                  nSimRunSave <- paste0("00", i)
            } else if (i >= 10 & i < 100) {
                  nSimRunSave <- paste0("0", i)
            } else {
                  nSimRunSave <- as.character(i)
            }

            sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

            logger::log_info("Load input file: {sDirInputLoad}")
            lData <- readRDS(sDirInputLoad)
            lDataStanInput <- lData$model_data

            ## estimate each model for a certain input data
            for (j in seq_along(lCMDmodels)) {

                  sModel <- names(lCMDmodels)[[j]]
                  logger::log_info("{sModel} (run {i})")

                  sModelSaveName <- tolower(gsub(" ", "_", sModel))

                  ## check if model output folder exist and create if necessary
                  sFolderModel <- here::here(sDir, "output", "01_stan_objects", sModelSaveName)
                  helper_createFolder(sFolderModel)

                  sDirStanModel <- here::here(sDir, "output", "01_stan_objects", sModelSaveName, glue::glue("{nSimRunSave}_{sModelSaveName}.rds"))

                  if (file.exists(sDirStanModel)) {

                     logger::log_info("File already exists - skipped to next...")

                  } else {

                     logger::log_info("Estimte Model...")
                     ## estimate and time model
                     start_time_sampling <- Sys.time()

                     md_stan <- lCMDmodels[[j]]

                     lStanObj <- md_stan$sample(
                           data = lDataStanInput,
                           seed = 123,
                           chains = 4,
                           parallel_chains = 4,
                           iter_sampling = 4000, iter_warmup = 4000,
                           save_warmup = FALSE,
                           adapt_delta = 0.99, max_treedepth = 10
                     )

                     end_time_sampling <- Sys.time()

                     lTimes[[j]] <- c(lTimes[[j]], end_time_sampling - start_time_sampling)

                     ## extract and save model output
                     stanfit <- rstan::read_stan_csv(lStanObj$output_files())
                     saveRDS(stanfit, sDirStanModel)

                     rm(lStanObj)
                     rm(stanfit)
                  }

                  #}

                  # lStanObj$save_object(file = sDirStanModel)

                  # lStanObj$save_object(file = paste0(here::here("inst/data/simulation_output/"),
                  #                                    tolower(gsub(" ", "_", sModel)), "/",
                  #                                    tolower(gsub(" ", "_", sModel)), "_", i, ".rds"))


            } # end of model estimation for one input

            ## remove input data to make sure no overlaps whatsoever

            rm(lData)
            rm(lDataStanInput)
            gc()
            Sys.sleep(5)
      }




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




      ## combine data from simulation for each model
      lCombined <- vector("list", length(sModelName))


      for (i in seq_along(sModelName)) {
            # for (i in seq_along(sModelName[1:2])) {
            # i <- 1
            sModel <- sModelName[[i]]

            sModelSaveName <- tolower(gsub(" ", "_", sModel))

            ## check if model output folder exist and create if necessary
            sFolderModel <- here::here(sDir, "output", "01_stan_objects", sModelSaveName)
            vAllSims <- sort(list.files(sFolderModel, full.names = TRUE))

            logger::log_info("Extract data for: {sModel}")

            ## create model folder for raw and summary data
            sFolderRaw <- here::here(sDir, "output", "02_raw_extracted", sModelSaveName)
            sFolderSummary <- here::here(sDir, "output", "03_summaries", sModelSaveName)
            helper_createFolder(sFolderRaw)
            helper_createFolder(sFolderSummary)

            for (j in seq_along(vAllSims)) {

               ## create names and folder
               nSimRunSave <- helper_create_number_name(j)

               sSimRaw <- here::here(sDir, "output", "02_raw_extracted", sModelSaveName, j)
               sSimSummary <- here::here(sDir, "output", "03_summaries", sModelSaveName, j)

               helper_createFolder(sSimRaw)
               helper_createFolder(sSimSummary)

               ## load stand object
               stanObj <- readRDS(vAllSims[[j]])

               list_of_draws <- rstan::extract(stanObj)
               vPars <- names(list_of_draws)

               #### prepare and save raw data
               tblRaw <-
                  purrr::map(vPars, . %>% moc_prepare_raw_data(stanObj = stanObj, sPars = .)) %>%
                  dplyr::bind_rows()

               ## calculate beta coefficient
               tblTheta <- tblRaw %>% dplyr::filter(parameter == "theta")
               tblBeta <- moc_create_beta_per_chain(tblTheta, nGroupSize = nGroupSize)

               tblRaw <- dplyr::bind_rows(tblRaw, tblBeta)

               arrow::write_parquet(tblRaw, glue::glue("{sSimRaw}/{nSimRunSave}_{sModelSaveName}_raw.parquet"))


               ## prepare and save summary file
               stanSummary <- rstan::summary(stanObj)
               stanSummary <- stanSummary$summary

               tblSummary <-
                  stanSummary %>%
                     as.data.frame() %>%
                     dplyr::mutate(parameter = rownames(.)) %>%
                     dplyr::mutate(key = as.numeric(regmatches(parameter, gregexpr("\\[\\K[^\\]]+(?=\\])", parameter, perl=TRUE)))) %>%
                     dplyr::mutate(parameter = gsub("\\[.*", "", parameter)) %>%
                     dplyr::select(parameter, key, dplyr::everything()) %>%
                     tidyr::pivot_longer(names_to = "summary", values_to = "value", -c(1, 2))


               arrow::write_parquet(tblRaw, glue::glue("{sSimSummary}/{nSimRunSave}_{sModelSaveName}_summary.parquet"))


            }

            # ## load all simulation files
            # comb_fit_model <- lapply(sort(list.files(sFolderModel, full.names = TRUE)), readRDS)
            #
            # # nSimRegex <- as.numeric(unlist(stringr::str_extract_all(sFileName, "\\(?[0-9]+\\)?")))
            # lYTrue_sort <- lYTrue_raw
            # lYTrue_sort <- lapply(lYTrue_sort, function(x, nSplit) {list("train" = x[1:nSplit], "test" = x[(nSplit+1): length(x)])},
            #                       nSplit = nTrain)
            # lYTrue_sort <- purrr::transpose(lYTrue_sort)
            #
            # ## create raw output for each simulation and model in parquet format
            # lCombined[[i]] <- moc_data_prep_model(lStanObj = comb_fit_model, lYTrue = lYTrue_sort,
            #                                       nGroupSize = nGroupSize, vLag = nLag, vP = nP,
            #                                       .sPolyMatrix = "beta", .bEndpoint = FALSE, bolTrue = bolTrue)
      }


      ## stopped here on 26.7.2022

}

