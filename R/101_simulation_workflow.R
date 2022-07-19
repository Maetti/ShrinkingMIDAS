#' Run Simulation
#'
#' This will create everything from directory to simulated data, running the model and saving the report.
#' This was chosen because {drake} seemed a little hard to work with for this specific purpose (I tried for hours on several days...).
#'
#' The inputs should be the simulation parameter, the runs of the simulation and an empty directory.
#'
#' @param sDirectory
#' @param nSimulation
#' @param nY
#' @param vBeta
#' @param vK
#' @param lBlock
#' @param vFreq
#' @param vLag
#' @param nBlockNames
#' @param nMu
#' @param lRho
#' @param nVar
#' @param nWithin
#' @param nBetween
#' @param nT1
#' @param nT2
#' @param nTrain
#' @param vP
#'
#' @return
#' @export
#'
#' @examples
#' sDirectory <- paste0(here::here("inst/data/"), "01_dgp")
#' nSimulation <- 10
#'
#'
#' vBeta <- c(0, 0.4, 1.2, 0, 0, 0.8, 1, 0, 0, 0)
#' bolTrue <- vBeta != 0
#' vTrueValue <- vBeta[bolTrue]
#' nY <- 200
#' nTrain <- 150
#' vK = c(10)
#' vP = 3
#' lBlock = list(c(10))
#' vFreq = 3
#' vLag = 6
#' nBlockNames = "m"
#' nMu = 0.1
#' lRho = list(c(0.5))
#' nVar = 1
#' nWithin = 0.5
#' nBetween = 0
#' nT1 = 0.0005
#' nT2 = -0.00007
#'
#'
run_simulation <- function(sDirectory, nSimulation,
                           nY = 200, nTrain = 150,
                           vBeta = c(0, 0.4, 1.2, 0, 0, 0.8, 1, 0, 0, 0),
                           vK = c(10), lBlock = list(c(10)), vFreq = 3, vLag = 6, vP = 3,
                           nBlockNames = "m",
                           nMu = 0.1, lRho = list(c(0.5)),
                           nVar = 1, nWithin = 0.5, nBetween = 0,
                           nT1 = 0.0005, nT2 = -0.00007) {



##  ............................................................................
##  Check Input                                                             ####

      stopifnot(dir.exists(sDirectory))

      dirSimulation <- sDirectory
      # dirSimulation <- paste0(sDirectory, "/simulation")
      dirLogging <- paste0(sDirectory, "/log")

      if (!dir.exists(dirSimulation)) dir.create(dirSimulation)
      if (!dir.exists(dirLogging)) dir.create(dirLogging)

      ## set logging
      sLogFile <- paste0(dirLogging,"/", gsub("-", "_", Sys.Date()), "_simulation_run")
      logger::log_appender(logger::appender_tee(sLogFile))

      nGroupSize <- rep(vP + 1, lBlock[[1]]) # group Index


#   ____________________________________________________________________________
#   Create Directories                                                      ####

      logger::log_info("Staring with directory structure")

      dirSim <- 1:nSimulation
      dirSim2 <- vector("character", length(dirSim))
      lDirSim <- vector("list", length(dirSim))
      for (i in seq_along(dirSim)) {
            if (dirSim[[i]] < 10) {
                  dirSim2[i] <- paste0("0", dirSim[[i]], "_simulation")
            } else {
                  dirSim2[i] <- paste0(dirSim[[i]], "_simulation")
            }

            dirNewSim <- paste0(dirSimulation, "/", dirSim2[i])
            lDirSim[[i]] <- dirNewSim

            if (!dir.exists(dirNewSim)) {
                  logger::log_info("Creating Directory Structure: {dirSim2[i]}")
                  dir.create(dirNewSim)

                  dirNew2 <- paste0(dirNewSim, "/", c("output", "input"))
                  lapply(dirNew2, dir.create)

                  dirNew3 <- paste0(dirNew2[[1]], "/stan_", c("extract", "object", "summary"))
                  lapply(dirNew3, dir.create)

                  dirNew4 <- paste0(dirNew3, "/", i)
                  lapply(dirNew4, dir.create)
            }
      }

      ## set seed
      nSeed <- 1001:(1001 + nSimulation)


#   ____________________________________________________________________________
#   Create Simulated Data                                                   ####

      logger::log_info("Creating Simulated Data")
      for (i in 1:nSimulation) {

            logger::log_info("Sim Data: {i}/{nSimulation}")
            x1 <- dgp_create_full(nY = nY, vBeta = vBeta,
                            vK = vK, lBlock = lBlock, vFreq = vFreq, vLag = vLag,
                            nBlockNames = nBlockNames,
                            nMu = nMu, lRho = lRho,
                            nVar = nVar, nWithin = nWithin, nBetween = nBetween,
                            nT1 = nT1, nT2 = nT2,
                            .sSeed = nSeed[[i]])

            x2 <- create_predictor_lag_poly(x1[["x_align"]],
                                            vLag = vLag, vP = vP,
                                            .sPolyMatrix = "almon", .bEndpoint = FALSE)

            lModelInput <- create_model_input(maY = x1[["y"]], dfDataX = x2[[1]],
                                              nTrain = nTrain, nG = lBlock[[1]],
                                              nGroupSize = nGroupSize)

            saveRDS(lModelInput,
                    paste0(lDirSim[[i]], "/input/simulatedData_", i, "_seed_", nSeed[[i]], ".rds"))
      }




#   ____________________________________________________________________________
#   Create Models (if necessary)                                            ####

      logger::log_info("Create Models")

      lCMDmodels <- vector("list", length(stanmodels))
      dirStan <- here::here("inst/stan")
      sStanFiles <- dir(dirStan)[grep("stan", dir(dirStan))]

      for (i in seq_along(sStanFiles)) {
            sStanModel <- glue::glue("{dirStan}/{sStanFiles[i]}")
            lCMDmodels[[i]] <- cmdstanr::cmdstan_model(stan_file = sStanModel)
      }

      sModelName <- stringr::str_to_title(gsub("_", " ", stringr::str_extract(sStanFiles, "[^\\.]+")))

      logger::log_info("Total Models: {length(names(sModelName))}- {paste0(names(sModelName), collapse = ', ')}")

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
      for (i in 1:nSimulation) {
            logger::log_info("----------------- Run: {i}/{length(nSim)} -----------------")

            dirSim_input <- paste0(lDirSim[[i]], "/raw_input")
            lData <- readRDS(paste0(dirSim_input, "/", dir(dirSim_input)))

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


}
