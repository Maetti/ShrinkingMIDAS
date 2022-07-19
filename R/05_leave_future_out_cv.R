#
# i <- 1
# j <- 1
# for (i in seq_along(nSim)) {
#       logger::log_info("----------------- Run: {i}/{length(nSim)} -----------------")
#       lData <- readRDS(here::here("inst/data/input", nSim[[i]]))
#
#       for (j in seq_along(lCMDmodels)) {
#             sModel <- names(lCMDmodels)[[j]]
#             logger::log_info("{sModel} (run {i})")
#
#             sDirStanModel <- paste0(here::here("inst/data/simulation_output/"),
#                                     tolower(gsub(" ", "_", sModel)), "/",
#                                     tolower(gsub(" ", "_", sModel)), "_", i, ".rds")
#
#             if (file.exists(sDirStanModel)) {
#
#                   logger::log_info("File already exists - skipped to next...")
#                   next()
#
#             } else {
#
#                   md_stan <- lCMDmodels[[j]]
#
#                   lStanObj <- md_stan$sample(
#                         data = lData,
#                         seed = 123,
#                         chains = 4,
#                         parallel_chains = 4,
#                         iter_sampling = 2000, iter_warmup = 2000,
#                         save_warmup = FALSE
#                   )
#
#                   stanfit <- rstan::read_stan_csv(lStanObj$output_files())
#                   saveRDS(stanfit, sDirStanModel)
#
#                   rm(lStanObj)
#                   rm(stanfit)
#
#             }
#
#             # lStanObj$save_object(file = paste0(here::here("inst/data/simulation_output/"),
#             #                                    tolower(gsub(" ", "_", sModel)), "/",
#             #                                    tolower(gsub(" ", "_", sModel)), "_", i, ".rds"))
#
#
#       }
#       rm(lData)
#       Sys.sleep(5)
# }
