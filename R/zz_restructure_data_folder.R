# ## create output folder
#
# s1 <- "/home/matthias/Schreibtisch/SoSe 2020/Empirical Finance and Financial Econometrics/assignment II/midasExperimentR/inst/data/01_dgp"
# dir(s1)
#
# for (i in dir(s1)) {
#       s2 <- paste0(s1, "/", i)
#       s3 <- paste0(s2, "/output")
#
#       # if (!dir.exists(s3))  {
#       #       dir.create(s3)
#       # }
#       #
#       s4_1 <- paste0(s3, "/object_stan")
#       s4_2 <- paste0(s3, "/stan_object")
#       file.rename(s4_1, s4_2)
#
#       s4_1 <- paste0(s3, "/extract_stan")
#       s4_2 <- paste0(s3, "/stan_extract")
#       file.rename(s4_1, s4_2)
#
#       s4_1 <- paste0(s3, "/summary")
#       s4_2 <- paste0(s3, "/stan_summary")
#       file.rename(s4_1, s4_2)
# }
#
#
# ## stan object to simulation folder
#
# rm(list = ls())
# s1 <- "/home/matthias/Schreibtisch/SoSe 2020/Empirical Finance and Financial Econometrics/assignment II/midasExperimentR/inst/data/simulation_output/"
# sModel <- dir(s1)
#
#
# sTo <- "/home/matthias/Schreibtisch/SoSe 2020/Empirical Finance and Financial Econometrics/assignment II/midasExperimentR/inst/data/01_dgp"
#
#
# for (i in sModel) {
#
#
#       dModel <- paste0(s1, i)
#
#       dirSim <- dir(dModel)
#
#       for (j in dirSim) {
#
#             nSimRegex <- as.numeric(unlist(stringr::str_extract_all(j, "\\(?[0-9]+\\)?")))
#             dirFrom <- paste0(dModel, "/", j)
#
#             if (length(nSimRegex) == 1) {
#                   nNumSim <- paste0(0, nSimRegex)
#             } else {
#                   nNumSim <- nSimRegex
#             }
#
#             dirTo <- paste0(sTo, "/", dir(sTo)[grep(nSimRegex, dir(sTo))[[1]]], "/output/object_stan/01_", nNumSim, "_object_stan_", i, ".rds")
#
#             file.copy(dirFrom, dirTo)
#
#       }
#
# }
#
#
#
# #   ____________________________________________________________________________
# #   Stan: Extract and Summary                                               ####
#
# ## create extract and summary files
# dirDGP <- here::here("inst/data/01_dgp/")
# for (i in dir(dirDGP)) {
#
#       print(i)
#
#       nSimRegex <- as.numeric(unlist(stringr::str_extract_all(i, "\\(?[0-9]+\\)?")))
#
#       sDirSimParq <- paste0(dirDGP, i, "/output/stan_extract/", nSimRegex)
#       if (!dir.exists(sDirSimParq)) {
#             dir.create(sDirSimParq)
#       }
#
#       sDirSummaryParq <- paste0(dirDGP, i, "/output/stan_summary/", nSimRegex)
#       if (!dir.exists(sDirSummaryParq)) {
#             dir.create(sDirSummaryParq)
#       }
#
#       ## extract stan object
#       sDirStanObj <- paste0(dirDGP, i,  "/output/stan_object")
#
#       lStanObj <- vector("list", length(dir(sDirStanObj)))
#       sExtractNams <- gsub(".*_stan_|\\.rds.*", "", dir(sDirStanObj))
#
#       for (j in 1:length(dir(sDirStanObj))) {
#             lStanObj[[j]] <- readRDS(paste0(sDirStanObj, "/", dir(sDirStanObj)[[j]]))
#       }
#
#       ## create Extract data
#       lExtract <- vector("list", length(lStanObj))
#       for (j in seq_along(lStanObj)) {
#
#             print(sExtractNams[[j]])
#
#             lE <- rstan::extract(lStanObj[[j]])
#
#             ## create Beta coefficient
#             lTheta <- moc_extract_theta_chains(lStanObj[[j]])
#             lBeta <- moc_create_beta_coef_parquet(lTheta, nGroupSize, vLag, vP = vP, .sPolyMatrix, .bEndpoint)
#
#             lE$beta <- lBeta
#
#             # lE <- lE[-which(names(lE) == "lp__")]
#             lE <- purrr::map2(lE, names(lE), function(x, y, z) {
#
#                   x <- as.data.frame(x)
#                   colnames(x) <- 1:length(x)
#                   x$iteration <- 1:nrow(x)
#
#
#
#                   x %>%
#                         tidyr::pivot_longer(names_to = "parameter_2", values_to = "value", -iteration) %>%
#                         dplyr::mutate(parameter_1 = y) %>%
#                         dplyr::select(iteration, parameter_1, parameter_2, value)
#
#             }, z = sExtractNams[[j]])
#
#             tblE <- purrr::reduce(lE, dplyr::bind_rows)
#
#             sDirStan_extract <- paste0(dirDGP, i, "/output/stan_extract/", nSimRegex, "/", sExtractNams[[j]])
#             if (!dir.exists(sDirStan_extract)) {
#                   dir.create(sDirStan_extract)
#             }
#
#             tf1 <- paste0(sDirStan_extract, "/01_", strsplit(i, "_")[[1]][[1]], "_stan_extract_", sExtractNams[[j]], ".parquet")
#             arrow::write_parquet(tblE, tf1)
#
#             ## create summary data
#             lS <- rstan::summary(lStanObj[[j]])$summary
#             dfS <- as.data.frame(lS)
#
#             tblS <-
#                   dfS %>%
#                   # dplyr::mutate(parameter = gsub("]", "", rownames(.))) %>%
#                   dplyr::mutate(parameter_raw = rownames(.)) %>%
#                   dplyr::mutate(parameter_1 = gsub("\\[.*", "", parameter_raw),
#                                 parameter_2 = stringr::str_extract_all(parameter_raw, "([0-9]+)", TRUE)) %>%
#                   dplyr::select(parameter_1, parameter_2, dplyr::everything(), -parameter_raw)
#
#             tblS$parameter_2 <- tblS$parameter_2[, 1]
#
#             tblS <- tblS %>% dplyr::as_tibble()
#
#             dfBetaSummary <-
#                   tblE[tblE$parameter_1 == "beta", ] %>%
#                         dplyr::group_by(parameter_2) %>%
#                         dplyr::summarise(mean = mean(value),
#                                          se_mean = mean(value) / n(),
#                                          sd = sd(value),
#                                          `2.5%` = quantile(value, 0.025),
#                                          `25%` = quantile(value, 0.25),
#                                          `50%` = quantile(value, 0.50),
#                                          `75%` = quantile(value, 0.75),
#                                          `97.5%` = quantile(value, 0.975),
#                                          n_eff = NA_real_,
#                                          Rhat = NA_real_)
#
#             dfBetaSummary$parameter_1 <- "beta"
#
#             tblS <- dplyr::bind_rows(tblS, dfBetaSummary)
#
#             sDirStan_summary <- paste0(dirDGP, i, "/output/stan_summary/", nSimRegex, "/", sExtractNams[[j]])
#             if (!dir.exists(sDirStan_summary)) {
#                   dir.create(sDirStan_summary)
#             }
#
#             tf2 <- paste0(sDirStan_summary, "/01_", strsplit(i, "_")[[1]][[1]], "_stan_summary_", sExtractNams[[j]], ".parquet")
#             arrow::write_parquet(tblS, tf2)
#       }
# }
