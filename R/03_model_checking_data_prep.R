#' Title
#'
#' @param lInpFit
#'
#' @return
#' @export
#'
#' @examples
moc_prep_summary <- function(lStanObj) {

      lOut <- vector("list", length(lStanObj))
      for (i in seq_along(lStanObj)) {
            lOut[[i]] <-
               rstan::summary(lStanObj[[i]])$summary %>%
               dplyr::mutate(simulation = i)
      }

      lOut %>% dplyr::bind_rows()
}


#' Title
#'
#' @param tblDiag
#'
#' @return
#' @export
#'
#' @examples
moc_diagnose_rename <- function(tblInput) {
   ## adjust names
   s1 <- gsub("\\[", "_", names(tblInput))
   s2 <- gsub("\\]", "", s1)
   colnames(tblInput) <- s2

   tblInput
}


#' Title
#'
#' @param lStanObj A fitted Stan model
#' @param sMatch The name of the variables
#' @param sExtract What to extracts (eg.: Rhat)
#'
#' @return
#' @export
#'
#' @examples
moc_diagnose_extract <- function(lStanObj, sMatch = "theta\\[", sExtract = "Rhat") {

   tblOut <-
      lStanObj %>%
         purrr::map(., function(x, y, z) {

            sObj <- shinystan::as.shinystan(x)
            sOut <- shinystan::retrieve(sObj, y)
            sOut[grep(z, names(sOut))]

         }, y = sExtract, z = sMatch) %>%
         dplyr::bind_rows() %>%
         dplyr::mutate(simulation = 1:nrow(.)) %>%
         dplyr::select(simulation, dplyr::everything())

   tblOut <- moc_diagnose_rename(tblOut)

   tblOut
}

#' Title
#'
#' @param lStanObj
#' @param maQ
#'
#' @return a tibble with run and chain as index and the theta coefficient as values
#' @export
#'
#' @examples
moc_extract_theta_chains <- function(stanObj) {

   ## create
   arStan <- as.array(stanObj, pars = c("theta"))
   lTheta <- lapply(seq_along(dimnames(arStan)[[3]]),
                    function(i) {
                       x <- arStan[, ,i]
                       x <- cbind(1:nrow(x), x)
                       colnames(x)[1] <- "run"
                       x %>%
                          as.data.frame() %>%
                          tidyr::pivot_longer(names_to = "chain", values_to = "value", cols = -1)
                  })
   names(lTheta) <- dimnames(arStan)[[3]]

   lTheta %>%
      purrr::map2(., names(lTheta), function(x, y) {
         colnames(x)[3] <- y
         x
      }) %>%
      purrr::reduce(., dplyr::left_join, by = c("run", "chain"))
}


#' Title
#'
#' @param dfTheta
#' @param nGroupSize
#' @param vLag
#' @param vP
#' @param .sPolyMatrix
#' @param .bEndpoint
#'
#' @return
#' @export
#'
#' @examples
moc_create_beta_coef <- function(dfTheta, nGroupSize,
                                 vLag, vP,
                                 .sPolyMatrix = "almon", .bEndpoint = FALSE) {

   ## Group Index
   vIndGroup <- c(0, cumsum(nGroupSize))

   ## Lag Matrix
   if (.sPolyMatrix == "almon") {
      mQ <- prep_model_lag_matrix(vP, vLag, .bEndpoint)
   } else if(.sPolyMatrix == "legendre") {
      mQ <- prep_model_lag_matrix_legendre(vP, 0, 1, vLag)
   }

   maTheta <- as.matrix(dfTheta[, -c(1, 2)])
   lBeta <- vector("list", length(nGroupSize))
   for (i in 1:(length(vIndGroup) - 1)) {
      lBeta[[i]] <- rowSums(maTheta[, (vIndGroup[[i]] + 1):vIndGroup[[i + 1]]] %*% mQ)
   }

   maBeta <- do.call("cbind", lBeta)

   dfBeta <- cbind(dfTheta[, c(1, 2)], maBeta)
   colnames(dfBeta)[-c(1, 2)] <- paste0("beta_", seq_along(nGroupSize))
   dfBeta
}

#' Title
#'
#' @param dfTheta
#' @param nGroupSize
#' @param vLag
#' @param vP
#' @param .sPolyMatrix
#' @param .bEndpoint
#'
#' @return
#' @export
#'
#' @examples
moc_create_beta_coef_parquet <- function(dfTheta, nGroupSize,
                                 vLag, vP,
                                 .sPolyMatrix = "almon", .bEndpoint = FALSE) {

   ## Group Index
   vIndGroup <- c(0, cumsum(nGroupSize))

   ## Lag Matrix
   if (.sPolyMatrix == "almon") {
      mQ <- prep_model_lag_matrix(vP, vLag, .bEndpoint)
   } else if(.sPolyMatrix == "legendre") {
      mQ <- prep_model_lag_matrix_legendre(vP, 0, 1, vLag)
   }

   maTheta <- as.matrix(dfTheta[, -c(1, 2)])
   lBeta <- vector("list", length(nGroupSize))
   for (i in 1:(length(vIndGroup) - 1)) {
      lBeta[[i]] <- rowSums(maTheta[, (vIndGroup[[i]] + 1):vIndGroup[[i + 1]]] %*% mQ)
   }

   maBeta <- do.call("cbind", lBeta)

   # colnames(dfBeta) <- seq_along(nGroupSize)
   maBeta
}



#' Title
#'
#' @param dfBeta
#' @param bolTrue logical vector to indicate which beta coeffiecients are truly non zero
#'
#' @return
#' @export
#'
#' @examples
moc_create_confusion_matrix <- function(maBeta, bolTrue) {

   bolSim <- apply(maBeta, 2, function(x) {findInterval(0, quantile(x, probs = c(0.05, 0.95))) == 0L })
   as.data.frame(table(bolSim, bolTrue))
}



#' Title
#'
#' @param lStanObj
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_summary <- function(lStanObj, lYTrue, .sY_extract = "y_pred") {

   lY <- vector("list", length(lStanObj))
   .sY_extract <- match.arg(.sY_extract, c("y_pred", "y_rep"))

   for (i in seq_along(lY)) {
      s1 <- rstan::summary(lStanObj[[i]], .sY_extract)$summary
      s1 <- cbind(s1[, c("mean", "2.5%", "97.5%")], "y_true" = lYTrue[[i]])

      lY[[i]] <-
         as.data.frame(s1) %>%
         dplyr::mutate(parameter = gsub("\\[", "_", gsub("]", "", rownames(.)))) %>%
         tidyr::pivot_longer(names_to = "stat", values_to = "value", -c("parameter")) %>%
         dplyr::mutate(simulation = i) %>%
         dplyr::select(simulation, dplyr::everything())
   }

   lY %>% dplyr::bind_rows()

}

#' Extract Data for all Simulation
#'
#' This function extracts the needed data (theta, beta) for a model for all simulations.
#'
#' @param lStanObj
#' @param nGroupSize
#' @param vLag
#' @param vP
#' @param .sPolyMatrix
#' @param .bEndpoint
#' @param bolTrue
#'
#' @return
#' @export
#'
#' @examples
moc_data_prep_model <- function(lStanObj, lYTrue,
                                nGroupSize, vLag, vP, .sPolyMatrix = "almon", .bEndpoint = FALSE,
                                nWarmUp = 2000,
                                bolTrue) {

   lTheta <- lBeta <- lConfusion <- vector("list", length(lStanObj))
   for (i in seq_along(lStanObj)) {
      print(i)
      lTheta[[i]] <- moc_extract_theta_chains(lStanObj[[i]])
      lBeta[[i]] <- moc_create_beta_coef(lTheta[[i]], nGroupSize, vLag, vP = vP, .sPolyMatrix, .bEndpoint)
      lConfusion[[i]] <- moc_create_confusion_matrix(lBeta[[i]][, -c(1, 2)], bolTrue)
   }

   foo_sim_bind <-
      . %>% purrr::map2(., seq_along(.), function(x, y) {
               x %>%
                  dplyr::mutate(simulation = y) %>%
                  dplyr::select(simulation, dplyr::everything())
            }
      ) %>%
      dplyr::bind_rows() %>%
      tibble::as_tibble()

   ## prep output
   tblConf <- lConfusion %>% foo_sim_bind()

   tblTheta <-
      lTheta %>%
         foo_sim_bind() %>%
         moc_diagnose_rename() %>%
         tidyr::pivot_longer(names_to = "parameter", values_to = "value", -c(1:3))

   tblBeta <-
      lBeta %>%
         foo_sim_bind() %>%
         tidyr::pivot_longer(names_to = "parameter", values_to = "value", -c(1:3))

   tblRhat <-
      lStanObj %>%
         moc_diagnose_extract() %>%
            tidyr::pivot_longer(names_to = "parameter", values_to = "value", -1)

   tblNeff <-
      lStanObj %>%
         moc_diagnose_extract(sExtract = "n_eff") %>%
         tidyr::pivot_longer(names_to = "parameter", values_to = "value", -1)

   lYrep <-
      lStanObj %>%
         moc_yrep_summary(., lYTrue$train, .sY_extract = "y_rep")

   lYPred <-
      lStanObj %>%
         moc_yrep_summary(., lYTrue$test, .sY_extract = "y_pred")

   ## output
   list(
      "theta" = tblTheta,
      "beta" =  tblBeta,
      "confusion" = tblConf,
      "rhat" = tblRhat,
      "yrep" = lYrep,
      "ypred" = lYPred
   )

}


#' Title
#'
#' @param lAll
#'
#' @return
#' @export
#'
#' @examples
#' lAll <- drake::readd(allModel)
#' glimpse(lAll)
# moc_data_prep_all <- function(lAll,
#                               nGroupSize, vLag, vP, .sPolyMatrix = "almon", .bEndpoint = FALSE,
#                               bolTrue) {
#
#    lAll <- unlist(lAll, recursive = FALSE)
#    sAllName <- names(lAll)
#    lOut <- vector("list", length(lAll))
#
#    for (i in seq_along(lAll)) {
#
#       sModelName <- sAllName[[i]]
#       lOut[[i]] <- moc_data_prep_model(lAll[[i]],
#                                        nGroupSize, vLag, vP, .sPolyMatrix = "almon", .bEndpoint = FALSE,
#                                        bolTrue)
#
#    }
#
#    names(lOut) <- sAllName
#    lOut <- purrr::transpose(lOut)
#    for (i in seq_along(lOut)) {
#
#       lOut[[i]] <-
#          lOut[[i]]  %>%
#             purrr::map2(., names(.), function(x, y) {
#                x %>% dplyr::mutate(model = y) %>% dplyr::select(model, dplyr::everything())
#             }) %>%
#             dplyr::bind_rows()
#    }
#
#    lOut
# }

#' Title
#'
#' @param lAll
#'
#' @return
#' @export
#'
#' @examples
moc_data_combine_all <- function(lAll) {

   # lAll <- unlist(lAll, recursive = FALSE)
   sAllName <- names(lAll)
   lAll <- purrr::transpose(lAll)

   for (i in seq_along(lAll)) {

      lAll[[i]] <-
         lAll[[i]]  %>%
            purrr::map2(., names(.), function(x, y) {
               x %>% dplyr::mutate(model = y) %>% dplyr::select(model, dplyr::everything())
            }) %>%
            dplyr::bind_rows()
   }

   lAll

}


