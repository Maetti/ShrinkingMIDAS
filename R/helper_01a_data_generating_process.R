#' Create the Correlation Matrix
#' This will create the correlation matrix with a special structure
#'
#' @param vBlock
#' @param nVar
#' @param nWithin
#' @param vBetween
#' @param nSeed
#'
#' @return
#' @export
#'
#' @examples
dgp_correl_matrix <- function(vBlock, nVar, nWithin, nBetween) {

      mOut <- matrix(data = nBetween, nrow = sum(vBlock), ncol = sum(vBlock))
      sStart <- c(0, cumsum(vBlock))

      for (i in seq_along(vBlock)) {
            s1 <- sStart[i] + 1
            e1 <- sStart[i+1]

            mOut[s1:e1, s1:e1] <- nWithin

      }

      diag(mOut) <- nVar

      mOut

}



#' Title
#'
#' @param nY
#' @param nMu
#' @param nK
#' @param nFreq
#' @param nLag
#' @param maError
#' @param vRho
#'
#' @return
#' @export
#'
#' @examples
dgp_pred_create <- function(nY, nK, nFreq, nLag, mCorr, nMu, vRho, .sSeed = 1245, .bRandom = FALSE) {

      nX <- nY * nFreq * 2 + nLag

      maX <- matrix(runif(nX), nrow = nX, ncol = nK, byrow = T)

      if (.bRandom) {

            maX <- matrix(rnorm(nX), nrow = nX, ncol = nK, byrow = T)

      } else {

            nAR <- length(vRho)

            set.seed(.sSeed)
            mEps <- MASS::mvrnorm(nX, rep(0, nK), mCorr)

            for (i in 1:(nX - nAR)) {
                  maX[i + nAR, ] <- nMu + colSums(vRho * maX[i:(i + nAR - 1), , drop = FALSE]) + mEps[i + nAR, ]
            }

            maX <-maX[(nX - (nY * nFreq + nLag)):nX, ]

      }

      maX

}


#' Title
#'
#' @param nY
#' @param mX
#' @param nFreq
#' @param nLag
#'
#' @return
#' @export
#'
#' @examples
dgp_pred_align <- function(nY, mX, nFreq, nLag) {

      nK <- ncol(mX)
      nX <- nrow(mX)

      lOut <- vector("list", nK)


      for (i in 1:nK) {

            vX <- mX[, i]
            vStart <- seq(nLag, nX, nFreq)[1:nY]
            mNew <- matrix(nrow = nY, ncol = nLag)

            for (j in seq_along(vStart)) {
                  mNew[j, ] <- rev(vX[(vStart[j] - nLag + 1):vStart[j]])
            }

            lOut[[i]] <- mNew

      }

      lOut

}


#' Title
#'
#' @param nLag
#' @param nT1
#' @param nT2
#'
#' @return
#' @export
#'
#' @examples
dgp_exp_almon_lag <- function(nLag, nT1, nT2) {
      vLag <- 1:nLag
      vPoly <- exp(nT1 * vLag + nT2 * vLag^2)
      vPoly / sum(vPoly)
}

#' Title
#'
#' @param nLag
#' @param nT1
#' @param nT2
#'
#' @return
#' @export
#'
#' @examples
dgp_pred_beta_lag <- function(nLag, nT1, nT2){

      checkmate::assertNumeric(nT1, lower = 0)
      checkmate::assertNumber(nT2, lower = 0)

      vPoly <- seq(0, 1, length.out = nLag)
      #vPoly <- vPoly^(nT1 - 1) * (1 - vPoly)^(nT2 - 1)

      vBeta <- dbeta(vPoly, nT1, nT2)
      vBeta / sum(vBeta)
}

#' Title
#'
#' @param mX
#' @param nLag
#' @param nT1
#' @param nT2
#'
#' @return
#' @export
#'
#' @examples
dgp_pred_exp_almon <- function(mX, nLag, nT1, nT2, .sPoly = "almon") {

      sPoly <- match.arg(.sPoly, choices = c("almon", "beta"))

      if (sPoly == "almon") {
            vLagPoly <- dgp_exp_almon_lag(nLag, nT1, nT2)
      } else if (sPoly == "beta") {
            vLagPoly <- dgp_pred_beta_lag(nLag, nT1, nT2)
      }

      mX %*% vLagPoly
}




#' Title
#'
#' @param lX
#' @param nLag
#' @param nT1
#' @param nT2
#'
#' @return
#' @export
#'
#' @examples
dgp_pred_all_weighted <- function(nY, lX, nLag, nT1, nT2, .sPoly = "almon") {

      mNew <- matrix(nrow = nY, ncol = length(lX))

      for (i in seq_along(lX)) {
            mNew[, i] <- dgp_pred_exp_almon(lX[[i]], nLag, nT1, nT2, .sPoly)
      }

      mNew
}

#' Title
#'
#' @param nY
#' @param vK
#' @param vFreq
#' @param vLag
#' @param nMu
#' @param vRho
#' @param nVar
#' @param nWithin
#' @param nBetween
#' @param nT1
#' @param nT2
#'
#' @return
#' @export
#'
#' @examples
dgp_create_full <- function(nY, nK, vBeta, nFreq = 3, nLag = 6,
                            nMu = 0.5, nRho = 0.5,
                            nVar = 1, nWithin = 0.5, nBetween = 0.2,
                            nT1 = 0.0005, nT2 = -0.00007, .sPoly = "beta",
                            .sSeed, .bRandom = FALSE) {


   ## corrleation matrix
   mCorr <- dgp_correl_matrix(nK, nVar, nWithin, nBetween)

   ## predictor variables
   mpredX_raw <- dgp_pred_create(nY, nK, nFreq, nLag, mCorr, nMu, nRho, .sSeed = .sSeed, .bRandom = .bRandom)

   ## aligning predictor variable
   mpredX_align <- dgp_pred_align(nY, mpredX_raw, nFreq, nLag)

   ## apply midas weighting scheme (e.g exponential almon or beta)
   mPredX <- dgp_pred_all_weighted(nY, mpredX_align, nLag, nT1, nT2, .sPoly)

   ## calculate Y
   set.seed(.sSeed)
   vEpsilon = rnorm(nY, 0, 1)
   nY <- mPredX %*% vBeta + vEpsilon

   list(
      "y" = nY,
      "x" = mPredX,
      "x_raw" = mpredX_raw,
      "x_align" = mpredX_align
   )
}


#' Title
#'
#' @param nY
#' @param vK
#' @param vFreq
#' @param vLag
#' @param nMu
#' @param vRho
#' @param nVar
#' @param nWithin
#' @param nBetween
#' @param nT1
#' @param nT2
#'
#' @return
#' @export
#'
#' @examples
OLD_dgp_create_full <- function(nY, vBeta,
                            vK = c(40, 30, 30),
                            lBlock = list(c(15, 15, 10), c(10, 10, 10), c(15, 15)),
                            vFreq = c(3, 12, 90), vLag = c(6, 12, 25),
      nBlockNames = c("m", "w", "d"),
      nMu = 0.5, lRho = list(c(0.5), c(0.7), c(0.9)),
      nVar = 1, nWithin = 0.5, nBetween = 0.2,
      nT1 = 0.0005, nT2 = -0.00007, .sPoly = "almon", .sSeed, .bRandom = FALSE) {


      ## create data
      lpredX_align <- lpredX_raw <- vector("list", length = length(vK))
      lpredX <- vector("list", length = length(vK))
      # lBlock <- list(c(15, 15, 10), c(10, 10, 10), c(15, 15))

      for (i in seq_along(vK)) {
            # lpredX[[i]] <- dgp_prep_wrapper(nY, nK = vK[i],
            #                                 nFreq = vFreq[i], nLag = vLag[i], nMu = 0.5, vRho = lRho[[i]],
            #                                 vBlock = lBlock[[i]], nVar = nVar, nWithin = nWithin, nBetween = nBetween,
            #                                 nT1, nT2, .sSeed)

            ## generate correlation matrix
            mCorr <- dgp_correl_matrix(lBlock[[i]], nVar, nWithin, nBetween)

            lpredX_raw[[i]] <- dgp_pred_create(nY, vK[i], vFreq[i], vLag[i], mCorr, nMu, lRho[[i]], .sSeed = .sSeed, .bRandom = .bRandom)

            lpredX_align[[i]] <- dgp_pred_align(nY, lpredX_raw[[i]], vFreq[i], vLag[i])

            lpredX[[i]] <- dgp_pred_all_weighted(nY, lpredX_align[[i]], vLag[i], nT1, nT2, .sPoly)

      }

      mPredX <- do.call("cbind", lpredX)

      names(lpredX_raw) <- nBlockNames

      nY <- mPredX %*% vBeta + rnorm(nY)

      list(
            "y" = nY,
            "x" = mPredX,
            "x_raw" = lpredX_raw,
            "x_align" = lpredX_align
      )
}
