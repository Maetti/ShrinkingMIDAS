#' Title
#'
#' @param variables
#'
#' @return
#' @export
#'
#' @examples
factor_midas <- function(lData) {


   ## prep training data
   nTrain <- lData$model_data$nY_train
   dfXTrain <- lData$x_raw[1:(6 + 3 * nTrain), ]

   ## PCA
   dfXAdj <- apply(dfXTrain, 2, function(x) x - mean(x))
   mdPCA <- prcomp(dfXAdj, center = FALSE, scale. = FALSE)
   mEigenVector <- mdPCA$rotation
   kInd <- 15

   maFeature <- dfXAdj %*% mEigenVector[, 1:kInd]

   ## align training features
   lFeature <- dgp_pred_align(nTrain, maFeature, nFreq = 3, nLag = 6)
   maFeatureAligned <- do.call("cbind", lFeature)

   dfFeature <- as.data.frame(maFeatureAligned)
   colnames(dfFeature) <- paste0("x_", 1:ncol(dfFeature))

   ## model

   mdFactor <- lm(y_train ~ . - 1, data = dfFeature)

   ## prep test data
   dfXAll <- lData$x_raw

   ## PCA
   dfXAdjAll <- apply(dfXAll, 2, function(x) x - mean(x))
   mdPCA_all <- prcomp(dfXAdjAll, center = FALSE, scale. = FALSE)
   mEigenVectorAll <- mdPCA_all$rotation
   maFeatureAll <- dfXAdjAll %*% mEigenVectorAll[, 1:kInd]

   ## align test data
   lFeatureAll <- dgp_pred_align(nY, maFeatureAll, nFreq = 3, nLag = 6)
   maTransAll <- do.call("cbind", lFeatureAll)


   maTest <- maTransAll[(nTrain + 1):nrow(maTransAll), ]
   dfTest <- as.data.frame(maTest)
   colnames(dfTest) <- paste0("x_", 1:ncol(dfTest))

   ## prediction
   predict.lm(mdFactor, dfTest)

}

#' Title
#'
#' @param vTheta
#'
#' @return
#' @export
#'
#' @examples
midas_coef <- function(vTheta) {

      sInd <- which(1:length(vTheta) %% 2 == 1)
      nLag <- 6

      lCoef <- vector("list", length(sInd))

      for (i in seq_along(sInd)) {

            nT1 <- vTheta[[sInd[[i]]]]
            nT2 <- vTheta[[sInd[[i]] + 1]]

            lCoef[[i]] <- dgp_exp_almon_lag(nLag, nT1, nT2)
      }

      unlist(lCoef)

}


#' Title
#'
#' @param vTheta
#'
#' @return
#' @export
#'
#' @examples
midas_optim <- function(vTheta) {

      coeffs <- midas_coef(vTheta)

      y_hat <- X %*% coeffs
      r <- y - y_hat
      sum(r^2, na.rm = TRUE)

}


#' Title
#'
#' @param vTheta
#'
#' @return
#' @export
#'
#' @examples
midas_rhs <- function(vTheta) {

      coeffs <- midas_coef(vTheta)
      X %*% coeffs

}


#' Title
#'
#' @param lInput
#'
#' @return
#' @export
#'
#' @examples
midas_coef_old <- function(lInput) {

      lCoef <- vector("list", length = length(lInput))
      for (i in seq_along(lInput)) {

            nLag <- lInput[[i]]$lag
            nT1 <- lInput[[i]]$theta[[1]]
            nT2 <- lInput[[i]]$theta[[2]]

            lCoef[[i]] <- dgp_exp_almon_lag(nLag, nT1, nT2)
      }

      unlist(lCoef)

}
