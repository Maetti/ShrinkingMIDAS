#' Title
#'
#' @param dfTheta_sim
#' @param gSize
#' @param nLag
#' @param nP
#' @param .bEndpoint
#'
#' @return
#' @export
#'
#' @examples
md_extract_beta_coef_test <- function(vBeta_est, gSize, nLag = c(6), nP = c(2), .bEndpoint = FALSE) {

      mQ <- dgp_lag_matrix(nP, nLag, .bEndpoint)

      maOut <- matrix(0, nrow = length(gSize), ncol = 1)
      vInd <- cumsum(gSize)

      for (i in seq_along(gSize)) {
            maOut[i, ] <- sum(vBeta_est[c((vInd[i] - (gSize[i] - 1)):vInd[i])] %*% mQ)
      }

      maOut

}
#' Title
#'
#' @return
#' @export
#'
#' @examples
md_extract_beta_coef <- function(dfTheta_sim, gSize, nLag = c(24), nP = c(2), .bEndpoint = FALSE) {

      mQ <- dgp_lag_matrix(nP, nLag, .bEndpoint)

      maOut <- matrix(0, nrow = nrow(dfTheta_sim), ncol = length(gSize))
      vInd <- cumsum(gSize)

      for (i in seq_along(gSize)) {
            maOut[, i] <- apply(dfTheta_sim[, c(vInd[i] - 1, vInd[i])], 1, function(x, mQ) {rowSums(x %*% mQ)}, mQ = mQ)
      }

}

#' Title
#'
#' @param maInput
#'
#' @return
#' @export
#'
#' @examples
check_summary = function(maInput) {

      sColNames <- factor(colnames(maInput), levels = colnames(maInput))

      maInput %>%
            as.data.frame() %>%
            tidyr::pivot_longer(values_to = "val", names_to = "input", cols = 1:ncol(.)) %>%
            dplyr::mutate(input = factor(input, levels = colnames(maInput))) %>%
            dplyr::group_by(input) %>%
            dplyr::summarise(
                  amount = dplyr::n(),
                  mean = mean(val),
                  se_mean = sd(val) / amount,
                  sd = sd(val),
                  p_0.025 = quantile(val, probs = 0.025),
                  p_0.25 = quantile(val, probs = 0.25),
                  p_0.5 = quantile(val, probs = 0.5),
                  p_0.75 = quantile(val, probs = 0.75),
                  p_0.975 = quantile(val, probs = 0.975)
            ) %>%
            dplyr::select(-amount)

}
#' Title
#'
#' @param dfInput
#'
#' @return
#' @export
#'
#' @examples
check_tpr <- function(dfInput, sTrueBeta, vTF = c(9, 21)) {

      dfBeta <- check_summary(dfInput)

      estBeta <-
            dfBeta %>%
            dplyr::filter(!c(p_0.025 < 0 & p_0.975 > 0)) %>%
            dplyr::pull(input)

      tpr <- sum(estBeta %in% sTrueBeta) / vTF[1]
      fpr <- sum(!(estBeta %in% sTrueBeta)) / vTF[2]

      data.frame(
            "tpr" = tpr,
            "fpr" = fpr
      )
}

#' Title
#'
#' @param dfInput
#' @param vTrueBeta
#'
#' @return
#' @export
#'
#' @examples
check_mse <- function(dfInput, vTrueBeta) {

      vMean <- apply(dfInput, 2, mean)
      nVar <- (1 / (nrow(dfInput) * ncol(dfInput))) * sum(t(t(dfInput) - vMean)^2)
      nBias <- (1/ncol(dfInput)) * sum((vMean - vBeta)^2)

      data.frame(
            "mse" = nVar + nBias,
            "var" = nVar,
            "bias_2" = nBias
      )

}

#' Title
#'
#' @param xTest_scale
#' @param md_stan_obj
#'
#' @return
#' @export
#'
#' @examples
pred_new_response_wrapper <- function(xTest_scale, md_stan_obj) {

      ## model specific
      mTheta <- md_stan_obj$list_of_draws$theta
      nSigma_2 <- md_stan_obj$list_of_draws$sigma_2

      ## select theta
      mTheta <- as.data.frame(mTheta)
      colnames(mTheta) <- paste0("theta_", 1:ncol(mTheta))
      sInd <-
            check_summary(mTheta) %>%
            dplyr::filter(!(p_0.025 < 0 & p_0.975 > 0)) %>%
            dplyr::pull(input)

      sInd <- match(sInd, colnames(mTheta))

      mTheta <- as.matrix(mTheta[, sInd])
      xTest_scale <- xTest_scale[, sInd]


      mY_pred <- pred_new_response_bayes(xTest_scale, mTheta, nSigma_2)

      mY_pred

}
#' Title
#'
#' @param dfX
#' @param mTheta
#' @param nSigma_2
#'
#' @return
#' @export
#'
#' @examples
pred_new_response_bayes <- function(dfX, mTheta, nSigma_2) {

      mYnew <- matrix(nrow = length(nSigma_2), ncol = nrow(dfX))

      for (i in 1:nrow(dfX)) {
            for (j in 1:length(nSigma_2)) {
                  mYnew[j, i] <- rnorm(1, mTheta[j, ] %*% dfX[i, ], nSigma_2[j])
            }
      }
      mYnew
}
#'
#' @param nTrue_y
#' @param ma_predY
#'
#' @return
#' @export
#'
#' @examples
pred_new_mse <- function(nTrue_y, ma_predY) {
      sum((apply(ma_predY, 2, mean) - nTrue_y)^2) / length(nTrue_y)
}
