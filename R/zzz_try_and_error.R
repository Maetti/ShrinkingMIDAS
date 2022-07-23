tryAndError_newLoad <- function() {

      ## check if data generating process works as intendend
      mX <- mpredX_align[[1]]
      vLagPoly <- dgp_pred_beta_lag(nLag, nT1, nT2)

      s1 <- mX %*% vLagPoly
      s1[2]
      sum(mX[2, ] * vLagPoly)

}
