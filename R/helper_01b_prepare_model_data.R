#' Title
#'
#' @param nP
#' @param nLag
#' @param .bEndpoint
#'
#' @return
#' @export
#'
#' @examples
prep_model_lag_matrix <- function(nP, nLag, .bEndpoint = FALSE) {

      if (.bEndpoint) {

            mQ <- matrix(1, nrow = 2, ncol = nLag)

            for (i in 1:2) {
                  # mQ[i, ] <- i * (nLag - 1)^(i + 1) - (i + 1) * (nLag - 1)^(i) + (0:(nLag-1))^(i + 1)
                  mQ[i, ] <- i * (nLag - 1)^(i + 1) - (i + 1) * (nLag - 1)^(i) + (0:(nLag-1))^(i + 1)

            }
            # mQ[1, 1] <- (nLag-1)^2
            # mQ[1, 2] <- 5^2 + 1 - 2 * 5
            # mQ[1, 3] <- 5^2 + 2^2 - 2 * 5
            # mQ[1, 4] <- 5^2 + 3^2 - 2 * 5
            # mQ[1, 5] <- 5^2 + 4^2 - 2 * 5
            # mQ[1, 6] <- 5^2 + 5^2 - 2 * 5^2

            mQ[, 1] <- c((nLag-1)^2, 2 * (nLag-1)^3)
            mQ[, nLag] <- 0

      } else {

            mQ <- matrix(1, nrow = nP + 1, ncol = nLag)

            for (i in 1:nP) {
                  mQ[i + 1, ] <- (0:(nLag-1))^i
            }

      }

      mQ

}


#' Title
#'
#' @param degree
#' @param a
#' @param b
#' @param jmax
#' @param X
#'
#' @return
#' @export
#'
#' @examples
prep_model_lag_matrix_legendre <- function(degree,a=0,b=1,jmax){
      # Notes:
      #   References:
      #   H. J. Weber, G. B. Arfken, Essential Mathematical Methods for Physicists,
      #   Elsevier Academic Press, San Diego, USA, 2004.
      #   Translated from matlab function lb to R. 2019-03-01, Jonas Striaukas (jonas.striaukas@gmail.com)
      X <- seq(0,1,length.out=jmax)
      n <- length(X)
      P <- matrix(1,nrow=n,ncol=degree+2)
      Psi <- matrix(1,nrow=n,ncol=degree+1)/ sqrt(b-a)
      P[, 2] <-  2*X/(b-a) - ((b+a)/(b-a))
      if(degree>0){
            for (i in 1:degree){
                  P[, i+2]   <- ((2*i+1)/(i+1)) * P[, 2]*P[, i+1] - i/(i+1) * P[, i]
                  Psi[, i+1] <- sqrt((2*i + 1) / (b-a)) %*% P[, i+1]
            }
      }
      t(Psi)
}

#' Title
#'
#' @param mX
#' @param nP
#' @param nLag
#'
#' @return
#' @export
#'
#' @examples
prep_model_input_data <- function(mX, nP, nLag, .LagMatrix = "almon", .bEndpoint = FALSE) {

      sLagMatrix <- match.arg(.LagMatrix, c("almon", "legendre"))

      if (sLagMatrix == "almon") {
            mQ <- prep_model_lag_matrix(nP, nLag, .bEndpoint)
      } else if(sLagMatrix == "legendre") {
            mQ <- prep_model_lag_matrix_legendre(nP, 0, 1, nLag)
      }

      mX %*% t(mQ)
      # t(mQ %*% t(mX))

}



#' Title
#'
#' @param lAlign
#' @param vLag
#'
#' @return
#' @export
#'
#' @examples
create_predictor_lag_poly <- function(lAlign, vLag = c(6, 12, 25), vP = c(2, 4, 6), .sPolyMatrix = "almon", .bEndpoint = FALSE) {

      lOut <- vector("list", length(lAlign))

      if (!(.bEndpoint)) {
            for (i in seq_along(lAlign)) {
                  lOut[[i]] <- lapply(lAlign[[i]], prep_model_input_data, vP[[i]], vLag[[i]], .sPolyMatrix)
                  lOut[[i]] <- do.call("cbind", lOut[[i]])
            }
      } else {
            for (i in seq_along(lAlign)) {
                  lOut[[i]] <- lapply(lAlign[[i]], prep_model_input_data, vP[[i]], vLag[[i]], .bEndpoint)
                  lOut[[i]] <- do.call("cbind", lOut[[i]])
            }
      }

      lOut
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
create_model_input <- function(maY, dfDataX, nG, nGroupSize, nTrain) {

   vY <- as.vector(maY)
   nY_full <- length(vY)

   nY_test <- nY - nTrain

   y_train <- vY[1:nTrain]
   y_test <- vY[(nTrain + 1):nY]

   x_train <- dfDataX[1:nTrain, ]
   x_test <- dfDataX[(nTrain + 1):nY, ]

   lData_stan_raw <- list(
      nY_train = nTrain,
      nY_test = nY_test,
      nZ = ncol(dfDataX),
      nG = nG,
      y_train = y_train,
      y_test = y_test,
      x_train = x_train,
      x_test = x_test,
      gSize = nGroupSize,
      gInd = rep(1:nG, nGroupSize),
      pr_sigma = c(2, 0.1),
      pr_lambda = c(2, 0.1)
   )

   return(lData_stan_raw)
}
