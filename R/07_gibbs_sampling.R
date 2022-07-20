#' Title
#'
#' @param vY
#' @param dfX
#'
#' @return
#' @export
#'
#' @examples
gibbs_horseshoe <- function(vY, mX, nG, nGroupSize, vY_new, mX_new, nRuns = 100000) {

      checkmate::assertVector(vY, len = nrow(mX))
      checkmate::assertMatrix(mX, nrows = len(vY))


      nN <- length(vY)
      nP <- ncol(mX)

      vIndGroup <- rep(1:length(nGroupSize), nGroupSize)

      #   ____________________________________________________________________________
      #   Variables                                                               ####


      ## posterior variables
      mBeta <- matrix(nrow = nRuns + 1, ncol = nP)
      mSigma <- matrix(nrow = nRuns + 1, ncol = 1)

      mLambda <- matrix(nrow = nRuns + 1, ncol = nP)
      mDelta <- matrix(nrow = nRuns + 1, ncol = nP)
      mTau <- matrix(nrow = nRuns + 1, ncol = 1)

      ## posterior prediction sampling
      mYRep <- matrix(nrow = nRuns + 1, ncol = nN)

      ## posterior predictive sampling
      mYRep <- matrix(nrow = nRuns + 1, ncol = length(vY_new))

      ## log likelihood
      mLogLik <- matrix(nrow = nRuns + 1, ncol = 1)

      #   ____________________________________________________________________________
      #   Inital values                                                           ####

      mBeta[1, ] <- rnorm(n = nP, mean = 0, sd = 1)
      mSigma[1, ] <- sqrt(t(vY - mX %*% initBeta) %*% (vY - mX %*% initBeta) / nN)

      mLambda[1, ] <- runif(nP, 0, 2)
      mDelta[1, ] <- runif(nG, 0, 2)[vIndGroup]
      mTau[1, ] <- runif(1, 0, 2)



      #   ____________________________________________________________________________
      #   Sampling                                                                ####

      start_time_sampling <- Sys.time()

      for (i in 2:nRuns) {

            print(i)

            ## previous run values
            vBeta <- mBeta[i - 1, ]
            nSigma <- mSigma[i - 1, ]

            diagLamba <- diag(mLambda[i - 1, ])
            diagDelta <- diag(mDelta[i - 1, ])
            nTau <- mTau[i - 1, ]

            invA <- gibbs_helper_inv_tau_delta_lambda(nTau, diagDelta, diagLamba)

            ## beta and sigma
            mBeta[i, ] <- gibbs_beta_coef(vY, mX, nSigma, invA)
            mSigma[i, ] <- gibbs_sigma(vY, mX, mBeta[i, ], invA)

            ## individual shrinkage
            for (j in 1:nP) {

                  nLambda <- mLambda[i - 1, j]
                  nDelta <- mDelta[i - 1, j]
                  mLambda[i, j] <- gibbs_horseshoe_lambda(nLambda, nBeta = mBeta[i, j], nSigma = mSigma[i, ], nTau = nTau, nDelta = nDelta)

            }

            ## group shrinkage
            for (j in 1:nP) {

                  ## needed var
                  nDelta <- mDelta[i - 1, j]

                  ## group level
                  nGroupBelong <- vIndGroup[j]
                  indGroup <- vIndGroup == nGroupBelong

                  vLambda <- mLambda[i, indGroup]
                  vBeta <- mBeta[i, indGroup]

                  mDelta[i, j] <- gibbs_horseshoe_delta(nDelta = nDelta, vBeta = vBeta, nSigma = mSigma[i, ], nTau = nTau, vLambda = vLambda, nGroupSize = sum(indGroup))
            }

            ## global shrinkage
            mTau[i, ] <- gibbs_tau(nTau = mTau[i - 1, ], mX, vBeta = mBeta[i, ], nSigma = mSigma[i, ],
                                   diagDelta = diag(mDelta[i, ]), diagLambda = diag(mLambda[i, ]))

      }

      end_time_sampling <- Sys.time()

      end_time_sampling - start_time_sampling


}


#' Title
#'
#' shrinkage on individual level
#'
#' @param nBeta
#' @param nSigma
#' @param nTau
#' @param nDelta
#'
#' @return
#' @export
#'
#' @examples
gibbs_horseshoe_lambda <- function(nLambda, nBeta, nSigma, nTau, nDelta) {

      nAlpha <- (2 * nSigma * nTau * nDelta) / nBeta^2
      nC <- invgamma::rinvgamma(1, shape = 1, rate = 1 / nLambda + 1)

      invgamma::rinvgamma(1, shape = 1, rate =  1 / nAlpha + 1 / nC)
}

#' Title
#'
#' shrinkage on group level
#'
#' @return
#' @export
#'
#' @examples
gibbs_horseshoe_delta <- function(nDelta, vBeta, nSigma, nTau, vLambda, nGroupSize) {

      nT <- invgamma::rinvgamma(1, 1 / nDelta + 1)

      shape <- (nGroupSize + 1) / 2
      rate <- 1 / (2 * nSigma * nTau) * sum(vBeta^2 / vLambda) + 1 / nT

      invgamma::rinvgamma(1, shape = shape, rate = rate)

}

#' Title
#'
#'
#' @param vY
#' @param mX
#' @param nTau
#' @param mDelta
#' @param mLambda
#'
#' @return
#' @export
#'
#' @examples
gibbs_beta_coef <- function(vY, mX, nSigma, invShrink) {

      A <- t(mX) %*% mX + invShrink
      A_inv <- solve(A)

      nMean <- A_inv %*% t(mX) %*% vY
      mSigma <- nSigma * A_inv

      MASS::mvrnorm(n = 1, nMean, mSigma)
}

#' Title
#'
#' @param n
#' @param p
#' @param vY
#' @param mX
#' @param vBeta
#' @param nTau
#' @param mDelta
#' @param mLambda
#'
#' @return
#' @export
#'
#' @examples
gibbs_sigma <- function(vY, mX, vBeta, invShrink) {

      n <- vY
      p <- nrow(mX)

      shape <- (n - 1 + p) / 2
      rate <- t(vY - mX %*% vBeta) %*% (vY - mX %*% vBeta) + t(vBeta) %*% invShrink %*% vBeta

      invgamma::rinvgamma(n = 1, shape = shape, rate = rate / 2)

}


#' Title
#'
#' @param mX
#' @param vBeta
#' @param nSigma
#' @param mDelta
#' @param mLambda
#' @param nV
#'
#' @return
#' @export
#'
#' @examples
gibbs_tau <- function(nTau, mX, vBeta, nSigma, diagDelta, diagLambda) {

      p <- ncol(mX)

      nV <- invgamma::rinvgamma(1, shape = 1, rate = 1/ nTau + 1)

      shape <- (p + 1) / 2
      rate <- (t(vBeta) %*% solve(diagDelta %*% diagLambda) %*% vBeta) / (2 * nSigma) + 1 / nV

      invgamma::rinvgamma(n = 1, shape = shape, rate = rate)
}

#' Title
#'
#' @param nTau
#'
#' @return
#' @export
#'
#' @examples
gibbs_helper_v <- function(nTau) {

      invgamma::rinvgamma(1, shape = 1, rate = 1 / nTau + 1)

}

#' Title
#'
#' @param nTau
#' @param mDelta
#' @param mLambda
#'
#' @return
#' @export
#'
#' @examples
gibbs_helper_inv_tau_delta_lambda <- function(nTau, mDelta, mLambda) {
      solve(nTau * mDelta %*% mLambda)
}
