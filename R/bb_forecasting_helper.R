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
