testthat::test_that("exponential almon", {


      exp_almon<-function(k,K,w1,w2){
            j<-1:K
            num <- exp(w1*k + w2*k^2)
            den <- sum( exp(w1*j + w2*j^2) )
            exp_func <- num/den
            return(exp_func)
      }


      expAlmon <- dgp_exp_almon_lag(6, 1, -0.5)
      expAlmon2 <- exp_almon(1:6, 6, 1, -0.5)

      testthat::expect_equal(expAlmon, expAlmon2)
})


testthat::test_check("midas optim", {


      library(midasr)

      ## fake data
      set.seed(1001)
      n <- 200
      trend <- 1:n
      x <- rnorm(3 * n)
      z <- rnorm(3 * n)
      fn_x <- midasr::nealmon(p = c(1, -0.5), d = 6)
      fn_z <- midasr::nealmon(p = c(0.5, -1), d = 6)
      y <- midasr::mls(x, 0:5, 3) %*% fn_x + midasr::mls(z, 0:5, 3) %*% fn_z + rnorm(n)

      ## midas package
      eq_r <- midas_r(y ~ mls(x, 0:5, 3, nealmon) + mls(z, 0:5, 3, nealmon) - 1, start = list(x = c(1, -0.5), z = c(0.5, -1)))

      eq_r2 <- midas_r(y ~ mls(x, 0:5, 3, nealmon) + mls(z, 0:5, 3, nealmon) - 1,
                      start = list(x = c(1, -0.5), z = c(0.5, -1)),
                      Ofunction = "optim", method = "Nelder-Mead")

      ## own midas optimization
      opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD",
                   "xtol_rel" = 1.0e-8)

      vTheta <- c(1, -0.5, 0.5, -1)
      X <- cbind(midasr::mls(x, 0:5, 3), mls(z, 0:5, 3, nealmon))
      res <- nloptr::nloptr(x0 = vTheta,
                            eval_f = midas_optim,
                            opts = opts)


      eval_grad_f2 <- function(x) {
         numDeriv::grad(func = midas_optim, x = x)
      }

      optim(fn = midas_optim, par = vTheta, method = "Nelder-Mead")



      opts <- list("algorithm" = "NLOPT_LN_BOBYQA",
                   "xtol_rel" = 1.0e-8)

      res2 <- nloptr::nloptr(x0 = vTheta,
                            eval_f = midas_optim,
                            opts = opts)

      formula(y ~ rhs(p))

      args <- list()
      nls_formula <- formula(y ~ midas_rhs(vTheta))
      args$formula <- nls_formula
      args$start$vTheta <- vTheta

      opt <- try(do.call("nls", args), silent = TRUE)


      y2 <- y[-1, , drop = FALSE]
      X2 <- X[-1, ]
      nls(formula = y2 ~ midas_rhs_new(vTheta) - 1, start = list(vTheta = vTheta))

      library("minpack.lm")
      curve.nlslrc = nlsLM(formula = y2 ~ midas_rhs_new(vTheta) - 1, start = list(vTheta = vTheta))


      midas_rhs_new <- function(vTheta) {
         coeffs <- midas_coef(vTheta)
         X2 %*% coeffs
      }


      vTheta_1_1 <- 1
      vTheta_1_2 <- 0.5

      vTheta_2_1 <- 1
      vTheta_2_2 <- 0.5

      nlsLM(formula = y2 ~ X2[, 1:6] %*% dgp_exp_almon_lag(6, vTheta_1_1, vTheta_1_2) +
                          X2[, 7:12] %*% dgp_exp_almon_lag(6, vTheta_2_1, vTheta_2_2) - 1,
            start = list(vTheta_1_1 = 3, vTheta_1_2 = 0.5, vTheta_2_1 = 2, vTheta_2_2 = 0.5))

})
