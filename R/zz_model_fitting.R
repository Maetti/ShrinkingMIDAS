fit_stan <-  function(sFile, sName, nIter = 2000, ...) {

                  ## fit the stan model
                  obFit <- rstan::stan(file = sFile, model_name = sName, data = self$lData, iter = nIter, ...)

                  ## extract usefull information
                  self$stan_model <- obFit
                  self$list_of_draws <- rstan::extract(obFit)
                  self$fit_summary <- summary(obFit)

                  ## compute beta coefficient
                  cat("Computing Beta Coefficient")

                  nX <- self$lData$nG
                  dfDataX <- self$lData$x
                  nP <- self$lData$gSize[[1]]
                  dfTheta_sim <- self$list_of_draws$theta

                  ## unstandardized theat
                  sdX <- attr(dfDataX, "scaled:scale")
                  dfTheta_sim_unstd <- t(t(dfTheta_sim) / sdX)

                  ## get weights and coefficient
                  lWeights <- vector("list", length = nX)
                  maWeigths <- matrix(0, nrow = nrow(dfTheta_sim), ncol = ncol(mQ))
                  maCoef <- matrix(0, nrow = nrow(dfTheta_sim), ncol = nX)
                  vInd <- c(seq(1, ncol(dfDataX), nP), ncol(dfDataX) + 1)

                  for (i in 1:nX) {
                        vB <- dfTheta_sim_unstd[, vInd[i]:(vInd[i+1] - 1)]
                        maWeigths <- t(apply(vB, 1, function(x, mQ) {x %*% mQ}, mQ = mQ))
                        maCoef[, i] <- apply(maWeigths, 1, sum)
                        lWeights[[i]] <- maWeigths / maCoef[, i]
                        colnames(lWeights[[i]]) <- paste0("weight_", 1:ncol(mQ))
                  }


                  colnames(maCoef) <- paste0("Beta_", 1:nX)
                  names(lWeights) <- paste0("covariate_", 1:nX)

                  self$beta_coef <- maCoef
                  self$lWeights <- lWeights
                  self$unstdr_theta <- dfTheta_sim_unstd


                  ## summaries
                  self$summary_beta <- foo_summary(maCoef)
                  self$summary_weight <- lapply(lWeights, foo_summary)
                  names(self$summary_weight) <- names(lWeights)
}

save_stan_fit = function(sPath) {
                  saveRDS(object = self$stan_model, file = sPath)
}

foo_summary = function(maInput) {

                  sColNames <- factor(colnames(maInput), levels = colnames(maInput))

                  maInput %>%
                        as.data.frame() %>%
                        tidyr::pivot_longer(values_to = "val", names_to = "input", cols = 1:ncol(.)) %>%
                        dplyr::mutate(input = factor(input, levels = colnames(maInput))) %>%
                        dplyr::group_by(input) %>%
                        dplyr::summarise(
                              # amount = dplyr::n(),
                              mean = mean(val),
                              # se_mean = sd(val) / amount,
                              sd = sd(val),
                              p_0.025 = quantile(val, probs = 0.025),
                              p_0.25 = quantile(val, probs = 0.25),
                              p_0.5 = quantile(val, probs = 0.5),
                              p_0.75 = quantile(val, probs = 0.75),
                              p_0.975 = quantile(val, probs = 0.975)
                        ) %>%
                        # dplyr::select(-amount) %>%
                        force()

}

pred_new = function() {

                  mY_new <- gen_pred_y(xTest_scale, mTheta, nSigma_2)
                  colnames(mY_new) <- paste0("y_", 1:50)

                  mY_new_summary <-
                        foo_summary(mY_new) %>%
                        dplyr::mutate(true_y = mY_test[, 1]) %>%
                        dplyr::select(input, true_y, dplyr::everything())

                  sum((mY_new_summary$true_y - mY_new_summary$mean)^2)/50

}

gen_pred_y = function(dfX, mTheta, nSigma_2) {

                  mYnew <- matrix(nrow = length(nSigma_2), ncol = nrow(dfX))

                  for (i in 1:nrow(dfX)) {
                        for (j in 1:length(nSigma_2)) {
                              mYnew[j, i] <- rnorm(1, mTheta[j, ] %*% dfX[i, ], nSigma_2[j])
                        }
                  }
                  mYnew
            }



create_plot = function() {
                  ## plots
                  self$plot_intervals <- bayesplot::mcmc_intervals(self$beta_coef)
                  self$plot_areas <- bayesplot::mcmc_areas(
                        self$beta_coef,
                        prob = 0.8, # 80% intervals
                        prob_outer = 0.99, # 99%
                        point_est = "mean"
                  )
                  self$plot_hist <- bayesplot::mcmc_hist(self$beta_coef)
                  self$plot_dens <- bayesplot::mcmc_dens(self$beta_coef)
            }

