# library("loo")
# library("brms")
# library("bayesplot")
# library("ggplot2")
# color_scheme_set("brightblue")
# theme_set(theme_default())
#
# CHAINS <- 4
# SEED <- 5838296
# set.seed(SEED)
#
#
# N <- length(LakeHuron)
# df <- data.frame(
#       y = as.numeric(LakeHuron),
#       year = as.numeric(time(LakeHuron)),
#       time = 1:N
# )
#
# # save plot labels to reuse them
# plot_labs <- labs(
#       y = "Water Level (ft)",
#       x = "Year",
#       title = "Water Level in Lake Huron (1875-1972)"
# )
#
# ggplot(df, aes(x = year, y = y)) +
#       geom_point(size = 1) +
#       plot_labs
#
# ## fit mode
# control <- list(adapt_delta = 0.99)
# fit <- brm(
#       y ~ ar(time, p = 4),
#       data = df,
#       prior = prior(normal(0, 0.5), class = "ar"),
#       control = control,
#       seed = SEED,
#       chains = CHAINS
# )
#
#
#
#
# ## approximate leave future out CV
# L <- 20
# k_thres <- 0.7
# approx_elpds_1sap <- rep(NA, N)
#
# # initialize the process for i = L
# past <- 1:L
# oos <- L + 1
# df_past <- df[past, , drop = FALSE]
# df_oos <- df[c(past, oos), , drop = FALSE]
# fit_past <- update(fit, newdata = df_past, recompile = FALSE)
# loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
# approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])
#
# # iterate over i > L
# i_refit <- L
# refits <- L
# ks <- NULL
# for (i in (L + 1):(N - 1)) {
#       past <- 1:i
#       oos <- i + 1
#       df_past <- df[past, , drop = FALSE]
#       df_oos <- df[c(past, oos), , drop = FALSE]
#       loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
#
#       logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
#       psis_obj <- suppressWarnings(psis(logratio))
#       k <- pareto_k_values(psis_obj)
#       ks <- c(ks, k)
#       if (k > k_thres) {
#             # refit the model based on the first i observations
#             i_refit <- i
#             refits <- c(refits, i)
#             fit_past <- update(fit_past, newdata = df_past, recompile = FALSE)
#             loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
#             approx_elpds_1sap[i + 1] <- log_mean_exp(loglik[, oos])
#       } else {
#             lw <- weights(psis_obj, normalize = TRUE)[, 1]
#             approx_elpds_1sap[i + 1] <- log_sum_exp(lw + loglik[, oos])
#       }
# }
#
#
# ## helper function
# # some helper functions we'll use throughout
#
# # more stable than log(sum(exp(x)))
# log_sum_exp <- function(x) {
#       max_x <- max(x)
#       max_x + log(sum(exp(x - max_x)))
# }
#
# # more stable than log(mean(exp(x)))
# log_mean_exp <- function(x) {
#       log_sum_exp(x) - log(length(x))
# }
#
# # compute log of raw importance ratios
# # sums over observations *not* over posterior samples
# sum_log_ratios <- function(loglik, ids = NULL) {
#       if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
#       rowSums(loglik)
# }
#
# # for printing comparisons later
# rbind_print <- function(...) {
#       round(rbind(...), digits = 2)
# }
