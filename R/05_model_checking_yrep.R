# ## connection to files
# dirDGP <- here::here("inst/data/01_dgp")
#
# dirExtract <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_extract/")
# l1 <- lapply(dirExtract, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
# arrExtract <- arrow::open_dataset(l1)
#
# sIteration <- seq(0, 16000, 8)
# sIteration[1] <- 1
#
# tblYRep <-
#       arrExtract %>%
#             dplyr::filter(parameter_1 == "y_rep",
#                           # simulation == 1,
#                           iteration %in% sIteration) %>%
#             dplyr::collect() %>%
#             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
#             dplyr::select(-parameter_1, parameter_2)
#
#
# tblYPred <-
#    arrExtract %>%
#       dplyr::filter(parameter_1 == "y_pred",
#                     iteration %in% sIteration) %>%
#       dplyr::collect() %>%
#       dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
#       dplyr::select(-parameter_1, parameter_2)
#
#
# tblLogLik <-
#    arrExtract %>%
#       dplyr::filter(parameter_1 == "log_lik",
#                     iteration %in% sIteration) %>%
#       dplyr::collect() %>%
#       dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
#       dplyr::select(-parameter_1, parameter_2)
#
#
# ## all true Y values
# lYTrue_raw <- vector("list", length(nSim))
# fSims <- dir(dirDGP)
# for (i in seq_along(nSim)) {
#    sFile <- dir(paste0(dirDGP, "/", fSims[[i]], "/raw_input"))
#    s1 <- readRDS(paste0(dirDGP, "/", fSims[[i]], "/raw_input/", sFile))
#    lYTrue_raw[[i]] <- data.frame(simulation = i,
#                                  key = seq(1, s1$nY_train +  s1$nY_test),
#                                  value = c(s1$y_train, s1$y_test))
# }
# tblYTrue <- dplyr::bind_rows(lYTrue_raw) %>% dplyr::as_tibble()
# tblYTrain <- tblYTrue %>% dplyr::filter(key <= 150)
# tblYTest <- tblYTrue %>% dplyr::filter(key > 150)
# tblYTest$key <- tblYTest$key - 150
# rm(tblYTrue)


#' Title
#'
#' @param tblInput
#'
#' @return
#' @export
#'
#' @examples
foo_check_yrep_stat <- function(tblInput) {
      tblInput %>%
            dplyr::summarise(
                  mean = mean(value),
                  sd = sd(value),
                  min = min(value),
                  max = max(value),
                  q75 = quantile(value, 0.75),
                  q25 = quantile(value, 0.25),
                  q975 = quantile(value, 0.975),
                  q025 = quantile(value, 0.025)
            )
}

#' Title
#'
#' @param tblYRep
#' @param lYtrue
#'
#' @return
#' @export
#'
#' @examples
md_check_yrep_stat <- function(tblYRep) {
      tblYrep_stat <-
            tblYRep %>%
                  dplyr::group_by(model, simulation, iteration) %>%
                  foo_check_yrep_stat() %>%
                  dplyr::ungroup()
}

#' Title
#'
#' @param lYtrueTrain
#'
#' @return
#' @export
#'
#' @examples
md_check_ytrain_stat <- function(tblYTrain) {

      tblYTrain %>%
            dplyr::group_by(simulation) %>%
            foo_check_yrep_stat() %>%
            dplyr::ungroup()


}

#' Title
#'
#' This functions gives a summary of p-values for mean, sd, max, and min of the true y value.
#' This is to check how good the models can capture the true underlying process
#'
#' @param tblYRep
#' @param lYtrueTrain
#'
#' @return
#' @export
#'
#' @examples
md_check_yrep_pval <- function(tblYRep, tblYTrain) {

      tblStat1 <- md_check_yrep_stat(tblYRep)
      tblStat2 <- md_check_ytrain_stat(tblYTrain)

      colnames(tblStat2)[-1] <- paste0("ytrue_", colnames(tblStat2)[-1])

      tblAll <- dplyr::left_join(tblStat1, tblStat2, by = "simulation")

      tblAll %>%
            dplyr::group_by(model, simulation) %>%
                  dplyr::summarise(
                        p_mean = sum(ytrue_mean <= mean) / n(),
                        p_sd = sum(ytrue_sd <= sd) / n(),
                        p_min = sum(ytrue_min <= min) / n(),
                        p_max = sum(ytrue_max <= max) / n()
                  ) %>%
            dplyr::ungroup() %>%
            tidyr::pivot_longer(names_to = "stat", values_to = "value", -c(1, 2))
}

#' Title
#'
#' @param tblInput
#'
#' @return
#' @export
#'
#' @examples
md_check_yrep_pval_plot <- function(tblInput, sTitle = "Basic Statistics") {

   tblInput %>%
      dplyr::mutate(simulation = factor(as.numeric(simulation),
                                        levels = sort(unique(as.numeric(simulation))))) %>%
      ggplot2::ggplot(., mapping = ggplot2::aes(x = simulation, y = value, color = model)) +
      ggplot2::geom_line(ggplot2::aes(group = model)) +
      ggplot2::ggtitle(sTitle) +
      ggplot2::theme_minimal() +
      ggplot2::facet_wrap(~stat, ncol = 2)

}

#' Title
#'
#' @param tblYrep
#'
#' @return
#' @export
#'
#' @examples
md_check_yrep_95_covered <- function(tblInput, tblYTrain) {

      tblStat <-
         tblInput %>%
            dplyr::group_by(model, simulation, parameter_2) %>%
            foo_check_yrep_stat() %>%
            dplyr::ungroup() %>%
            dplyr::select(model, simulation, parameter_2, q975, q025) %>%
            dplyr::mutate(parameter_2 = as.numeric(parameter_2))

      tblYPrep <-
         dplyr::left_join(tblStat, tblYTrain, by = c("simulation", "parameter_2" = "key")) %>%
         md_foo_yrep_covered()

      tblYPrep %>%
            dplyr::group_by(model, simulation) %>%
            dplyr::summarise(covered = sum(covered) / n()) %>%
            dplyr::mutate(simulation = as.numeric(simulation))





}

#' Title
#'
#' @param tblYCovered
#'
#' @return
#' @export
#'
#' @examples
md_check_yrep_95_covered_plot <- function(tblYCovered) {

   tblYCovered %>%
      ggplot2::ggplot(ggplot2::aes(x = simulation, y = covered, color = model)) +
         ggplot2::geom_point() +
         ggplot2::geom_line() +
         ggplot2::scale_x_continuous(
            breaks = seq(min(tblYCovered$simulation), max(tblYCovered$simulation), by = 1)) +
         ggplot2::ylab("Percentage of true Y within 95% intervall") +
         ggplot2::labs(
            title = "Percentage of true Y observation in 95% prediction Intervall"
         )

   # ggplot2::ggplot(ggplot2::aes(x = simulation, y = Freq, fill = bolSim)) +
   # ggplot2::geom_bar(position="fill", stat="identity") +
   # ggplot2::facet_wrap(~model)
}

#' Title
#'
#' @param tblInput
#'
#' @return
#' @export
#'
#' @examples
md_foo_yrep_covered <- function(tblInput) {
   tblInput %>%
      dplyr::mutate(
         covered = dplyr::case_when(
            value >= q025 & value <= q975 ~ TRUE,
            TRUE ~ FALSE
         )
      )
}


#' Title
#'
#' @param tblYRep
#' @param tblYTrain
#' @param sSim
#' @param sStat
#'
#' @return
#' @export
#'
#' @examples
md_yrep_stat_simulation_plot <- function(tblInput, tblYTrain,
                                    sSim = 1, sStat = "mean",
                                    .fInput = "stat_group",
                                    sTitle = "Basic stat") {

   .fInput <- match.arg(.fInput, choices = c("stat_group", "ribbon", "intervals", "violin", "freqpoly", "scatter"))

   lYrep1 <-
      tblInput %>%
         dplyr::filter(simulation == sSim) %>%
         dplyr::select(iteration, value, model, parameter)

   lYrep1 <-
      lYrep1 %>%
         split(., lYrep1$model) %>%
         purrr::map(., function(x) {
            x %>%
               dplyr::select(-model) %>%
               tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
               dplyr::select(-iteration) %>%
               as.matrix() %>%
               unname()
         })


   vGroup <- rep(names(lYrep1), purrr::map_dbl(lYrep1, ncol))

   tbl1 <- lYrep1 %>% purrr::reduce(., dplyr::bind_cols) %>% as.matrix() %>% unname()

   s1 <- tblYTrain %>% dplyr::filter(simulation == sSim)
   s2 <- rep(s1$value, length(unique(tblInput$model)))

   .f <-switch (.fInput,
      "stat_group" = bayesplot::ppc_stat_grouped,
      "ribbon" = bayesplot::ppc_ribbon_grouped,
      "intervals" = bayesplot::ppc_intervals_grouped,
      "scatter" = bayesplot::ppc_scatter_avg_grouped
   )

   if (.fInput == "stat_group") {
      pOut <- .f(s2, tbl1, group = vGroup, stat = sStat)
   } else {
      pOut <- .f(s2, tbl1, group = vGroup)
   }

   pOut +
      ggplot2::ggtitle(sTitle)

}


#' Title
#'
#' @param tblInput
#' @param tblYTrain
#' @param nSim
#'
#' @return
#' @export
#'
#' @examples
md_yrep_stat_simulation_foo <- function(tblInput, tblYTrain, nSim = 1) {


   lYrep1 <-
      tblInput %>%
         dplyr::filter(simulation == sSim) %>%
         dplyr::select(iteration, value, model, parameter)

   lYrep1 <-
      lYrep1 %>%
            split(., lYrep1$model) %>%
            purrr::map(., function(x) {
               x %>%
                  dplyr::select(-model) %>%
                  tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
                  dplyr::select(-iteration) %>%
                  as.matrix() %>%
                  unname()
            })

   s1 <- tblYTrain %>% dplyr::filter(simulation == sSim) %>% dplyr::pull(value)


   lOut <- vector("list", length(lYrep1))

   for (i in seq_along(lOut)) {
      lOut[[i]] <- bayesplot::ppc_dens_overlay(s1, lYrep1[[i]])
   }

   do.call("grid.arrange", lOut, ncol = 3)

}



#' Title
#'
#' @param tblYPred
#' @param tblYTest
#'
#' @return
#' @export
#'
#' @examples
md_ypred_mse <- function(tblYPred, tblYTest) {

   tbl1 <- tblYPred %>% dplyr::mutate(parameter_2 = as.numeric(parameter_2))
   tbl2 <- tblYTest %>% dplyr::rename(y_test = value)

   dplyr::left_join(tbl1, tbl2, by = c("simulation", "parameter_2" = "key")) %>%
         dplyr::group_by(model, simulation) %>%
         dplyr::summarise(
            value = sum((value - y_test)^2)
         ) %>%
         dplyr::ungroup()

}


#' Title
#'
#' @param tblInput
#'
#' @return
#' @export
#'
#' @examples
md_ypre_loo_prep_foo1 <- function(tblInput) {

   split(tblInput, tblInput$simulation) %>%
      purrr::map(., function(x) {
         split(x, x$model) %>%
            purrr::map(., function(y) {
               y %>%
                  dplyr::arrange(iteration) %>%
                  dplyr::select(iteration, parameter_2, value) %>%
                  tidyr::pivot_wider(names_from = parameter_2, values_from = value) %>%
                  dplyr::select(-iteration) %>%
                  as.matrix()
         })
      })

}


#' Title
#'
#' @param tblInput
#'
#' @return
#' @export
#'
#' @examples
md_ypred_loo_compute_foo2 <- function(tblInput) {

   lLogLik <- md_ypre_loo_prep_foo1(tblInput)

   lLoo <- lLogLik

   for (i in seq_along(lLogLik)) {
      print(paste0("Current run: ", i, "/", length(lLogLik)))
      for (j in seq_along(lLogLik[[i]])) {
         lLoo[[i]][[j]] <- loo(lLogLik[[i]][[j]], cores = 4)
      }
   }

   lLoo

}


#' Title
#'
#' @param lLoo
#'
#' @return
#' @export
#'
#' @examples
md_ypred_loo_clean <- function(tblInput) {

   lOut <- lLoo <- md_ypred_loo_compute_foo2(tblInput)


   for (i in seq_along(lLoo)) {
      sNames <- names(lLoo[[i]])
      for (j in seq_along(lLoo[[i]])) {
         x <- as.data.frame(lLoo[[i]][[j]]$estimates)
         x$key <- rownames(x)
         x$model <- sNames[[j]]
         x$simulation <- i
         lOut[[i]][[j]] <- dplyr::as_tibble(x)
      }

      lOut[[i]] <- dplyr::bind_rows(lOut[[i]])
   }

   tblOut <- dplyr::bind_rows(lOut)
   colnames(tblOut) <- tolower(colnames(tblOut))

   tblOut %>% dplyr::select(model, simulation, key, value = estimate, dplyr::everything())

}




#' Title
#'
#' @param tblInput
#' @param sTitle
#'
#' @return
#' @export
#'
#' @examples
md_check_general_sim_plot <- function(tblInput,
                                      sTitle = "Basic Statistics",
                                      sYtitle, sXtitle) {

   gpOut <-
         tblInput %>%
            dplyr::mutate(simulation = factor(as.numeric(simulation),
                                              levels = sort(unique(as.numeric(simulation))))) %>%
            ggplot2::ggplot(., mapping = ggplot2::aes(x = simulation, y = value, color = model)) +
            # ggplot2::geom_line(ggplot2::aes(group = model)) +
            ggplot2::geom_point() +
            ggplot2::ggtitle(sTitle) +
            ggplot2::theme_minimal() +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 22),
                           axis.title = ggplot2::element_text(size = 16),
                           axis.text = ggplot2::element_text(size = 14),
                           legend.title = ggplot2::element_text(size = 16),
                           legend.text = ggplot2::element_text(size = 14))

   if (!missing(sYtitle)) {
      gpOut <- gpOut + ggplot2::ylab(sYtitle)
   }

   if (!missing(sXtitle)) {
      gpOut <- gpOut + ggplot2::xlab(sXtitle)
   }

   gpOut
}


#' Title
#'
#' @param tblInput
#'
#' @return
#' @export
#'
#' @examples
md_ypred_stat_foo <- function(tblInput) {
   tblInput %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(
         mean = mean(value),
         min = min(value),
         max = max(value)
      ) %>%
      dplyr::ungroup()
}


#' Title
#'
#' @param tblYPred
#' @param tblYTest
#' @param sSim
#'
#' @return
#' @export
#'
#' @examples
md_ypred_pred_interval_plot <- function(tblYPred, tblYTest, nSim = 1) {


   tb1 <-
      tblYPred %>%
         dplyr::filter(simulation == nSim) %>%
         dplyr::group_by(model, parameter_2) %>%
         dplyr::summarise(
            mean = mean(value),
            q005 = quantile(value, 0.05),
            q95 = quantile(value, 0.95)
         ) %>%
         dplyr::ungroup() %>%
         dplyr::mutate(parameter_2 = as.numeric(parameter_2)) %>%
         dplyr::arrange(model, parameter_2)

   tb2 <- tblYTest %>% dplyr::filter(simulation == nSim) %>% dplyr::rename(ytrue = value)

   tbOut <- dplyr::left_join(tb1, tb2, by = c("parameter_2" = "key"))


   ggplot2::ggplot(tbOut, ggplot2::aes(x = parameter_2, y = mean)) +
      ggplot2::geom_smooth(ggplot2::aes(ymin = q005, ymax = q95), stat = "identity", size = 0.5) +
      ggplot2::geom_point(ggplot2::aes(y = ytrue)) +
      ggplot2::facet_wrap(~model) +
      ggplot2::ggtitle("90% Predicition Interval vs Observered Data") +
      ggplot2::labs(subtitle = "Mean (blue) and 90% predictive intervals (gray) vs. observed data (black)") +
      ggplot2::xlab("Time") +
      ggplot2::ylab("Obsereved Y from Train data") +
      ggplot2_theme_custom

}



#' Title
#'
#' @param gp
#'
#' @return
#' @export
#'
#' @examples
ggplot2_theme_custom <-
   ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 22),
                     axis.title = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size = 14),
                     legend.title = ggplot2::element_text(size = 16),
                     legend.text = ggplot2::element_text(size = 14))
