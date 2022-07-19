#' Create an Rhat Plot
#'
#' @param sValTo
#' @param sTitle
#' @param tblRhat
#'
#' @return
#' @export
#'
#' @examples
moc_diagnose_plot <- function(tblRhat, sValTo = "Rhat", sTitle = "Theta Parameter") {

      colnames(tblRhat) <- gsub(".*_", "", colnames(tblRhat))
      nYmin <- min(lOut$rhat$value, 0.9)
      nYmax <- max(lOut$rhat$value + 0.5)

      tblRhat %>%
            dplyr::mutate(parameter = factor(as.numeric(gsub(".*_", "", parameter)),
                                             levels = sort(unique(as.numeric(gsub(".*_", "", parameter)))))) %>%
            dplyr::mutate(simulation = factor(as.numeric(simulation),
                                              levels = sort(unique(as.numeric(simulation))))) %>%
            ggplot2::ggplot(., mapping = ggplot2::aes(x = parameter, y = value, color = simulation)) +
            ggplot2::geom_point() +
            ggplot2::ylim(c(nYmin, nYmax)) +
            ggplot2::ggtitle(sTitle) +
            ggplot2::theme_minimal() +
            ggplot2::facet_wrap(~model, ncol = 1)

}



#' Title
#'
#' @param tblConf
#' @param .sType
#'
#' @return
#' @export
#'
#' @examples
moc_tp_fp_plot <- function(tblConf, .sType = c("TP", "FN")) {

   .sType <- match.arg(.sType, choices = c("TP", "FN"))

   if (.sType == "TP") {

      pOut <-
         tblConf %>%
            dplyr::filter(bolTrue == TRUE) %>%
            dplyr::mutate(simulation = factor(simulation, levels = sort(unique(simulation)))) %>%
            ggplot2::ggplot(ggplot2::aes(x = simulation, y = Freq, fill = bolSim)) +
            ggplot2::geom_bar(position="fill", stat="identity") +
            ggplot2::facet_wrap(~model)

   } else {

      pOut <-
         tblConf %>%
            dplyr::filter(bolTrue == FALSE) %>%
            dplyr::mutate(bolTrue = !as.logical(bolTrue)) %>%
            dplyr::mutate(bolSim = !as.logical(bolSim)) %>%
            dplyr::mutate(simulation = factor(simulation, levels = sort(unique(simulation)))) %>%
            ggplot2::ggplot(ggplot2::aes(x = simulation, y = Freq, fill = bolSim)) +
            ggplot2::geom_bar(position="fill", stat="identity") +
            ggplot2::facet_wrap(~model)

   }

   pOut
}


#' Title
#'
#' @param tblBeta
#'
#' @return
#' @export
#'
#' @examples
#' tblBeta <- lOut$beta
#' vBeta <- vBeta
# moc_mse_table <- function(tblBeta, vBeta) {
#
#    tblBetaTrue <- data.frame(
#       parameter = paste0("beta_", seq_along(vBeta)),
#       beta_true = vBeta
#    )
#
#    tblBetaPrep <-
#       tblBeta %>%
#          dplyr::left_join(., tblBetaTrue, by = "parameter") %>%
#          dplyr::group_by(model, simulation, parameter) %>%
#             dplyr::mutate(beta_mean = mean(value)) %>%
#                           # mse_var = (value - beta_mean)^2,
#                           # mse_bias = (beta_mean - beta_true)^2) %>%
#          dplyr::ungroup()
#
#
#    tblVar <-
#       tblBetaPrep %>%
#          dplyr::group_by(model) %>%
#             dplyr::summarise(mse_var = sum((value - beta_mean)^2) / n())
#
#    tblBias <-
#       tblBetaPrep %>%
#          dplyr::distinct(model, simulation, parameter, beta_true, beta_mean) %>%
#          dplyr::group_by(model) %>%
#          dplyr::summarise(mse_bias = sum((beta_mean - beta_true)^2) / n())
#
#    dplyr::left_join(tblVar, tblBias, by = "model") %>%
#       dplyr::mutate(mse = mse_var + mse_bias) %>%
#       dplyr::select(model, mse, mse_var, mse_bias)
#
#
# }


#' Title
#'
#' @param tblConf
#' @param bSplit
#'
#' @return
#' @export
#'
#' @examples
moc_tpr_fpr_tbl <- function(tblConf, bSplit = "sim") {


   bSplit <- match.arg(bSplit, c("sim", "model"))

   tblTPR <-
      tblConf %>%
         dplyr::filter(bolTrue == TRUE) %>%
            dplyr::select(-bolTrue) %>%
            dplyr::group_by(model, simulation) %>%
            dplyr::mutate(tpr = Freq / sum(Freq)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(bolSim == TRUE) %>%
            dplyr::select(model, simulation, tpr)


   tblFPR <-
      tblConf %>%
         dplyr::filter(bolTrue == FALSE) %>%
         dplyr::select(-bolTrue) %>%
         dplyr::group_by(model, simulation) %>%
         dplyr::mutate(fpr = Freq / sum(Freq)) %>%
         dplyr::ungroup() %>%
         dplyr::filter(bolSim == TRUE) %>%
         dplyr::select(model, simulation, fpr)


   tblOut <- dplyr::left_join(tblTPR, tblFPR, by = c("model", "simulation"))

   if (bSplit == "model") {
      tblOut <-
         tblOut %>%
            dplyr::group_by(model) %>%
            dplyr::summarise(tpr = mean(tpr),
                             fpr = mean(fpr)) %>%
            dplyr::ungroup()
   }

   tblOut

}


#' Title
#'
#' @return
#' @export
#'
#' @examples
moc_matthes_corr <- function(tblConf) {

   tblConf %>%
         dplyr::mutate(
            bolSim = as.logical(bolSim),
            bolTrue = as.logical(bolTrue)) %>%
         dplyr::mutate(conf_matrix = dplyr::case_when(
            bolSim & bolTrue ~ "tp",
            !bolSim & !bolTrue ~ "tn",
            bolSim & !bolTrue ~ "fp",
            !bolSim & bolTrue ~ "fn")) %>%
         dplyr::group_by(model, conf_matrix) %>%
         dplyr::summarise(model_freq = sum(Freq)) %>%
         dplyr::ungroup() %>%
         tidyr::pivot_wider(names_from = conf_matrix, values_from = model_freq) %>%
         dplyr::mutate(
            tpr = tp / (tp + fn),
            fpr = fp / (tn + fp),
            mcc = (tp * tn - fp * fn) / (sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
         )

}
#' Title
#'
#' @param tblBeta
#' @param bolTrue
#'
#' @return
#' @export
#'
#' @examples
moc_tp_fp_table <- function(tblBeta, bolTrue) {

   ## count the selection in simulation of the true beta per model
   tblBetaSum <- foo_confusion_per_beta(tblBeta,bolTrue)

   tblBetaSum %>%
      dplyr::mutate(parameter = gsub("_", " ", stringr::str_to_title(parameter))) %>%
      dplyr::mutate(avg = round(avg, 4)) %>%
      tidyr::pivot_wider(names_from = "parameter", values_from = "avg") %>%
      kableExtra::kbl() %>%
      kableExtra::kable_material(c("striped", "hover"))

}


#' Title
#'
#' @param tblBeta
#' @param bolTrue
#' @param vTrueValue
#'
#' @return
#' @export
#'
#' @examples
moc_beta_distr_plot <- function(tblBeta, bolTrue, vTrueValue) {

      sSelectBeta <- paste0("beta_", which(bolTrue))

      tblBeta <-
            tblBeta %>%
                  dplyr::filter(parameter %in% sSelectBeta) %>%
                  foo_simulation_as_factor()

      lPlot <- vector("list", length(sSelectBeta))

      for (i in seq_along(lPlot)) {

            sTitle <- paste0("Denisty ", gsub("_", " ", stringr::str_to_title(sSelectBeta[[i]])))
            sSubtitle <- paste0("True value: ", vTrueValue[[i]])

            lPlot[[i]] <-
                  tblBeta %>%
                        dplyr::filter(parameter == sSelectBeta[[i]]) %>%
                        ggplot2::ggplot(mapping = ggplot2::aes(value, color = simulation)) +
                        ggplot2::geom_density() +
                        ggplot2::geom_vline(xintercept = vTrueValue[[i]], color = "red", size=0.5) +
                        ggplot2::ggtitle(sTitle, subtitle = sSubtitle) +
                        ggplot2::facet_wrap(~model) +
                        ggplot2::theme_minimal()
      }

      lPlot

}

#' Title
#'
#' @param sType
#' @param tblYrep
#' @param nY
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_plot <- function(tblYrep, sType = c("mean", "95ci"), nY) {

   ## plot true y vs replicated y
   ## use gghighlight package
   ## does not make sense -> each y is different for each simulation...
   ## can compare for each simulation each model for 9 random simulations

   tblYrep$parameter <- as.numeric(gsub(".*_", "", tblYrep$parameter))

   tblYTrue <-
      tblYrep %>%
         dplyr::select(-model) %>%
         dplyr::filter(stat == "y_true") %>%
         dplyr::distinct(simulation, parameter, stat, value)

   tblYTrue$stat <- sType
   tblYTrue$model <- "True Y"

   tblYPlot <- dplyr::bind_rows(tblYrep, tblYTrue)

   tblYPlot %>%
      dplyr::filter(stat == sType) %>%
      ggplot2::ggplot(ggplot2::aes(x = parameter, y = value, color = model)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~simulation)

}



#' Title
#'
#' @param tblYrep
#' @param nSimulation
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_plot_per_simulation <- function(tblYrep, nSimulation) {

   tblYSim <- tblYrep %>% dplyr::filter(simulation == nSimulation)
   tblYSim$parameter <- as.numeric(gsub(".*_", "", tblYSim$parameter))

   tblYPlot <-
      tblYSim %>%
         tidyr::pivot_wider(names_from = stat, values_from = value) %>%
         dplyr::rename(Q97 = "97.5%", Q025 = "2.5%")


   tblYPlot %>%
      ggplot2::ggplot(ggplot2::aes(x = parameter, y = mean)) +
         ggplot2::geom_smooth(ggplot2::aes(ymin = Q025, ymax = Q97), stat = "identity", size = 0.5) +
         ggplot2::geom_point(ggplot2::aes(y = y_true)) +
         ggplot2::labs(
            title = glue::glue("Simulation {nSimulation} per Model"),
            subtitle = "Mean (blue) and 95% predictive intervals (gray) vs. observed data (black)") +
         ggplot2::facet_wrap(~model)
}



#' Title
#'
#' @param tblYrep
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_95_covered <- function(tblYrep) {

   tblYPrep <- foo_yrep_covered(tblYrep)

   tblYPrep %>%
      dplyr::group_by(model, simulation) %>%
         dplyr::summarise(covered = sum(covered) / n()) %>%
         dplyr::mutate(simulation = as.numeric(simulation)) %>%
         ggplot2::ggplot(ggplot2::aes(x = simulation, y = covered, color = model)) +
         ggplot2::geom_point() +
         ggplot2::geom_line() +
         ggplot2::scale_x_continuous(
            breaks = seq(min(tblYPrep$simulation), max(tblYPrep$simulation), by = 1)) +
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
foo_simulation_as_factor <- function(tblInput) {

      tblInput %>%
            dplyr::mutate(simulation = factor(as.integer(simulation),
                                              levels = sort(unique(as.integer(simulation)))))

}


#' Title
#'
#' @param tblBeta
#' @param bolTrue
#'
#' @return
#' @export
#'
#' @examples
# foo_confusion_per_beta <- function(tblBeta, bolTrue) {
#
#    sSelectBeta <- paste0("beta_", which(bolTrue))
#
#    tblBeta %>%
#       dplyr::filter(parameter %in% sSelectBeta) %>%
#       dplyr::group_by(model, simulation, parameter) %>%
#       dplyr::summarise(donna = findInterval(0, quantile(value, probs = c(0.05, 0.95))) == 0L) %>%
#       dplyr::group_by(model, parameter) %>%
#       dplyr::arrange(model, parameter) %>%
#       dplyr::summarise(avg = sum(donna) / n()) %>%
#       dplyr::ungroup()
#
# }


#' Title
#'
#' @param tblYrep
#'
#' @return
#' @export
#'
#' @examples
foo_yrep_covered <- function(tblYrep) {
   tblYrep %>%
      tidyr::pivot_wider(names_from = stat, values_from = value) %>%
      dplyr::rename(Q97 = "97.5%", Q025 = "2.5%") %>%
      dplyr::mutate(
         covered = dplyr::case_when(
            y_true >= Q025 & y_true <= Q97 ~ TRUE,
            TRUE ~ FALSE
         )
      )
}
