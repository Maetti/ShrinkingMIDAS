# library(dplyr)
# library(ggplot2)
# ## connection to files
# dirDGP <- here::here("inst/data/01_dgp")
#
# dirExtract <- paste0(dirDGP, "/", dir(dirDGP), "/output/stan_extract/")
# l1 <- lapply(dirExtract, function(x) arrow::dataset_factory(x, partitioning = c("simulation", "model")))
# arrData <- arrow::open_dataset(l1)
#
# sIteration <- seq(0, 16000, 8)
# sIteration[1] <- 1
#
# tblBeta <-
#       arrData %>%
#             dplyr::filter(parameter_1 == "beta", iteration %in% sIteration) %>%
#             dplyr::collect() %>%
#             dplyr::mutate(parameter = paste0(parameter_1, "_", parameter_2)) %>%
#             dplyr::select(-parameter_1, -parameter_2)
#
# vBeta <- c(0, 0.4, 1.2, 0, 0, 0.8, 1, 0, 0, 0)
#
# tblBetaTrue <- data.frame(
#       parameter = paste0("beta_", seq_along(vBeta)),
#       beta_true = vBeta,
#       bolTrue = bolTrue
# )
# bolTrue <- vBeta != 0
# vTrueValue <- vBeta[bolTrue]
#
# tblConf <-
#       tblBeta %>%
#             dplyr::group_by(model, simulation) %>%
#             tidyr::nest() %>%
#             dplyr::mutate(data = purrr::map(data, function(x, bolTrue) {
#                   x %>%
#                         tidyr::pivot_wider(names_from = "parameter", values_from = "value") %>%
#                         dplyr::select(-iteration) %>%
#                         as.matrix() %>%
#                         md_check_beta_create_confusion_matrix_foo(., bolTrue)
#             }, bolTrue)) %>%
#             tidyr::unnest(cols = c(data)) %>%
#             dplyr::ungroup()

## beta coeffsys
# moc_tp_fp_plot(tblConf, .sType = "TP")
# moc_tp_fp_table(lOut$beta, bolTrue)
# moc_beta_distr_plot(lOut$beta, bolTrue, vTrueValue)
# moc_mse_table(lOut$beta, vBeta)
# moc_matthes_corr(tblConf)




#' Title
#'
#' @param tblBeta
#' @param vBeta
#'
#' @return
#' @export
#'
#' @examples moc_mse_table
md_check_beta_mse <- function(tblBeta, tblBetaTrue) {

      tblBetaPrep <-
            tblBeta %>%
                  dplyr::group_by(model, simulation, key) %>%
                  dplyr::mutate(beta_mean = mean(value)) %>%
                  # mse_var = (value - beta_mean)^2,
                  # mse_bias = (beta_mean - beta_true)^2) %>%
                  dplyr::ungroup() %>%
                  dplyr::left_join(., tblBetaTrue, by = c("parameter", "key"))

      tblVar <-
            tblBetaPrep %>%
                  dplyr::group_by(model, simulation) %>%
                  dplyr::summarise(mse_var = sum((value - beta_mean)^2) / n())

      tblBias <-
            tblBetaPrep %>%
                  dplyr::distinct(model, simulation, parameter, beta_true, beta_mean) %>%
                  dplyr::group_by(model, simulation) %>%
                  dplyr::summarise(mse_bias = sum((beta_mean - beta_true)^2) / n())

      dplyr::left_join(tblVar, tblBias, by = c("model", "simulation")) %>%
            dplyr::mutate(mse = mse_var + mse_bias) %>%
            dplyr::select(model, simulation, mse, mse_var, mse_bias)

}


#' Title
#'
#' @param tblBeta
#' @param bolTrue
#'
#' @return
#' @export
#'
#' @examples moc_tp_fp_table(lOut$beta, bolTrue)
md_check_beta_tpfp_table <- function(tblBeta, tblBetaTrue) {

      ## count the selection in simulation of the true beta per model
      tblBetaSum <- md_check_foo_confusion_per_beta(tblBeta, tblBetaTrue)

      tblBetaSum %>%
            dplyr::arrange(model, key) %>%
            dplyr::mutate(avg = round(avg, 4))

}


#' Title
#'
#' @param tblBeta
#' @param bolTrue
#'
#' @return
#' @export
#'
#' @examples foo_confusion_per_beta
md_check_foo_confusion_per_beta <- function(tblBeta, tblBetaTrue) {

      # sSelectBeta <- paste0("beta_", which(bolTrue))

      tbl1 <-
            tblBeta %>%
               # dplyr::filter(parameter %in% sSelectBeta) %>%
               dplyr::group_by(model, simulation, key) %>%
               dplyr::summarise(donna = findInterval(0, quantile(value, probs = c(0.05, 0.95))) == 0L) # 0 if 0 is not in it

      tbl1 %>%
            dplyr::left_join(., tblBetaTrue, by = c("key")) %>%
            dplyr::mutate(donna = dplyr::case_when(
                  bolTrue ~ donna,
                  !bolTrue ~ !donna
            )) %>%
            dplyr::group_by(model, key) %>%
            dplyr::arrange(model, key) %>%
            dplyr::summarise(avg = sum(donna) / n()) %>%
            dplyr::ungroup()

}

#' Title
#'
#' @param tblBeta
#' @param bolTrue
#' @param vTrueValue
#' @param sWarp
#'
#' @return
#' @export
#'
#' @examples moc_beta_distr_plot
md_check_beta_distribution <- function(tblInput, tblBetaTrue,
                                       bolTrueOnly = TRUE, sWrap = "simulation") {

      # bolTrueOnly <- match.arg(bolTrueOnly, choices = c(TRUE, FALSE))
      # sWrap <- match.arg(sWrap, choices = c("simulation", "model"))


      if (bolTrueOnly) {
            tblBTrue <-
                  tblBetaTrue %>%
                  dplyr::filter(bolTrue)

            sSelectBeta <- tblBTrue$key
            vTrueValue <- tblBTrue$beta_true

      } else {
            sSelectBeta <- tblBetaTrue$key
            vTrueValue <- tblBetaTrue$beta_true
      }


      tbl1 <-
            tblInput %>%
            dplyr::filter(key %in% sSelectBeta) %>%
            foo_simulation_as_factor()

      lPlot <- vector("list", length(sSelectBeta))

      simWrap <- rlang::sym(sWrap)

      for (i in seq_along(lPlot)) {

            sTitle <- paste0("Denisty ", gsub("_", " ", stringr::str_to_title(sSelectBeta[[i]])))
            sSubtitle <- paste0("True value: ", vTrueValue[[i]])

            gp <-
                  tbl1 %>%
                  dplyr::filter(key == sSelectBeta[[i]]) %>%
                  ggplot2::ggplot(mapping = ggplot2::aes(value, color = !!simWrap)) +
                  ggplot2::geom_density() +
                  ggplot2::geom_vline(xintercept = vTrueValue[[i]], color = "red", size=0.5) +
                  ggplot2::ggtitle(sTitle, subtitle = sSubtitle)

            if (sWrap == "simulation") {
                  lPlot[[i]] <-
                        gp +
                        ggplot2::facet_wrap(~model) +
                        ggplot2::theme_minimal()

            } else {
                  lPlot[[i]] <-
                        gp +
                        ggplot2::facet_wrap(~simulation) +
                        ggplot2::theme_minimal()
            }
      }

      lPlot

}



#' Title
#'
#' @param tblInput
#' @param tblBetaTrue
#' @param bolTrueOnly
#' @param sWrap
#'
#' @return
#' @export
#'
#' @examples
md_check_beta_distribution_per_coef <- function(tblInput, tblBetaTrue, bolTrueOnly = TRUE) {

   # bolTrueOnly <- match.arg(bolTrueOnly, choices = c(TRUE, FALSE))
   # sWrap <- match.arg(sWrap, choices = c("simulation", "model"))


   if (bolTrueOnly) {
      tblBTrue <-
         tblBetaTrue %>%
         dplyr::filter(bolTrue)

      sSelectBeta <- tblBTrue$key
      vTrueValue <- tblBTrue$beta_true

   } else {
      sSelectBeta <- tblBetaTrue$key
      vTrueValue <- tblBetaTrue$beta_true
   }


   tbl1 <-
      tblInput %>%
         dplyr::filter(key %in% sSelectBeta) %>%
         foo_simulation_as_factor()

   sModel <- unique(tblInput$model)

   lPlot <- vector("list", length(sModel))

   dfTrueValue <- data.frame(key = 1:length(vTrueValue), "value" = vTrueValue)

   for (i in seq_along(sModel)) {

      sTitle <- gsub("_", " ", stringr::str_to_title(sModel[[i]]))

      lPlot[[i]] <-
         tbl1 %>%
            dplyr::filter(model == sModel[[i]]) %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = value, color = simulation)) +
            ggplot2::geom_density() +
            ggplot2::geom_vline(data = dfTrueValue, mapping = ggplot2::aes(xintercept = vTrueValue), color = "grey", size = 0.5) +
            ggplot2::ggtitle(sTitle) +
            ggplot2::facet_wrap(~key) +
            ggplot2::theme_minimal()
   }

   lPlot

}


#' Title
#'
#' @param tblInput
#' @param sStat
#'
#' @return
#' @export
#'
#' @examples
md_check_boxplot_tp_fp_mcc <- function(tblInput, sStat = "mcc") {

   sStat <- checkmate::matchArg(sStat, choices = c("fn", "fp", "tn", "tp", "tpr", "tnr", "fnr", "fpr", "mcc"))

   tblPlot <- tblInput[, c("model", sStat)]
   colnames(tblPlot) <- c("model", "value")

   tblPlot <-
      tblPlot %>%
         dplyr::mutate(model = stringr::str_to_title(gsub("_", " ", model)))

   sTitle <- switch (sStat,
      "fn" = "False Negative",
      "fp" = "False Positive",
      "tn" = "True Negative",
      "tp" = "True Positive",
      "fnr" = "False Negative Rate",
      "fpr" = "False Positive Rate",
      "tnr" = "True Negative Rate",
      "tpr" = "True Positive Rate",
      "mcc" = "Mathews Correlation Coefficient"
   )

   ggplot2::ggplot(data = tblPlot, ggplot2::aes(x = model, y = value)) +
      ggplot2::geom_violin() +
      ggplot2::ggtitle(sTitle) +
      ggplot2::xlab("Models")

}
#' Title
#'
#' @param tblConf
#' @param .sType
#'
#' @return
#' @export
#'
#' @examples moc_tp_fp_plot
md_check_beta_tp_fp_plot <- function(tblConf, .sType = c("TP", "FN")) {

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
#' @param tblConf
#'
#' @return
#' @export
#'
#' @examples
md_check_beta_true_false_mcc <- function(tblConf) {

      tblConf %>%
            dplyr::mutate(
                  bolSim = as.logical(bolSim),
                  bolTrue = as.logical(bolTrue)) %>%
            dplyr::mutate(conf_matrix = dplyr::case_when(
                  bolSim & bolTrue ~ "tp",
                  !bolSim & !bolTrue ~ "tn",
                  bolSim & !bolTrue ~ "fp",
                  !bolSim & bolTrue ~ "fn")) %>%
            dplyr::group_by(model, simulation, conf_matrix) %>%
            dplyr::summarise(model_freq = sum(Freq)) %>%
            dplyr::ungroup() %>%
            tidyr::pivot_wider(names_from = conf_matrix, values_from = model_freq) %>%
            dplyr::mutate(
                  total_pos = tp + fn,
                  total_neg = tn + fp,
                  tpr = tp / total_pos,
                  tnr = tn / total_neg,
                  fnr = fn / total_pos,
                  fpr = fp / total_neg,
                  mcc = (tp * tn - fp * fn) / (sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
            )

}


#' Title
#'
#' @param maBeta
#' @param bolTrue
#'
#' @return
#' @export
#'
#' @examples
md_check_beta_create_confusion_matrix_foo <- function(maBeta, bolTrue) {

      bolSim <- apply(maBeta, 2, function(x) {findInterval(0, quantile(x, probs = c(0.05, 0.95))) == 0L })
      as.data.frame(table(bolSim, bolTrue))
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
