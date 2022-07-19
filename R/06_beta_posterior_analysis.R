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
            dplyr::left_join(., tblBetaTrue, by = "parameter") %>%
            dplyr::group_by(model, simulation, parameter) %>%
            dplyr::mutate(beta_mean = mean(value)) %>%
            # mse_var = (value - beta_mean)^2,
            # mse_bias = (beta_mean - beta_true)^2) %>%
            dplyr::ungroup()


      tblVar <-
            tblBetaPrep %>%
            dplyr::group_by(model) %>%
            dplyr::summarise(mse_var = sum((value - beta_mean)^2) / n())

      tblBias <-
            tblBetaPrep %>%
            dplyr::distinct(model, simulation, parameter, beta_true, beta_mean) %>%
            dplyr::group_by(model) %>%
            dplyr::summarise(mse_bias = sum((beta_mean - beta_true)^2) / n())

      dplyr::left_join(tblVar, tblBias, by = "model") %>%
            dplyr::mutate(mse = mse_var + mse_bias) %>%
            dplyr::select(model, mse, mse_var, mse_bias)


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
            dplyr::mutate(parameter = gsub("_", " ", stringr::str_to_title(parameter))) %>%
            dplyr::mutate(parameter_sort = as.numeric(stringr::str_extract_all(parameter, "[0-9]+"))) %>%
            dplyr::arrange(model, parameter_sort) %>%
            dplyr::select(-parameter_sort) %>%
            dplyr::mutate(avg = round(avg, 4)) %>%
            tidyr::pivot_wider(names_from = "parameter", values_from = "avg")

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
            dplyr::group_by(model, simulation, parameter) %>%
            dplyr::summarise(donna = findInterval(0, quantile(value, probs = c(0.05, 0.95))) == 0L)

      tbl1 %>%
         dplyr::left_join(., tblBetaTrue, by = "parameter") %>%
         dplyr::mutate(donna = dplyr::case_when(
            bolTrue ~ donna,
            !bolTrue ~ !donna
         )) %>%
         dplyr::group_by(model, parameter) %>%
         dplyr::arrange(model, parameter) %>%
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

      sSelectBeta <- tblBTrue$parameter
      vTrueValue <- tblBTrue$beta_true

   } else {
      sSelectBeta <- tblBetaTrue$parameter
      vTrueValue <- tblBetaTrue$beta_true
   }


   tbl1 <-
      tblInput %>%
         dplyr::filter(parameter %in% sSelectBeta) %>%
         foo_simulation_as_factor()

   lPlot <- vector("list", length(sSelectBeta))

   simWrap <- rlang::sym(sWrap)

   for (i in seq_along(lPlot)) {

         sTitle <- paste0("Denisty ", gsub("_", " ", stringr::str_to_title(sSelectBeta[[i]])))
         sSubtitle <- paste0("True value: ", vTrueValue[[i]])

         gp <-
               tbl1 %>%
               dplyr::filter(parameter == sSelectBeta[[i]]) %>%
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
md_check_beta_matthes_cor <- function(tblConf) {

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
