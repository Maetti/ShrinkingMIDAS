#' Title
#'
#' @param tblStat
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_stat_plot <- function(tblStat, sStat = "mean") {

      sStat <- checkmate::matchArg(sStat, choices = unique(tblStat$stat))

      tblStat %>%
            dplyr::filter(stat == sStat) %>%
                  dplyr::mutate(simulation = as.factor(simulation)) %>%
                  ggplot(aes(x = simulation, y = value)) +
                  geom_violin() +
                  geom_point(aes(x = simulation, y = true))

}



#' Title
#'
#' @param tblStat
#' @param sStat
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_stat_violin <- function(tblStat, sStat, sTitle = "Mean",
                                 sXlab = "Simulation", sYlab = "Value", sLegend = "In/Out"
                                 ) {

      sStat <- checkmate::matchArg(sStat, choices = c("mean", "median", "min", "max", "sd"))

      # tblPlot <- tblStat %>% dplyr::filter(stat == sStat)
      tblPlot <-
            tblStat %>%
                  dplyr::filter(stat == sStat) %>%
                  dplyr::mutate(pInOut = dplyr::case_when(
                        p_value == 0 ~ "Out",
                        p_value == 1 ~ "Out",
                        TRUE ~ "In")
                  ) %>%
                  dplyr::mutate(pInOut = factor(pInOut, levels = c("In", "Out")))

      tblTrueStat <- tblPlot %>% dplyr::select(simulation, stat, true) %>% dplyr::distinct()

      ggplot2::ggplot(data = tblPlot, ggplot2::aes(x = factor(simulation), y = value, color = pInOut)) +
            ggplot2::geom_violin() +
            ggplot2::geom_point(data = tblTrueStat,
                                ggplot2::aes(x = factor(simulation), y = true),
                                color = "red", fill = "red", shape = 23) +
            ggplot2::ggtitle(sTitle) +
            ggplot2::xlab(sXlab) + ggplot2::ylab(sYlab) +
            ggplot2::labs(color = sLegend) +
            ggplot2::scale_color_manual(values = c("In" = "#67a9cf", "Out" = "#ef8a62"))
            # ggplot2::scale_color_manual(values = c("TRUE" = "#80cdc1", "FALSE" = "#dfc27d"))
            # ggplot2::theme_classic() +
            # my_thesis_scale_colour() +
            # theme_custom_thesis()

}



#' Title
#'
#' @param tblStat
#' @param sStat
#' @param sTitle
#' @param xlab
#' @param ylab
#' @param sLegend
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_stat_violin_facet <- function(tblStat, sStat, sTitle = "Mean",
                                 sXlab = "Simulation", sYlab = "Value", sLegend = "In/Out") {

   sStat <- checkmate::matchArg(sStat, choices = c("mean", "median", "min", "max", "sd"))

   # tblPlot <- tblStat %>% dplyr::filter(stat == sStat)
   tblPlot <-
      tblStat %>%
      dplyr::filter(stat == sStat) %>%
      dplyr::mutate(pInOut = dplyr::case_when(
         p_value == 0 ~ "Out",
         p_value == 1 ~ "Out",
         TRUE ~ "In")
      ) %>%
      dplyr::mutate(pInOut = factor(pInOut, levels = c("In", "Out")))

   tblTrueStat <- tblPlot %>% dplyr::select(simulation, stat, true) %>% dplyr::distinct()

   ggplot2::ggplot(data = tblPlot, ggplot2::aes(x = factor(simulation), y = value, color = pInOut)) +
      ggplot2::geom_violin() +
      ggplot2::geom_point(data = tblTrueStat,
                          ggplot2::aes(x = factor(simulation), y = true),
                          color = "red", fill = "red", shape = 23) +
      ggplot2::ggtitle(sTitle) +
      ggplot2::xlab(sXlab) + ggplot2::ylab(sYlab) +
      ggplot2::labs(color = sLegend) +
      ggplot2::scale_color_manual(values = c("In" = "#67a9cf", "Out" = "#ef8a62")) +
      ggplot2::facet_wrap(~model, nrow = 3)
   # ggplot2::scale_color_manual(values = c("TRUE" = "#80cdc1", "FALSE" = "#dfc27d"))
   # ggplot2::theme_classic() +
   # my_thesis_scale_colour() +
   # theme_custom_thesis()

}


#' Title
#'
#' @param tblStat
#' @param sStat
#' @param sTitle
#' @param sXlab
#' @param sYlab
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_ecdf_per_y <- function(tblStat, sStat, sTitle = "F(X)",
                                sXlab = "Simulation", sYlab = "F(x)") {

   ggplot2::ggplot(data = tblStat, ggplot2::aes(x = factor(simulation), y = ecdf)) +
      ggplot2::geom_boxplot() +
      # ggplot2::geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
      ggplot2::ggtitle(sTitle) +
      ggplot2::xlab(sXlab) + ggplot2::ylab(sYlab) +
      ggplot2::labs(color = sLegend) +
      ggplot2::facet_wrap(~model, nrow = 3)

}



#' Title
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_intervall_cover <- function(tblInput, nIntervall = 0.9,
                                     sTitle = "F(X)", sXlab = "Simulation", sYlab = "F(x)", sLegend = "within range") {

   nLower <- 0.5 - nIntervall / 2
   nUpper <- 0.5 + nIntervall / 2
   tblPlot <-
      tblInput %>%
      dplyr::mutate(intervall = dplyr::case_when(
         ecdf > nUpper ~ FALSE,
         ecdf < nLower ~ FALSE,
         TRUE ~ TRUE
      )) %>%
      dplyr::group_by(model, simulation) %>%
      dplyr::summarise("outOfIntervall" = 1 - sum(intervall) / n()) %>%
      dplyr::mutate(withinRange = dplyr::case_when(
         outOfIntervall <= (1 - nIntervall) ~ TRUE,
         TRUE ~ FALSE
      )) %>%
      dplyr::mutate(withinRange = factor(withinRange, levels = c("TRUE", "FALSE")))

   tblPlot %>%
      ggplot2::ggplot(aes(x = simulation, y = outOfIntervall, fill = withinRange)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#67a9cf", "FALSE" = "#ef8a62")) +
      # ggplot2::ylim(0, 1) +
      # ggplot2::geom_hline(yintercept = (1 - nIntervall)) +
      ggplot2::ggtitle(sTitle) +
      ggplot2::xlab(sXlab) + ggplot2::ylab(sYlab) +
      ggplot2::labs(fill = sLegend) +
      ggplot2::facet_wrap(~model, nrow = 3) +
      theme_custom_thesis()


}


#' Title
#'
#' @param tblInput
#' @param sModel
#'
#' @return
#' @export
#'
#' @examples
moc_yrep_best_worst_lines <- function(tblInput, sModel,
                                      sTitle = "F(X)", sXlab = "Simulation", sYlab = "F(x)", sLegend = "within range") {

   sModel <- checkmate::matchArg(sModel, choices = c("Group Lasso", "Horseshoe Plus", "Horseshoe"))

   tblRank <-
      tblInput %>%
         dplyr::filter(model == sModel) %>%
         dplyr::mutate(mse_rank = dplyr::dense_rank(mse))

   gpltBest <-
      tblRank %>%
         dplyr::filter(mse_rank %in% c(1, 2, 3)) %>%
            foo_moc_yrep_best_worst_plot(tblRank = ., sTitle = paste0(sModel, ": Replicated Data with lowest MSE"))

   gpltWorst <-
      tblRank %>%
         dplyr::filter(mse_rank %in% rev(sort(unique(tblRank$mse_rank)))[1:3]) %>%
            foo_moc_yrep_best_worst_plot(., sTitle = paste0(sModel, ": Replicated Data with highest MSE"))


   ## output
   list(
      "best" = gpltBest + theme_custom_thesis(base_size = 12),
      "worst" = gpltWorst + theme_custom_thesis(base_size = 12)
   )

}


#' Title
#'
#' @param tblRank
#' @param sTitle
#' @param sSubtitle
#' @param sXLab
#' @param sYlab
#'
#' @return
#' @export
#'
#' @examples
foo_moc_yrep_best_worst_plot <- function(tblRank, sTitle = "Replicated Data",
                                         sSubtitle ="Mean (blue) and 90% predictive intervals (gray) vs. observed data (black)",
                                         sXLab = "Time", sYlab = "") {

   tblRank %>%
      ggplot2::ggplot(., aes(x = key, y = mean)) +
         ggplot2::geom_smooth(aes(ymin = q05, ymax = q95), stat = "identity", size = 0.5) +
         ggplot2::geom_point(aes(y = ytrue)) +
         ggplot2::labs(
            y = sYlab,
            x = sXLab,
            title = sTitle,
            subtitle = sSubtitle
         ) +
         ggplot2::facet_wrap(~simulation, nrow = 3)

}

## TODO: create custom theme

# from https://designpieces.com/2012/12/facebook-colour-palette/
# and darkblue from http://www.color-hex.com/color/3b5998
my_thesis_colors <- c(blue = "#3b5998", medblue = "#6d84b4",
                      lightblue = "#afbdd4", lightestblue = "#d8dfea",
                      white = "#ffffff", black = "#000000",
                      darkblue = "#111a2d")

# facebook color palette (see ggthemes)
my_thesis_pal <- function () {
      function(n) {
            colors <- my_thesis_colors[rev(c("darkblue", "blue", "lightblue"))]
            unname(colors[seq_len(n)])
      }
}

# discrete scale colors for facebook
my_thesis_scale_colour <- function (...) {
      ggplot2::discrete_scale("colour", "myThesis", my_thesis_pal(),...)
}

# theme_custom_thesis <- function(base_size = 11,
#                                 base_family = "Lucida Grande",
#                                 customColor = my_thesis_colors){
#       half_line <- base_size/2
#       theme_bw() ggplot2::`%+replace%`
#             theme(text = element_text(family = base_family, face = "plain",
#                                       color = customColor["medblue"], size = base_size,
#                                       hjust = .5, vjust = .5, angle = 0, lineheight = 1.1,
#                                       margin = margin(), debug = FALSE ),
#                   panel.grid.major = element_line(colour = my_thesis_colors["black"]),
#                   panel.grid.minor = element_line(colour = my_thesis_colors["black"],
#                                                   size = .25),
#                   plot.title = element_text(size = rel(1.2), hjust = 0,
#                                             vjust = 1,
#                                             margin = margin(b = half_line * 1.2),
#                                             face = "bold", color = my_thesis_colors["blue"]),
#                   panel.border = element_rect(fill = NA,
#                                               colour = my_thesis_colors["black"]),
#                   strip.background = element_rect(fill = my_thesis_colors["black"],
#                                                   colour = my_thesis_colors["black"]),
#                   strip.text = element_text(size = rel(0.8),
#                                             colour = my_thesis_colors["white"]),
#                   axis.text = element_text(size = rel(0.8),
#                                            colour = my_thesis_colors["black"]),
#                   axis.ticks = element_line(colour = my_thesis_colors["black"])
#             )
# }


theme_custom_thesis <- function(base_size = 11, base_family = "Lucida Grande") {

   half_line <- base_size / 2

   ggplot2::theme_bw() +
   ggplot2::theme(
      text = ggplot2::element_text(# family = base_family,
                          face = "plain",
                          size = base_size,
                          hjust = .5, vjust = .5, angle = 0, lineheight = 1.1,
                          margin = ggplot2::margin(), debug = FALSE ),
      plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0,
                                vjust = 0.5,
                                margin = ggplot2::margin(b = half_line * 1.2),
                                # face = "bold"
                                ),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.2))
   )

}
