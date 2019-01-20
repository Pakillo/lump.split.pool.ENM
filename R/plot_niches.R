

#' Plot niches
#'
#' @param df Data frame
#' @param suitab.column Character. Name of the column containing suitability values
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#'

plot_niches <- function(df, suitab.column) {

  plot.niches <- ggplot(df, aes(x = env, colour = taxa)) +
    + aes_string(suitab.column) +
    ylim(0, 1) +
    labs(x = "Temperature", y = "Suitability") +
    geom_line(size = 2) +
    viridis::scale_color_viridis(discrete = TRUE) +
    theme(plot.title = element_text(size = 18)) +
    theme(legend.position = "none")

  plot.niches

}
