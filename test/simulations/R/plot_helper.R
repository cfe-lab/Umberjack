library(grid) # for ggplot in grid formation  (unit.pmax)
library(gridExtra)  # for arranging ggplot in grid formation
library(gtable)  # for arranging in grid formation
library(ggplot2)

# Stolen from http://stackoverflow.com/questions/15016995/how-to-align-multiple-ggplot2-plots-and-add-shadows-over-all-of-them
# http://stackoverflow.com/questions/26159495/align-multiple-ggplot-graphs-with-and-without-legends
#+ fig.width=20, fig.height=12

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}