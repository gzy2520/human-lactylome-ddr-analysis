heatmap_publication_font <- "Arial Unicode MS"
heatmap_figure_width_in <- 19.5
heatmap_figure_height_in <- 13.5
heatmap_left_axis_width_mm <- 112
heatmap_right_strip_width_mm <- 48

align_heatmap_gtable <- function(plot) {
  metric_device <- tempfile(fileext = ".png")
  grDevices::png(
    metric_device,
    width = 2600,
    height = 1800,
    res = 144,
    type = "cairo"
  )
  gtable <- tryCatch(
    ggplot2::ggplotGrob(plot),
    finally = {
      grDevices::dev.off()
      unlink(metric_device)
    }
  )

  axis_left_columns <- unique(
    gtable$layout$l[grepl("^axis-l", gtable$layout$name)]
  )
  if (length(axis_left_columns)) {
    gtable$widths[axis_left_columns] <- grid::unit(
      heatmap_left_axis_width_mm,
      "mm"
    )
  }

  strip_right_columns <- unique(
    gtable$layout$l[grepl("^strip-r", gtable$layout$name)]
  )
  if (length(strip_right_columns)) {
    gtable$widths[strip_right_columns] <- grid::unit(
      heatmap_right_strip_width_mm,
      "mm"
    )
  }

  gtable
}

save_aligned_heatmap <- function(plot, figure_dir, stem) {
  aligned <- align_heatmap_gtable(plot)
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    aligned,
    width = heatmap_figure_width_in,
    height = heatmap_figure_height_in,
    dpi = 320,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    aligned,
    width = heatmap_figure_width_in,
    height = heatmap_figure_height_in,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}
