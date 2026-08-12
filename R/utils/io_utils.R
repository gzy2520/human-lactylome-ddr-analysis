# io_utils.R - Shared I/O utilities for kla reanalysis
lib_loaded <- TRUE

relative_path <- function(path, root) {
  sub(
    paste0("^", root, "/"),
    "",
    normalizePath(path)
  )
}

valid_maxquant_rows <- function(data) {
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) {
    keep <- keep & (is.na(data$Reverse) | data$Reverse != "+")
  }
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep &
      (is.na(data$`Potential contaminant`) |
         data$`Potential contaminant` != "+")
  }
  if ("Only identified by site" %in% names(data)) {
    keep <- keep &
      (is.na(data$`Only identified by site`) |
         data$`Only identified by site` != "+")
  }
  keep[is.na(keep)] <- FALSE
  keep
}

read_delimited <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  }
}

write_csv_std <- function(frame, path) {
  write.csv(frame, path, row.names = FALSE, na = "")
}

save_figure <- function(path_stem, plot, width, height, dpi = 350, ...) {
  png_path <- paste0(path_stem, ".png")
  pdf_path <- paste0(path_stem, ".pdf")
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    ...
  )
  grDevices::cairo_pdf(
    file = pdf_path,
    width = width,
    height = height,
    ...
  )
  print(plot)
  grDevices::dev.off()
}
