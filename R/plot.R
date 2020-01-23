#' Plot barcode distribution
#'
#' @param empty_drops object produced by DropletUtils::emptyDrops
#' @param barcode_ranks object produced by DropletUtils::barcodeRanks
#'
#' @importFrom S4Vectors metadata
#' @export
plot_bc <- function(empty_drops,
                    barcode_ranks,
                    fdr = 0.01){

  plt_dat <- cbind(as.data.frame(empty_drops),
                   as.data.frame(barcode_ranks)) %>%
    mutate(is_cell = ifelse(!is.na(FDR),
                            FDR < fdr,
                            FALSE))

  knee <- S4Vectors::metadata(barcode_ranks)$knee
  inflection <- S4Vectors::metadata(barcode_ranks)$inflection
  n_cells <- sum(plt_dat$is_cell, na.rm = TRUE)

  p <- ggplot(plt_dat, aes(rank, Total)) +
    geom_line(aes(color = is_cell)) +
    geom_vline(xintercept = knee,
               color = "blue", linetype = 2) +
    geom_vline(xintercept = inflection,
               color = "green", linetype = 2) +
    annotate("text", y = 1000, x = 1.5 * c(knee,
                                           inflection),
             label = c("knee", "inflection"),
             color = c("blue", "green")) +
    annotate("text",
             x = 5,
             y = 10,
             label = paste0("n = ", n_cells),
             color = "black") +
    scale_color_brewer(palette = "Set1",
                       name = "Is Cell?") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Barcode rank",
         y = "Total UMI count") +
    cowplot::theme_cowplot()

  p
}
