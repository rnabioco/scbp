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

#' Plot cell proportions across each sample
#'
#' @param obj Seurat object
#' @param sample_id column name containing the per sample id
#' @param group_id column name with cluster id
#'
#' @export
plot_cell_proportions <- function(obj,
                           sample_id = "orig.ident",
                           group_id = "coarse_clusters",
                           facet_by = NULL,
                           cols = NULL){

  mdata <- get_metadata(obj, embedding = NULL)

  to_keep <- c(sample_id, group_id, facet_by)
  cell_summary <- group_by(mdata, !!sym(sample_id)) %>%
    mutate(n_cells = n()) %>%
    ungroup() %>%
    select(all_of(to_keep), n_cells) %>%
    mutate(n_cells = str_c("n = ", n_cells),
           n_cells = str_pad(n_cells, max(nchar(n_cells)), "right")) %>%
    unique()

  if(!is.null(facet_by)){
    per_patient <- group_by(mdata, !!sym(sample_id)) %>%
      mutate(n_cells = n()) %>%
      group_by(!!sym(sample_id), !!sym(group_id), !!sym(facet_by)) %>%
      summarize(n = n(),
                prop_cell_type = n / unique(n_cells))
  } else {
    per_patient <- group_by(mdata, !!sym(sample_id)) %>%
      mutate(n_cells = n()) %>%
      group_by(!!sym(sample_id), !!sym(group_id)) %>%
      summarize(n = n(),
                prop_cell_type = n / unique(n_cells))
  }

  p <- ggplot(per_patient,
              aes_string(sample_id, "prop_cell_type")) +
    geom_col(aes_string(fill = group_id)) +
    labs(x = "Sample ID",
         y = "Proportion of each cell type")

  if(!is.null(cols)){
    p <- p + scale_fill_manual(values = cols)
  }


  if(!is.null(facet_by)){
    p <- p + facet_grid(as.formula(paste0("~", facet_by)), scales = "free_x", space = "free_x")
  }

  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top",
          strip.background = element_rect(fill = "white")) +
    geom_text(data = cell_summary,
              aes_string(x = sample_id, y = 0.15,
                  label = "n_cells"),
              angle = 90)
  p

}
