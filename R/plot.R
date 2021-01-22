#' Plot cells in reduced dimensionality 2D space
#'
#' @description Cells can be colored by gene or feature in meta.data dataframe
#'
#' @param seurat_obj object of class Seurat
#' @param feature feature to plot, either gene name or column in seurat_obj@meta.data
#' @param plot_dat supplemental data.frame containing feature to plot.
#' Must have a column named cell that contains matching colnames in seurat_obj@data
#' @param pt_size size of points produced by geom_point
#' @param pt_alpha alpha value for points plotted by geom_point
#' @param label_text if TRUE display feature labels on plot
#' @param label_size size of label text
#' @param label_color color of label text, defaults to same colors as feature
#' @param .cols vector of colors to use for plot.
#' @param cell_filter character vector of cell names to include in plot
#' @param palette_type color palette type to use (either viridis, brewer, or cloupe)
#' defaults to using cellranger loupe-like colors
#' @param col_pal palette name to use if palette_type is brewer
#' @param max_y maximum feature value to set scale to. Defaults to max of the feature
#' @param legend_title string to supply for title for the legend
#' @param embedding dimensionality reduction to extract from seurat_obj. Can be any
#' reduction present in seurat_obj@reductions (e.g. umap, pca, tsne). defaults to tsne
#' @param show_negative By default the legend value for continuous features will be clipped at zero.
#' If false, then the minumum value for the plotted feature will be used.
#' @param minimal_theme plot bare minimum
#' @param group grouping variable to split plots via faceting
#' @param group_rows number of rows of plots when faceting by group (defaults to NULL)
#' @param group_cols  number of cols of plots when faceting by group
#' @param dims which dims to plot from embedding, defaults to first and second, i.e. c(1,2).
#' @param sorted should the plotting be determined by sorting in ascending order? Default
#' is sorted by_feature (one of "by_feature", "none", "random")
#' @param transform adrgument o be passed to scale_color_gradientn for continuous data. defaults
#' to no transformation (i.e. "identity") See ?continous_scale for available transforms.
#' @param na_col Color for NA values (default = "grey")
#' @param ggrepel_opts named list of options to pass to ggrepel geom_text_repel.
#' @export
plot_feature <- function(seurat_obj,
                         feature = NULL,
                         plot_dat = NULL,
                         pt_size = 0.001,
                         pt_alpha = 1,
                         label_text = FALSE,
                         label_size = 6,
                         label_color = NULL,
                         .cols = NULL,
                         cell_filter = NULL,
                         palette_type = "cloupe",
                         col_pal = "Reds",
                         max_y = NULL,
                         legend_title = NULL,
                         embedding = "tsne",
                         show_negative = FALSE,
                         minimal_theme = FALSE,
                         group = NULL,
                         group_rows = NULL,
                         group_cols = NULL,
                         dims = c(1, 2),
                         sorted = c("by_feature", "none", "random"),
                         transform = "identity",
                         na_col = "grey",
                         ggrepel_opts = list()){

  if(length(feature) > 1){
    args <- as.list(match.call(expand.dots = TRUE)[-1])

    plts <- map(feature,function(x)  {
      args$feature <- x
      do.call(plot_feature, args)
    })

    return(plts)
  }

  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")

  if(!embedding %in% names(seurat_obj@reductions)){
    stop(paste0(embedding, " not found in seurat object"))
  }

  embed_dat <- seurat_obj@reductions[[embedding]]@cell.embeddings %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  embed_cols <- colnames(embed_dat)
  dims_to_plot <- dims + 1
  xcol <- embed_cols[dims_to_plot[1]]
  ycol <- embed_cols[dims_to_plot[2]]

  embed_dat <- left_join(mdata, embed_dat, by = "cell")

  if (!is.null(cell_filter)){
    embed_dat <- dplyr::filter(embed_dat,
                               cell %in% cell_filter)
  }

  meta_data_col <- feature %in% colnames(embed_dat)

  if (!is.null(feature) & !meta_data_col) {
    feature_dat <- FetchData(seurat_obj, feature) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

    # if data is pulled from another assay seurat prefixes the assay name
    # therefore fix column name
    if(colnames(feature_dat)[2] != feature){
      warning("renaming ", colnames(feature_dat)[2], " to ", feature, call. = FALSE)
      colnames(feature_dat)[2] <- feature
    }

    embed_dat <- left_join(embed_dat, feature_dat, by = "cell")
  }

  if (!is.null(plot_dat)){
    embed_dat <- left_join(embed_dat, plot_dat, by = "cell")
  }

  color_aes_str <- feature

  color_aes_str_q <- quo(color_aes_str)

  if(sorted[1] == "by_feature"){
    embed_dat <- embed_dat %>% arrange_at(.vars = color_aes_str)
  } else if (sorted[1] == "random"){
    set.seed(42)
    idx <- sample(1:nrow(embed_dat), nrow(embed_dat), replace = FALSE)
    embed_dat <- embed_dat[idx, ]
  }

  p <- ggplot(embed_dat,
              aes_string(xcol, ycol)) +
    geom_point(aes_string(color = paste0("`", color_aes_str, "`")),
               size = pt_size,
               alpha = pt_alpha)

  p <- p + labs(x = str_replace(xcol, "_", " "),
                y = str_replace(ycol, "_", " "))

  ## discrete or continuous data?
  if (typeof(embed_dat[[feature]]) %in% c(
    "character",
    "logical"
  ) | is.factor(embed_dat[[feature]])) {
    discrete <- TRUE
  } else {
    discrete <- FALSE
  }

  ## increase legend size
  if (discrete) {
    p <- p + guides(colour = guide_legend(override.aes = list(size = 4))) +
      theme(legend.title = element_blank())
  }

  if (label_text) {
    ggrepel_def_opts <- list(
      seed = 42,
      force = 1,
      max.overlaps = 1e6,
      segment.color = NA,
      size = label_size
    )

    ggrepel_opts <- modifyList(ggrepel_def_opts, ggrepel_opts)

    if(discrete) {
      embed_med_dat <- embed_dat %>%
        group_by_at(vars(one_of(feature))) %>%
        mutate(median_x = median(get(xcol)),
               median_y = median(get(ycol))) %>%
        mutate(new_id = ifelse(closest_to_point(data.frame(.data[[xcol]],
                                                           .data[[ycol]]),
                                                c(unique(median_x),
                                                  unique(median_y))),
                               as.character(!!sym(feature)),
                               ""))

      # use same colors as each feature
      if(is.null(label_color)){
        p <- p +
          do.call(ggrepel::geom_text_repel,
                  c(list(data = embed_med_dat,
                          aes_string(x = "median_x",
                                     y = "median_y",
                                     label = "new_id",
                                     color = feature)),
                    ggrepel_opts))
      } else {
        ggrepel_opts$color <- label_color[1]
        p <- p +
          do.call(ggrepel::geom_text_repel,
                  c(list(data = embed_med_dat,
                         aes_string(x = "median_x",
                                    y = "median_y",
                                    label = "new_id")),
                    ggrepel_opts))
      }

    } else {
      warning("label_text not compatible with continuous features")
    }
  }

  ## handle legend limit
  if (is.null(max_y) & !discrete) {
    min_value <- ifelse(show_negative, min(embed_dat[[color_aes_str]], na.rm = TRUE), 0L)
    max_y <- c(min_value, max(embed_dat[[color_aes_str]], na.rm = TRUE))
  } else if (discrete & is.null(max_y)){
    max_y <- c(NA, NA)
  }

  # loupe-like colors
  cols <- rev(brewer.pal(11, "RdGy")[c(1:5, 7)])

  #handle legend name
  if(is.null(legend_title)) legend_title <- color_aes_str

  ## handle zero expression
  if (!all(is.na(max_y)) && all(max_y == c(0, 0))){
    p <- p + scale_color_gradient(low = cols[1], high = cols[1], name = legend_title)
    return(p)
  }

  ## handle colors
  if (is.null(.cols) && !discrete){
    if (palette_type == "viridis") {
      p <- p + scale_color_viridis(discrete = F,
                                   direction = -1,
                                   option = col_pal,
                                   limits = max_y,
                                   name = legend_title,
                                   trans = transform,
                                   na.value = na_col)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1,
                                     name = legend_title,
                                     trans = transform,
                                     na.value = na_col)
    } else if (palette_type == "cloupe") {
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = cols,
                                     name = legend_title,
                                     trans = transform,
                                     na.value = na_col)
    }
  } else if (!is.null(.cols) && !discrete){
    p <- p + scale_color_gradientn(limits = max_y,
                                   colors = .cols,
                                   name = legend_title,
                                   trans = transform,
                                   na.value = na_col)
  } else {

    if(!is.null(.cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = .cols,
        name = legend_title,
        na.value = na_col
      )
    } else {
      p <- p + scale_color_manual(
        values = discrete_palette_default,
        name = legend_title,
        na.value = na_col
      )
    }
  }

  # drop axes, labels, and legend, just plot feature title and projection
  if(minimal_theme){
    p <- p +
      labs(title = feature) +
      theme_void() +
      theme(legend.position="none",
            plot.title = element_text(hjust = 0.5))
  } else {
    p <- p + cowplot::theme_cowplot()
  }

  if(!is.null(group)){
    p <- p +
      facet_wrap(as.formula(paste0("~", group)),
                 nrow = group_rows,
                 ncol = group_cols) +
      theme(strip.background = element_rect(fill = "white"))
  }
  p
}

#' Plot cells in UMAP space
#'
#' @param seurat_obj seurat object
#' @param \dots Additional parameters to pass to plot_feature
#'
#' @rdname plot_feature
#' @importFrom Gmisc fastDoCall
#' @export
plot_umap <- function(seurat_obj, ...){
  cmd_args <- list(seurat_obj = seurat_obj,
                   embedding = "umap",
                   ...)
  Gmisc::fastDoCall(plot_feature, cmd_args)
}

#' Plot cells in tSNE space
#'
#' @param seurat_obj seurat object
#' @param \dots Additional parameters to pass to plot_feature
#'
#' @rdname plot_feature
#' @export
plot_tsne <- function(seurat_obj, ...){
  cmd_args <- list(seurat_obj = seurat_obj,
                   embedding = "tsne",
                   ...)
  # do.call is really terrible with large objects and
  # will not propagate errors quickly, use other function
  Gmisc::fastDoCall(plot_feature, cmd_args)
}

#' Plot cells in PCA space
#'
#' @param seurat_obj seurat object
#' @param \dots Additional parameters to pass to plot_feature
#'
#' @rdname plot_feature
#' @export
plot_pca <- function(seurat_obj, ...){
  cmd_args <- list(seurat_obj = seurat_obj,
                   embedding = "pca",
                   ...)
  Gmisc::fastDoCall(plot_feature, cmd_args)
}

#' Plot cells in Harmony space
#'
#' @param seurat_obj seurat object
#' @param \dots Additional parameters to pass to plot_feature
#'
#' @rdname plot_feature
#' @export
plot_harmony <- function(seurat_obj, ...){
  cmd_args <- list(seurat_obj = seurat_obj,
                   embedding = "harmony_umap",
                   ...)
  Gmisc::fastDoCall(plot_feature, cmd_args)
}

#' plot feature across multiple panels split by group
#'
#' @description See also plot_feature group argument
#' @param seurat_obj seurat object
#' @param feature feature to plot
#' @param group grouping varible to split plots
#' @param embedding dimensionality reduction to use for plotting
#' @param cols vector of cols to identity class, used to keep consistent colors
#' between plots
#' @param add_title want a title?
#' @param ... additional args passed to plot_feature
#'
#' @importFrom stats as.formula median na.omit
#' @importFrom utils write.table
#'
#' @rdname plot_feature
#' @export
plot_features_split <- function(seurat_obj, feature, group = "orig.ident",
                                embedding = "umap", cols = NULL,
                                ...) {

  plot_feature(seurat_obj,
               feature = feature,
               .cols = cols,
               embedding = embedding,
               ...) +
    facet_wrap(as.formula(paste0("~", group))) +
    theme(strip.background = element_rect(fill = "white"))

}


#' @export
plot_violin <- function(df, .x, .y,
                        .fill = NULL,
                        .size = 0.50,
                        .width = 1,
                        .scale = "width",
                        .alpha = 1,
                        cols = ggplot2::scale_fill_viridis_d(),
                        single_col = NULL,
                        jitter = F,
                        rotate_x_text = TRUE,
                        arrange_by_fill = TRUE){

  if (arrange_by_fill && !is.null(.fill)){
    tmp <- sym(.fill)
    df <- arrange(df, !!tmp)
    df[[.x]] <- factor(df[[.x]], levels = unique(df[[.x]]))
  }

  p <- ggplot(df, aes_string(x = .x, y = .y))

  if (jitter){
    p <- p  + geom_jitter(size = 0.1, alpha = 0.2, color = "black")
  }

  if (!is.null(single_col)){
    p <- p +
      geom_violin(size = .size,
                  scale = .scale,
                  fill = single_col,
                  alpha = .alpha)
  } else {
    p <- p +
      geom_violin(aes_string(fill = .fill),
                  size = .size,
                  scale = .scale,
                  alpha = .alpha) +
      cols
  }

  if(rotate_x_text){
    p <- p + theme(axis.text.x = element_text(angle = 90,
                                              hjust = 1,
                                              vjust = 0.5))
  }
  p <- p + theme(legend.title = element_blank())
  p
}

#' Make summary violin plots
#' @param seurat_obj seurat object
#' @param group x axis grouping variable
#' @param features features to plot
#' @param split_by additional metadata columns to split violins per group
#' @param .size violin size
#' @param .width violin width
#' @param .scale violin scaling parameter
#' @param .alpha violin alpha
#' @param cols color palette
#' @param rotate_x_text flip x and y
#' @param arrange_by_fill not sure
#' @param order_by_input maintain order of features supplied
#'
#' @export
plot_violins <- function(seurat_obj, group, features,
                         split_by = NULL,
                        .size = 0.50,
                        .width = 1,
                        .scale = "width",
                        .alpha = 1,
                        cols = discrete_palette_default,
                        rotate_x_text = TRUE,
                        arrange_by_fill = TRUE,
                        order_by_input = TRUE){

  if(length(features) > 1){
    multiple_features <- TRUE
    df <- get_metadata(seurat_obj, features, embedding = NULL) %>%
      tidyr::pivot_longer(cols = one_of(features),
                          names_to = "feature",
                          values_to = "expr")
  } else {
    multiple_features <- FALSE
    df <- get_metadata(seurat_obj, features, embedding = NULL)
  }

  if(order_by_input && multiple_features){
    df <- mutate(df, feature = factor(feature, levels = features))
  }

  if(!is.null(split_by)){
    fill_value <- split_by
  } else {
    fill_value <- group
  }

  if (arrange_by_fill){
    tmp <- rlang::sym(group)
    df <- dplyr::arrange(df, !!tmp)
    if(!is.factor(df[[group]])){
      df[[group]] <- factor(df[[group]], levels = unique(df[[group]]))
    }
  }

  if(multiple_features){
    p <- ggplot(df, aes_string(x = group, y = "expr")) +
      geom_violin(aes_string(fill = fill_value),
                  size = .size,
                  scale = .scale,
                  alpha = .alpha) +
      facet_grid(as.formula("feature ~ ."),
                 scales = "free_y", switch = "y")
  } else {
    p <- ggplot(df, aes_string(x = group,
                               y = str_c("`", features, "`"))) +
      geom_violin(aes_string(fill = fill_value),
                  size = .size,
                  scale = .scale,
                  alpha = .alpha)
  }

  if(rotate_x_text){
    p <- p + theme(axis.text.x = element_text(angle =90,
                                              hjust = 1,
                                              vjust = 0.5))
  }
  p <- p + scale_fill_manual(values = cols)

  p
}

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
                           cols = discrete_palette_default){

  mdata <- get_metadata(obj, embedding = NULL)

  to_keep <- c(sample_id, group_id, facet_by)


  if(!is.null(facet_by)){
    per_patient <- group_by(mdata, !!sym(sample_id)) %>%
      mutate(n_cells = n()) %>%
      group_by(!!sym(sample_id), !!sym(group_id), !!sym(facet_by)) %>%
      summarize(n = n(),
                prop_cell_type = n / unique(n_cells))

    cell_summary <- group_by(mdata, !!sym(sample_id), !!sym(facet_by)) %>%
      mutate(n_cells = n()) %>%
      ungroup() %>%
      select(all_of(to_keep), n_cells) %>%
      mutate(n_cells = str_c("n = ", n_cells),
             n_cells = str_pad(n_cells, max(nchar(n_cells)), "right")) %>%
      unique()

  } else {
    per_patient <- group_by(mdata, !!sym(sample_id)) %>%
      mutate(n_cells = n()) %>%
      group_by(!!sym(sample_id), !!sym(group_id)) %>%
      summarize(n = n(),
                prop_cell_type = n / unique(n_cells))

    cell_summary <- group_by(mdata, !!sym(sample_id)) %>%
      mutate(n_cells = n()) %>%
      ungroup() %>%
      select(all_of(to_keep), n_cells) %>%
      mutate(n_cells = str_c("n = ", n_cells),
             n_cells = str_pad(n_cells, max(nchar(n_cells)), "right")) %>%
      unique()
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

#' Plot heatmap per cluster/cell-type
#'
#'@param obj Seurat object
#'@param features features to plot as rows
#'@param group column from meta data to use to group cells across columns
#'@param annotations additional meta.data columns to add as column annotations
#'supplied as a character vector, defaults to display just group
#'@param average if TRUE, return average values per group rather than display all cells
#'@param normalize_cell_counts if TRUE, downsample cell counts to minimum cell count per group
#'@param slot data slot to retrive values defaults to scale.data
#'@param max_disp max/min value to display (2.5), only applied for scale.data
#'@param col_palettes list of alternative palettes to use for each annotation + group
#'variable supplied.
#'@param default_discrete_pal defatul discrete palette, defaults to scbp::discrete_palette_default,
#'@param default_continuous_pal default continuous color palette, defaults to viridis::inferno(256)
#'@param hmap_options named list of options that are passed to ComplexHeatmap::Heatmap
#'@examples
#'
#' plot_heatmap(Seurat::pbmc_small,
#'              features = rownames(Seurat::pbmc_small@assays$RNA@scale.data),
#'              group = "letter.idents",
#'              annotations = colnames(Seurat::pbmc_small@meta.data))
#'
#'@importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#'@importFrom circlize colorRamp2
#'@importFrom viridis inferno viridis
#'@import Seurat
#'@export
plot_heatmap <- function(obj,
                         features,
                         group,
                         annotations = group,
                         average = FALSE,
                         slot = "scale.data",
                         max_disp = 2.5,
                         col_palettes = NULL,
                         default_discrete_pal = discrete_palette_default,
                         default_continuous_pal_fxn = viridis::inferno(256),
                         hmap_options = NULL) {

  Idents(obj) <- group
  assay <- DefaultAssay(obj)
  annotations <- union(group, annotations)
  check_in_metadata(obj, annotations)
  features <- unique(features)

  if(is.null(col_palettes)){
    col_palettes <- map(seq_along(annotations), ~default_discrete_pal)
  } else {
     if(!(length(col_palettes) == length(annotations))){
      stop("col_palettes should be a list of col_palettes the same length as the # of annotations")
     }
  }

  if(average){
    mat <- AverageExpression(obj,
                             slot = slot,
                             assays = assay,
                             features = features)[[assay]][features, ] %>%
      as.matrix()

    numeric_cols <- annotations[which(sapply(annotations,
                                             function(x){
                                               is.numeric(obj@meta.data[[x]])
                                               }))]

    numeric_cols <- setdiff(numeric_cols, group)
    discrete_cols <- setdiff(annotations, c(numeric_cols, group))

    annot_df <- obj@meta.data[, annotations, drop = FALSE] %>%
      group_by(!!sym(group)) %>%
      mutate_at(numeric_cols, mean, na.rm = TRUE) %>%
      mutate_at(discrete_cols, first) %>%
      distinct() %>%
      arrange(!!sym(group)) %>%
      as.data.frame()
    show_cols <- TRUE
  } else {
    mat <- GetAssayData(obj, slot = slot)[features, ] %>%
      as.matrix()
    annot_df <- obj@meta.data[, annotations, drop = FALSE]
    annot_df <- annot_df[order(annot_df[group]), , drop = FALSE]
    mat <- mat[, rownames(annot_df)]
    show_cols <- FALSE
  }

  annot_cols <- map2(annotations,
                     col_palettes,
                     function(x, cp){
                       to_map <- annot_df[[x]]
                       if(is.numeric(to_map)){
                         res <- circlize::colorRamp2(seq(min(to_map, na.rm = TRUE),
                                        max(to_map, na.rm = TRUE),
                                        length = length(default_continuous_pal_fxn)),
                                    default_continuous_pal_fxn)

                       } else {
                       vals <- to_map %>%
                         unique() %>%
                         sort()
                       res <- cp[1:length(vals)]
                       names(res) <- vals
                       }
                       res
                       })

  names(annot_cols) <- annotations

  if(slot == "scale.data"){
    mat[mat > max_disp] <- max_disp
    mat[mat < -max_disp] <- -max_disp
  }
  ha <- HeatmapAnnotation(df = annot_df,
                          col = annot_cols)

  hmap_option_defaults <- list(
    name = "Log Normalized\nAverage Z-scores",
    col = viridis::viridis(256),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 0,
    show_column_names = show_cols)

  if(is.null(hmap_options)){
    hmap_options <- hmap_option_defaults
  } else {
    opts_to_add <- hmap_option_defaults[which(!names(hmap_option_defaults) %in% names(hmap_options))]
    hmap_options <- c(hmap_options, opts_to_add)
  }

  hmap_options$top_annotation <- ha

  hmat <- do.call(
    function(...){Heatmap(mat, ...)},
    hmap_options)


  hmat

}

#' Get data  closest to a point
#' @returns Logical vector for each row in input matrix, TRUE if closest point
#' @importFrom RANN nn2
closest_to_point <- function(mat, point){

  if(length(point) != ncol(mat)){
    stop("incompatible dimensions for finding closest row to point")
  }

  idx <- RANN::nn2(mat, matrix(point, nrow = 1),k = 1)$nn.idx[1, 1]
  is_closest <- 1:nrow(mat) == idx
  is_closest
}

