
ggplot2::theme_set(cowplot::theme_cowplot())

#'@export
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#'@export
tableu_classic_palatte <-
  c("#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#7f7f7f",
    "#c7c7c7",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5")

#'@export
discrete_palette_default <- c(tableu_classic_palatte,
                             brewer.pal(8, "Dark2"),
                             palette_OkabeIto)

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
#' @param label_color color of label text
#' @param .cols vector of colors to use for plot.
#' @param cell_filter character vector of cell names to include in plot
#' @param palette_type color palette type to use (either viridis, brewer, or cloupe)
#' defaults to using cellranger loupe-like colors
#' @param col_pal palette name to use if palette_type is brewer
#' @param max_y maximum feature value to set scale to. Defaults to max of the feature
#' @param legend_title string to supply for title for the legend
#' @param embedding dimensionality reduction to extract from seurat_obj. Can be any
#' dr method present in seurat_obj@dr (e.g. umap, pca, tsne). defaults to tsne
#' @param show_negative By default the legend value for continuous features will be clipped at zero.
#' If false, then the minumum value for the plotted feature will be used.
#' @param minimal_theme plot bare minimum
#' @param dims which dims to plot from embedding, defaults to first and second, i.e. c(1,2).
#' @param sorted should the plotting be determined by sorting in ascending order? Default
#' is sorted by_feature (one of "by_feature", "none", "random")
#' @export
plot_feature <- function(seurat_obj,
                         feature = NULL,
                         plot_dat = NULL,
                         pt_size = 0.001,
                         pt_alpha = 1,
                         label_text = FALSE,
                         label_size = 6,
                         label_color = "black",
                         .cols = NULL,
                         cell_filter = NULL,
                         palette_type = "cloupe",
                         col_pal = "Reds",
                         max_y = NULL,
                         legend_title = NULL,
                         embedding = "tsne",
                         show_negative = FALSE,
                         minimal_theme = FALSE,
                         dims = c(1, 2),
                         sorted = c("by_feature", "none", "random")){

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
    if(discrete) {
      embed_mean_dat <- embed_dat %>%
        group_by_at(vars(one_of(feature))) %>%
        summarize(med_dim_1 = median(get(xcol)),
                  med_dim_2 = median(get(ycol)))

      p <- p +
        geom_text(data = embed_mean_dat,
                  aes_string(x = "med_dim_1",
                             y = "med_dim_2",
                             label = feature),
                  size = label_size,
                  color = label_color)
    } else {
      warning("label_text not compatible with continuous features")
    }
  }

  ## handle legend limit
  if (is.null(max_y) & !discrete) {
    min_value <- ifelse(show_negative, min(embed_dat[[color_aes_str]]), 0L)
    max_y <- c(min_value, max(embed_dat[[color_aes_str]]))
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
                                   limits = max_y, name = legend_title)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1, name = legend_title)
    } else if (palette_type == "cloupe") {
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = cols, name = legend_title)
    }
  } else if (!is.null(.cols) && !discrete){
    p <- p + scale_color_gradientn(limits = max_y,
                                   colors = .cols, name = legend_title)
  } else {

    if(!is.null(.cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = .cols,
        name = legend_title
      )
    } else {
      p <- p + scale_color_manual(
        values = discrete_palette_default,
        name = legend_title
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

  p
}


#' @export
plot_umap <- function(...){
  plot_feature(..., embedding = "umap")
}

#' @export
plot_tsne <- function(...){
  plot_feature(..., embedding = "tsne")
}

#' @export
plot_pca <- function(...){
  plot_feature(..., embedding = "pca")
}

#' @export
plot_harmony <- function(...){
  plot_feature(..., embedding = "harmony_umap")
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


#' @export
ilovehue_pal <- c(
   "#e03d6e",
   "#e27c8b",
   "#a64753",
   "#da453e",
   "#db8364",
   "#a54423",
   "#dc652e",
   "#de8c31",
   "#d2a46c",
   "#8f672b",
   "#cea339",
   "#b2b939",
   "#717822",
   "#627037",
   "#a3b46c",
   "#7ba338",
   "#67c042",
   "#3d8829",
   "#35773e",
   "#55c267",
   "#5ca76a",
   "#277257",
   "#5fcea4",
   "#399d82",
   "#40c2d1",
   "#5099cf",
   "#7490df",
   "#615ea5",
   "#716bdf",
   "#c291d6",
   "#984db6",
   "#d558c2",
   "#e17fc0",
   "#995580",
   "#bd3c80"
)

#' @export
get_distinct_cols <- function(vec, seed = 42) {

  seq_col_pals <- c("Blues", "Greens", "Oranges", "Purples", "Reds", "Greys")
  #seq_cols <- map(seq_col_pals, ~brewer.pal(9, .x) %>% .[1:9] %>% rev(.))

  vec <- sort(vec)
  n_needed <- rle(as.character(vec))$lengths
  n_groups <- length(levels(factor(vec)))

  if(n_groups > 6){
    stop("not enough palettes for ", n_groups, " groups", call. = FALSE)
  }

  seq_col_pals <- seq_col_pals[order(n_needed, decreasing = T)]

  vals <- list()
  for (i in 1:n_groups){
    n <- n_needed[i]
    cols <- suppressWarnings(brewer.pal(n, seq_col_pals[i]))
    if (n < 3){
      cols <- cols[n:3]
    }
    vals[[i]] <- cols
  }
  unlist(vals)
}


#' @export
set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}

#' @export
seurat_marker_description <- function(description = ""){
  tibble(
  Columns = c(
    description,
    "",
    "Columns",
    "pval",
    "avg_logFC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "cluster",
    "gene"
  ), Description = c(
    "",
    "",
    "",
    "p-value from wilcox test of indicated cluster compared to other clusters",
    "average fold change expressed in natural log",
    "percent of cells expressing gene (UMI > 0) in cluster",
    "percent of cell expressing gene (UMI > 0) in all other clusters",
    "Bonferroni corrected p-value",
    "cluster name",
    "gene name"
  ))
}

#' @export
presto_marker_description <- function(description = ""){
  tibble(
  Columns = c(
    description,
    "",
    "Columns",
    "feature",
    "group",
    "avgExpr",
    "logFC",
    "statistic",
    "auc",
    "pval",
    "padj",
    "pct_in",
    "pct_out"
  ), Description = c(
    "",
    "",
    "",
    "gene/feature name",
    "group/cluster id",
    "averge normalized expression",
    "average fold change expressed in natural log",
    "test statistic from wilcox test",
    "AUC stat using this single gene as a classifier for the cluster",
    "unadjusted p-value from wilcox test of indicated cluster compared to other clusters",
    "Benjamini-Hochberge corrected p-value from wilcox test of indicated cluster compared to other clusters",
    "percent of cells expressing gene (UMI > 0) in cluster",
    "percent of cell expressing gene (UMI > 0) in all other clusters"
  ))
}

marker_type <- function(df){
  if(is.character(df)){
    df <- suppressMessages(read_tsv(df))
  }
  seurat_cols <- c("cluster",
                   "gene",
                   "p_val",
                   "avg_logFC",
                   "pct.1",
                   "pct.2",
                   "p_val_adj")

  presto_cols <- c("feature",
                   "group",
                   "avgExpr",
                   "logFC",
                   "statistic",
                   "auc",
                   "pval",
                   "padj",
                   "pct_in",
                   "pct_out")
  if(all(colnames(df) %in% seurat_cols)){
    marker_type <- "Seurat"
  } else if (all(colnames(df) %in% presto_cols)) {
    marker_type <- "presto"
  } else {
    stop("unknown marker file")
  }
  marker_type
}

#' Write marker table to excel for easier perusal by collaborators
#' @param mrkr_list list of marker data.frames from Seurat or presto
#' @param path path to output xlsx
#' @param description_string String to add to first excel sheet to describe data
#' . Text is included with scbp::presto_marker_description or scbp::seurat_marker_description
#' @export
write_markers_xlsx <- function(mrkr_list,
                               path,
                               description_string = "Genes differentially expressed between each cluster and all other cells"){

  if(marker_type(mrkr_list[[1]]) == "Seurat"){
    readme_sheet <- seurat_marker_description(description_string)
    gene_col <- "gene"
  } else {
    readme_sheet <- presto_marker_description(description_string)
    gene_col <- "feature"
  }

  readme_sheet <- list(README = readme_sheet)
  names(readme_sheet) <- "README"

  xcel_out <- map(mrkr_list,
                  ~set_xlsx_class(.x, gene_col, "Text"))

  xcel_out <- c(readme_sheet,  xcel_out)

  # santize for spreadsheet tab names
  names(xcel_out) <- str_replace_all(names(xcel_out), "[[:punct:]]", " ")
  openxlsx::write.xlsx(xcel_out,
                       path)

}

#' Write average expression matrix to a file
#' @param sobj seurat object
#' @param col metadata column for averaging
#' @param path output path
#' @param assay assay to write out, defaults to RNA
#' @export
write_avg_expr <- function(sobj, col, path, assay = "RNA") {
  Idents(sobj) <- col
  expr <- AverageExpression(sobj, return.seurat = FALSE)
  expr <- as.data.frame(expr[[assay]]) %>%
    tibble::rownames_to_column("gene")
  write_tsv(expr, path)
}

#' Extract out reduced dimensions and cell metadata to tibble
#'
#' @param obj Seurat Object
#' @param ... additional features to extract into output via FetchData
#' @param embedding dr slot to extract (defaults to all embeddings (2D)),
#' if NULL, will not extract embeddings
#' @export
get_metadata <- function(obj, ..., embedding = names(obj@reductions)) {

  res <- as_tibble(obj@meta.data, rownames = "cell")

  if (!is.null(embedding)) {
    if (any(!embedding %in% names(obj@reductions))) {
      stop(paste0(embedding, " not found in seurat object\n"), call. = FALSE)
    }
    embed_dat <- map(names(obj@reductions),
                     ~obj@reductions[[.x]]@cell.embeddings[, 1:2]) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

    res <- left_join(res,
                     embed_dat,
                     by = "cell")

  }
  vargs <-  list(...)

  if(length(vargs) > 0){
    res <- Seurat::FetchData(obj, vars = ...) %>%
      rownames_to_column("cell") %>%
      left_join(res,
                .,
                by = "cell")
  }
  res
}

#' Get cell counts from seurat object
#' @param sobj seurat object
#' @param ... grouping columns (quoted)
#'
#' @export
get_counts <- function(obj, ...) {
  group_cols <- c(...)
  mdata <- get_metadata(obj, embedding = NULL)
  res <- group_by_at(mdata, vars(all_of(group_cols)))
  res <- res %>%
    summarize(cell_counts = n()) %>%
    ungroup()
  res
}

#' Get cell counts from seurat object
#' @param sobj seurat object
#' @param row_var meta.data column to group for counts. will
#' be rows in the output matrix
#' @param col_var meta.data column to group for counts. will
#' be columns in the output matrix
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @export
get_cell_count_matrix <- function(obj, row_var, col_var){
  res <- get_counts(obj, row_var, col_var)
  res <- tidyr::pivot_wider(res,
              names_from = col_var,
              values_from = "cell_counts")
  res <- tibble::column_to_rownames(res, var = row_var)
  res
}

#'
#' plot feature across multiple panels split by group
#' @param sobj seurat object
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
#' @export
plot_features_split <- function(sobj, feature, group = "orig.ident",
                                embedding = "umap", cols = NULL,
                                ...) {

  plot_feature(sobj,
               feature = feature,
               .cols = cols,
               embedding = embedding,
               ...) +
    facet_wrap(as.formula(paste0("~", group))) +
    theme(strip.background = element_rect(fill = "white"))

}


#' Run fGSEA on gene lists
#' @export
run_fgsea <- function(ranked_gene_list,
                      species = c("human", "mouse"),
                      database = reactome.db::reactome.db,
                      convert_ids = TRUE,
                      min_size=15,
                      max_size=500,
                      n_perm=10000,
                      ...){

  stopifnot(requireNamespace("reactome.db"))
  stopifnot(requireNamespace("fgsea"))
  stopifnot(requireNamespace("AnnotationDbi"))
  stopifnot(requireNamespace("org.Hs.eg.db"))
  stopifnot(requireNamespace("org.Mm.eg.db"))

  if(convert_ids){
    message("converting gene symbols to entrez ids")
    if(species[1] == "human") {
      gs_db <- org.Hs.eg.db
    } else if (species[1] == "mouse") {
      gs_db <- org.Mm.eg.db
    } else {
      gs_db <- species
    }

    e_ids <- mapIds(gs_db, names(ranked_gene_list), 'ENTREZID', 'SYMBOL')
    e_ids <- e_ids[names(ranked_gene_list)]

    new_gene_list <- ranked_gene_list
    names(new_gene_list) <- e_ids
    new_gene_list <- new_gene_list[!is.na(names(new_gene_list))]

  } else {
    new_gene_list <- gene_list
  }

  pathways <- fgsea::reactomePathways(names(new_gene_list))

  res <- list()
  res$fgsea <- fgsea(pathways = pathways,
                    stats = new_gene_list,
               minSize=min_size,
               maxSize=max_size,
               nperm=n_perm,
                ...)

  res$pathways <- pathways
  res$ids <- new_gene_list
  res

}

is_discrete <- function(x) {
  is.character(x) | is.logical(x) | is.factor(x)
}

#' Compute similarities between vectors
#' @export
jaccard_index <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

#' Compute similarities between lists of character vectors
#' @param x input list of character vectors
#' @param y input list of character vectors
#' @return a matrix of pairwise jaccard coefficents
#' @export
jaccard_lists <- function(x, y){
  mat <- matrix(nrow = length(x),
                ncol = length(y),
                dimnames = list(names(x), names(y)))

  to_comp <- expand.grid(names(x), names(y), stringsAsFactors = FALSE)

  for (i in 1:nrow(to_comp)){
    x_id <- to_comp[i, 1]
    y_id <- to_comp[i, 2]

    mat[x_id, y_id] <- jaccard_index(x[[x_id]], y[[y_id]])

  }
  mat
}

#' Compare two seurat objects from different species and return
#' matrices with one-to-one orthologs.
#' @param so1 seurat object
#' @param so2 seurat object
#' @param orthologs data.frame with ortholog info
#' @param col_species_1 column name with species 1 in data.frame
#' @param col_species_2 column name with species 2 in data.frame
#' @export
set_shared_orthologs <- function(so1, so2, orthologs,
                                 col_species_1 = "external_gene_name",
                                 col_species_2 = "mmusculus_homolog_associated_gene_name"){
  mat1 <- so1@assays$RNA@data
  mat2 <- so2@assays$RNA@data

  mat1_orthos <- left_join(tibble(id = rownames(mat1)), orthologs, by = c("id" = col_species_1))
  mat2_orthos <- left_join(tibble(id = rownames(mat2)), orthologs, by = c("id" = col_species_2))

  if(!all(mat2_orthos$id == rownames(mat2))){
    stop("check ortholog table, not 1 to 1 mapping")
  }

  if(!all(mat1_orthos$id == rownames(mat1))){
    stop("check ortholog table, not 1 to 1 mapping")
  }
  shared_orthos <- intersect(mat1_orthos$ortho_id, mat2_orthos$ortho_id) %>% na.omit()

  mat1_orthos <- filter(mat1_orthos, ortho_id %in% shared_orthos)
  mat2_orthos <- filter(mat2_orthos, ortho_id %in% shared_orthos)
  mat1 <- mat1[mat1_orthos$id, ]
  mat2 <- mat2[mat2_orthos$id, ]

  rownames(mat1) <- mat1_orthos$ortho_id
  rownames(mat2) <- mat2_orthos$ortho_id

  #reorder mat1 and mat2 to have same row orders
  mat2 <- mat2[rownames(mat1), ]

  list(mat1 = mat1, mat2 = mat2)
}
#' Cell browser wrapper
#' Builds cell browser with better defaults
#' @param so seurat object
#' @param column_list columns to keep in browser
#' @param primary_color_palette palette for catagorical variables
#' @param secondary_color_palette secondary palette for catagorical variables
#' @param secondary_cols columns to color by secondary_color_palette
#' @param outdir output directory for cellbrowser
#' @param project project string
#' @param marker_file seurat marker file path
#' @param marker_file seurat marker file path
#' @param ident default variable for labeling
#' @param embeddings embeddings to show in browser
#' @param color_skip no idea
#' @param skip_expr_matrix dont overwrite expression matrix
#' @param skip_markers dont overwrite markers
#' @param overwrite_cb_config overwritec cellbrowser.conf file
#' @param annotate_markers annotate markers using cellbrowser functionality
#' @param cellbrowser_dir directory where cellbrowser binary lives
#' @param color_skip no idea
#' @param default_assay assay to export (RNA)
#' @param assays character vector of additional assay expression matrices to export.
#'  If supplied then the data matrix with be concatenated via rowbinding to the default assay matrix.
#'  If a named character vector is supplied the name will be prefixed to the rows (i.e features).
#'
#' @export
make_cellbrowser <- function(so,
                             column_list = NULL,
                             primary_color_palette = discrete_palette_default,
                             secondary_color_palette = palette_OkabeIto,
                             secondary_cols = NULL,
                             outdir = "cellbrowser",
                             project = "seurat",
                             marker_file = NULL,
                             ident = "clusters",
                             embeddings = names(so@reductions),
                             color_skip = NULL,
                             skip_expr_matrix = FALSE,
                             skip_markers = FALSE,
                             overwrite_cb_config = TRUE,
                             annotate_markers = TRUE,
                             default_assay = "RNA",
                             assays = NULL,
                             cellbrowser_dir = "/miniconda3/bin/"
                             ) {

  dir.create(file.path(outdir, "markers"),
             recursive = TRUE, showWarnings = FALSE)

  col_file <- file.path(outdir, paste0(project, "_colorMap.csv"))
  cbmarker_file <- file.path(outdir, "markers", paste0(project, "_markers.tsv"))
  cb_config_file <- file.path(outdir, project, "cellbrowser.conf")

  if(overwrite_cb_config){
    unlink(cb_config_file)
  }

  cols <- colnames(so@meta.data)

  if(!is.null(embeddings)){
    embeddings <- intersect(names(so@reductions), embeddings)
  } else {
    embeddings <- names(so@reductions)
  }

  if(!(all(column_list %in% cols))){
    stop("columns in column_list not found in object",
         call. = FALSE)
  }

  column_list <- c(column_list,
                  structure(c("nCount_RNA", "nFeature_RNA"),
                            names = c("nCount_RNA", "nFeature_RNA")))

  so@meta.data <- so@meta.data[, column_list]
  colnames(so@meta.data) <- names(column_list)
  Idents(so) <- ident

  ## Set colors
  col_palette <- primary_color_palette
  short_col_palette <- secondary_color_palette

  ## assign colors per cluster annotations for discrete types
  if(is.null(secondary_cols)){
    to_primary_cols <- names(column_list)
  } else {
    to_primary_cols <- setdiff(names(column_list), secondary_cols)
  }

  to_map <- to_primary_cols[map_lgl(to_primary_cols,
                                   ~is_discrete(so@meta.data[[.x]]))]

  if(!is.null(color_skip)){
    to_map <- setdiff(to_map, color_skip)
  }

  col_map <- as.list(so@meta.data[, to_map, drop = FALSE]) %>%
    map(~as.character(unique(.x)))

  col_res <- map(col_map,
       function(x) {
         cols = col_palette[1:length(x)]
         structure(cols, names = x)
       })

  if(!is.null(secondary_cols)){
    to_map <- secondary_cols[map_lgl(secondary_cols,
                                      ~is_discrete(so@meta.data[[.x]]))]

    col_map <- as.list(so@meta.data[, to_map, drop = FALSE]) %>%
      map(~as.character(unique(.x)))

    col_res_secondary <- map(col_map,
                   function(x) {
                     cols = short_col_palette[1:length(x)]
                     structure(cols, names = x)
                   })
    col_res <- c(col_res, col_res_secondary)
  }

  map_dfr(col_res,
          ~tibble(clusterName = names(.x),
                  color = .x)) %>%
    as.data.frame() %>%
    readr::write_csv(col_file, col_names = F,  quote_escape = "none")

  cols <- colnames(so@meta.data)
  names(cols) <- colnames(so@meta.data)

  # resave marker file in custom format
  # not that the third column will sorted in rev order
  # only if it has these strings
  #["p_val", "p-val", "p.val", "pval", "fdr"]
  # otherwise ascending order
  # see https://github.com/maximilianh/cellBrowser/blob/c643946d160c9729833a47d1bc44cd49fface6f6/src/cbPyLib/cellbrowser/cellbrowser.py#L2260
  if(!is.null(marker_file)){
    if(marker_type(marker_file) == "Seurat"){
      mkrs <- read_tsv(marker_file) %>%
        select(cluster, gene, fdr = p_val_adj, avg_logFC, everything(), -p_val)
    } else {
      mkrs <- read_tsv(marker_file) %>%
        select(cluster = group,
               gene = feature,
               fdr = padj,
               logFC,
               everything(),
               -pval,
               -statistic)
    }
    print(cbmarker_file)
    write_tsv(mkrs, cbmarker_file)

    if(annotate_markers){
      system2(file.path(cellbrowser_dir, "cbMarkerAnnotate"),
              args = c(cbmarker_file, paste0(cbmarker_file, ".tmp")))
      file.rename(paste0(cbmarker_file, ".tmp"),
                  cbmarker_file)
    }
  } else {
    cbmarker_file <- NULL
  }

  # by default ExportToCellBrowser will not overwrite the markers.tsv.gz file
  # if it exists. This can cause some unexpected issues
  if(!skip_markers){
    unlink(file.path(outdir, project, "markers.tsv"))
  }

  if(!is.null(assays)){
    for(i in seq_along(assays)){
      assay <- assays[i]
      so <- add_features(so, so@assays[[assay]]@data, names(assay), default_assay)
    }
  }

  DefaultAssay(so) <- default_assay

  do.call(function(...) {ExportToCellbrowserFast(so,
                                             dir = file.path(outdir, project),
                                             dataset.name = project,
                                             reductions = embeddings,
                                             markers.file = cbmarker_file,
                                             cluster.field = ident,
                                             skip.expr.matrix = skip_expr_matrix,
                                             ...)},
          as.list(cols))

  # add color line to config

  outline <- paste0("\ncolors=", '"', normalizePath(col_file), '"')
  readr::write_lines(outline, cb_config_file, append = TRUE)

}

#' Build cellbrowser
#' @param dataset_paths vector of paths to cellbrowser directories
#' @param outdir output path
#' @param cbBuild_path path to cbBuild binary
#' @param command print command run
#' @export
build_cellbrowser <- function(dataset_paths,
                              outdir = "cellbrowser_build",
                              cbBuild_path = "/miniconda3/envs/py37/bin/cbBuild",
                              command = TRUE){

  cb_args <- unlist(map(dataset_paths, ~c("-i", .x)))
  out_args <- c("-o", outdir)
  cb_args <- c(cb_args, out_args)

  system2(cbBuild_path,
          args = cb_args)

  if(command){
    message(paste(c(cbBuild_path, cb_args), collapse = " "))
  }
}


#'  get clonotype information per cell
#'
#' @param vdj_dir path to clonotypes file
#' @param prefix prefix to add to cell_id
#' @export
get_clonotypes <- function(vdj_dir, prefix = "") {
  ctypes <- read_csv(file.path(vdj_dir, "filtered_contig_annotations.csv"))

  ctypes <- mutate(ctypes,
                   cell_id = str_remove(barcode, "-1$"),
                   cell_id = str_c(prefix, cell_id)) %>%
    select(cell_id,
           raw_clonotype_id) %>%
    unique()

  consensus_ctypes <- read_csv(file.path(vdj_dir,
                                         "clonotypes.csv"))

  ctypes <- left_join(ctypes,
                      consensus_ctypes,
                      by = c("raw_clonotype_id" = "clonotype_id"))

  ctypes
}

#'  add clonotype information per cell
#'
#' @param sobj Seurat object
#' @param vdj_dirs path to clonotypes file
#' @param prefixes prefixes to add to cell_id
#' @export
add_clonotypes <- function(sobj, vdj_dirs, prefixes = NULL) {

  if(is.null(prefixes)){
    prefixes <- rep("", length(vdj_dirs))
  } else {
    if(length(prefixes) != length(vdj_dirs)){
      stop("prefixes must be the same length as the # of vdj_dirs")
    }
  }
  out <- map2_dfr(vdj_dirs, prefixes, get_clonotypes)

  cells <- tibble(cell_id = Cells(sobj))

  res <- left_join(cells, out, by = "cell_id") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("cell_id")

  sobj <- AddMetaData(sobj, metadata = res)

  sobj

}




#' Export Seurat object for UCSC cell browser
#' optionally uses Readr instead of base R for writing to disk
#'
#' @param object Seurat object
#' @param dir path to directory where to save exported files. These are:
#' exprMatrix.tsv, tsne.coords.tsv, meta.tsv, markers.tsv and a default cellbrowser.conf
#' @param dataset.name name of the dataset. Defaults to Seurat project name
#' @param reductions vector of reduction names to export
#' @param markers.file path to file with marker genes
#' @param cluster.field name of the metadata field containing cell cluster
#' @param cb.dir path to directory where to create UCSC cellbrowser static
#' website content root, e.g. an index.html, .json files, etc. These files
#' can be copied to any webserver. If this is specified, the cellbrowser
#' package has to be accessible from R via reticulate.
#' @param port on which port to run UCSC cellbrowser webserver after export
#' @param skip.expr.matrix whether to skip exporting expression matrix
#' @param skip.metadata whether to skip exporting metadata
#' @param skip.reductions whether to skip exporting reductions
#' @param ... specifies the metadata fields to export. To supply field with
#' human readable name, pass name as \code{field="name"} parameter.
#'
#' @return This function exports Seurat object as a set of tsv files
#' to \code{dir} directory, copying the \code{markers.file} if it is
#' passed. It also creates the default \code{cellbrowser.conf} in the
#' directory. This directory could be read by \code{cbBuild} to
#' create a static website viewer for the dataset. If \code{cb.dir}
#' parameter is passed, the function runs \code{cbBuild} (if it is
#' installed) to create this static website in \code{cb.dir} directory.
#' If \code{port} parameter is passed, it also runs the webserver for
#' that directory and opens a browser.
#'
#' @author Maximilian Haeussler, Nikolay Markov
#'
#' @importFrom utils browseURL
#' @importFrom reticulate py_module_available import
#' @importFrom tools file_ext
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ExportToCellbrowserFast(object = pbmc_small, dataset.name = "PBMC", dir = "out")
#' }
#'
ExportToCellbrowserFast <- function(
  object,
  dir,
  dataset.name = Project(object = object),
  reductions = "tsne",
  markers.file = NULL,
  cluster.field = "Cluster",
  cb.dir = NULL,
  port = NULL,
  skip.expr.matrix = FALSE,
  skip.metadata = FALSE,
  skip.reductions = FALSE,
  ...
) {
  vars <- c(...)

  use_readr <- Seurat:::PackageCheck("readr", error = FALSE)
  use_datatable <- Seurat:::PackageCheck("data.table", error = FALSE)
  if (is.null(x = vars)) {
    vars <- c("nCount_RNA", "nFeature_RNA")
    if (length(x = levels(x = Seurat::Idents(object = object))) > 1) {
      vars <- c(vars, cluster.field)
      names(x = vars) <- c("", "", "ident")
    }
  }
  names(x = vars) <- names(x = vars) %||% vars
  names(x = vars) <- sapply(
    X = 1:length(x = vars),
    FUN = function(i) {
      return(ifelse(
        test = nchar(x = names(x = vars)[i]) > 0,
        yes = names(x = vars[i]),
        no = vars[i]
      ))
    }
  )
  if (!is.null(x = port) && is.null(x = cb.dir)) {
    stop("cb.dir parameter is needed when port is set")
  }
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir)
  }
  if (!dir.exists(paths = dir)) {
    stop("Output directory ", dir, " cannot be created or is a file")
  }
  if (dataset.name == "SeuratProject") {
    warning("Using default project name means that you may overwrite project with the same name in the cellbrowser html output folder")
  }
  order <- colnames(x = object)
  enum.fields <- c()
  # Export expression matrix:
  if (!skip.expr.matrix) {
    # Relatively memory inefficient - maybe better to convert to sparse-row and write in a loop, row-by-row?
    df <- as.data.frame(x = as.matrix(x = Seurat::GetAssayData(object = object)))
    df <- data.frame(gene = rownames(x = object), df, check.names = FALSE)
    gzPath <- file.path(dir, "exprMatrix.tsv.gz")

    message("Writing expression matrix to ", gzPath)

    if(!use_readr){
      z <- gzfile(gzPath, "w")
      write.table(x = df, sep = "\t", file = z, quote = FALSE, row.names = FALSE)
      close(con = z)
    } else if (use_datatable) {
      data.table::fwrite(x = df, file = gzPath,
                         sep = "\t", compress = "gzip")
    } else {
      readr::write_tsv(x = df, path = gzPath)
    }

  }
  # Export cell embeddings
  embeddings.conf <- c()
  for (reduction in reductions) {
    if (!skip.reductions) {
      df <- Seurat::Embeddings(object = object, reduction = reduction)
      if (ncol(x = df) > 2) {
        warning(
          'Embedding ',
          reduction,
          ' has more than 2 coordinates, taking only the first 2'
        )
        df <- df[, 1:2]
      }
      colnames(x = df) <- c("x", "y")
      df <- data.frame(cellId = rownames(x = df), df)
      fname <- file.path(dir, paste0(reduction, '.coords.tsv'))
      message("Writing embeddings to ", fname)
      if(!use_readr){
        write.table(
          x = df[order, ],
          sep = "\t",
          file = fname,
          quote = FALSE,
          row.names = FALSE
        )
      } else if (use_datatable) {
        data.table::fwrite(x = df[order, ], file = fname,
                           sep = "\t")
      } else {
        readr::write_tsv(x = df[order, ], path = fname)
      }

    }
    conf <- sprintf(
      '{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
      reduction
    )
    embeddings.conf <- c(embeddings.conf, conf)
  }
  # Export metadata
  df <- data.frame(row.names = rownames(x = object[[]]))
  df <- Seurat::FetchData(object = object, vars = names(x = vars))
  colnames(x = df) <- vars
  enum.fields <- Filter(
    f = function(name) {!is.numeric(x = df[[name]])},
    x = vars
  )
  if (!skip.metadata) {
    fname <- file.path(dir, "meta.tsv")
    message("Writing meta data to ", fname)
    df <- data.frame(Cell = rownames(x = df), df, check.names = FALSE)

    if(!use_readr){
      write.table(
        x = df[order, ],
        sep = "\t",
        file = fname,
        quote = FALSE,
        row.names = FALSE
      )
    } else if (use_datatable) {
      data.table::fwrite(x = df[order, ],
                        file = fname,
                        sep = "\t")
    } else {
      readr::write_tsv(x = df[order, ],
                       path = fname)
    }

  }
  # Export markers
  markers.string <- ''
  if (!is.null(x = markers.file)) {
    ext <- tools::file_ext(x = markers.file)
    fname <- paste0("markers.", ext)
    file.copy(from = markers.file, to = file.path(dir, fname))
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      fname
    )
  }
  config <- c(
    'name="%s"',
    'shortLabel="%1$s"',
    'exprMatrix="exprMatrix.tsv.gz"',
    '#tags = ["10x", "smartseq2"]',
    'meta="meta.tsv"',
    '# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"',
    'geneIdType="auto"',
    'clusterField="%s"',
    'labelField="%2$s"',
    'enumFields=%s',
    '%s',
    'coords=%s'
  )
  config <- paste(config, collapse = '\n')
  enum.string <- paste0(
    "[",
    paste(paste0('"', enum.fields, '"'), collapse = ", "),
    "]"
  )
  coords.string <- paste0(
    "[",
    paste(embeddings.conf, collapse = ",\n"),
    "]"
  )
  config <- sprintf(
    config,
    dataset.name,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  fname <- file.path(dir, "cellbrowser.conf")
  if (file.exists(fname)) {
    message(
      "`cellbrowser.conf` already exists in target directory, refusing to ",
      "overwrite it"
    )
  } else {
    cat(config, file = fname)
  }
  message("Prepared cellbrowser directory ", dir)
  if (!is.null(x = cb.dir)) {
    if (!reticulate::py_module_available(module = "cellbrowser")) {
      stop(
        "The Python package `cellbrowser` is required to prepare and run ",
        "Cellbrowser. Please install it ",
        "on the Unix command line with `sudo pip install cellbrowser` (if root) ",
        "or `pip install cellbrowser --user` (as a non-root user). ",
        "To adapt the Python that is used, you can either set the env. variable RETICULATE_PYTHON ",
        "or do `require(reticulate) and use one of these functions: use_python(), use_virtualenv(), use_condaenv(). ",
        "See https://rstudio.github.io/reticulate/articles/versions.html; ",
        "at the moment, R's reticulate is using this Python: ",
        reticulate::import(module = 'sys')$executable,
        ". "
      )
    }
    if (!is.null(x = port)) {
      port <- as.integer(x = port)
    }
    message("Converting cellbrowser directory to html/json files")
    cb <- reticulate::import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
    if (!is.null(port)) {
      message("Starting http server")
      cb$cellbrowser$stop()
      cb$cellbrowser$serve(cb.dir, port)
      Sys.sleep(time = 0.4)
      utils::browseURL(url = paste0("http://localhost:", port))
    }
  }
}


#' Preprocessing for Bustools output
#'
#' @param paths paths to directories with mtx files. If a
#' named vector is supplied then the cellbarcode will be
#' prefixed with "name_".
#' @param fdr fdr cutoff for cell calling
#' @param cell_bc list of character vectors of cell ids to keep in output,
#' if supplied cell calling will be skipped. Order must match paths argument
#' @param ... additional arguments passed onto emptyDrops
#'
#' @return A list of named lists with slots:
#'   mtx a sparseMatrix with only cell associated barcodes
#'   empty_drops a data.frame with information from emptyDrops
#'   barcode_ranks a data.frame with information from barcodeRanks
#' @importFrom BUSpaRse read_count_output
#' @importFrom DropletUtils emptyDrops barcodeRanks
#' @export
preprocess_bus <- function(paths,
                           fdr = 0.01,
                           cell_bc = NULL,
                           ...){

  res <- list()

  for (i in seq_along(paths)){

    mtx_path <- paths[i]
    message("loading matrix: ", mtx_path)

    mtx_prefix <- stringr::str_subset(dir(mtx_path), ".mtx$") %>%
      stringr::str_remove(".mtx")

    mat <- BUSpaRse::read_count_output(mtx_path,
                                       name = mtx_prefix,
                                       tcc = FALSE)
    if(is.null(cell_bc)){
      barcode_ranks <- DropletUtils::barcodeRanks(mat)

      empty_drops <- DropletUtils::emptyDrops(mat, ...)

      is_cell <- empty_drops$FDR <= fdr
      cell_bcs <- rownames(empty_drops)[which(is_cell)]
    } else {
      cell_bcs <- cell_bc[[i]]
      barcode_ranks <- data.frame()
      empty_drops <- data.frame()
    }

    cell_mat <- mat[, cell_bcs, drop = FALSE]

    prefix <- names(paths[i])

    if(!is.null(prefix)){
      colnames(cell_mat) <- stringr::str_c(
        prefix,
        "_",
        colnames(cell_mat))

      if(is.null(cell_bc)){
        rownames(empty_drops) <- stringr::str_c(
          prefix,
          "_",
          rownames(empty_drops))

        rownames(barcode_ranks) <- stringr::str_c(
          prefix,
          "_",
          rownames(barcode_ranks))
      }
    }
    res[[i]] <- list(mtx = cell_mat,
                     barcode_ranks = barcode_ranks,
                     empty_drops = empty_drops)
  }

  res
}

#' Calculate cluster sample diversity using Shannon Entropy
#'
#' @param obj Seurat object
#' @param sample_id column name containing the per sample id
#' @param group_id column name with cluster id
#'
#' @export
calc_diversity <- function(obj,
                           sample_id = "orig.ident",
                           group_id = "coarse_clusters"){
  mdata <- get_metadata(obj, embedding = NULL)
  props <- group_by(mdata, !!sym(group_id)) %>%
    mutate(n_cells = n()) %>%
    group_by(!!sym(group_id), !!sym(sample_id)) %>%
    summarize(n = n(),
              n_cells = unique(n_cells),
              prop = n / n_cells)

  props <- split(props, props[[group_id]]) %>%
    purrr::map(~pull(.x, prop))

  etp <- map_dfr(props, shannon_entropy) %>%
    imap_dfr(~tibble(!!sym(group_id) := .y,
                entropy = .x))

  if(is.factor(mdata[[group_id]])){
    etp <- mutate(etp,
                  !!sym(group_id) := factor(!!sym(group_id),
                                            levels = levels(mdata[[group_id]])))
  }

  res <- left_join(mdata, etp, by = group_id) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("cell")

  Seurat::AddMetaData(obj, res)

}

shannon_entropy <- function(x){
  # H' = sum of pi * ln(pi) for each proportion
  if(any(x > 1 | x < 0)){
    stop("values must be between 0 and 1")
  }

  if(any(is.na(x))){
    warning("NA values found, ignoring")
  }

  if(!isTRUE(all.equal(sum(x), 1))){
    stop("proportions should sum to 1")
  }

  se <- -sum(x * log(x), na.rm = TRUE)

  # scale by log(length(x)) to be between 0 and 1
  if(length(x) == 1){
    res <- 0
  } else {
    res <- se / log(length(x))
  }

  res
}


#' add features (i.e. rows) to a Seurat object
#' used for exporting various cell level measurements
#' for labeled in ucsc cell browser. Note that column names much be
#' identical to input seurat obj. Only the 'data' matrix will be altered
#' @param so seurat object
#' @param mat data to add to expression matrix
#' @param prefix optional prefix to add to rownames of mat before
#' @param assay assay to which data will be appended, defaults to RNA
#' merging
#' @importFrom stringr str_c
add_features <- function(so,
                         mat,
                         prefix = NULL,
                         assay = "RNA") {

  mat <- mat[, colnames(so), drop = FALSE]
  mat <- as(mat, "dgCMatrix")
  if(!is.null(prefix)){
    rownames(mat) <- stringr::str_c(prefix, rownames(mat))
  }
  so@assays[[assay]]@data <- rbind(so@assays[[assay]]@data, mat)
  so
}




