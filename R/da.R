#' Calculate differential abundance of cell count labels using edgeR
#'
#' @param sobj seurat object
#' @param cluster_col cluster or cell type column in meta.data
#' @param sample_col sample column in meta.data
#' @param condition_col  column in meta.data indicating condition, defaults to sample_col if not
#' specified
#' @param custom_design design matrix for testing, defaults to `~condition_col`
#'
#' @details See [Orchestraing Single Cell Analysis](https://osca.bioconductor.org/multi-sample-comparisons.html#differential-abundance)
#' @import edgeR
calc_da <- function(sobj,
                    cluster_col = NULL,
                    sample_col = NULL,
                    condition_col = NULL,
                    custom_design = NULL) {

  if(is.null(cluster_col) || is.null(sample_col)){
    stop("please specify a cluster column and sample column")
  }

  abundances <- table(sobj@meta.data[[cluster_col]], sobj@meta.data[[sample_col]])
  abundances <- unclass(abundances)

  extra.info <- sobj@meta.data[match(colnames(abundances), sobj@meta.data[[sample_col]]),]
  y.ab <- DGEList(abundances, samples=extra.info)

  if(is.null(custom_design)){
    if(is.null(condition_col)){
      condition_col <- sample_col
    }
    design <- model.matrix(as.formula(paste0("~", condition_col)), y.ab$samples)
  }

  y.ab <- estimateDisp(y.ab, design, trend="none")
  fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
  res <- glmQLFTest(fit.ab, coef=ncol(design))

  as.data.frame(topTags(res)) %>%
    rownames_to_column(cluster_col) %>%
    select(matches(cluster_col), logFC, PValue, FDR)
}
