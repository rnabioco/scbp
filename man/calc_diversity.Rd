% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_diversity}
\alias{calc_diversity}
\title{Calculate cluster sample diversity using Shannon Entropy}
\usage{
calc_diversity(obj, sample_id = "orig.ident", group_id = "coarse_clusters")
}
\arguments{
\item{obj}{Seurat object}

\item{sample_id}{column name containing the per sample id}

\item{group_id}{column name with cluster id}
}
\description{
Calculate cluster sample diversity using Shannon Entropy
}
