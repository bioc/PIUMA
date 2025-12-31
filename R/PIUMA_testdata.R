#' Example datasets for PIUMA package
#'
#' We tested PIUMA on a subset of the single-cell RNA Sequencing dataset
#' (GSE:GSE193346 generated and published by Feng et al. (2022)
#' on Nature Communication to demonstrate that distinct
#' transcriptional profiles are present in specific cell types of each heart
#' chambers, which were attributed to have roles in cardiac development.
#' In this tutorial, our aim will be to exploit PIUMA for identifying
#' sub-population of vascular endothelial cells, which can be associated
#' with specific heart developmental stages. The original dataset consisted
#' of three layers of heterogeneity: cell type, stage and zone
#' (i.e., heart chamber). Our testing dataset was obtained by
#' subsetting vascular endothelial cells (cell type) by Seurat object,
#' extracting raw counts and metadata. Thus, we filtered low expressed genes
#' and normalized data by DaMiRseq
#'
#' @format A data frame with 1180 rows (cells) and 2 columns (outcomes).
#' @usage data(vascEC_meta)
#' @name vascEC_meta
#' @docType data
#' @keywords datasets
NULL

#' We tested PIUMA on a subset of the single-cell RNA Sequencing dataset
#' (GSE:GSE193346 generated and published by Feng et al. (2022)
#' on Nature Communication to demonstrate that distinct
#' transcriptional profiles are present in specific cell types of each heart
#' chambers, which were attributed to have roles in cardiac development.
#' In this tutorial, our aim will be to exploit PIUMA for identifying
#' sub-population of vascular endothelial cells, which can be associated
#' with specific heart developmental stages. The original dataset consisted
#' of three layers of heterogeneity: cell type, stage and zone
#' (i.e., heart chamber). Our testing dataset was obtained by
#' subsetting vascular endothelial cells (cell type) by Seurat object,
#' extracting raw counts and metadata. Thus, we filtered low expressed genes
#' and normalized data by DaMiRseq
#'
#' @format A matrix with 1180 rows (cells) and 838 columns (genes).
#' @usage data(vascEC_norm)
#' @name vascEC_norm
#' @docType data
#' @keywords datasets
NULL

#' A dataset to test the \code{\link{dfToProjection}} and
#' \code{\link{dfToDistance}} funtions of \code{PIUMA} package.
#'
#' @format A data frame with 15 rows (cells) and 15 columns (genes).
#' @usage data(df_test_proj)
#' @name df_test_proj
#' @docType data
#' @keywords datasets
NULL

#' A TDAobj to test the \code{PIUMA} package.
#'
#' A \code{TDAobj} with data in all slots for testing.
#'
#' @format A \code{TDAobj}.
#' @usage data(tda_test_data)
#' @name tda_test_data
#' @docType data
#' @keywords datasets
NULL
