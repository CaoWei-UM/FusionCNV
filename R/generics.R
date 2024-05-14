#' Infer CNR from expression matrix
#'
#' Run InferCNR to expression matrix
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns Seurat object or a list with InferCNR matrix and its metadata
#'
#' @export
#'
#' @rdname InferCNR
#' @export InferCNR
#'
InferCNR <- function(object, ...) {
  UseMethod(generic = 'InferCNR', object = object)
}

#' Plot CNR
#'
#' Plot DNA CNR or RNA CNR
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns pheatmap object or NA if specify file to draw
#'
#' @export
#'
#' @rdname PlotCNR
#' @export PlotCNR
#'
PlotCNR <- function(object, ...) {
  UseMethod(generic = 'PlotCNR', object = object)
}

#' Plot CNV
#'
#' Plot DNA CNV or RNA CNV
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns ggplot2 object or NULL if specify file to draw
#'
#' @export
#'
#' @rdname PlotCNV
#' @export PlotCNV
#'
PlotCNV <- function(object, ...) {
  UseMethod(generic = 'PlotCNV', object = object)
}

#' Plot MST of CNV sub clusters
#'
#' Plot minimal spanning tree of CNV sub clusters
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns ggplot2 object or development pair and weight
#'
#' @export
#'
#' @rdname PlotMST
#' @export PlotMST
#'
PlotMST <- function(object, ...) {
  UseMethod(generic = 'PlotMST', object = object)
}

#' Normalize CNR to CNS
#'
#' Segment CNR matrix to segment matrix
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns CNS matrix or CNR object or Seurat object
#'
#' @export
#'
#' @rdname NormalizeCNR
#' @export NormalizeCNR
#'
NormalizeCNR <- function(object, ...) {
  UseMethod(generic = 'NormalizeCNR', object = object)
}

#' Ploidy normalized CNR data
#'
#' Ploidy normalized CNR matrix to CNV matrix
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns CNV matrix to CNR object or Seurat object
#'
#' @export
#'
#' @rdname PloidyCNV
#' @export PloidyCNV
#'
PloidyCNV <- function(object, ...) {
  UseMethod(generic = 'PloidyCNV', object = object)
}

#' Do PCA to weighted segment matrix
#'
#' Use PCA to reduce dimension segments with weight
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a CNV matrix or Seurat object
#'
#' @export
#'
#' @rdname RunCNVPCA
#' @export RunCNVPCA
#'
RunCNVPCA <- function(object, ...) {
  UseMethod(generic = 'RunCNVPCA', object = object)
}

#' Clustering CNV to clusters
#'
#' Use UMAP or PCA to cluster CNV segments with weight to cluster cells
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns cluster result or Seurat object
#'
#' @export
#'
#' @rdname ClusterCNV
#' @export ClusterCNV
#'
ClusterCNV <- function(object, ...) {
  UseMethod(generic = 'ClusterCNV', object = object)
}

#' Use UMAP to reduce dimension CNV matrix
#'
#' Use UMAP to reduce dimension CNV matrix
#'
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a UMAP reduction or Seurat object
#'
#' @export
#'
#' @rdname RunCNVUMAP
#' @export RunCNVUMAP
#'
RunCNVUMAP <- function(object, ...) {
  UseMethod(generic = 'RunCNVUMAP', object = object)
}
