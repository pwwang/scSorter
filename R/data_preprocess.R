#' Preprocess Data
#'
#' This function validates and preprocesses the input data for the downstream analysis.
#' @keywords internal
#'
#' @param expr A matrix of input data. Each row represents a gene and each column represents a cell.
#' @param anno_processed A list of processed annotation information that consists of the design matrix and the weight matrix for marker genes.
#'
#' @return A list contains processed expression matrix, design matrix, and weight matrix.
#'

data_preprocess <- function(expr, anno_processed) {
  designmat <- anno_processed[[1]]
  weightmat <- anno_processed[[2]]

  rownames(expr) <- make.unique(toupper(rownames(expr)))
  rownames(designmat) <- make.unique(toupper(rownames(designmat)))

  markers <- rownames(designmat)

  markers_avail <- markers %in% rownames(expr)

  if (sum(markers_avail) < length(markers)) {
    warning(paste(
      'The following specified marker genes are not found from the expression data: ',
      paste(markers[!markers_avail], collapse = ', '),
      '.',
      sep = ''
    ))

    designmat <- designmat[markers_avail, ]
    csdmat <- colSums(designmat)

    if (all(csdmat == 0)) {
      stop(
        'No markers in the annotation are found from the expression data. Please check your input data and annotation information.'
      )
    }
    # Remove cell types with no marker genes available
    zero_cols <- which(csdmat == 0)
    if (length(zero_cols) > 0) {
      warning(paste(
        'The following cell types have no marker genes available and will be removed: ',
        paste(colnames(designmat)[zero_cols], collapse = ', '),
        '.',
        sep = ''
      ))
      designmat <- designmat[, -zero_cols, drop = FALSE]
      weightmat <- weightmat[, -zero_cols, drop = FALSE]
    }
  }

  dmat_gene_names <- rownames(designmat)
  dat_gene_names <- rownames(expr)

  picker <- dat_gene_names %in% dmat_gene_names

  expr_mk <- expr[picker, ]
  expr_rt <- expr[!picker, ]

  #reorder genes so that the order of expr mat and design mat matches
  ror <- rep(0, nrow(designmat))
  for (i in 1:nrow(designmat)) {
    ror[i] <- which(rownames(expr_mk) == dmat_gene_names[i])
  }
  expr_mk <- expr_mk[ror, ]

  expr_cb <- rbind(expr_mk, expr_rt)

  return(list(dat = expr_cb, designmat = designmat, weightmat = weightmat))
}
