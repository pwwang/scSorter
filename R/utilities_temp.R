#' Normalize scRNA-seq Data
#'
#' Normalize scRNA-seq data. Please only use this function when you do not have access to Seurat package. More details are available in the vignette of this package.
#'
#' @param expr A matrix of input scRNA-seq data. Rows correspond to genes and columns correpond to cells.
#'
#' @return A matrix of normalized expression data.
#'
#' @export
#'
xnormalize_scData <- function(expr) {
  expr <- as.matrix(expr)

  cs <- colSums(expr)
  normexpr <- t(log(1 + (t(expr) / cs * 10000)))

  rownames(normexpr) <- rownames(expr)
  colnames(normexpr) <- colnames(expr)

  return(normexpr)
}

#' Select Highly Variable Genes
#'
#' Select Highly Variable Genes following the vst approach. Please only use this function when you do not have access to Seurat package. More details are available in the vignette of this package.
#'
#' @param expr A matrix of input scRNA-seq data. Rows correspond to genes and columns correpond to cells.
#' @param ngenes The number of most variable genes to be selected.
#'
#' @return A vector of top highly variable genes with the total number determined by @ngenes option.
#' @importFrom stats loess var
#' @export
#'
xfindvariable_genes <- function(expr, ngenes = 2000) {
  expr <- as.matrix(expr)
  hvginfo <- data.frame(mean = rowMeans(expr))
  hvginfo$var <- apply(expr, 1, var)
  hvginfo$varest <- 0
  hvginfo$varstd <- 0
  varpks <- hvginfo$var > 0

  rfit <- loess(
    formula = log10(var) ~ log10(mean),
    data = hvginfo[varpks, ],
    span = 0.3
  )
  hvginfo$varest[varpks] <- 10^rfit$fitted

  ct <- sqrt(ncol(expr))
  for (i in 1:nrow(expr)) {
    if (hvginfo$var[i] == 0) {
      next
    }
    lvec <- expr[i, ]
    lpks <- lvec > 0
    hvginfo$varstd[i] <- (sum(!lpks) *
      (hvginfo$mean[i] / sqrt(hvginfo$varest[i]))^2 +
      sum(sapply(
        (lvec[lpks] - hvginfo$mean[i]) / sqrt(hvginfo$varest[i]),
        function(x) return(min(x, ct)^2)
      ))) /
      (ncol(expr) - 1)
  }

  rownames(hvginfo) <- rownames(expr)

  hvginfo <- hvginfo[which(hvginfo[, 1] != 0), ]
  hvginfo <- hvginfo[order(hvginfo$varstd, decreasing = T), ]
  topgenes <- utils::head(rownames(hvginfo), ngenes)

  return(topgenes)
}
