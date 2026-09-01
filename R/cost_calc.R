#' Cost Function
#'
#' Calculates the cost.
#' @keywords internal
#'
#' @param dat A matrix of input data.
#' @param clus A vector of predicted cell types.
#' @param mu Parameter estimates from \code{update_mu}.
#' @param designmat An indicator variable matrix records specified marker genes of each cell type.
#' @param rt_sq Per-cell squared sums of the non-marker block (computed once in
#' \code{update_func}); recomputed here when not supplied.
#'

cost_func <- function(dat, clus, mu, designmat, rt_sq = NULL) {
  nmk <- nrow(designmat)
  ngenes <- nrow(dat)
  nclus <- ncol(designmat)
  designmat <- as.matrix(designmat)

  base_mu_vec <- rep(0, nmk)
  for (bv in 1:nmk) {
    zero_cols <- which(designmat[bv, ] == 0)
    if (length(zero_cols) > 0) {
      base_mu_vec[bv] <- mu[bv, zero_cols[1]]
    }
  }

  delta <- mu[1:nmk, ] * designmat
  diff <- dat[1:nmk, ] - base_mu_vec

  ## Marker block: for cell c in cluster clus_c, gene g contributes
  ## min(diff^2, (diff - delta_g,clus_c)^2) = diff^2 + c*ws with
  ## c = delta*(delta - 2*diff) and ws = (c < -eps2). The diff^2 part is
  ## constant across cells (sum(colSums(diff^2))); non-marker genes have
  ## delta = 0 so their correction is 0.
  eps2 <- .Machine$double.eps^0.5
  dcl <- delta[, clus]
  c <- dcl * (dcl - 2 * diff)
  ws <- c < -eps2
  cost <- sum(colSums(diff^2)) + sum(c * ws)

  ## HVG block: sum over cells and genes of (X - mu)^2 splits as
  ## total(X^2) - 2*sum(mu * per-cluster gene sums) + n * sum(mu^2).
  ind <- matrix(0, ncol(dat), nclus)
  ind[cbind(seq_len(ncol(dat)), clus)] <- 1
  cnt <- colSums(ind)
  if (is.null(rt_sq)) {
    hvg_sq_tot <- sum(dat[(nmk + 1):ngenes, , drop = FALSE]^2)
  } else {
    hvg_sq_tot <- sum(rt_sq)
  }
  gene_sums <- dat %*% ind
  gene_sums <- gene_sums[(nmk + 1):ngenes, , drop = FALSE]
  mu_hvg <- mu[(nmk + 1):ngenes, , drop = FALSE]
  cost <- cost + hvg_sq_tot - 2 * sum(gene_sums * mu_hvg) +
    sum(cnt * colSums(mu_hvg^2))

  return(cost)
}
