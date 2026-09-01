#' Update Cluster
#'
#' Updates cluster assignments based on center estimates from \code{update_mu}
#' @keywords internal
#'
#' @param dat A matrix of input data.
#' @param mu_mat Center estimates from \code{update_mu}
#' @param designmat An indicator variable matrix records specified marker genes of each cell type.
#'

update_C <- function(dat, mu_mat, designmat, dat_rt_t = NULL, rt_sq = NULL) {
  # dat_rt_t: cells x HVG-genes transpose of the (constant) non-marker block,
  # rt_sq: its per-cell squared sums. Precompute both once in update_func to
  # avoid re-copying the block every iteration; fall back to computing them
  # here when not supplied.
  ncmk <- nrow(designmat)
  ngenes <- nrow(dat)
  ncells <- ncol(dat)
  nclusters <- ncol(mu_mat)
  designmat <- as.matrix(designmat)  # row/column indexing a data.frame is slow
  if (is.null(dat_rt_t)) {
    dat_rt_t <- t(dat[(ncmk + 1):ngenes, , drop = FALSE])
    rt_sq <- rowSums(dat_rt_t^2)
  }

  dat_dist_mat_mk <- dat_dist_mat_mk_cache <- matrix(0, ncells, nclusters)

  base_mu_vec <- rep(0, ncmk)
  for (bv in 1:ncmk) {
    zero_cols <- which(designmat[bv, ] == 0)
    if (length(zero_cols) > 0) {
      base_mu_vec[bv] <- mu_mat[bv, zero_cols[1]]
    }
  }

  ## Marker block: non-marker genes contribute the same diff^2 to every
  ## cluster, so sum them once (base) and only loop over each cluster's own
  ## marker genes. Per gene/cell, min(diff^2, (diff - delta)^2) = diff^2 +
  ## c*ws with c = delta*(delta - 2*diff) and ws = (c < -eps2); the diff^2
  ## part sums to base, so only the c*ws correction is needed. Non-marker
  ## genes have delta = 0, so their correction is 0 and self-cancels.
  delta <- mu_mat[1:ncmk, ] * designmat
  diff <- dat[1:ncmk, ] - base_mu_vec
  base <- colSums(diff^2)
  eps2 <- .Machine$double.eps^0.5
  for (j in 1:nclusters) {
    gj <- which(designmat[, j] == 1)
    if (length(gj) == 0) {
      dat_dist_mat_mk[, j] <- base
      next
    }
    diffj <- diff[gj, , drop = FALSE]
    dd <- delta[gj, j]
    c <- dd * (dd - 2 * diffj)
    ws <- c < -eps2
    dat_dist_mat_mk[, j] <- base + colSums(c * ws)
    dat_dist_mat_mk_cache[, j] <- colSums(ws)
  }

  ## Additional block: all cluster distances in one BLAS call instead of
  ## one colSums pass per cluster. dat_rt_t %*% mu_ad == crossprod(dat_ad,
  ## mu_ad) but avoids re-copying the constant HVG block every iteration.
  mu_ad <- mu_mat[(ncmk + 1):ngenes, , drop = FALSE]
  dat_dist_mat_ad <- sweep(
    sweep(-2 * (dat_rt_t %*% mu_ad), 1, rt_sq, "+"),
    2,
    colSums(mu_ad^2),
    "+"
  )

  # calculate the distance matrix
  dat_dist_mat <- dat_dist_mat_mk + dat_dist_mat_ad

  # I also re-wrote the rest to make it quicker and more concise
  # typically tie is not a big problem as it happens rarely; but if you indeed
  # worry about it, my code gives an easy way to get around it.
  dat_dist_mat_rand <- dat_dist_mat +
    stats::rnorm(length(dat_dist_mat)) * .Machine$double.eps^0.5
  # max.col is much faster than apply(which.min); the eps jitter makes exact
  # ties impossible, and ties.method = "first" picks the first minimum anyway.
  clus <- max.col(-dat_dist_mat_rand, ties.method = "first")

  return(list(clus, dat_dist_mat_mk_cache))
}
