#' Update Cluster
#'
#' Updates cluster assignments based on center estimates from \code{update_mu}
#' @keywords internal
#'
#' @param dat A matrix of input data.
#' @param mu_mat Center estimates from \code{update_mu}
#' @param designmat An indicator variable matrix records specified marker genes of each cell type.
#'

update_C <- function(dat, mu_mat, designmat) {
  ncmk <- nrow(designmat)
  ngenes <- nrow(dat)
  ncells <- ncol(dat)
  nclusters <- ncol(mu_mat)

  dat_dist_mat_mk <- dat_dist_mat_mk_cache <- matrix(0, ncells, nclusters)

  base_mu_vec <- rep(0, ncmk)
  for (bv in 1:ncmk) {
    zero_cols <- which(designmat[bv, ] == 0)
    if (length(zero_cols) > 0) {
      base_mu_vec[bv] <- mu_mat[bv, zero_cols[1]]
    }
  }

  ## Marker block: non-marker genes contribute the same diff^2 to every
  ## cluster, so sum them once (base) and only loop over each cluster's
  ## own marker genes.
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
    mat1 <- (diffj - delta[gj, j])^2
    mat2 <- diffj^2
    which.smaller <- (mat1 < mat2 - eps2)
    # On the above, when you compare two float numbers, it is always a good idea to specify
    # whether you want to include or exclude the equal case
    dat_dist_mat_mk[, j] <- base - colSums(mat2) + colSums(
      mat1 * which.smaller + mat2 * (!which.smaller)
    )
    dat_dist_mat_mk_cache[, j] <- colSums(which.smaller)
  }

  ## Additional block: all cluster distances in one BLAS call instead of
  ## one colSums pass per cluster.
  dat_ad <- dat[(ncmk + 1):ngenes, , drop = FALSE]
  mu_ad <- mu_mat[(ncmk + 1):ngenes, , drop = FALSE]
  dat_dist_mat_ad <- sweep(
    sweep(-2 * crossprod(dat_ad, mu_ad), 1, colSums(dat_ad^2), "+"),
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
  clus <- apply(dat_dist_mat_rand, 1, which.min)

  return(list(clus, dat_dist_mat_mk_cache))
}
