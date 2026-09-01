#' Mu Update
#'
#' Solves mu and delta given sample cluster assignment.
#' @keywords internal
#'
#' @param dat A matrix of input data.
#' @param designmat An indicator variable matrix records marker genes of each pre-specified cell type.
#' @param clus A vector of cluster assignment.
#'
#' @return A matrix of parameter estimates.
#'

update_mu <- function(dat, designmat, clus) {
  nclus <- ncol(designmat)
  nsample <- ncol(dat)

  # row indexing a data.frame is much slower than a matrix; convert once
  designmat <- as.matrix(designmat)
  mu_mat <- designmat
  rownames(mu_mat) <- NULL
  colnames(mu_mat) <- NULL

  # Per-cluster cell indices and per-gene marker/zero-column structure are
  # constant within a call; precompute once instead of re-deriving the
  # `clus == k` masks in every inner loop.
  cells <- split(seq_len(nsample), factor(clus, levels = seq_len(nclus)))
  cnt <- lengths(cells)
  mk_of <- apply(designmat, 1, function(r) which(r == 1), simplify = FALSE)
  zc_of <- apply(designmat, 1, function(r) which(r == 0), simplify = FALSE)

  # One BLAS call gives every per-gene per-cluster sum: the marker rows feed
  # the baseline estimate and the HVG rows the p2 means below. Row sums over
  # all cells follow as rowSums(sums).
  ind <- matrix(0, nsample, nclus)
  ind[cbind(seq_len(nsample), clus)] <- 1
  sums <- dat %*% ind
  row_sums <- rowSums(sums[1:nrow(designmat), , drop = FALSE])

  n_marker_gene <- nrow(designmat)
  # Rows of mu_mat are independent: baseline and the delta/basemu fixed point
  # of one gene never touch another gene's row, so each gene is solved alone.
  for (g in 1:n_marker_gene) {
    mk_all <- mk_of[[g]]
    zc <- zc_of[[g]]
    if (length(zc) == 0) {
      mu_mat[g, ] <- 0
      next
    }

    # pre determine baseline level
    n_zero <- nsample - sum(cnt[mk_all])
    if (n_zero == 0) {
      mu_mat[g, ] <- 0
    } else {
      s_zero <- row_sums[g] - sum(sums[g, mk_all])
      mu_mat[g, ] <- s_zero / n_zero
    }

    mk <- mk_all[cnt[mk_all] > 0]
    if (length(mk) == 0) {
      mu_mat[g, zc] <- row_sums[g] / nsample
      next
    }

    for (z in 1:20) {
      mu_old <- mu_mat[g, ]

      # delta estimation: pick out the samples satisfy delta < 2(x-mu)
      # and calculate delta by an iterative approach
      base_mu <- mu_mat[g, zc[1]]
      for (k in mk) {
        obs <- dat[g, cells[[k]]] - base_mu
        obs <- obs[obs > 0]
        if (length(obs) == 0) {
          mu_mat[g, k] <- 0
          next
        }
        for (i in 1:20) {
          delta <- sum(obs) / length(obs)
          picker <- delta < 2 * obs
          if (all(picker)) {
            mu_mat[g, k] <- delta
            break
          } else {
            obs <- obs[picker]
          }
        }
      }

      # baseline nu update: samples not used in delta estimation
      # are used to estimate the baseline nu
      base_mu <- mu_mat[g, zc[1]]
      mu_l <- row_sums[g] / nsample
      for (k in mk) {
        J <- dat[g, cells[[k]]] > base_mu + .5 * mu_mat[g, k]
        mu_l <- mu_l - sum(J) * mu_mat[g, k] / nsample
      }
      mu_mat[g, zc] <- mu_l

      if (mean(abs(mu_old - mu_mat[g, ])) < 10^-6) break
    }
  }

  #this part estimates mu for other highly variable genes which follows kmeans approach.
  n_total_gene <- nrow(dat)
  n_hvg <- n_total_gene - n_marker_gene
  mu_mat_p2 <- matrix(0, n_hvg, nclus)
  if (n_hvg > 0) {
    # sums already holds per-gene per-cluster totals; divide by cluster size
    mu_mat_p2 <- sweep(
      sums[(n_marker_gene + 1):n_total_gene, , drop = FALSE], 2, cnt, "/"
    )
    mu_mat_p2[, cnt == 0] <- 0
  }

  return(rbind(mu_mat, mu_mat_p2))
}
