#' Update Function
#'
#' Implements the scSorter method by iteratively running \code{update_mu} and \code{update_C}.
#' @keywords internal
#'
#' @param dat A matrix of input data.
#' @param design_mat An indicator variable matrix records specified marker genes of each cell type.
#' @param weightmat A matrix of weights assigned to each marker gene.
#' @param unknown_threshold1 The parameter determines undecided cells cutoff. The default value is 0.
#' @param unknown_threshold2 The parameter determines whether undecided cells are further processed. The default value is 0.05.
#' @param max_iter The maximum number of iterations for the algorithm to update parameters. The default value is 100.
#' @importFrom stats sd pchisq
#' @return A list contains parameter estimates, type assignments, and the corresponding cost.
#'

update_func <- function(
  dat,
  design_mat,
  weightmat,
  unknown_threshold1 = 0,
  unknown_threshold2 = 0.05,
  max_iter = 100
) {
  cluster <- sample(1:ncol(design_mat), ncol(dat), replace = T)

  rgwmat <- range(weightmat)
  if (rgwmat[1] == rgwmat[2]) {
    ldt <- nrow(dat)
    lmk <- nrow(design_mat)
    marker_w <- rep(weightmat[1, 1], lmk)
    rest_w <- rep(min(lmk / (ldt - lmk), 1), ldt - lmk)
    w <- c(marker_w, rest_w)
    dat <- dat * sqrt(w)

    # HVG block is constant across iterations: precompute its transpose and
    # per-cell squared sums once instead of re-copying 160MB every iteration.
    dat_rt_t <- t(dat[(lmk + 1):ldt, , drop = FALSE])
    rt_sq <- rowSums(dat_rt_t^2)

    for (i in 1:max_iter) {
      mu <- update_mu(dat, design_mat, cluster)
      cluster_old <- cluster
      cluster_ot <- update_C(dat, mu, design_mat, dat_rt_t, rt_sq)
      cluster <- cluster_ot[[1]]
      if (sum(cluster != cluster_old) == 0) {
        break
      }
    }
  } else {
    ldt <- nrow(dat)
    lmk <- nrow(design_mat)
    rest_w <- rep(min(lmk / (ldt - lmk), 1), ldt - lmk)

    dat_mk <- dat[1:lmk, ]
    dat_rt <- dat[(lmk + 1):ldt, ]

    dat_rt <- dat_rt * sqrt(rest_w)
    wmat_sqrt <- sqrt(weightmat)

    # Build the per-iteration matrix once; only the small marker block is
    # re-weighted each iteration. rbind-ing the full matrix every iteration
    # is much slower than rewriting one block.
    dat2 <- matrix(0, ldt, ncol(dat))
    dat2[(lmk + 1):ldt, ] <- dat_rt

    # HVG block is constant across iterations: precompute its transpose and
    # per-cell squared sums once instead of re-copying 160MB every iteration.
    dat_rt_t <- t(dat_rt)
    rt_sq <- rowSums(dat_rt_t^2)

    for (i in 1:max_iter) {
      # weightmat[, cluster] indexes each cell's per-cluster marker weights
      # in one vectorized call, replacing the per-cluster loop
      dat2[1:lmk, ] <- dat_mk * wmat_sqrt[, cluster]
      mu <- update_mu(dat2, design_mat, cluster)
      cluster_old <- cluster
      cluster_ot <- update_C(dat2, mu, design_mat, dat_rt_t, rt_sq)
      cluster <- cluster_ot[[1]]
      if (sum(cluster != cluster_old) == 0) {
        break
      }
    }

    dat <- dat2
  }

  cache_mat <- cluster_ot[[2]]
  numofmarkergenes <- colSums(design_mat)
  numofhighexpmkgenes <- cache_mat[cbind(1:nrow(cache_mat), cluster)] /
    numofmarkergenes[cluster]
  pks <- numofhighexpmkgenes <= unknown_threshold1

  cluster_ukn_helper <- cluster[pks]
  cluster_ukn <- cluster[pks]

  cluster_kn <- cluster[!pks]

  #detect unknown cells from a given cluster.
  #Now we need to determine whether to put those potential unknown cells back to the corresponding cluster.
  #The total cost of two groups of cells (cells of that cluster and unknown cells) should reach a minimum when true unknown cells are picked out while the rest are put back into the cluster.
  #Under chi-square distribution, a cutoff is set to pick out true unknown cells. Since unknown cells are expected to be far away from the samples of that cluster. The search starts from one standard deviation.
  #50 cutoffs are selected and the one leads to the minimum total cost are used as final cutoff.
  unknown_detector <- function(dat_kn, dat_ukn) {
    pknz <- rowSums(dat_kn != 0) != 0
    df <- sum(pknz)

    dat_kn <- dat_kn[pknz, , drop = FALSE]
    dat_ukn <- dat_ukn[pknz, , drop = FALSE]

    rm <- rowMeans(dat_kn)
    sdr <- apply(dat_kn, 1, sd)

    dat_kn2 <- (dat_kn - rm) / sdr
    dat_ukn2 <- (dat_ukn - rm) / sdr

    nkn <- ncol(dat_kn)
    nuk <- ncol(dat_ukn)

    # cost_calc2(x) = sum_g [colsum2_g - (colsum_g)^2 / n]; the "known" set
    # for cutoff p is a prefix of the unknown cells ordered by their
    # chi-square p-value, so all 50 costs come from one pair of per-gene
    # cumulative sums instead of 50 full re-computations.
    s_kn <- rowSums(dat_kn2)
    q_kn <- rowSums(dat_kn2^2)
    s_uk <- rowSums(dat_ukn2)
    q_uk <- rowSums(dat_ukn2^2)

    pv <- stats::pchisq(colSums(dat_ukn2^2), df)
    ord <- order(pv)
    cum_s <- apply(dat_ukn2[, ord, drop = FALSE], 1, cumsum)
    cum_q <- apply(dat_ukn2[, ord, drop = FALSE]^2, 1, cumsum)

    pl <- c(
      seq(stats::pchisq(df + sqrt(2 * df), df), 0.9, length.out = 10),
      seq(0.9, 0.95, length.out = 11)[2:11],
      seq(0.95, 0.99, length.out = 11)[2:11],
      seq(0.99, 1, length.out = 21)[2:21]
    )
    ms <- findInterval(pl, pv[ord])

    calc_cost <- function(m) {
      if (m == 0) {
        c1 <- sum(q_kn) - sum(s_kn^2) / nkn
      } else {
        s <- s_kn + cum_s[m, ]
        q <- q_kn + cum_q[m, ]
        n <- nkn + m
        c1 <- sum(q) - sum(s^2) / n
      }
      if (m == nuk) {
        c2 <- 0
      } else {
        s <- s_uk - cum_s[m, ]
        q <- q_uk - cum_q[m, ]
        n2 <- nuk - m
        c2 <- sum(q) - sum(s^2) / n2
      }
      return(c1 + c2)
    }
    ot <- vapply(ms, calc_cost, numeric(1))

    pl_chosen <- pl[which.min(ot)[1]]

    pkfinal <- pv <= pl_chosen
    return(pkfinal)
  }

  uni_clus_ukn <- unique(cluster_ukn_helper)

  nc <- ncol(design_mat)
  cells_kn <- which(!pks)
  cells_pks <- which(pks)

  for (ucukn in uni_clus_ukn) {
    #this part could be further modified to account for bad quality clusters
    if (
      sum(cluster_kn == ucukn) <=
        max(
          1,
          round(
            unknown_threshold2 *
              (sum(cluster_kn == ucukn) + sum(cluster_ukn_helper == ucukn))
          )
        )
    ) {
      cluster_ukn[cluster_ukn_helper == ucukn] <- nc + ucukn
    } else if (sum(cluster_ukn_helper == ucukn) == 1) {
      cluster_ukn[cluster_ukn_helper == ucukn] <- nc + ucukn
    } else {
      uknrt <- unknown_detector(
        dat[(lmk + 1):ldt, cells_kn[cluster_kn == ucukn], drop = F],
        dat[(lmk + 1):ldt, cells_pks[cluster_ukn_helper == ucukn], drop = F]
      )
      cluster_ukn[cluster_ukn_helper == ucukn][!uknrt] <- nc + ucukn
    }
  }

  cluster[pks] <- cluster_ukn

  uniclus <- unique(cluster)

  for (allclus in 1:max(uniclus)) {
    if (allclus > nc) {
      design_mat[, paste('Unknown_', allclus - nc, sep = '')] <- 0
      if (allclus %in% uniclus) {
        mu <- cbind(
          mu,
          c(
            rep(0, lmk),
            rowMeans(dat[(lmk + 1):ldt, cells_pks[cluster_ukn == allclus], drop = F])
          )
        )
      } else {
        mu <- cbind(mu, rep(0, nrow(mu)))
      }
    } else {
      if (allclus %in% uniclus) {
        mu[(lmk + 1):ldt, allclus] <- rowMeans(dat[
          (lmk + 1):ldt,
          cluster == allclus,
          drop = F
        ])
      } else {
        mu[(lmk + 1):ldt, allclus] <- rep(0, ldt - lmk)
      }
    }
  }

  cost <- cost_func(dat, cluster, mu, design_mat, rt_sq)
  return(list(mu, cluster, cost))
}
