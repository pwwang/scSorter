#' scSorter
#'
#' This is the main function that implements the scSorter method.
#'
#' @param expr A matrix of the input expression data. Each row represents a gene and each column represents a cell. Each row of this matrix should be named by the gene name it represents.
#' @param anno A matrix or data frame that contains marker genes specified for cell types of interest.
#' It should contain three columns named "Type", "Marker", and "Weight" that records the name and weight of marker genes specified for each cell type.
#' "Weight" column is optional. If it is not specified, the \code{default_weight} will be applied to all marker genes.
#' @param default_weight The default weight assigned to marker genes. The default value is 2.
#' @param n_start The number of possible cluster initializations. The default value is 10.
#' @param alpha The parameter determines the cutoff whether the cell type of a cell should be considered as undecided during unknown cell calling. The default value is 0.
#' @param u The parameter determines whether undecided cells are further processed. The default value is 0.05.
#' @param max_iter The maximum number of iterations for the algorithm to update parameters. The default value is 100.
#' @param setseed Random seed for cluster initialization. The default value is 0.
#' @param mc.cores The number of cores to use for parallel processing. If NULL, the function will use the maximum number of available cores. The default value is NULL.
#'
#' @return A list contains the elements:
#'  \code{Pred_Type}: The predicted cell types.
#'  \code{Pred_param}: The parameter estimates of \code{mu} and \code{delta}.
#'
#' @export
scSorter <- function(
  expr,
  anno,
  default_weight = 2,
  n_start = 10,
  alpha = 0,
  u = 0.05,
  max_iter = 100,
  setseed = 0,
  mc.cores = NULL
) {
  #this is a wrapper function that implements the whole method based on the rest functions.
  #Rfast package is needed to run this method.
  message("[scSorter] Building the design matrix based on the provided marker genes...")
  anno_processed <- design_matrix_builder(anno, default_weight)

  message("[scSorter] Preprocessing the expression data...")
  dt <- data_preprocess(expr, anno_processed)

  dat <- dt[[1]]
  designmat <- dt[[2]]
  weightmat <- dt[[3]]

  # The n_start initializations are independent and each is seeded, so they
  # run in parallel when multiple cores are available. Set
  # options(mc.cores = 1) for a serial run. Results are identical either way.
  if (is.null(mc.cores)) {
    mc.cores <- max(1, min(n_start, getOption("mc.cores", parallel::detectCores())))
  }
  message(paste0("[scSorter] Running ", n_start, " initializations on ", mc.cores, " core(s)..."))

  if (mc.cores == 1) {
    pred_ots <- lapply(
      1:n_start,
      .scsorter_run_one,
      dat = dat, designmat = designmat, weightmat = weightmat,
      alpha = alpha, u = u, max_iter = max_iter, setseed = setseed
    )
  } else {
    # PSOCK workers are fresh R processes that read the BLAS thread-count env
    # vars at startup, so set the cap in the parent before they launch. Forked
    # workers would inherit the parent's already-initialized multithreaded
    # BLAS pool and oversubscribe the machine (n_start workers x 32 BLAS
    # threads), which runs slower than serial.
    Sys.setenv(
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      OMP_NUM_THREADS = "1"
    )
    cl <- parallel::makePSOCKcluster(mc.cores, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(scSorter))
    # Export the big matrices once per worker; the per-task function below
    # stays small so parLapply does not re-serialize 160MB for every task.
    parallel::clusterExport(cl, c("dat", "designmat", "weightmat"), envir = environment())
    pred_ots <- parallel::parLapply(
      cl,
      1:n_start,
      .scsorter_par_one,
      alpha = alpha, u = u, max_iter = max_iter, setseed = setseed
    )
  }

  c_cost <- vapply(pred_ots, function(x) x[[3]], numeric(1))
  c_mu <- lapply(pred_ots, function(x) x[[1]])
  c_clus <- lapply(pred_ots, function(x) x[[2]])

  pk <- which.min(c_cost)

  pred_clus <- c_clus[[pk]]
  pred_clus <- c(colnames(designmat), rep('Unknown', ncol(designmat)))[
    pred_clus
  ]
  pred_mu <- c_mu[[pk]]

  return(list(Pred_Type = pred_clus, Pred_param = pred_mu))
}


# Run one seeded initialization; defined at package level so that it
# serializes small (no captured 160MB expression matrix).
.scsorter_run_one <- function(i, dat, designmat, weightmat, alpha, u, max_iter, setseed) {
  set.seed(i + setseed)
  t1 <- Sys.time()
  pred_ot <- update_func(
    as.matrix(dat),
    designmat,
    weightmat,
    unknown_threshold1 = alpha,
    unknown_threshold2 = u,
    max_iter = max_iter
  )
  t2 <- Sys.time()
  message(paste0("[scSorter] Initialization ", i, " completed in ", round(difftime(t2, t1, units = "secs"), 2), " seconds."))
  return(pred_ot)
}

# Worker entry for the PSOCK cluster: the expression data is read from the
# worker's global environment where clusterExport put it.
.scsorter_par_one <- function(i, alpha, u, max_iter, setseed) {
  e <- .GlobalEnv
  return(.scsorter_run_one(i, e$dat, e$designmat, e$weightmat, alpha, u, max_iter, setseed))
}

#' Run scSorter on a Seurat object
#'
#' This function runs the scSorter method on a Seurat object.
#'
#' @inheritParams scSorter
#' @param object A Seurat object containing the expression data.
#' `FindVariableFeatures` should be run on the Seurat object before using this function.
#' @param layer The name of the layer in the Seurat object to use for expression data. If NULL, the default data will be used. The default value is NULL.
#' @param assay The name of the assay in the Seurat object to use for expression data. If NULL, the default assay will be used. The default value is NULL.
#' @param top_vf The number of top variable features to use for scSorter. If NULL, all variable features will be used. The default value is NULL.
#' @param min_pct The minimum fraction of cells that must have non-zero expression for a gene to be retained. The default value is 0.1.
#' @param set_ident A logical value indicating whether to set the predicted cell types as the active identity class in the Seurat object. The default value is TRUE.
#' @param name The name of the metadata column to store the predicted cell types. The default value is "scSorter_celltype".
#' @param mc.cores The number of cores to use for parallel processing. If NULL, the function will use the maximum number of available cores. The default value is NULL.
#' @param ... Additional arguments to pass to the scSorter function.
#' @return The Seurat object with the predicted cell types added to the metadata.
#' @importFrom SeuratObject VariableFeatures GetAssayData
#' @export
RunScSorter <- function(
  object,
  anno,
  layer = NULL,
  assay = NULL,
  top_vf = 2000,
  min_pct = 0.1,
  set_ident = TRUE,
  name = "scSorter_celltype",
  mc.cores = NULL,
  ...
) {
  if (!inherits(object, "Seurat")) {
    stop("The input object must be a Seurat object.")
  }

  message("[RunScSorter] Fetching variable features from the Seurat object...")
  if (is.null(top_vf)) {
    hvf <- VariableFeatures(object)
  } else {
    hvf <- utils::head(VariableFeatures(object), top_vf)
  }

  message("[RunScSorter] Fetching expression data from the Seurat object...")
  expr <- GetAssayData(object, layer = layer, assay = assay)

  message("[RunScSorter] Filtering genes based on minimum expression percentage...")
  hvf_filter <- rowSums(as.matrix(expr)[hvf, ] != 0) > ncol(expr) * min_pct
  hvf <- hvf[hvf_filter]

  message("[RunScSorter] Running scSorter on the filtered expression data...")
  picked_genes <- unique(c(anno$Marker, hvf))
  expr <- expr[intersect(rownames(expr), picked_genes), , drop = FALSE]
  result <- scSorter(expr, anno, mc.cores = mc.cores, ...)

  message("[RunScSorter] Adding predicted cell types to the Seurat object metadata...")
  object@meta.data[[name]] <- factor(
    result$Pred_Type,
    levels = sort(unique(result$Pred_Type))
  )

  if (set_ident) {
    message("[RunScSorter] Setting the predicted cell types as the active identity class in the Seurat object...")
    SeuratObject::Idents(object) <- name
  }
  return(object)
}
