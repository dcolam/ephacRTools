#' @importFrom magrittr %>%
NULL

# ============================================================================
# Multi-modal unsupervised clustering pipeline for HT-APC data
# ============================================================================
#
# Pipeline:
#   prepareClusterFeatures()   — scaled, modality-weighted feature matrix
#   optimalClusters()          — silhouette / WSS / gap to suggest k
#   clusterSE()                — kmeans / hierarchical / Louvain / GMM
#   reducedDimMultimodal()     — PCA + UMAP + t-SNE with modality weighting
#   plotDimRed()               — plot reducedDims with colData colouring
#   clusterHeatmap()           — ComplexHeatmap: features × wells, split by cluster
#   clusterSummary()           — per-cluster boxplots + summary table
#   seMOFAFeatures()           — build scaled view matrices for MOFA2
#   fitMOFASE()                — fit MOFA2 model, return mofa_fit_se object
#   predictMOFASE()            — write factor scores to colData / reducedDims
#   mofaWeights()              — feature weights per view/factor (→ ldaLoadings)
#   mofaVariance()             — R² per view per factor  (→ ldaPCAVariance)
#   clusterMOFA()              — deprecated one-shot wrapper (back-compat)
#   clusterPipeline()          — all-in-one wrapper

# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

.scale_zscore_col <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - mean(x, na.rm = TRUE)) / s
}

.scale_minmax_col <- function(x) {
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - r[1]) / diff(r)
}

.scale_block <- function(mat, method) {
  switch(method,
    zscore = apply(mat, 2, .scale_zscore_col),
    minmax = apply(mat, 2, .scale_minmax_col),
    none   = mat
  )
}

# Impute a matrix BEFORE scaling so imputed values participate in the
# distribution calculation. "impute_zero" puts true zeros into the raw
# distribution (correct for morphology channels where 0 = no signal).
# "impute_median" uses the raw-space median (robust center before outliers
# are compressed by scaling). "omit" is a no-op here — handled post-combine.
.impute_block_pre_scale <- function(mat, action) {
  if (action == "impute_zero") {
    mat[is.na(mat)] <- 0
  } else if (action == "impute_median") {
    for (j in seq_len(ncol(mat))) {
      na_j <- is.na(mat[, j])
      if (any(na_j)) {
        med <- median(mat[!na_j, j], na.rm = TRUE)
        if (is.na(med)) {
          warning(".impute_block_pre_scale: column '", colnames(mat)[j],
                  "' is entirely NA — setting to 0.")
          med <- 0
        }
        mat[na_j, j] <- med
      }
    }
  }
  mat
}

# Discrete colour palette for cluster labels (Tableau-10 inspired)
.cluster_palette <- function(n) {
  base_pal <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                "#9C755F", "#BAB0AC")
  if (n <= length(base_pal)) return(base_pal[seq_len(n)])
  grDevices::colorRampPalette(base_pal)(n)
}


# ============================================================================
# prepareClusterFeatures
# ============================================================================

#' Build a scaled, modality-weighted feature matrix from a SingleCellExperiment
#'
#' Extracts electrophysiology features (from \code{colData} scalar columns
#' and/or raw assay matrices) and morphology features (from \code{colData}),
#' scales each modality independently, then rescales the two blocks so their
#' total variance contributions match the requested weights before
#' concatenating into a single matrix for clustering.
#'
#' @param se A \code{SingleCellExperiment}.
#' @param ephys_cols Character vector of \code{colData} column names for ephys
#'   scalar features (e.g. \code{"Imax.minima"}, \code{"Vhalf.minima"},
#'   \code{"Capacitance_mean"}).  \code{NULL} skips this block.
#' @param morpho_cols Character vector of \code{colData} column names for
#'   morphology / imaging features (e.g. \code{"GFP_Mean_z"},
#'   \code{"GFP_normArea"}).  \code{NULL} skips this block.
#' @param assay_names Character vector of assay names whose full sweep × well
#'   matrices are used as additional ephys features (one feature per sweep).
#'   \code{NULL} (default) skips assay views.
#' @param ephys_weight Relative variance contribution of the ephys block.
#'   Default \code{1}.
#' @param morpho_weight Relative variance contribution of the morphology block.
#'   Default \code{1}.  Setting \code{morpho_weight = 2} doubles the morphology
#'   block's influence relative to ephys.
#' @param scale_method Scaling applied \emph{within} each modality block before
#'   weighting.  \code{"zscore"} (default), \code{"minmax"} (0–1), or
#'   \code{"none"}.
#' @param na_action Default NA strategy applied to all blocks that do not have
#'   a modality-specific override (see \code{ephys_na_action} /
#'   \code{morpho_na_action}).  \code{"impute_median"} (default) replaces each
#'   \code{NA} with the column median; \code{"impute_zero"} replaces with 0;
#'   \code{"omit"} removes wells that have any \code{NA} (the returned matrix
#'   may have fewer rows than \code{ncol(se)}).
#'   \strong{Important:} imputation is applied \emph{before} scaling so that
#'   imputed values participate correctly in the distribution calculation.
#'   With z-score scaling and \code{"impute_zero"}, a zero is placed in the raw
#'   distribution and will end up below the mean — which is the correct
#'   biological position for a morphology channel where zero means no signal.
#' @param ephys_na_action NA strategy for the ephys block (colData scalars +
#'   assay matrices).  \code{NULL} (default) falls back to \code{na_action}.
#'   Typical choice: \code{"impute_median"} (a failed patch has no
#'   interpretable zero) or \code{"omit"} after pre-filtering with
#'   \code{subsetSE()}.
#' @param morpho_na_action NA strategy for the morphology block.  \code{NULL}
#'   (default) falls back to \code{na_action}.  Typical choice:
#'   \code{"impute_zero"} — a well with no detected particles genuinely has
#'   zero fluorescence signal, and that zero should sit below the distribution
#'   mean after z-scoring.
#'
#' @return A numeric matrix with wells as rows and features as columns.
#'   Ephys columns are prefixed \code{"e__"}, morphology columns \code{"m__"}.
#'   The attribute \code{"modality"} is a named character vector mapping each
#'   column to \code{"ephys"} or \code{"morpho"}.
#'
#' @examples
#' \dontrun{
#' se <- get_metric(se, assay = "Minima")
#' se <- colAG(se, c("Capacitance", "Seal_Resistance"))
#' mat <- prepareClusterFeatures(
#'   se,
#'   ephys_cols    = c("Imax.minima", "Vhalf.minima", "Capacitance_mean"),
#'   morpho_cols   = c("GFP_Mean_z", "GFP_normArea",
#'                     "mCherry_Mean_z", "mCherry_normArea"),
#'   ephys_weight  = 1,
#'   morpho_weight = 1
#' )
#' dim(mat)   # wells × features
#' }
#'
#' @importFrom SummarizedExperiment colData assay assayNames
#' @export
prepareClusterFeatures <- function(se,
                                    ephys_cols       = NULL,
                                    morpho_cols      = NULL,
                                    assay_names      = NULL,
                                    ephys_weight     = 1,
                                    morpho_weight    = 1,
                                    scale_method     = c("zscore", "minmax", "none"),
                                    na_action        = c("impute_median",
                                                         "impute_zero", "omit"),
                                    ephys_na_action  = NULL,
                                    morpho_na_action = NULL) {
  scale_method <- match.arg(scale_method)
  na_action    <- match.arg(na_action)

  valid_na_actions <- c("impute_median", "impute_zero", "omit")
  if (!is.null(ephys_na_action))
    ephys_na_action  <- match.arg(ephys_na_action,  valid_na_actions)
  if (!is.null(morpho_na_action))
    morpho_na_action <- match.arg(morpho_na_action, valid_na_actions)

  # Resolve effective per-block NA action (modality-specific overrides na_action)
  eff_ephys_na  <- if (!is.null(ephys_na_action))  ephys_na_action  else na_action
  eff_morpho_na <- if (!is.null(morpho_na_action)) morpho_na_action else na_action

  cd <- as.data.frame(SummarizedExperiment::colData(se))

  # Keep only atomic (non-nested) columns to avoid silent coercion failures
  atomic_cols <- names(cd)[vapply(cd, function(x) is.atomic(x) || is.factor(x), logical(1))]

  blocks <- list()

  # --- ephys block: colData scalars -----------------------------------------
  if (!is.null(ephys_cols) && length(ephys_cols) > 0) {
    valid <- intersect(ephys_cols, atomic_cols)
    missing <- setdiff(ephys_cols, atomic_cols)
    if (length(missing) > 0)
      warning("prepareClusterFeatures: skipping ephys_cols not found or nested: ",
              paste(missing, collapse = ", "))
    if (length(valid) > 0) {
      mat_e <- matrix(as.numeric(as.matrix(cd[, valid, drop = FALSE])),
                      nrow = nrow(cd), ncol = length(valid),
                      dimnames = list(rownames(cd), valid))
      if (eff_ephys_na != "omit")
        mat_e <- .impute_block_pre_scale(mat_e, eff_ephys_na)
      mat_e <- .scale_block(mat_e, scale_method)
      colnames(mat_e) <- paste0("e__", valid)
      blocks[["ephys_cd"]] <- list(mat = mat_e, modality = "ephys")
    }
  }

  # --- ephys block: full assay matrices -------------------------------------
  if (!is.null(assay_names) && length(assay_names) > 0) {
    valid_a <- intersect(assay_names, SummarizedExperiment::assayNames(se))
    if (length(valid_a) == 0)
      warning("prepareClusterFeatures: none of assay_names found in SE")
    for (a in valid_a) {
      mat_a <- t(as.matrix(SummarizedExperiment::assay(se, a)))  # wells × sweeps
      storage.mode(mat_a) <- "double"
      if (eff_ephys_na != "omit")
        mat_a <- .impute_block_pre_scale(mat_a, eff_ephys_na)
      mat_a <- .scale_block(mat_a, scale_method)
      colnames(mat_a) <- paste0("e__", a, "_s", seq_len(ncol(mat_a)))
      blocks[[paste0("assay_", a)]] <- list(mat = mat_a, modality = "ephys")
    }
  }

  # --- morpho block ---------------------------------------------------------
  if (!is.null(morpho_cols) && length(morpho_cols) > 0) {
    valid_m <- intersect(morpho_cols, atomic_cols)
    missing_m <- setdiff(morpho_cols, atomic_cols)
    if (length(missing_m) > 0)
      warning("prepareClusterFeatures: skipping morpho_cols not found or nested: ",
              paste(missing_m, collapse = ", "))
    if (length(valid_m) > 0) {
      mat_m <- matrix(as.numeric(as.matrix(cd[, valid_m, drop = FALSE])),
                      nrow = nrow(cd), ncol = length(valid_m),
                      dimnames = list(rownames(cd), valid_m))
      if (eff_morpho_na != "omit")
        mat_m <- .impute_block_pre_scale(mat_m, eff_morpho_na)
      mat_m <- .scale_block(mat_m, scale_method)
      colnames(mat_m) <- paste0("m__", valid_m)
      blocks[["morpho"]] <- list(mat = mat_m, modality = "morpho")
    }
  }

  if (length(blocks) == 0)
    stop("prepareClusterFeatures: no valid features found.",
         " Check ephys_cols, morpho_cols, and assay_names.")

  # --- modality weighting ---------------------------------------------------
  # Rescale each block so its total variance contribution = weight × n_cols,
  # which makes blocks with equal weight contribute equally per feature on average.
  only_one_modality <- length(unique(vapply(blocks, `[[`, character(1), "modality"))) == 1
  if (!only_one_modality && scale_method != "none") {
    for (nm in names(blocks)) {
      mod <- blocks[[nm]]$modality
      w   <- if (mod == "ephys") ephys_weight else morpho_weight
      m   <- blocks[[nm]]$mat
      tv  <- sum(apply(m, 2, var, na.rm = TRUE))
      if (!is.na(tv) && tv > 0) {
        blocks[[nm]]$mat <- m * sqrt(w * ncol(m) / tv)
      }
    }
  }

  mat_list <- lapply(blocks, `[[`, "mat")
  modality_vec <- unlist(lapply(blocks, function(b) {
    setNames(rep(b$modality, ncol(b$mat)), colnames(b$mat))
  }))

  combined <- do.call(cbind, mat_list)
  rownames(combined) <- colnames(se)

  # --- Post-combine NA check and omit handling --------------------------------
  # Imputation already happened per-block before scaling.
  # Any remaining NAs here come from blocks using na_action="omit" (pre-scale
  # imputation was skipped for those blocks), or from degenerate all-NA columns.
  na_per_col <- colSums(is.na(combined))
  if (any(na_per_col > 0)) {
    bad_cols <- names(na_per_col[na_per_col > 0])
    message("prepareClusterFeatures: NAs remain in ",
            length(bad_cols), " feature(s) after pre-scale imputation:\n",
            paste0("  ", bad_cols, ": ", na_per_col[bad_cols],
                   "/", nrow(combined), " wells NA", collapse = "\n"))
  }

  # "omit" strategy: drop any well with a remaining NA across all blocks
  if (anyNA(combined)) {
    n_before <- nrow(combined)
    combined  <- combined[stats::complete.cases(combined), , drop = FALSE]
    n_dropped <- n_before - nrow(combined)
    if (n_dropped > 0)
      message("prepareClusterFeatures: dropped ", n_dropped, "/", n_before,
              " wells with remaining NA features.\n",
              "  These wells will have NA cluster labels.\n",
              "  Use ephys_na_action/morpho_na_action='impute_median' or\n",
              "  pre-filter with subsetSE() to keep all wells.")
  }

  # Final guard — should not be reachable
  if (anyNA(combined)) {
    stop("prepareClusterFeatures: NA values remain after all imputation steps. ",
         "Please inspect your features for fully-NA columns.")
  }

  attr(combined, "modality") <- modality_vec
  combined
}


# ============================================================================
# optimalClusters
# ============================================================================

#' Evaluate cluster quality across a range of k
#'
#' Runs k-means for each \code{k} in \code{k_range} and computes one or more
#' quality metrics to help select the number of clusters:
#' \describe{
#'   \item{silhouette}{Average silhouette width — higher is better
#'     (well-separated clusters).  The k with the highest value is suggested.}
#'   \item{wss}{Total within-cluster sum of squares — look for an elbow
#'     (inflection point).  The elbow is detected via the second derivative.}
#'   \item{gap}{Gap statistic (Tibshirani 2001) vs. a uniform reference
#'     distribution — higher is better.  Uses the \dQuote{firstSEmax} rule.
#'     \strong{Note:} requires \code{n_boot} k-means runs and can be slow for
#'     large datasets or large \code{k_range}.}
#' }
#'
#' @param feature_mat Numeric matrix (wells × features), typically from
#'   \code{\link{prepareClusterFeatures}}.
#' @param k_range Integer vector of k values to evaluate.  Default \code{2:8}.
#' @param methods Character vector of quality metrics.  Default
#'   \code{c("silhouette", "wss")}; add \code{"gap"} explicitly (slow).
#' @param n_boot Number of bootstrap replicates for the gap statistic.
#'   Default \code{50}.
#' @param seed Random seed.  Default \code{42}.
#'
#' @return A named list:
#'   \describe{
#'     \item{scores}{data.frame with columns \code{k}, \code{metric},
#'       \code{value}.}
#'     \item{plot}{A \code{ggplot} showing all metrics vs. k with vertical
#'       dashed lines at the suggested k for each metric.}
#'     \item{suggested_k}{Named integer vector — one suggested k per method.}
#'   }
#'
#' @examples
#' \dontrun{
#' mat <- prepareClusterFeatures(se,
#'   ephys_cols  = c("Imax.minima", "Vhalf.minima"),
#'   morpho_cols = c("GFP_Mean_z", "GFP_normArea"))
#' opt <- optimalClusters(mat, k_range = 2:8)
#' opt$plot
#' opt$suggested_k
#' }
#'
#' @importFrom cluster silhouette clusGap maxSE
#' @export
optimalClusters <- function(feature_mat,
                             k_range = 2:8,
                             methods = c("silhouette", "wss"),
                             n_boot  = 50,
                             seed    = 42) {
  set.seed(seed)

  .km <- function(mat, k) {
    stats::kmeans(mat, centers = k, nstart = 25, iter.max = 200)
  }

  results  <- list()
  suggested <- integer(0)

  # ---- Silhouette ----------------------------------------------------------
  if ("silhouette" %in% methods) {
    d <- stats::dist(feature_mat)
    sil_vals <- vapply(k_range, function(k) {
      km  <- .km(feature_mat, k)
      sil <- cluster::silhouette(km$cluster, d)
      mean(sil[, "sil_width"])
    }, numeric(1))
    results[["silhouette"]] <- data.frame(
      k = k_range, metric = "Silhouette width", value = sil_vals,
      stringsAsFactors = FALSE
    )
    suggested["silhouette"] <- k_range[which.max(sil_vals)]
  }

  # ---- WSS (elbow) ---------------------------------------------------------
  if ("wss" %in% methods) {
    wss_raw <- vapply(k_range, function(k) .km(feature_mat, k)$tot.withinss,
                      numeric(1))
    # Normalise to [0, 1] so the y-axis is meaningful alongside other metrics
    r        <- range(wss_raw)
    wss_norm <- if (diff(r) > 0) (wss_raw - r[1]) / (r[2] - r[1]) else
      rep(0, length(wss_raw))
    results[["wss"]] <- data.frame(
      k = k_range, metric = "WSS (normalised)", value = wss_norm,
      stringsAsFactors = FALSE
    )
    # Elbow: largest second-derivative magnitude
    if (length(k_range) >= 3) {
      d2  <- diff(diff(wss_raw))
      suggested["wss"] <- k_range[which.max(abs(d2)) + 1L]
    } else {
      suggested["wss"] <- k_range[1L]
    }
  }

  # ---- Gap statistic -------------------------------------------------------
  if ("gap" %in% methods) {
    message("optimalClusters: computing gap statistic (n_boot = ", n_boot,
            ") — this may take a moment...")
    gap_obj <- cluster::clusGap(
      feature_mat,
      FUNcluster = function(x, k) list(cluster = .km(x, k)$cluster),
      K.max   = max(k_range),
      B       = n_boot,
      verbose = FALSE
    )
    gap_tab  <- as.data.frame(gap_obj$Tab)
    gap_vals <- gap_tab$gap[k_range]
    results[["gap"]] <- data.frame(
      k = k_range, metric = "Gap statistic", value = gap_vals,
      stringsAsFactors = FALSE
    )
    opt_k <- tryCatch(
      cluster::maxSE(gap_obj$Tab[, "gap"],
                     gap_obj$Tab[, "SE.sim"],
                     method = "firstSEmax"),
      error = function(e) k_range[which.max(gap_vals)]
    )
    suggested["gap"] <- as.integer(opt_k)
  }

  scores_df <- do.call(rbind, results)
  rownames(scores_df) <- NULL

  # ---- Plot ----------------------------------------------------------------
  vline_df <- data.frame(
    metric = names(suggested),
    xint   = as.integer(suggested),
    stringsAsFactors = FALSE
  )
  # Map internal names to display names
  metric_map <- c(silhouette = "Silhouette width",
                  wss        = "WSS (normalised)",
                  gap        = "Gap statistic")
  vline_df$metric <- metric_map[vline_df$metric]

  p <- ggplot2::ggplot(scores_df,
                       ggplot2::aes(x = .data$k, y = .data$value)) +
    ggplot2::geom_line(colour = "#4E79A7", linewidth = 0.8) +
    ggplot2::geom_point(colour = "#4E79A7", size = 2.5) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
    ggplot2::scale_x_continuous(breaks = k_range) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = "Number of clusters (k)", y = NULL,
                  title = "Cluster quality metrics",
                  subtitle = "Dashed line = suggested k") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (nrow(vline_df) > 0) {
    p <- p + ggplot2::geom_vline(
      data     = vline_df,
      mapping  = ggplot2::aes(xintercept = .data$xint),
      linetype = "dashed", colour = "grey40", linewidth = 0.5
    )
  }

  list(scores = scores_df, plot = p, suggested_k = suggested)
}


# ============================================================================
# clusterSE
# ============================================================================

#' Cluster wells and store labels in colData
#'
#' Applies a clustering algorithm to a feature matrix and writes the resulting
#' integer cluster labels as a factor into \code{colData(se)[[col_name]]}.
#'
#' @param se A \code{SingleCellExperiment}.
#' @param feature_mat Numeric matrix (wells × features) from
#'   \code{\link{prepareClusterFeatures}}.  Must have rownames matching
#'   \code{colnames(se)}.
#' @param method Clustering algorithm.  One of:
#'   \describe{
#'     \item{\code{"kmeans"}}{k-means (\code{nstart = 25}).  Requires \code{k}.}
#'     \item{\code{"hierarchical"}}{Ward-linkage hierarchical clustering +
#'       \code{\link[stats]{cutree}}.  Requires \code{k}.}
#'     \item{\code{"louvain"}}{Louvain community detection on a k-nearest-
#'       neighbour graph.  Does \emph{not} require \code{k}; the resolution
#'       is controlled implicitly by the graph topology.  Requires
#'       \pkg{igraph}.}
#'     \item{\code{"gmm"}}{Gaussian Mixture Model via \code{mclust::Mclust()}.
#'       \code{k = NULL} lets mclust choose via BIC.  Requires \pkg{mclust}.}
#'   }
#' @param k Integer number of clusters.  Required for \code{"kmeans"} and
#'   \code{"hierarchical"}; optional for \code{"gmm"}; ignored for
#'   \code{"louvain"}.
#' @param col_name Name of the colData column to write labels into.
#'   Default \code{"cluster"}.
#' @param hclust_method Linkage method for \code{stats::hclust()}.
#'   Default \code{"ward.D2"}.
#' @param louvain_k Number of nearest neighbours for the KNN graph used by
#'   Louvain.  Default \code{10}.
#' @param seed Integer random seed.  Default \code{42}.
#' @param ... Additional arguments passed to the underlying clustering function.
#'
#' @return The input \code{se} with cluster labels added to \code{colData}.
#'
#' @examples
#' \dontrun{
#' mat <- prepareClusterFeatures(se,
#'   ephys_cols  = c("Imax.minima", "Vhalf.minima"),
#'   morpho_cols = c("GFP_Mean_z", "GFP_normArea"))
#' se  <- clusterSE(se, mat, method = "kmeans", k = 4)
#' table(se$cluster)
#' }
#'
#' @export
clusterSE <- function(se, feature_mat,
                       method        = c("kmeans", "hierarchical",
                                         "louvain", "gmm"),
                       k             = NULL,
                       col_name      = "cluster",
                       hclust_method = "ward.D2",
                       louvain_k     = 10,
                       seed          = 42,
                       ...) {
  method <- match.arg(method)
  set.seed(seed)

  if (method %in% c("kmeans", "hierarchical") && is.null(k))
    stop("clusterSE: 'k' must be specified for method = '", method, "'")

  labels <- switch(method,

    kmeans = {
      km <- stats::kmeans(feature_mat, centers = k,
                          nstart = 25, iter.max = 200, ...)
      km$cluster
    },

    hierarchical = {
      d  <- stats::dist(feature_mat)
      hc <- stats::hclust(d, method = hclust_method)
      stats::cutree(hc, k = k)
    },

    louvain = {
      if (!requireNamespace("igraph", quietly = TRUE))
        stop("clusterSE: package 'igraph' is required for method = 'louvain'.\n",
             "  Install with: install.packages('igraph')")
      n   <- nrow(feature_mat)
      knn <- min(louvain_k, n - 1L)
      d   <- as.matrix(stats::dist(feature_mat))
      # Build directed KNN edge list, then symmetrise
      edges <- do.call(rbind, lapply(seq_len(n), function(i) {
        nbrs <- order(d[i, ])[2:(knn + 1L)]
        cbind(from = i, to = nbrs)
      }))
      g  <- igraph::graph_from_edgelist(edges, directed = FALSE)
      g  <- igraph::simplify(g)
      cl <- igraph::cluster_louvain(g, ...)
      as.integer(igraph::membership(cl))
    },

    gmm = {
      if (!requireNamespace("mclust", quietly = TRUE))
        stop("clusterSE: package 'mclust' is required for method = 'gmm'.\n",
             "  Install with: install.packages('mclust')")
      fit <- mclust::Mclust(feature_mat, G = k, verbose = FALSE, ...)
      fit$classification
    }
  )

  # Align labels to SE column order (handles na_action = "omit" case)
  all_labels <- rep(NA_integer_, ncol(se))
  names(all_labels) <- colnames(se)
  well_names <- rownames(feature_mat)
  if (!is.null(well_names)) {
    idx <- match(well_names, colnames(se))
    valid <- !is.na(idx)
    all_labels[idx[valid]] <- as.integer(labels[valid])
  } else {
    all_labels[seq_along(labels)] <- as.integer(labels)
  }

  SummarizedExperiment::colData(se)[[col_name]] <- as.factor(all_labels)
  se
}


# ============================================================================
# reducedDimMultimodal
# ============================================================================

#' Multi-modal dimensionality reduction with modality weighting
#'
#' Builds a modality-weighted feature matrix via
#' \code{\link{prepareClusterFeatures}}, runs PCA, and optionally t-SNE and
#' UMAP, storing results in \code{reducedDims(se)}.  Cluster labels are
#' computed on the PCA coordinates and stored in
#' \code{colData(se)$cluster}.
#'
#' This function supersedes the legacy \code{reducedDim.Cellwise()}.
#'
#' @param se A \code{SingleCellExperiment}.
#' @param ephys_cols colData columns for ephys scalar features.
#' @param morpho_cols colData columns for morphology features.
#' @param assay_names Assay names to include as full IV-curve feature blocks.
#' @param ephys_weight Relative variance contribution of the ephys block.
#'   Default \code{1}.
#' @param morpho_weight Relative variance contribution of the morphology block.
#'   Default \code{1}.
#' @param scale_method Feature scaling within each block.  Default
#'   \code{"zscore"}.
#' @param na_action NA handling.  Default \code{"impute_median"}.
#' @param n_pcs Maximum number of PCs to retain.  Default \code{50}.
#' @param method Character vector of reductions to compute.  Any combination
#'   of \code{"pca"}, \code{"tsne"}, \code{"umap"}.  Default all three.
#' @param k_clusters Number of clusters for the clustering step.  Default
#'   \code{3}.
#' @param cluster_method Clustering algorithm passed to \code{\link{clusterSE}}.
#'   Default \code{"kmeans"}.
#' @param seed Random seed.  Default \code{42}.
#'
#' @return \code{se} with \code{reducedDims} populated (\code{PCA}, \code{TSNE},
#'   \code{UMAP} as requested) and \code{colData(se)$cluster} set.
#'
#' @examples
#' \dontrun{
#' se <- reducedDimMultimodal(
#'   se,
#'   ephys_cols    = c("Imax.minima", "Vhalf.minima", "Capacitance_mean"),
#'   morpho_cols   = c("GFP_Mean_z", "GFP_normArea"),
#'   ephys_weight  = 1,
#'   morpho_weight = 1,
#'   k_clusters    = 4
#' )
#' }
#'
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom S4Vectors DataFrame
#' @export
reducedDimMultimodal <- function(se,
                                  ephys_cols       = NULL,
                                  morpho_cols      = NULL,
                                  assay_names      = NULL,
                                  ephys_weight     = 1,
                                  morpho_weight    = 1,
                                  scale_method     = c("zscore", "minmax", "none"),
                                  na_action        = c("impute_median",
                                                       "impute_zero", "omit"),
                                  ephys_na_action  = NULL,
                                  morpho_na_action = NULL,
                                  n_pcs            = 50,
                                  method           = c("pca", "tsne", "umap"),
                                  k_clusters       = 3,
                                  cluster_method   = c("kmeans", "hierarchical",
                                                       "louvain", "gmm"),
                                  seed             = 42) {
  scale_method   <- match.arg(scale_method)
  na_action      <- match.arg(na_action)
  cluster_method <- match.arg(cluster_method)
  set.seed(seed)

  mat <- prepareClusterFeatures(
    se,
    ephys_cols       = ephys_cols,
    morpho_cols      = morpho_cols,
    assay_names      = assay_names,
    ephys_weight     = ephys_weight,
    morpho_weight    = morpho_weight,
    scale_method     = scale_method,
    na_action        = na_action,
    ephys_na_action  = ephys_na_action,
    morpho_na_action = morpho_na_action
  )

  n_pcs_actual <- min(n_pcs, ncol(mat), nrow(mat) - 1L)
  pca_result   <- stats::prcomp(mat, rank. = n_pcs_actual,
                                 center = FALSE, scale. = FALSE)

  rd <- list()

  if ("pca" %in% method) {
    rd[["PCA"]] <- pca_result$x
    attr(rd[["PCA"]], "percentVar") <- (pca_result$sdev^2 /
                                          sum(pca_result$sdev^2) * 100)
  }

  n_input_pcs <- min(30L, ncol(pca_result$x))
  pca_top <- pca_result$x[, seq_len(n_input_pcs), drop = FALSE]

  if ("tsne" %in% method) {
    tsne_out <- Rtsne::Rtsne(pca_top, pca = FALSE,
                              check_duplicates = FALSE, verbose = FALSE)
    rd[["TSNE"]] <- S4Vectors::DataFrame(
      tsne1 = tsne_out$Y[, 1],
      tsne2 = tsne_out$Y[, 2]
    )
  }

  if ("umap" %in% method) {
    umap_out <- umap::umap(pca_top)
    rd[["UMAP"]] <- S4Vectors::DataFrame(
      UMAP1 = umap_out$layout[, 1],
      UMAP2 = umap_out$layout[, 2]
    )
  }

  # mat_wells / se_wells — determined from names, with a safe NULL guard
  mat_wells <- rownames(mat)     # NULL when colnames(se) was NULL
  se_wells  <- colnames(se)

  # Bug guard: if mat_wells is NULL we cannot do name-based alignment.
  # Fall back to positional alignment (requires mat to have nrow == ncol(se)).
  if (is.null(mat_wells)) {
    if (nrow(mat) != length(se_wells))
      stop("reducedDimMultimodal: feature matrix has no rownames and its row ",
           "count (", nrow(mat), ") does not match ncol(se) (", length(se_wells),
           "). Set colnames(se) before calling this function.")
    mat_wells <- se_wells   # treat as positionally aligned
    rownames(mat) <- mat_wells
  }

  # Fill reducedDims to full SE size when na_action = "omit" dropped some wells
  omit_wells <- setdiff(se_wells, mat_wells)
  if (length(omit_wells) > 0) {
    idx <- match(mat_wells, se_wells)
    for (nm in names(rd)) {
      m    <- as.matrix(rd[[nm]])
      full <- matrix(NA_real_, nrow = length(se_wells), ncol = ncol(m),
                     dimnames = list(se_wells, colnames(m)))
      full[idx, ] <- m
      rd[[nm]] <- S4Vectors::DataFrame(full)
    }
  }

  SingleCellExperiment::reducedDims(se) <- rd

  # Store PCA model so plotFeatureImportance() can access loadings later
  S4Vectors::metadata(se)[["PCA_model"]]       <- pca_result
  S4Vectors::metadata(se)[["feature_modality"]] <- attr(mat, "modality")

  # Cluster on the wells that are in mat (i.e. had complete features).
  # Extract their PCA coordinates directly from pca_result$x (never from the
  # filled DataFrame, which may have NA rows for omitted wells).
  pca_for_clust <- pca_result$x   # already has rownames = mat_wells

  # Sanity check: warn if any NaN/Inf crept in from degenerate features
  if (anyNA(pca_for_clust) || any(!is.finite(pca_for_clust))) {
    n_bad <- sum(!is.finite(pca_for_clust))
    warning("reducedDimMultimodal: PCA produced ", n_bad,
            " non-finite value(s). This usually means one or more feature ",
            "columns had zero variance (constant after scaling). ",
            "Those wells will receive NA cluster labels.\n",
            "  Tip: remove constant features before clustering.")
    pca_for_clust <- pca_for_clust[stats::complete.cases(pca_for_clust), , drop = FALSE]
  }

  se <- clusterSE(se, pca_for_clust,
                   method   = cluster_method,
                   k        = if (cluster_method == "louvain") NULL else k_clusters,
                   col_name = "cluster",
                   seed     = seed)

  se
}


# ============================================================================
# plotDimRed  (moved from utilityTools.R, improved)
# ============================================================================

#' Plot dimensionality reduction results
#'
#' Renders scatter plots of \code{reducedDims} results, coloured by one or
#' more \code{colData} columns.  Each \code{colorColumns} entry produces a
#' separate panel; numeric columns use a viridis colour scale, categorical
#' columns use a discrete palette.
#'
#' @param se A \code{SingleCellExperiment} with \code{reducedDims} populated.
#' @param redDim.method One of \code{"PCA"}, \code{"TSNE"}, \code{"UMAP"}.
#' @param colorColumns Character vector of \code{colData} column names to use
#'   as colour variables.  When empty (default), the function looks for a
#'   column whose name contains the reduction name (e.g. \code{"cluster"}).
#' @param point_size Geom point size.  Default \code{2}.
#' @param alpha Point transparency.  Default \code{0.7}.
#'
#' @return A \code{ggplot} (single colour variable) or
#'   \code{ggpubr::ggarrange} panel grid (multiple colour variables).
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @export
plotDimRed <- function(se,
                        redDim.method,
                        colorColumns = character(),
                        point_size   = 2,
                        alpha        = 0.7) {

  cd    <- as.data.frame(SummarizedExperiment::colData(se))
  red   <- as.data.frame(SingleCellExperiment::reducedDim(se, redDim.method))
  ax1   <- colnames(red)[1]
  ax2   <- colnames(red)[2]

  # Default colour variable: column whose name contains the reduction name
  if (length(colorColumns) == 0) {
    auto <- grep(tolower(redDim.method), colnames(cd), value = TRUE, ignore.case = TRUE)
    colorColumns <- if (length(auto) > 0) auto[1] else colnames(cd)[1]
  }

  colorColumns <- intersect(colorColumns, colnames(cd))
  if (length(colorColumns) == 0) {
    warning("plotDimRed: none of colorColumns found in colData")
    return(invisible(NULL))
  }

  .one_plot <- function(col_var) {
    vals <- cd[[col_var]]
    df   <- cbind(red, colour_var = vals)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]],
                                           colour = .data$colour_var)) +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      ggplot2::ggtitle(col_var) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(colour = col_var)

    if (is.numeric(vals)) {
      p <- p + viridis::scale_colour_viridis(option = "plasma")
    }
    p
  }

  plots <- lapply(colorColumns, .one_plot)
  if (length(plots) == 1) plots[[1]] else ggpubr::ggarrange(plotlist = plots)
}


# ============================================================================
# clusterHeatmap
# ============================================================================

#' Heatmap of feature profiles per cluster
#'
#' Renders a \code{ComplexHeatmap} with features as rows and wells as columns,
#' split by cluster assignment.  A row annotation distinguishes ephys from
#' morphology features when the feature matrix carries a \code{"modality"}
#' attribute (as returned by \code{\link{prepareClusterFeatures}}).
#'
#' @param se A \code{SingleCellExperiment}.
#' @param feature_mat Scaled feature matrix (wells × features) from
#'   \code{\link{prepareClusterFeatures}}.
#' @param cluster_col Name of the \code{colData} column containing cluster
#'   labels.  Default \code{"cluster"}.
#' @param scale Logical.  If \code{TRUE} (default), z-score each feature row
#'   before plotting so that the colour scale reflects relative variation
#'   within each feature rather than absolute magnitude.
#' @param show_well_names Logical.  Show column labels (well IDs).  Default
#'   \code{FALSE} (too dense for full plates).
#' @param ... Additional arguments forwarded to \code{ComplexHeatmap::Heatmap}.
#'
#' @return A \code{Heatmap} object.  Print or call
#'   \code{ComplexHeatmap::draw()} to render.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @importFrom SummarizedExperiment colData
#' @export
clusterHeatmap <- function(se,
                            feature_mat,
                            cluster_col     = "cluster",
                            scale           = TRUE,
                            show_well_names = FALSE,
                            ...) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))

  # Align feature_mat rows to SE column order
  mat_wells <- rownames(feature_mat)
  se_wells  <- colnames(se)
  if (!is.null(mat_wells)) {
    idx     <- match(se_wells, mat_wells)
    valid   <- !is.na(idx)
    mat_ord <- feature_mat[idx[valid], , drop = FALSE]
    cd_ord  <- cd[valid, , drop = FALSE]
  } else {
    mat_ord <- feature_mat
    cd_ord  <- cd
  }

  mat_t <- t(mat_ord)  # features × wells

  if (scale) {
    mat_t <- t(apply(mat_t, 1, .scale_zscore_col))
    col_name_hm <- "z-score"
  } else {
    col_name_hm <- "value"
  }

  cluster_labels <- as.character(cd_ord[[cluster_col]])
  n_clust        <- length(unique(stats::na.omit(cluster_labels)))
  pal            <- .cluster_palette(n_clust)
  clust_colors   <- setNames(pal[seq_len(n_clust)],
                             sort(unique(stats::na.omit(cluster_labels))))

  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = cluster_labels,
    col     = list(Cluster = clust_colors),
    show_annotation_name = FALSE
  )

  # Optional row annotation for modality
  mod_attr   <- attr(feature_mat, "modality")
  row_annot  <- NULL
  feat_names <- rownames(mat_t)
  if (!is.null(mod_attr) && !is.null(feat_names)) {
    mod_vals <- mod_attr[feat_names]
    mod_vals[is.na(mod_vals)] <- "unknown"
    row_annot <- ComplexHeatmap::rowAnnotation(
      Modality = mod_vals,
      col      = list(Modality = c(ephys   = "#4E79A7",
                                   morpho  = "#F28E2B",
                                   unknown = "#BAB0AC")),
      show_annotation_name = FALSE,
      width = grid::unit(3, "mm")
    )
  }

  # Strip e__ / m__ prefixes for display
  display_names <- gsub("^[em]__", "", feat_names)

  ComplexHeatmap::Heatmap(
    mat_t,
    name              = col_name_hm,
    column_split      = cluster_labels,
    top_annotation    = col_annot,
    left_annotation   = row_annot,
    row_labels        = display_names,
    show_column_names = show_well_names,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    column_title      = paste0("Cluster (", cluster_col, ")"),
    col               = circlize::colorRamp2(
      c(-2, 0, 2),
      c("#4E79A7", "white", "#E15759")
    ),
    ...
  )
}


# ============================================================================
# clusterSummary
# ============================================================================

#' Per-cluster summary plots and statistics
#'
#' Produces a panel of box or violin plots (one facet per feature) coloured by
#' cluster, plus a tidy summary table of per-cluster mean ± SD for each
#' feature.
#'
#' @param se A \code{SingleCellExperiment}.
#' @param feature_cols Character vector of \code{colData} column names to
#'   summarise.
#' @param cluster_col Name of the \code{colData} column containing cluster
#'   labels.  Default \code{"cluster"}.
#' @param plot_type \code{"boxplot"} (default) or \code{"violin"}.
#' @param ncol Number of columns in the faceted plot grid.  Default \code{3}.
#'
#' @return A named list:
#'   \describe{
#'     \item{plot}{A \code{ggplot} faceted by feature.}
#'     \item{table}{A data.frame with per-cluster \code{n}, \code{mean}, and
#'       \code{sd} per feature.}
#'   }
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise n
#' @export
clusterSummary <- function(se,
                            feature_cols,
                            cluster_col = "cluster",
                            plot_type   = c("boxplot", "violin"),
                            ncol        = 3) {
  plot_type <- match.arg(plot_type)

  cd    <- as.data.frame(SummarizedExperiment::colData(se))
  valid <- intersect(feature_cols,
                     names(cd)[vapply(cd, function(x) is.atomic(x) || is.factor(x),
                                      logical(1))])
  if (length(valid) == 0)
    stop("clusterSummary: no valid feature_cols found in colData (check for nested DataFrames)")

  plot_df           <- cd[, c(cluster_col, valid), drop = FALSE]
  names(plot_df)[1] <- "cluster"
  plot_df$cluster   <- as.factor(plot_df$cluster)

  long_df        <- tidyr::pivot_longer(plot_df, cols = -cluster,
                                         names_to  = "feature",
                                         values_to = "value")
  long_df$value  <- as.numeric(long_df$value)

  n_clust <- nlevels(plot_df$cluster)
  pal     <- .cluster_palette(n_clust)

  geom_layer <- if (plot_type == "violin") {
    list(
      ggplot2::geom_violin(ggplot2::aes(fill = .data$cluster),
                           alpha = 0.6, colour = NA),
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA,
                            fill = "white", alpha = 0.8)
    )
  } else {
    list(ggplot2::geom_boxplot(ggplot2::aes(fill = .data$cluster),
                               alpha = 0.7, outlier.size = 0.8))
  }

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = .data$cluster, y = .data$value,
                                    fill = .data$cluster)) +
    geom_layer +
    ggplot2::geom_jitter(width = 0.18, size = 0.5,
                          alpha = 0.35, colour = "grey20") +
    ggplot2::facet_wrap(~feature, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Cluster", y = "Value") +
    ggplot2::theme(legend.position = "none")

  tbl <- long_df %>%
    dplyr::group_by(.data$cluster, .data$feature) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      mean = mean(.data$value,  na.rm = TRUE),
      sd   = sd(.data$value,    na.rm = TRUE),
      .groups = "drop"
    )

  list(plot = p, table = as.data.frame(tbl))
}


# ============================================================================
# MOFA2 multi-modal factor analysis — SE-level API
#
# Mirrors the ldaTools.R pattern:
#   seMOFAFeatures()  →  fitMOFASE()  →  predictMOFASE()
#                                     →  mofaWeights()
#                                     →  mofaVariance()
#
# Key differences vs the LDA pipeline:
#   • Views stay separate (features × wells per modality) rather than being
#     concatenated, so MOFA can handle cells missing from one view.
#   • "assay" scaling is omitted — MOFA's scale_views handles cross-view
#     variance; plate + feature steps handle batch drift and unit differences.
#   • NA values are kept by default ("pass") so MOFA fills them via the
#     generative model — the main advantage over LDA for partially-imaged data.
#   • predictMOFASE extracts already-fitted factor scores rather than
#     projecting new data (MOFA does not support out-of-sample projection).
# ============================================================================

#' MOFA2-based multi-modal factor analysis for SummarizedExperiment objects
#'
#' A modular pipeline for multi-modal factor analysis using MOFA2.  Mirrors
#' the \code{\link{ldaTools}} API:
#' \enumerate{
#'   \item \code{\link{seMOFAFeatures}} — build scaled view matrices
#'   \item \code{\link{fitMOFASE}} — fit the MOFA2 model
#'   \item \code{\link{predictMOFASE}} — write factor scores to \code{colData}
#'   \item \code{\link{mofaWeights}} — feature weights per view/factor
#'   \item \code{\link{mofaVariance}} — R² per view per factor
#' }
#'
#' @name mofaTools
NULL

.mofa_check <- function() {
  if (!requireNamespace("MOFA2", quietly = TRUE))
    stop("Package 'MOFA2' is required. ",
         "Install with: BiocManager::install('MOFA2')")
}

.scale_mofa_block <- function(mat, scale, plates) {
  # mat: wells × features
  plate_centers <- NULL
  if ("plate" %in% scale && !is.null(plates)) {
    upl <- unique(as.character(plates))
    plate_centers <- matrix(NA_real_, length(upl), ncol(mat),
                             dimnames = list(upl, colnames(mat)))
    for (p in upl) {
      idx <- which(as.character(plates) == p)
      plate_centers[p, ] <- colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      mat[idx, ] <- sweep(mat[idx, , drop = FALSE], 2, plate_centers[p, ], "-")
    }
  }
  feat_center <- feat_scale <- NULL
  if ("feature" %in% scale) {
    feat_center <- colMeans(mat, na.rm = TRUE)
    feat_scale  <- apply(mat, 2, sd, na.rm = TRUE)
    feat_scale[is.na(feat_scale) | feat_scale == 0] <- 1
    mat <- sweep(sweep(mat, 2, feat_center, "-"), 2, feat_scale, "/")
  }
  list(mat = mat, plate_centers = plate_centers,
       feat_center = feat_center, feat_scale = feat_scale)
}

.impute_mofa_block <- function(mat, action) {
  if (action == "impute_zero") {
    mat[is.na(mat)] <- 0
  } else if (action == "impute_median") {
    for (j in seq_len(ncol(mat))) {
      na_j <- is.na(mat[, j])
      if (any(na_j)) {
        med <- median(mat[!na_j, j], na.rm = TRUE)
        mat[na_j, j] <- if (is.na(med)) 0 else med
      }
    }
  }
  mat  # "pass" returns as-is; MOFA handles the NAs
}

# ── seMOFAFeatures ────────────────────────────────────────────────────────────

#' Build MOFA2 view matrices from a SummarizedExperiment
#'
#' Constructs separate \emph{ephys} and \emph{morphology} view matrices for
#' \code{MOFA2}, applying the same plate-centering and feature z-scoring as
#' \code{\link{seLDAFeatures}}.  Unlike the LDA pipeline, \code{NA} values
#' are preserved by default so MOFA2 can handle them natively — cells that
#' are missing from one view (e.g. border wells without imaging) still receive
#' factor scores driven by the other view.
#'
#' @param se A \code{SummarizedExperiment} or \code{SingleCellExperiment}.
#' @param ephys_cols Character vector of numeric \code{colData} columns for the
#'   ephys view.  \code{NULL} skips this view.
#' @param morpho_cols Character vector of numeric \code{colData} columns for
#'   the morphology view.  \code{NULL} skips this view.
#' @param assay_names Assay names to include as individual views (one view per
#'   assay, features = sweeps).  \code{NULL} (default) skips.
#' @param scale Character vector of scaling steps to apply \emph{within each
#'   view}: \code{"plate"} (subtract per-plate means, removes batch drift) and
#'   \code{"feature"} (z-score per feature, normalises units).  Default both.
#'   Cross-view variance balancing is handled by \code{scale_views} in
#'   \code{\link{fitMOFASE}}, not here.
#' @param plate_col \code{colData} column for plate IDs.  Default
#'   \code{"Plate_ID"}.  Matched to the default used by \code{\link{seLDAFeatures}}.
#' @param na_action \code{"pass"} (default — preserve \code{NA}s for MOFA),
#'   \code{"impute_median"}, or \code{"impute_zero"}.
#' @param ephys_na_action Override NA strategy for the ephys block.
#'   \code{NULL} falls back to \code{na_action}.
#' @param morpho_na_action Override NA strategy for the morphology block.
#'   \code{NULL} falls back to \code{na_action}.  Typical choice:
#'   \code{"impute_zero"} when missing particles mean zero signal, though
#'   \code{"pass"} lets MOFA learn this from context.
#' @param wells Optional subset of well names (\code{colnames(se)}).
#'   \code{NULL} uses all wells.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{views}}{Named list of \code{[features × wells]} matrices
#'     suitable for \code{MOFA2::create_mofa()}.}
#'   \item{\code{scale_params}}{Per-view scaling parameters (plate centers,
#'     feature means/SDs) for diagnostics.}
#'   \item{\code{well_names}}{Character vector of well IDs.}
#'   \item{\code{ephys_cols}, \code{morpho_cols}, \code{assay_names}}{Feature
#'     specification forwarded to the fit object and \code{mofaWeights()}.}
#' }
#'
#' @seealso \code{\link{fitMOFASE}}, \code{\link{predictMOFASE}},
#'   \code{\link{mofaWeights}}, \code{\link{mofaVariance}}
#'
#' @importFrom SummarizedExperiment colData assay assayNames
#' @export
seMOFAFeatures <- function(se,
                             ephys_cols       = NULL,
                             morpho_cols      = NULL,
                             assay_names      = NULL,
                             scale            = c("plate", "feature"),
                             plate_col        = "Plate_ID",
                             na_action        = c("pass", "impute_median",
                                                   "impute_zero"),
                             ephys_na_action  = NULL,
                             morpho_na_action = NULL,
                             wells            = NULL) {
  na_action <- match.arg(na_action)
  valid_na  <- c("pass", "impute_median", "impute_zero")
  if (!is.null(ephys_na_action))
    ephys_na_action  <- match.arg(ephys_na_action,  valid_na)
  if (!is.null(morpho_na_action))
    morpho_na_action <- match.arg(morpho_na_action, valid_na)
  eff_ephys_na  <- if (!is.null(ephys_na_action))  ephys_na_action  else na_action
  eff_morpho_na <- if (!is.null(morpho_na_action)) morpho_na_action else na_action

  if (is.null(ephys_cols) && is.null(morpho_cols) && is.null(assay_names))
    stop("Provide at least one of 'ephys_cols', 'morpho_cols', or 'assay_names'.")

  if (!is.null(wells)) {
    wells <- intersect(wells, colnames(se))
    if (!length(wells)) stop("No matching wells found.")
    se <- se[, wells, drop = FALSE]
  }

  cd          <- as.data.frame(SummarizedExperiment::colData(se))
  atomic_cols <- names(cd)[vapply(cd, function(x) is.atomic(x) || is.factor(x),
                                   logical(1))]
  plates      <- cd[[plate_col]]
  if (is.null(plates))
    warning("plate_col '", plate_col, "' not found in colData; skipping plate centering.")
  well_names  <- colnames(se)

  views     <- list()
  scale_lst <- list()

  .build_view <- function(cols, na_act, view_nm) {
    valid <- intersect(cols, atomic_cols)
    miss  <- setdiff(cols, atomic_cols)
    if (length(miss))
      warning("seMOFAFeatures: skipping ", view_nm, "_cols not found or nested: ",
              paste(miss, collapse = ", "))
    if (!length(valid)) return(invisible(NULL))
    mat <- matrix(as.numeric(as.matrix(cd[, valid, drop = FALSE])),
                  nrow = nrow(cd), ncol = length(valid),
                  dimnames = list(well_names, valid))
    if (na_act != "pass") mat <- .impute_mofa_block(mat, na_act)
    res <- .scale_mofa_block(mat, scale, plates)
    views[[view_nm]]     <<- t(res$mat)   # features × wells for MOFA
    scale_lst[[view_nm]] <<- res[c("plate_centers", "feat_center", "feat_scale")]
  }

  if (!is.null(ephys_cols)  && length(ephys_cols))
    .build_view(ephys_cols,  eff_ephys_na,  "ephys")
  if (!is.null(morpho_cols) && length(morpho_cols))
    .build_view(morpho_cols, eff_morpho_na, "morpho")

  # Assay views
  if (!is.null(assay_names) && length(assay_names)) {
    valid_a <- intersect(assay_names, SummarizedExperiment::assayNames(se))
    bad_a   <- setdiff(assay_names, valid_a)
    if (length(bad_a))
      warning("seMOFAFeatures: assays not found (skipped): ",
              paste(bad_a, collapse = ", "))
    for (a in valid_a) {
      mat_a <- as.matrix(SummarizedExperiment::assay(se, a))   # sweeps × wells
      storage.mode(mat_a) <- "double"
      colnames(mat_a) <- well_names
      mat_t <- t(mat_a)   # wells × sweeps for scaling
      if (eff_ephys_na != "pass") mat_t <- .impute_mofa_block(mat_t, eff_ephys_na)
      res_a <- .scale_mofa_block(mat_t, scale, plates)
      v_nm  <- paste0("ephys_", a)
      views[[v_nm]]     <- t(res_a$mat)   # sweeps × wells for MOFA
      scale_lst[[v_nm]] <- res_a[c("plate_centers", "feat_center", "feat_scale")]
    }
  }

  if (length(views) == 0)
    stop("seMOFAFeatures: no valid views built. ",
         "Check ephys_cols, morpho_cols, and assay_names.")

  list(
    views        = views,
    scale_params = list(per_view  = scale_lst,
                        steps     = scale,
                        plate_col = plate_col),
    well_names   = well_names,
    ephys_cols   = ephys_cols,
    morpho_cols  = morpho_cols,
    assay_names  = assay_names
  )
}


# ── fitMOFASE ─────────────────────────────────────────────────────────────────

#' Fit a MOFA2 multi-modal factor model on a SummarizedExperiment
#'
#' Takes the output of \code{\link{seMOFAFeatures}}, trains a MOFA2 model,
#' and returns a \code{mofa_fit_se} S3 object analogous to the
#' \code{lda_fit_se} returned by \code{\link{fitLDASE}}.
#'
#' @param se A \code{SummarizedExperiment} (used for well alignment only;
#'   the actual feature data is in \code{mofa_features}).
#' @param mofa_features Output of \code{\link{seMOFAFeatures}}.
#' @param n_factors Number of latent factors requested.  ARD priors will
#'   shrink inactive factors toward zero; the actual number retained may be
#'   lower.  Default \code{5}.
#' @param scale_views Logical.  Rescale views to equal total variance before
#'   fitting (MOFA's cross-view balancing, applied on top of the per-view
#'   scaling done in \code{\link{seMOFAFeatures}}).  Default \code{TRUE}.
#' @param seed Integer random seed.  Default \code{42}.
#' @param use_basilisk Logical.  Use \pkg{basilisk} for Python management.
#'   Default \code{FALSE}.
#' @param condaenv Name of a conda environment with \code{mofapy2} installed.
#'   When provided, \code{use_basilisk} is forced to \code{FALSE}.
#' @param ... Additional arguments passed to \code{MOFA2::run_mofa()}.
#'
#' @return An S3 object of class \code{mofa_fit_se} containing the fitted
#'   MOFA2 object, variance explained, factor scores, and all parameters
#'   needed by \code{\link{predictMOFASE}}, \code{\link{mofaWeights}}, and
#'   \code{\link{mofaVariance}}.
#'
#' @seealso \code{\link{seMOFAFeatures}}, \code{\link{predictMOFASE}},
#'   \code{\link{mofaWeights}}, \code{\link{mofaVariance}}
#'
#' @export
fitMOFASE <- function(se,
                       mofa_features,
                       n_factors    = 5,
                       scale_views  = TRUE,
                       seed         = 42,
                       use_basilisk = FALSE,
                       condaenv     = NULL,
                       ...) {
  .mofa_check()
  stopifnot(is.list(mofa_features), !is.null(mofa_features$views))

  if (!is.null(condaenv)) {
    use_basilisk <- FALSE
    reticulate::use_condaenv(condaenv, required = TRUE)
  }

  mofa_obj <- MOFA2::create_mofa(mofa_features$views)

  data_opts             <- MOFA2::get_default_data_options(mofa_obj)
  data_opts$scale_views <- scale_views

  model_opts             <- MOFA2::get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors

  train_opts         <- MOFA2::get_default_training_options(mofa_obj)
  train_opts$seed    <- seed
  train_opts$verbose <- FALSE

  mofa_obj <- MOFA2::prepare_mofa(mofa_obj,
                                   data_options     = data_opts,
                                   model_options    = model_opts,
                                   training_options = train_opts)
  mofa_obj <- MOFA2::run_mofa(mofa_obj, use_basilisk = use_basilisk, ...)

  var_exp    <- tryCatch(MOFA2::get_variance_explained(mofa_obj),
                          error = function(e) NULL)
  factor_mat <- do.call(rbind, MOFA2::get_factors(mofa_obj))   # samples × factors

  structure(
    list(
      mofa_obj     = mofa_obj,
      views        = names(mofa_features$views),
      n_factors    = ncol(factor_mat),
      scale_params = mofa_features$scale_params,
      well_names   = mofa_features$well_names,
      ephys_cols   = mofa_features$ephys_cols,
      morpho_cols  = mofa_features$morpho_cols,
      assay_names  = mofa_features$assay_names,
      scale_views  = scale_views,
      var_exp      = var_exp,
      factor_mat   = factor_mat,
      seed         = seed
    ),
    class = "mofa_fit_se"
  )
}

#' @export
print.mofa_fit_se <- function(x, ...) {
  cat("MOFA fit (ephacRTools)\n")
  cat("  Views:      ", paste(x$views, collapse = ", "), "\n")
  cat("  Factors:    ", x$n_factors, "\n")
  cat("  Wells:      ", length(x$well_names), "\n")
  cat("  scale_views:", x$scale_views, "\n")
  if (!is.null(x$var_exp)) {
    r2 <- x$var_exp$r2_per_factor
    if (!is.null(r2)) {
      cat("  Variance explained (R²) per view:\n")
      for (v in names(r2)) {
        vals <- colMeans(as.matrix(r2[[v]]), na.rm = TRUE)
        cat(sprintf("    %-12s: %s\n", v,
                    paste(sprintf("F%d=%.1f%%", seq_along(vals), vals),
                          collapse = "  ")))
      }
    }
  }
  invisible(x)
}


# ── predictMOFASE ─────────────────────────────────────────────────────────────

#' Write MOFA2 factor scores to a SummarizedExperiment
#'
#' Extracts factor scores from a fitted \code{mofa_fit_se} object and writes
#' them to \code{colData(se)} as \code{\{prefix\}_Factor1},
#' \code{\{prefix\}_Factor2}, \ldots — the same \code{out_prefix} convention
#' used by \code{\link{predictLDASE}} for LD scores.  Factor scores are also
#' stored in \code{reducedDims(se)[["MOFA"]]}.
#'
#' @param se A \code{SummarizedExperiment} compatible with the SE used for
#'   fitting (same well names).
#' @param mofa_fit Output of \code{\link{fitMOFASE}}.
#' @param out_prefix Prefix for new colData columns.  Default \code{"mofa"}.
#'
#' @return The \code{se} with factor score columns appended to
#'   \code{colData} and factor matrix stored in
#'   \code{reducedDims(se)[["MOFA"]]}.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom S4Vectors DataFrame
#' @export
predictMOFASE <- function(se, mofa_fit, out_prefix = "mofa") {
  stopifnot(inherits(mofa_fit, "mofa_fit_se"))

  factor_mat <- mofa_fit$factor_mat   # samples × factors (rownames = well names)
  n_factors  <- mofa_fit$n_factors
  se_wells   <- colnames(se)
  factor_nms <- paste0(out_prefix, "_Factor", seq_len(n_factors))

  # Align factor scores to SE column order
  aligned <- matrix(NA_real_, nrow = length(se_wells), ncol = n_factors,
                    dimnames = list(se_wells, factor_nms))
  idx   <- match(rownames(factor_mat), se_wells)
  valid <- !is.na(idx)
  aligned[idx[valid], ] <- factor_mat[valid, ]

  # Write to colData
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  for (k in seq_len(n_factors))
    cd[[factor_nms[k]]] <- aligned[, k]
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)

  # Write to reducedDims
  SingleCellExperiment::reducedDims(se)[["MOFA"]] <- aligned

  se
}


# ── mofaWeights ───────────────────────────────────────────────────────────────

#' Feature weights from a fitted MOFA2 model
#'
#' Returns the per-feature weights for each view and factor, analogous to
#' \code{\link{ldaLoadings}} for LDA models.  Importance is defined as the
#' mean absolute weight across selected factors.
#'
#' @param mofa_fit Output of \code{\link{fitMOFASE}}.
#' @param view Character vector.  One or more view names to include.
#'   \code{NULL} (default) returns all views.
#' @param factor Integer vector.  Which factors to include.
#'   \code{NULL} (default) returns all factors.
#' @param n_top Return only the top \code{n_top} features per view by
#'   importance.  \code{NULL} returns all.
#'
#' @return A \code{data.frame} with columns \code{feature}, \code{view},
#'   one column per selected factor (\code{Factor1}, \code{Factor2}, \ldots),
#'   \code{importance} (mean |weight| across factors), and \code{rank}.
#'
#' @export
mofaWeights <- function(mofa_fit, view = NULL, factor = NULL, n_top = NULL) {
  .mofa_check()
  stopifnot(inherits(mofa_fit, "mofa_fit_se"))

  w_list    <- MOFA2::get_weights(mofa_fit$mofa_obj)
  views_use <- if (!is.null(view)) intersect(view, names(w_list)) else names(w_list)
  if (!length(views_use))
    stop("No matching views. Available: ", paste(names(w_list), collapse = ", "))

  rows <- lapply(views_use, function(v) {
    mat <- w_list[[v]]   # features × factors
    if (!is.null(factor)) {
      keep <- intersect(paste0("Factor", factor), colnames(mat))
      mat  <- mat[, keep, drop = FALSE]
    }
    factor_cols <- colnames(mat)
    df              <- as.data.frame(mat)
    df$feature      <- rownames(mat)
    df$view         <- v
    df$importance   <- rowMeans(abs(mat), na.rm = TRUE)
    df              <- df[order(-df$importance), , drop = FALSE]
    df$rank         <- seq_len(nrow(df))
    if (!is.null(n_top)) df <- head(df, n_top)
    lead <- c("feature", "view", factor_cols, "importance", "rank")
    df[, intersect(lead, colnames(df)), drop = FALSE]
  })

  do.call(rbind, rows)
}


# ── mofaVariance ──────────────────────────────────────────────────────────────

#' Variance explained by MOFA2 factors per view
#'
#' Returns the \eqn{R^2} variance explained by each factor in each view,
#' analogous to \code{\link{ldaPCAVariance}} for LDA models.
#'
#' @param mofa_fit Output of \code{\link{fitMOFASE}}.
#'
#' @return A \code{data.frame} with columns \code{view}, \code{factor}, and
#'   \code{r2} (fraction of variance explained, averaged across groups when
#'   the model was fit with multiple sample groups).
#'
#' @export
mofaVariance <- function(mofa_fit) {
  .mofa_check()
  stopifnot(inherits(mofa_fit, "mofa_fit_se"))

  var_exp <- mofa_fit$var_exp
  if (is.null(var_exp))
    var_exp <- tryCatch(
      MOFA2::get_variance_explained(mofa_fit$mofa_obj),
      error = function(e)
        stop("mofaVariance: could not retrieve variance explained. ",
             "Did the model converge?")
    )

  r2_list <- var_exp$r2_per_factor
  rows <- lapply(names(r2_list), function(v) {
    mat <- as.matrix(r2_list[[v]])   # groups × factors
    r2v <- colMeans(mat, na.rm = TRUE)
    data.frame(view   = v,
               factor = paste0("Factor", seq_along(r2v)),
               r2     = as.numeric(r2v),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}


# ── clusterMOFA (deprecated wrapper) ─────────────────────────────────────────

#' Multi-modal factor analysis using MOFA2 (deprecated wrapper)
#'
#' Calls the modular \code{\link{seMOFAFeatures}} →
#' \code{\link{fitMOFASE}} → \code{\link{predictMOFASE}} pipeline.
#' For new code, use those functions directly for full control over scaling,
#' NA handling, and output prefix.
#'
#' @inheritParams seMOFAFeatures
#' @inheritParams fitMOFASE
#' @param out_prefix Prefix for colData factor columns.  Default \code{"mofa"}.
#'
#' @return \code{se} updated by \code{\link{predictMOFASE}}.
#'
#' @export
clusterMOFA <- function(se,
                         ephys_cols   = NULL,
                         morpho_cols  = NULL,
                         assay_names  = NULL,
                         n_factors    = 5,
                         scale_views  = TRUE,
                         seed         = 42,
                         use_basilisk = FALSE,
                         condaenv     = NULL,
                         out_prefix   = "mofa",
                         ...) {
  .Deprecated("fitMOFASE",
              msg = paste("clusterMOFA() is a compatibility wrapper.",
                          "Use seMOFAFeatures() + fitMOFASE() + predictMOFASE()",
                          "for the full pipeline."))
  mofa_feats <- seMOFAFeatures(se,
                                ephys_cols  = ephys_cols,
                                morpho_cols = morpho_cols,
                                assay_names = assay_names)
  mofa_fit   <- fitMOFASE(se, mofa_feats,
                           n_factors    = n_factors,
                           scale_views  = scale_views,
                           seed         = seed,
                           use_basilisk = use_basilisk,
                           condaenv     = condaenv,
                           ...)
  predictMOFASE(se, mofa_fit, out_prefix = out_prefix)
}


# ============================================================================
# clusterPipeline  (all-in-one wrapper)
# ============================================================================

#' All-in-one unsupervised multi-modal clustering pipeline
#'
#' Convenience wrapper that runs the full clustering pipeline in sequence:
#' \enumerate{
#'   \item \code{\link{prepareClusterFeatures}} — build weighted feature matrix
#'   \item \code{\link{optimalClusters}} — evaluate k (skipped if \code{k}
#'     is supplied)
#'   \item \code{\link{reducedDimMultimodal}} — PCA + UMAP + t-SNE
#'   \item \code{\link{clusterSE}} — assign cluster labels via
#'     \code{cluster_method}
#'   \item \code{\link{clusterHeatmap}} — feature profile heatmap
#'   \item \code{\link{clusterSummary}} — per-cluster boxplots + table
#' }
#'
#' @param se A \code{SingleCellExperiment}.
#' @param ephys_cols colData columns for ephys features.
#' @param morpho_cols colData columns for morphology features.
#' @param assay_names Assay names to include as full IV-curve feature blocks.
#' @param ephys_weight Relative variance contribution of the ephys block.
#'   Default \code{1}.
#' @param morpho_weight Relative variance contribution of the morphology block.
#'   Default \code{1}.
#' @param scale_method Feature scaling.  Default \code{"zscore"}.
#' @param na_action Default NA handling for both modality blocks.  Default
#'   \code{"impute_median"}.  See \code{\link{prepareClusterFeatures}} for
#'   details.
#' @param ephys_na_action Override NA strategy for the ephys block.
#'   \code{NULL} (default) falls back to \code{na_action}.
#' @param morpho_na_action Override NA strategy for the morphology block.
#'   \code{NULL} (default) falls back to \code{na_action}.  Typical choice:
#'   \code{"impute_zero"} when missing particles mean zero signal.
#' @param k Integer number of clusters.  If \code{NULL} (default),
#'   \code{optimalClusters()} is run and the silhouette-suggested k is used.
#' @param k_range Integer vector for \code{optimalClusters()}.  Default
#'   \code{2:8}.
#' @param opt_methods Quality metrics passed to \code{optimalClusters()}.
#'   Default \code{c("silhouette", "wss")}.
#' @param cluster_method Clustering algorithm.  Default \code{"kmeans"}.
#' @param dim_methods Reductions to compute.  Default
#'   \code{c("pca", "tsne", "umap")}.
#' @param summary_cols colData columns for \code{clusterSummary()}.  If
#'   \code{NULL}, uses all of \code{ephys_cols} and \code{morpho_cols}.
#' @param seed Random seed.  Default \code{42}.
#' @param verbose Print progress messages.  Default \code{TRUE}.
#'
#' @return The input \code{SingleCellExperiment} updated with:
#' \describe{
#'   \item{colData(se)$cluster}{Factor of cluster labels.}
#'   \item{reducedDims(se)}{PCA / UMAP / t-SNE matrices.}
#'   \item{metadata(se)[["PCA_model"]]}{Fitted \code{prcomp} object for
#'     \code{\link{plotFeatureImportance}}.}
#'   \item{metadata(se)[["clustering"]]}{Named list with:
#'     \code{$feature_mat} (scaled matrix),
#'     \code{$optimal_k} (output of \code{optimalClusters()}, or \code{NULL}),
#'     \code{$heatmap} (a \code{ComplexHeatmap::Heatmap} — print to render),
#'     \code{$summary} (list with \code{$plot} and \code{$table}).}
#' }
#'
#' @examples
#' \dontrun{
#' data(se_iN)
#' se_iN <- get_metric(se_iN, assay = "Minima")
#' se_iN <- colAG(se_iN, c("Capacitance", "Seal_Resistance"))
#'
#' se_iN <- clusterPipeline(
#'   se_iN,
#'   ephys_cols       = c("Imax.minima", "Vhalf.minima", "Capacitance_mean"),
#'   morpho_cols      = c("GFP_Mean_z", "GFP_normArea"),
#'   ephys_weight     = 1,
#'   morpho_weight    = 1,
#'   morpho_na_action = "impute_zero",
#'   k_range          = 2:6
#' )
#'
#' # All results live inside the SE
#' metadata(se_iN)$clustering$optimal_k$plot   # elbow / silhouette plot
#' print(metadata(se_iN)$clustering$heatmap)   # feature profile heatmap
#' metadata(se_iN)$clustering$summary$plot     # per-cluster boxplots
#' metadata(se_iN)$clustering$summary$table    # per-cluster mean / SD table
#' plotDimRed(se_iN, "UMAP", colorColumns = c("cluster", "Imax.minima"))
#' }
#'
#' @importFrom S4Vectors metadata
#' @export
clusterPipeline <- function(se,
                             ephys_cols       = NULL,
                             morpho_cols      = NULL,
                             assay_names      = NULL,
                             ephys_weight     = 1,
                             morpho_weight    = 1,
                             scale_method     = c("zscore", "minmax", "none"),
                             na_action        = c("impute_median",
                                                  "impute_zero", "omit"),
                             ephys_na_action  = NULL,
                             morpho_na_action = NULL,
                             k                = NULL,
                             k_range          = 2:8,
                             opt_methods      = c("silhouette", "wss"),
                             cluster_method   = c("kmeans", "hierarchical",
                                                  "louvain", "gmm"),
                             dim_methods      = c("pca", "tsne", "umap"),
                             summary_cols     = NULL,
                             seed             = 42,
                             verbose          = TRUE) {
  scale_method   <- match.arg(scale_method)
  na_action      <- match.arg(na_action)
  cluster_method <- match.arg(cluster_method)

  .msg <- function(...) if (verbose) message("[clusterPipeline] ", ...)

  # Step 1 — feature matrix
  .msg("(1/5) Preparing feature matrix...")
  mat <- prepareClusterFeatures(
    se,
    ephys_cols       = ephys_cols,
    morpho_cols      = morpho_cols,
    assay_names      = assay_names,
    ephys_weight     = ephys_weight,
    morpho_weight    = morpho_weight,
    scale_method     = scale_method,
    na_action        = na_action,
    ephys_na_action  = ephys_na_action,
    morpho_na_action = morpho_na_action
  )
  .msg(sprintf("       %d wells × %d features  (ephys_weight=%.1f, morpho_weight=%.1f)",
               nrow(mat), ncol(mat), ephys_weight, morpho_weight))

  # Step 2 — optimal k
  opt_k_result <- NULL
  if (is.null(k)) {
    .msg("(2/5) Evaluating optimal k (k_range = ", min(k_range), "–",
         max(k_range), ")...")
    opt_k_result <- optimalClusters(mat, k_range = k_range,
                                     methods = opt_methods, seed = seed)
    k_chosen <- opt_k_result$suggested_k["silhouette"]
    if (is.na(k_chosen) || length(k_chosen) == 0)
      k_chosen <- opt_k_result$suggested_k[1]
    k <- as.integer(k_chosen)
    .msg(sprintf("       → using k = %d (silhouette suggestion)", k))
  } else {
    .msg("(2/5) Using supplied k = ", k, " (skipping optimalClusters)")
  }

  # Step 3 — dimensionality reduction + clustering
  .msg("(3/5) Running reducedDimMultimodal (method = c(",
       paste(dim_methods, collapse = ", "), "))...")
  se <- reducedDimMultimodal(
    se,
    ephys_cols       = ephys_cols,
    morpho_cols      = morpho_cols,
    assay_names      = assay_names,
    ephys_weight     = ephys_weight,
    morpho_weight    = morpho_weight,
    scale_method     = scale_method,
    na_action        = na_action,
    ephys_na_action  = ephys_na_action,
    morpho_na_action = morpho_na_action,
    method           = dim_methods,
    k_clusters       = k,
    cluster_method   = cluster_method,
    seed             = seed
  )
  .msg("       Cluster distribution: ",
       paste(names(table(se$cluster)), collapse = " | "),
       "  (n = ",
       paste(as.integer(table(se$cluster)), collapse = " / "), ")")

  # Step 4 — heatmap
  .msg("(4/5) Building cluster feature heatmap...")
  ht <- clusterHeatmap(se, feature_mat = mat, cluster_col = "cluster")

  # Step 5 — summary
  .msg("(5/5) Computing per-cluster summary...")
  s_cols <- if (!is.null(summary_cols)) summary_cols else
    c(ephys_cols, morpho_cols)
  s_cols <- intersect(s_cols,
                       names(as.data.frame(SummarizedExperiment::colData(se))))
  smry <- if (length(s_cols) > 0)
    clusterSummary(se, feature_cols = s_cols, cluster_col = "cluster")
  else NULL

  # Store all pipeline outputs inside the SE — nothing falls outside the object
  S4Vectors::metadata(se)[["clustering"]] <- list(
    feature_mat = mat,
    optimal_k   = opt_k_result,
    heatmap     = ht,
    summary     = smry
  )

  .msg("Done. Results stored in metadata(se)$clustering")
  .msg("  Access: metadata(se)$clustering$heatmap  |  $optimal_k  |  $summary  |  $feature_mat")
  se
}


# ============================================================================
# plotFeatureImportance
# ============================================================================

#' Plot feature importance: PCA loadings and cluster discriminability
#'
#' Two complementary views of which features matter:
#' \describe{
#'   \item{PCA loadings}{Which features drive the principal components —
#'     i.e. which features explain most of the \emph{global} variance in the
#'     data, regardless of the clustering result.  A high absolute loading on
#'     PC1 means the feature strongly separates wells along the main axis of
#'     variation.}
#'   \item{Cluster discriminability (η²)}{Which features best separate the
#'     identified clusters — i.e. the fraction of each feature's total variance
#'     that falls \emph{between} clusters rather than within them.  η² = 1
#'     means the feature perfectly predicts cluster membership; η² = 0 means
#'     it is identical across clusters.  This is computed as the ratio of
#'     between-cluster SS to total SS (one-way ANOVA decomposition, no
#'     parametric assumptions assumed on interpretation).}
#' }
#'
#' The two views often diverge: a feature can dominate PC1 but show low η²
#' (high global variance but not cluster-specific), or have low loadings but
#' high η² (a quiet feature that perfectly sorts cells into groups).  The
#' combined table lets you identify features that are informative in both senses.
#'
#' @param se A \code{SingleCellExperiment} with a PCA model stored in
#'   \code{S4Vectors::metadata(se)[["PCA_model"]]} (set automatically by
#'   \code{\link{reducedDimMultimodal}} or \code{\link{clusterPipeline}}).
#' @param feature_mat Scaled feature matrix (wells × features) from
#'   \code{\link{prepareClusterFeatures}}.  Needed for η² computation.
#' @param cluster_col Name of the colData column with cluster labels.
#'   Default \code{"cluster"}.
#' @param pcs Integer vector of PCs to show loadings for.  Default \code{1:3}.
#' @param n_top Maximum number of features shown per PC loading plot
#'   (those with the largest absolute loading).  Default \code{15}.
#'
#' @return A named list:
#'   \describe{
#'     \item{loadings_plot}{Bar chart of signed PCA loadings, faceted by PC.
#'       Blue = ephys features, orange = morphology features.}
#'     \item{discriminability_plot}{Horizontal bar chart of η² per feature,
#'       sorted descending.  Features at the top discriminate clusters best.}
#'     \item{combined_plot}{\code{ggpubr::ggarrange} of both plots side by
#'       side.}
#'     \item{table}{data.frame with one row per feature containing loadings
#'       for each requested PC and the η² value.}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- clusterPipeline(se, ephys_cols = c("Imax.minima", "Vhalf.minima"),
#'                            morpho_cols = c("GFP_Mean_z", "GFP_normArea"))
#' imp <- plotFeatureImportance(result$se, result$feature_mat)
#' imp$combined_plot
#' imp$table
#' }
#'
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @export
plotFeatureImportance <- function(se,
                                   feature_mat,
                                   cluster_col = "cluster",
                                   pcs         = 1:3,
                                   n_top       = 15) {

  # ---- PCA model -----------------------------------------------------------
  pca_model <- S4Vectors::metadata(se)[["PCA_model"]]
  if (is.null(pca_model))
    stop("plotFeatureImportance: PCA model not found.\n",
         "  Run reducedDimMultimodal() or clusterPipeline() first.")

  rotation <- pca_model$rotation   # features × PCs
  pvar     <- pca_model$sdev^2 / sum(pca_model$sdev^2) * 100

  pcs      <- intersect(pcs, seq_len(ncol(rotation)))
  pc_names <- colnames(rotation)[pcs]

  # Modality annotation — prefer stored version, fall back to feature_mat attr
  mod_stored <- S4Vectors::metadata(se)[["feature_modality"]]
  mod_attr   <- if (!is.null(mod_stored)) mod_stored else attr(feature_mat, "modality")
  mod_colors <- c(ephys = "#4E79A7", morpho = "#F28E2B", unknown = "#BAB0AC")

  # ---- Build long loading data frame ---------------------------------------
  load_list <- lapply(pcs, function(pc_i) {
    pc_nm  <- pc_names[pc_i]
    loads  <- rotation[, pc_i]
    # Keep only top n_top by absolute value
    top_idx <- order(abs(loads), decreasing = TRUE)[seq_len(min(n_top, length(loads)))]
    data.frame(
      feature  = names(loads)[top_idx],
      loading  = loads[top_idx],
      pc       = pc_nm,
      pct_var  = round(pvar[pc_i], 1),
      stringsAsFactors = FALSE
    )
  })
  load_df <- do.call(rbind, load_list)

  load_df$display  <- gsub("^[em]__", "", load_df$feature)
  load_df$modality <- mod_attr[load_df$feature]
  load_df$modality[is.na(load_df$modality)] <- "unknown"

  # pc label with % variance
  load_df$pc_label <- paste0(load_df$pc, "\n(",
                              load_df$pct_var, "% var)")

  # Within each PC, sort features by signed loading so plot reads cleanly
  load_df$display <- factor(
    load_df$display,
    levels = unique(load_df$display[order(load_df$loading)])
  )

  load_plot <- ggplot2::ggplot(load_df,
    ggplot2::aes(x = .data$display,
                 y = .data$loading,
                 fill = .data$modality)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~pc_label, scales = "free", nrow = 1) +
    ggplot2::scale_fill_manual(values = mod_colors, name = "Modality") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(x = NULL, y = "PCA loading",
                  title = paste0("PCA loadings  (top ", n_top,
                                 " features per PC)")) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold")
    )

  # ---- Cluster discriminability (η²) per feature ---------------------------
  cd          <- as.data.frame(SummarizedExperiment::colData(se))
  cluster_vec <- cd[[cluster_col]]

  # Align feature_mat rows to SE column order
  mat_wells <- rownames(feature_mat)
  se_wells  <- colnames(se)
  if (!is.null(mat_wells)) {
    idx     <- match(se_wells, mat_wells)
    valid   <- !is.na(idx)
    mat_ord <- feature_mat[idx[valid], , drop = FALSE]
    cl_ord  <- cluster_vec[valid]
  } else {
    mat_ord <- feature_mat
    cl_ord  <- cluster_vec
  }
  cl_ord <- as.character(cl_ord)

  eta2_vals <- vapply(seq_len(ncol(mat_ord)), function(j) {
    y   <- as.numeric(mat_ord[, j])
    grp <- cl_ord
    ok  <- !is.na(y) & !is.na(grp)
    y   <- y[ok];  grp <- grp[ok]
    if (length(unique(grp)) < 2L || length(y) == 0L) return(NA_real_)
    grand_mean <- mean(y)
    ss_total   <- sum((y - grand_mean)^2)
    if (ss_total == 0) return(0)
    grp_means  <- tapply(y, grp, mean)
    grp_n      <- as.integer(table(grp))
    ss_between <- sum(grp_n * (grp_means - grand_mean)^2)
    ss_between / ss_total
  }, numeric(1))

  disc_df <- data.frame(
    feature  = colnames(mat_ord),
    display  = gsub("^[em]__", "", colnames(mat_ord)),
    eta2     = eta2_vals,
    modality = mod_attr[colnames(mat_ord)],
    stringsAsFactors = FALSE
  )
  disc_df$modality[is.na(disc_df$modality)] <- "unknown"
  disc_df <- disc_df[!is.na(disc_df$eta2), ]
  disc_df <- disc_df[order(disc_df$eta2, decreasing = TRUE), ]
  disc_df$display <- factor(disc_df$display, levels = rev(disc_df$display))

  disc_plot <- ggplot2::ggplot(disc_df,
    ggplot2::aes(x = .data$display,
                 y = .data$eta2 * 100,   # show as percentage
                 fill = .data$modality)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = mod_colors, name = "Modality") +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      x        = NULL,
      y        = expression(eta^2 ~ "(% between-cluster variance)"),
      title    = "Cluster discriminability per feature",
      subtitle = paste0("Cluster: '", cluster_col, "'  |  ",
                        length(unique(stats::na.omit(cl_ord))), " clusters  |  ",
                        sum(!is.na(cl_ord)), " wells")
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold")
    )

  # ---- Combined side-by-side -----------------------------------------------
  combined <- ggpubr::ggarrange(load_plot, disc_plot,
                                 common.legend = TRUE,
                                 legend        = "bottom",
                                 widths        = c(1.5, 1))

  # ---- Summary table -------------------------------------------------------
  # Wide loadings
  load_wide <- stats::reshape(
    load_df[, c("feature", "pc", "loading")],
    idvar     = "feature",
    timevar   = "pc",
    direction = "wide"
  )
  names(load_wide) <- gsub("^loading\\.", "", names(load_wide))

  tbl <- merge(
    disc_df[, c("feature", "display", "modality", "eta2")],
    load_wide,
    by = "feature", all.x = TRUE
  )
  tbl$eta2_pct <- round(tbl$eta2 * 100, 1)
  tbl$eta2     <- round(tbl$eta2,       4)
  tbl          <- tbl[order(tbl$eta2, decreasing = TRUE), ]
  rownames(tbl) <- NULL

  list(
    loadings_plot          = load_plot,
    discriminability_plot  = disc_plot,
    combined_plot          = combined,
    table                  = tbl
  )
}
