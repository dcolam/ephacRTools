#' @importFrom magrittr %>%
NULL
#' @import SingleCellExperiment
NULL

#' Subset a SummarizedExperiment by a colData expression
#'
#' A convenience wrapper around \code{se[, keep]} that works like
#' \code{base::subset()} — you write a plain expression using colData column
#' names directly, and \code{NA}s in the result are automatically treated as
#' \code{FALSE} so they never cause a subscript error.
#'
#' @param se A \code{SummarizedExperiment} (or \code{SingleCellExperiment}).
#' @param subset A logical expression evaluated in the context of
#'   \code{colData(se)}.  Column names can be used directly without the
#'   \code{se$} prefix.
#'
#' @return A subsetted \code{SummarizedExperiment} containing only the columns
#'   (wells) for which \code{subset} is \code{TRUE}.
#'
#' @examples
#' \dontrun{
#' se <- subsetSE(se, Condition != "Empty" & Induction != "empty")
#' se <- subsetSE(se, plate_id %in% c("P1", "P2") & qc == "Pass")
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @export
subsetSE <- function(se, subset) {
  cd   <- as.data.frame(SummarizedExperiment::colData(se))
  expr <- substitute(subset)
  keep <- eval(expr, cd, parent.frame())
  keep[is.na(keep)] <- FALSE
  se[, keep]
}
#' Add column-wise aggregation such as mean of any given assay and store it into colData
#' @param assayName list of assay names to check
#' @param assayList list of assays in the SE
#' @return updated assayList
#' @export
checkAssay <- function(assayList, assayList.se){
  ## check whether assayList contains assays
  for(assayName in assayList){
    if(!(assayName %in% assayList.se)){
      print(paste(assayName, "is not an assay, skipped"))
      assayList <- assayList[!(assayName==assayList)]
    }
  }
  return(assayList)
}

#' Add column-wise aggregation of assays and store results into colData
#'
#' Applies a summary function column-wise (i.e. across sweeps) for each
#' requested assay and stores the result as a new \code{colData} column named
#' \code{<Assay>_<fun>[suffix]}, e.g. \code{Minima_mean} or
#' \code{Minima_mean.exp1}.
#'
#' @param se A \code{SummarizedExperiment} (or \code{SingleCellExperiment}).
#' @param assayList Character vector of assay names to aggregate.
#' @param fun Function to apply column-wise.  Must accept a numeric vector and
#'   return a single value (e.g. \code{mean}, \code{median}, \code{sd},
#'   \code{max}).  Default \code{mean}.  The function name is inferred
#'   automatically and used in the column name.
#' @param sweeps Row names (sweeps) to include before aggregating.  Default all
#'   rows.
#' @param suffix Optional string appended verbatim after the function name in
#'   the column name, e.g. \code{".exp1"} produces \code{Minima_mean.exp1}.
#'   Default \code{""} (no suffix).
#' @param ... Additional arguments passed to \code{fun} (e.g. \code{trim = 0.1}
#'   for \code{mean}).  Note: \code{NA} and non-finite values are already
#'   removed before \code{fun} is called, so \code{na.rm = TRUE} is redundant
#'   but harmless.
#' @return The input \code{se} with new \code{colData} columns added.
#' @examples
#' \dontrun{
#' se <- colAG(se, c("Minima", "Capacitance"))                        # -> Minima_mean, Capacitance_mean
#' se <- colAG(se, "Minima", suffix = ".exp1")                        # -> Minima_mean.exp1
#' se <- colAG(se, "Minima", fun = sd,     suffix = ".exp1")          # -> Minima_sd.exp1
#' se <- colAG(se, "Minima", fun = median, suffix = ".baseline")      # -> Minima_median.baseline
#' se <- colAG(se, "Minima", fun = mean,   trim = 0.1)                # -> Minima_mean (trimmed)
#' }
#' @export
colAG <- function(se, assayList, fun = mean, sweeps = row.names(se), suffix = "", ...) {
  fun_name  <- as.character(substitute(fun))
  assayList <- assayList[assayList %in% assayNames(se)]

  for (assayName in assayList) {
    colName <- paste0(assayName, "_", fun_name, suffix)
    subse   <- se[sweeps, ]
    mat     <- assay(subse, assayName)
    se[[colName]] <- apply(mat, 2, function(x) fun(x[is.finite(x)], ...))
  }
  return(se)
}
#' Perform cell-wise dimensionality reduction based on assays and colData
#'
#' Computes PCA, t-SNE, and UMAP on a feature matrix built from assay matrices
#' and/or flat \code{colData} columns, stores the results in
#' \code{reducedDims(se)}, and adds k-means cluster assignments as new
#' \code{colData} columns.
#'
#' @param se A \code{SummarizedExperiment} (or \code{SingleCellExperiment}).
#' @param assayList Character vector of assay names to include as features.
#'   Each assay is transposed (sweeps become features) before scaling.
#' @param colNames Character vector of flat \code{colData} column names to
#'   include as additional features.
#' @param scaling Scaling strategy: \code{"within"} (default) scales each
#'   assay/column independently using \code{sechm::safescale};
#'   \code{"global"} scales the combined feature matrix; \code{"none"} skips
#'   scaling.
#' @param byRow Logical. Scale by row instead of by column after transposing.
#'   Default \code{FALSE}.
#' @param center Logical. Center the data during scaling. Default \code{FALSE}.
#' @param method Character vector listing which reductions to compute.
#'   Any subset of \code{c("pca", "tsne", "umap")}. All three are computed by
#'   default.
#' @param k_clusters Integer. Number of k-means clusters to assign for each
#'   reduction. Default \code{3}.
#' @return The input \code{se} with \code{reducedDims} populated (\code{PCA},
#'   \code{TSNE}, \code{UMAP}) and three new \code{colData} columns:
#'   \code{cluster.pca}, \code{cluster.tsne}, \code{cluster.umap}.
#' @export
reducedDim.Cellwise <- function(se, assayList = c(), colNames = c(),
                                scaling = "within", byRow = FALSE,
                                center = FALSE,
                                method = c("pca", "tsne", "umap"),
                                k_clusters = 3) {

  if (length(assayList) != 0) {
    assayList <- assayList[assayList %in% assayNames(se)]
    pca_data <- lapply(assayList, function(x) {
      if (scaling == "within") {
        temp <- sechm::safescale(t(assay(se, x)), center = center, byRow = byRow)
      } else {
        temp <- t(assay(se, x))
      }
      temp <- as.data.frame(temp)
      colnames(temp) <- paste(x, colnames(temp))
      temp
    })
    names(pca_data) <- assayList
  } else {
    pca_data <- list()
  }

  if (length(colNames) != 0) {
    flattened.df <- as.data.frame(colData(se))
    colNames <- colNames[colNames %in% colnames(flattened.df)]
    col_Data <- lapply(colNames, function(x) {
      if (scaling == "within") {
        temp <- sechm::safescale(flattened.df[[x]], center = FALSE, byRow = byRow)
      } else {
        temp <- flattened.df[[x]]
      }
      temp <- as.data.frame(temp)
      colnames(temp) <- paste(x, colnames(temp))
      temp
    })
    names(col_Data) <- colNames
    pca_data <- list(pca_data, col_Data)
  }

  pca_data <- dplyr::bind_cols(pca_data)

  if (scaling == "global") {
    pca_data <- sechm::safescale(as.matrix(pca_data), byRow = byRow)
  }

  pca_data <- as.data.frame(pca_data)
  pca_data <- pca_data[, colSums(is.na(pca_data)) != nrow(pca_data)]
  pca_data[is.na(pca_data)] <- 0

  pca_result <- prcomp(pca_data, rank = 50)
  tsne_data  <- Rtsne::Rtsne(pca_data, pca = TRUE, check_duplicates = FALSE)
  tsne_data  <- tsne_data$Y %>%
    as.data.frame() %>%
    dplyr::rename(tsne1 = "V1", tsne2 = "V2")

  umap_data <- umap::umap(pca_data)
  umap_df   <- umap_data$layout %>%
    as.data.frame() %>%
    dplyr::rename(UMAP1 = "V1", UMAP2 = "V2")

  SingleCellExperiment::reducedDims(se) <- list(
    PCA  = pca_result$x,
    TSNE = S4Vectors::DataFrame(tsne_data),
    UMAP = S4Vectors::DataFrame(umap_df)
  )

  se$cluster.umap <- as.factor(kmeans(SingleCellExperiment::reducedDim(se, "UMAP")[, 1:2], k_clusters, iter.max = 100)$cluster)
  se$cluster.tsne <- as.factor(kmeans(SingleCellExperiment::reducedDim(se, "TSNE")[, 1:2], k_clusters, iter.max = 100)$cluster)
  se$cluster.pca  <- as.factor(kmeans(SingleCellExperiment::reducedDim(se, "PCA")[,  1:2], k_clusters, iter.max = 100)$cluster)

  return(se)
}

# reducedDim.Cellwise() above is the original simpler pipeline.
# A more feature-rich replacement (reducedDimMultimodal) lives in R/clusteringTools.R.

#' Wrapper function for plotting assays vs sweeps curves
#' @param se SummarizedExperiment Object with reducedDim data
#' @param assayList assays to be plotted as y-value
#' @param rowCol numeric row Column to plot as x-value
#' @param colorGroup color grouping
#' @param wrapFormula formula to be given to facet_wrap
#' @param grouped boolean to have grouped value (mean of the grouping) or single wells
#' @return a ggplot
#' @export
plotAssayVSSweeps <- function(se, assayList, rowCol, colorGroup=c(), wrapFormula=NULL, grouped=TRUE){
  assayList <- assayList[assayList %in% assayNames(se)]
  rowCol <- rowCol[rowCol %in% colnames(rowData(se))]
  melted.se <- sechm::meltSE(se, features=row.names(se), assayName=assayList, rowDat.columns = rowCol)

  melted.se <-reshape2::melt(melted.se, measure.vars = assayList)
  if(!grouped){
  ggplot2::ggplot(melted.se, aes(x=melted.se[[rowCol]], y=value, color=well_id)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::guides(color="none")
}else{

  p <- ggplot2::ggplot(melted.se, aes(x=.data[[rowCol]], y=value, color=.data[[colorGroup]])) +
    ggplot2::stat_summary(geom='errorbar',fun.data=mean_se, size=1, alpha=0.6) +
    ggplot2::stat_summary(geom='line', fun = "mean", size=1, alpha=1) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::ylab("Current (nA)") +
    ggplot2::xlab("Holding Potential (mV)") +
    ggplot2::geom_hline(yintercept=0, linetype="dashed")

  if (!is.null(wrapFormula)) {
    # Use facet_grid or facet_wrap depending on the formula
    p <- p + ggplot2::facet_wrap(wrapFormula, scales="free")
    # Or: p <- p + ggplot2::facet_wrap(group_formula)
  }
  p
  }
}
#' Extract IV-curve metrics from a step-wise voltage clamp assay
#'
#' Analyzes current-voltage (IV) curves from a step-wise voltage clamp recording
#' stored in a \code{SingleCellExperiment} object. For each well (identified by
#' \code{well_id} + \code{plate_id}), the function extracts three metrics and stores
#' them as new columns in \code{colData}:
#'
#' \describe{
#'   \item{\code{Imax.<assay>}}{Peak current amplitude -- minimum for inward
#'     currents (e.g. sodium, \code{inward = TRUE}) or maximum for outward
#'     currents (e.g. potassium, \code{inward = FALSE}).}
#'   \item{\code{Vmax.<assay>}}{Holding potential (mV) at which \code{Imax}
#'     occurs.}
#'   \item{\code{Vhalf.<assay>}}{Half-activation voltage (mV), estimated by
#'     fitting a smoothing spline to the activation limb of the IV curve
#'     (voltages \eqn{\leq} \code{Vmax}) and locating the voltage where current
#'     equals \code{Imax / 2}. Returns \code{NA} when spline fitting fails.}
#' }
#'
#' The assay matrix must have sweeps as rows and wells as columns. \code{v_clamp_mV}
#' must be present in \code{rowData(se)} and contain the holding potential for
#' each sweep, as produced by \code{prepareSE()}.
#'
#' @param se A \code{SingleCellExperiment} (or \code{SummarizedExperiment}) with
#'   step-wise voltage clamp data. Must contain \code{v_clamp_mV} in
#'   \code{rowData} and \code{well_id} / \code{plate_id} in \code{colData}.
#' @param assay Name of the assay to analyze. Default \code{"Minima"} -- the
#'   per-sweep minimum current output by DataControl Online Analysis, which
#'   captures peak inward current for step protocols.
#' @param inward Logical. \code{TRUE} (default) for inward currents: peak is
#'   the minimum value. \code{FALSE} for outward currents: peak is the maximum.
#'
#' @return The input \code{se} with three new \code{colData} columns:
#'   \code{Imax.<assay>}, \code{Vmax.<assay>}, and \code{Vhalf.<assay>}
#'   (all lowercase assay suffix, e.g. \code{Imax.minima}).
#'
#' @examples
#' \dontrun{
#' data(se_iN)
#' se_iN <- get_metric(se_iN, assay = "Minima", inward = TRUE)
#' head(colData(se_iN)[, c("Imax.minima", "Vmax.minima", "Vhalf.minima")])
#' }
#' @export
get_metric <- function(se, assay = "Minima", inward = TRUE) {

  suffix <- tolower(assay)
  imax_col <- paste0("Imax.", suffix)
  vhalf_col <- paste0("Vhalf.", suffix)
  vmax_col <- paste0("Vmax.", suffix)

  # Melt once into long format, pulling V_Clamp from rowData
  dat <- sechm::meltSE(
    se,
    rownames(se),
    assayName = assay,
    rowDat.columns = "v_clamp_mV"
  )

  # Split once by well_id + plate_id
  groups <- split(dat, interaction(dat$well_id, dat$plate_id, drop = TRUE))

  res <- lapply(groups, function(subdat) {

    x_vals <- subdat$v_clamp_mV
    y_vals <- subdat[[assay]]

    if (all(is.na(y_vals))) {
      return(c(Imax = NA_real_, Vhalf = NA_real_, Vmax = NA_real_))
    }

    if (inward) {
      Imax <- min(y_vals, na.rm = TRUE)
      Vmax1 <- min(x_vals[y_vals == Imax])
    } else {
      Imax <- max(y_vals, na.rm = TRUE)
      Vmax1 <- max(x_vals[y_vals == Imax])
    }

    y_vals[!is.finite(y_vals)] <- 1

    spl <- tryCatch(
      smooth.spline(x_vals, y_vals),
      error = function(e) NULL
    )

    if (is.null(spl)) {
      return(c(Imax = Imax, Vhalf = NA_real_, Vmax = Vmax1))
    }

    fit <- predict(spl, seq(min(x_vals), max(x_vals), length.out = 100))

    idx <- which(fit$x < Vmax1)
    Vhalf <- NA_real_

    if (length(idx)) {
      i2 <- which.min(abs(fit$y[idx] - Imax / 2))
      Vhalf <- fit$x[idx][i2]
    }

    c(Imax = Imax, Vhalf = Vhalf, Vmax = Vmax1)
  })

  res <- do.call(rbind, res)
  res <- as.data.frame(res)

  keys <- do.call(rbind, strsplit(names(groups), "\\."))
  colnames(keys) <- c("well_id", "plate_id")

  res$well_id <- keys[, 1]
  res$plate_id <- keys[, 2]

  idx <- match(
    interaction(se$well_id, se$plate_id),
    interaction(res$well_id, res$plate_id)
  )

  colData(se)[[imax_col]] <- res$Imax[idx]
  colData(se)[[vhalf_col]] <- res$Vhalf[idx]
  colData(se)[[vmax_col]] <- res$Vmax[idx]

  se
}

fit_boltzmann_se <- function(
    se,
    assay = "Gmin",
    id_cols = c("well_id", "plate_id"),
    vclamp_col = "v_clamp_mV",
    vmax_col = "Vmax.minima",              # if NULL, uses paste0("Vmax.", tolower(assay))
    use_activation_limb = TRUE,
    min_points = 8,
    r2_min = 0.85,
    k_abs_max = 50,
    require_span = 0.4,           # require max(y)-min(y) on limb
    strict_norm = TRUE            # for Gmin, enforce y in [0,1] with tolerance
) {

  suffix <- tolower(assay)
  if (is.null(vmax_col)) vmax_col <- paste0("Vmax.", suffix)

  out_vhalf <- paste0("Vhalf_boltz.", suffix)
  out_k     <- paste0("k_boltz.", suffix)
  out_y0    <- paste0("y0_boltz.", suffix)
  out_ymax  <- paste0("ymax_boltz.", suffix)
  out_conv  <- paste0("boltz_converged.", suffix)
  out_r2    <- paste0("boltz_r2.", suffix)
  out_ok    <- paste0("boltz_ok.", suffix)
  out_msg   <- paste0("boltz_msg.", suffix)

  # Melt once
  dat <- sechm::meltSE(
    se,
    features = rownames(se),
    assayName = assay,
    rowDat.columns = vclamp_col
  )

  # Create the group key in the melted data (well_id.plate_id)
  key <- interaction(dat[[id_cols[1]]], dat[[id_cols[2]]], drop = TRUE)
  dat$.key <- key

  # Lookup Vmax per cell from colData(se) (computed in pass 1)
  cd_key <- interaction(colData(se)[[id_cols[1]]], colData(se)[[id_cols[2]]], drop = TRUE)
  if (vmax_col %in% colnames(colData(se))) {

    vmax_lookup <- setNames(colData(se)[[vmax_col]], cd_key)

  } else {

    message("Vmax column not found — computing Vmax from ", assay)

    # Melt only once for this assay
    dat_vmax <- sechm::meltSE(
      se,
      features = rownames(se),
      assayName = assay,
      rowDat.columns = vclamp_col
    )

    # Create key
    dat_vmax$.key <- interaction(
      dat_vmax[[id_cols[1]]],
      dat_vmax[[id_cols[2]]],
      drop = TRUE
    )

    # Split per cell
    groups_vmax <- split(dat_vmax, dat_vmax$.key)

    vmax_lookup <- sapply(groups_vmax, function(subdat) {

      x <- subdat[[vclamp_col]]
      y <- subdat[[assay]]

      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]; y <- y[ok]

      if (length(y) == 0) return(NA_real_)

      ymax <- max(y, na.rm = TRUE)

      if (!is.finite(ymax)) return(NA_real_)

      # If multiple identical maxima, choose the first occurrence
      # (more stable for activation curves)
      idx <- which(y == ymax)[1]

      x[idx]
    })

  }

  groups <- split(dat, dat$.key)

  res <- lapply(names(groups), function(k) {
    subdat <- groups[[k]]

    x <- subdat[[vclamp_col]]
    y <- subdat[[assay]]

    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]

    if (length(y) < min_points) {
      return(c(Vhalf=NA, k=NA, y0=NA, ymax=NA, conv=0, r2=NA, ok=0, msg="too_few_points"))
    }

    # Optionally fit only activation limb (<= Vmax from pass 1)
    if (use_activation_limb && !is.na(vmax_lookup[[k]])) {
      keep <- x <= vmax_lookup[[k]]
      x <- x[keep]; y <- y[keep]
    }

    if (length(y) < min_points) {
      return(c(Vhalf=NA, k=NA, y0=NA, ymax=NA, conv=0, r2=NA, ok=0, msg="too_few_points_limb"))
    }

    # Sanity for Gmin-like normalization
    if (strict_norm) {
      # tolerate small numerical overshoot
      if (min(y, na.rm=TRUE) < -0.1 || max(y, na.rm=TRUE) > 1.1) {
        return(c(Vhalf=NA, k=NA, y0=NA, ymax=NA, conv=0, r2=NA, ok=0, msg="not_normalized"))
      }
    }

    span <- max(y) - min(y)
    if (!is.finite(span) || span < require_span) {
      return(c(Vhalf=NA, k=NA, y0=NA, ymax=NA, conv=0, r2=NA, ok=0, msg="insufficient_span"))
    }

    # For stability, fit y scaled 0..1 (even if already normalized)
    y01 <- (y - min(y)) / (max(y) - min(y))

    # Starting guesses: Vhalf at ~0.5, k ~ 8-12
    v0 <- x[which.min(abs(y01 - 0.5))]
    start <- list(V_half = v0, k = 10, y0 = 0, ymax = 1)

    fit <- tryCatch(
      nls(
        y01 ~ y0 + (ymax - y0) / (1 + exp((V_half - x) / k)),
        start = start,
        control = nls.control(maxiter = 200, warnOnly = TRUE)
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(c(Vhalf=NA, k=NA, y0=NA, ymax=NA, conv=0, r2=NA, ok=0, msg="nls_error"))
    }

    conv <- isTRUE(fit$convInfo$isConv)

    co <- coef(fit)
    Vhalf <- unname(co["V_half"])
    kk    <- unname(co["k"])
    y0    <- unname(co["y0"])
    ymax  <- unname(co["ymax"])

    yhat <- predict(fit)
    rss <- sum((y01 - yhat)^2)
    tss <- sum((y01 - mean(y01))^2)
    r2 <- if (tss > 0) 1 - rss/tss else NA_real_

    # QC: convergence + sane slope + good fit + monotonic-ish amplitude
    ok_fit <- conv &&
      is.finite(Vhalf) && is.finite(kk) &&
      abs(kk) <= k_abs_max &&
      is.finite(r2) && r2 >= r2_min &&
      (ymax - y0) >= 0.3

    if (!ok_fit) Vhalf <- NA_real_

    c(Vhalf=Vhalf, k=kk, y0=y0, ymax=ymax,
      conv=as.numeric(conv), r2=r2, ok=as.numeric(ok_fit),
      msg=ifelse(ok_fit, "ok", "qc_fail"))
  })

  res <- do.call(rbind, res)
  res <- as.data.frame(res)
  res$.key <- names(groups)

  # Map back to colData
  idx <- match(cd_key, res$.key)

  colData(se)[[out_vhalf]] <- as.numeric(res$Vhalf[idx])
  colData(se)[[out_k]]     <- as.numeric(res$k[idx])
  colData(se)[[out_y0]]    <- as.numeric(res$y0[idx])
  colData(se)[[out_ymax]]  <- as.numeric(res$ymax[idx])
  colData(se)[[out_conv]]  <- as.numeric(res$conv[idx])
  colData(se)[[out_r2]]    <- as.numeric(res$r2[idx])
  colData(se)[[out_ok]]    <- as.numeric(res$ok[idx])
  colData(se)[[out_msg]]   <- as.character(res$msg[idx])
  colData(se)
  se
}

# ----------------------------------------------------------------------------
# Internal helpers for aggregateSE
# ----------------------------------------------------------------------------

# Build a character group-key vector from one or more data.frame columns.
.make_key <- function(df, cols, sep = ".") {
  if (length(cols) == 1) return(as.character(df[[cols]]))
  do.call(paste, c(lapply(cols, function(c) as.character(df[[c]])), list(sep = sep)))
}

# Summarise a single data.frame into one row per unique key.
# - grouping columns  → taken from first occurrence (they are the key)
# - numeric columns   → aggregated with agg_fun (NA/Inf removed first)
# - character/factor  → kept if constant within group, otherwise dropped
# - anything else     → dropped silently
.summarise_metadata <- function(df, by_cols, agg_fun, key_vec, unique_keys) {

  first_idx <- match(unique_keys, key_vec)

  out <- df[first_idx, by_cols, drop = FALSE]
  rownames(out) <- unique_keys

  other_cols <- setdiff(colnames(df), by_cols)

  keep_cols <- list()
  for (col in other_cols) {
    x <- df[[col]]

    if (is.numeric(x)) {
      vals <- sapply(unique_keys, function(k) {
        v <- x[key_vec == k]
        v <- v[is.finite(v)]
        if (length(v) > 0) agg_fun(v) else NA_real_
      })
      keep_cols[[col]] <- vals

    } else if (is.character(x) || is.factor(x)) {
      vals <- sapply(unique_keys, function(k) {
        v <- unique(as.character(x[key_vec == k]))
        v <- v[!is.na(v)]
        if (length(v) == 1) v else NA_character_
      })
      # drop if all NA (i.e. never constant within any group)
      if (!all(is.na(vals))) keep_cols[[col]] <- vals
    }
    # other types silently dropped
  }

  if (length(keep_cols) > 0) {
    out <- cbind(out, as.data.frame(keep_cols, stringsAsFactors = FALSE))
  }
  out
}


#' Aggregate a SummarizedExperiment by row and column grouping factors
#'
#' Produces a new \code{SingleCellExperiment} where:
#' \itemize{
#'   \item Each \strong{row} corresponds to a unique combination of
#'     \code{row_by} columns in \code{rowData} (e.g. liquid period).
#'   \item Each \strong{column} corresponds to a unique combination of
#'     \code{col_by} columns in \code{colData} (e.g. Condition + Plate_ID).
#'   \item Each \strong{assay} is a summary matrix named
#'     \code{<Assay>_<funName>} (e.g. \code{Minima_mean}, \code{Minima_sd}).
#'     All values in the input sub-matrix (row-group × col-group block) are
#'     pooled and passed to the function after removing non-finite values.
#'   \item \strong{rowData} retains the grouping columns plus any
#'     \code{rowData} column that is numeric (aggregated with
#'     \code{row_agg_fun}) or constant within each row group.
#'   \item \strong{colData} retains the grouping columns plus any
#'     \code{colData} column that is numeric (aggregated with
#'     \code{col_agg_fun}) or constant within each column group.
#'     Non-summarisable columns are silently dropped.
#' }
#'
#' @param se A \code{SummarizedExperiment} or \code{SingleCellExperiment}.
#' @param row_by Character vector of \code{rowData} column names to group rows
#'   by (e.g. \code{"LP"}).  \code{NULL} collapses all sweeps into one row.
#' @param col_by Character vector of \code{colData} column names to group
#'   columns by (e.g. \code{c("Condition", "plate_id")}).
#' @param assayList Character vector of assay names to aggregate.  Defaults to
#'   all assays in \code{se}.
#' @param funs Named list of functions.  Each function is applied to the
#'   pooled numeric values of each row-group × col-group block.  Names become
#'   part of the output assay names (e.g. \code{list(mean = mean, sd = sd)}).
#' @param col_agg_fun Function used to summarise numeric \code{colData}
#'   columns within each column group.  Default \code{mean}.
#' @param row_agg_fun Function used to summarise numeric \code{rowData}
#'   columns within each row group.  Default \code{mean}.
#' @param sep Separator used to paste multi-column keys into row/column names.
#'   Default \code{"."}.
#' @param suffix Optional string appended to every output assay name after the
#'   function name (e.g. \code{".exp1"} → \code{Minima_mean.exp1}).
#'   Default \code{""}.
#' @param ... Additional arguments forwarded to every function in \code{funs}.
#'
#' @return A new \code{SingleCellExperiment} with aggregated assays, rowData,
#'   and colData.
#'
#' @examples
#' \dontrun{
#' # 3 liquid periods x (Condition x Plate_ID) columns, mean + sd assays
#' se_agg <- aggregateSE(
#'   se,
#'   row_by   = "LP",
#'   col_by   = c("Condition", "plate_id"),
#'   funs     = list(mean = mean, sd = sd)
#' )
#' # assayNames(se_agg) -> "Minima_mean", "Minima_sd", "Maxima_mean", ...
#' # nrow(se_agg)       -> 3   (one per LP)
#' # ncol(se_agg)       -> n unique Condition.Plate_ID combinations
#' }
#' @importFrom SummarizedExperiment assay assayNames rowData colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @export
aggregateSE <- function(se,
                        row_by       = NULL,
                        col_by,
                        assayList    = assayNames(se),
                        funs         = list(mean = mean),
                        col_agg_fun  = mean,
                        row_agg_fun  = mean,
                        sep          = ".",
                        suffix       = "",
                        ...) {

  assayList <- assayList[assayList %in% assayNames(se)]
  if (length(assayList) == 0) stop("No valid assays found in 'assayList'.")
  if (length(funs) == 0)      stop("'funs' must be a non-empty named list.")
  if (is.null(names(funs)) || any(names(funs) == ""))
    stop("Every element of 'funs' must be named.")

  cd <- as.data.frame(colData(se))
  rd <- as.data.frame(rowData(se))

  # --- column groups -------------------------------------------------------
  col_by  <- col_by[col_by %in% colnames(cd)]
  if (length(col_by) == 0) stop("None of 'col_by' columns found in colData.")
  col_key  <- .make_key(cd, col_by, sep)
  ucol_key <- unique(col_key)
  n_cols   <- length(ucol_key)

  # --- row groups ----------------------------------------------------------
  if (!is.null(row_by)) {
    row_by   <- row_by[row_by %in% colnames(rd)]
    if (length(row_by) == 0) stop("None of 'row_by' columns found in rowData.")
    row_key  <- .make_key(rd, row_by, sep)
    urow_key <- unique(row_key)
  } else {
    row_key  <- rep("all", nrow(se))
    urow_key <- "all"
  }
  n_rows <- length(urow_key)

  # --- build aggregated assays ---------------------------------------------
  new_assays <- list()

  for (assayName in assayList) {
    mat <- assay(se, assayName)

    # pre-split columns and rows by group for efficiency
    col_groups <- lapply(ucol_key, function(k) which(col_key == k))
    row_groups <- lapply(urow_key, function(k) which(row_key == k))

    for (fun_name in names(funs)) {
      fun      <- funs[[fun_name]]
      new_mat  <- matrix(NA_real_, nrow = n_rows, ncol = n_cols,
                         dimnames = list(urow_key, ucol_key))

      for (i in seq_len(n_rows)) {
        for (j in seq_len(n_cols)) {
          vals <- as.vector(mat[row_groups[[i]], col_groups[[j]], drop = FALSE])
          vals <- vals[is.finite(vals)]
          if (length(vals) > 0)
            new_mat[i, j] <- fun(vals, ...)
        }
      }

      out_name <- paste0(assayName, "_", fun_name, suffix)
      new_assays[[out_name]] <- new_mat
    }
  }

  # --- summarise rowData ---------------------------------------------------
  new_rd <- if (!is.null(row_by) && nrow(rd) > 0) {
    .summarise_metadata(rd, row_by, row_agg_fun, row_key, urow_key)
  } else {
    data.frame(row.names = urow_key)
  }

  # --- summarise colData ---------------------------------------------------
  new_cd <- .summarise_metadata(cd, col_by, col_agg_fun, col_key, ucol_key)

  # --- assemble new SCE ----------------------------------------------------
  SingleCellExperiment::SingleCellExperiment(
    assays  = new_assays,
    rowData = S4Vectors::DataFrame(new_rd),
    colData = S4Vectors::DataFrame(new_cd)
  )
}


