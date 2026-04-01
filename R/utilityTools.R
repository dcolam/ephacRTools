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
#' se <- subsetSE(se, Plate_ID %in% c("P1", "P2") & QC == "Pass")
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

#' Add column-wise aggregation such as mean of any given assay and store it into colData
#' @param se SummarizedExperiment Object with the Ephys-Data
#' @param assayNames list of assays to aggregate column-wise
#' @param fun to be changed, function used to aggregate, for now only mean available
#' @param sweeps which sweeps to take, tbd, default is all
#' @return se with updated colData
#' @export
colAG <- function(se, assayList, fun=mean, sweeps=row.names(se)){
  ## check whether assayList contains assays
  assayList <- assayList[assayList %in% assayNames(se)]

  for (assayName in assayList) {
    colName <- paste(assayName, "mean", sep = "_")
    subse <- se[sweeps,]
    se[[colName]] <- colMeans(assay(subse, assayName), na.rm = TRUE)
  }
  return(se)
}
#' Perform cell-wise dimensionality reduction based on assays and colData.
#' Results are stored in reducedDims
#' @param se SummarizedExperiment Object with the Ephys-Data
#' @param assayNames list of assays to include
#' @param colNames list of columns from colData
#' @param scaling option of which scaling to apply, within assay or all features together (scaling = "global"), default within assay
#' @param byRow scaling by row instead of by column (after transforming), default FALSE
#' @param method list of types of reductionality methods, default all (pca, tsne, umap)
#' @param k_clusters number of clusters
#' @return se with updated results
#' @export
reducedDim.Cellwise <- function(se, assayList=c(), colNames=c(), scaling = "within",
                                byRow=FALSE,center = FALSE, method=c("pca", "tsne", "umap"),
                                k_clusters=3){

  if(length(assayList) != 0){
  assayList <- assayList[assayList %in% assayNames(se)]
  pca_data <- lapply(assayList, function(x){
    if(scaling == "within"){
      #scale(t(assay(se, x)))
      temp <- sechm::safescale(t(assay(se, x)),center = center, byRow = byRow)
    }else{
      temp <-t(assay(se, x))
    }
    temp <- as.data.frame(temp)
    colnames(temp) <- paste(x, colnames(temp))
    temp
    }

    )
  names(pca_data) <- assayList
  }else{pca_data <- list()}


  if(length(colNames) != 0){
  flattened.df <- as.data.frame(colData(se))
  colNames <- colNames[colNames %in% colnames(flattened.df)]

  col_Data <- lapply(colNames, function(x){
                     if(scaling == "within"){

                       temp <- sechm::safescale(flattened.df[[x]], center = F, byRow = byRow)
                     }else{
                       temp <-flattened.df[[x]]
                     }
    temp <- as.data.frame(temp)
    colnames(temp) <- paste(x, colnames(temp))
    temp
                     })

  names(col_Data) <- colNames
  pca_data <- list(pca_data, col_Data)
  }
  #names(pca_data) <- assayList
  pca_data <- dplyr::bind_cols(pca_data)
  #print(class(pca_data))
  if(scaling == "global"){
    print(pca_data)
    pca_data <- sechm::safescale(as.matrix(pca_data), byRow = byRow)
  }
  ## handling missing values
  pca_data <-as.data.frame(pca_data)
  #return(pca_data)
  pca_data <- pca_data[, colSums(is.na(pca_data)) != nrow(pca_data)]
  pca_data[is.na(pca_data)] <- 0

  pca_result <- prcomp(pca_data, rank=50)

  tsne_data <- Rtsne::Rtsne(pca_data, pca = TRUE,  check_duplicates = FALSE)

  tsne_data <- tsne_data$Y %>%
    as.data.frame()%>%
    dplyr::rename(tsne1="V1",
           tsne2="V2")

  umap_data <- umap::umap(pca_data)
  umap_df <- umap_data$layout %>%
    as.data.frame()%>%
    dplyr::rename(UMAP1="V1",
           UMAP2="V2")

  SingleCellExperiment::reducedDims(se) <- list(PCA=pca_result$x, TSNE=S4Vectors::DataFrame(tsne_data), UMAP=S4Vectors::DataFrame(umap_df))

  se$cluster.umap <- as.factor(kmeans(SingleCellExperiment::reducedDim(se, "UMAP")[,1:2], k_clusters, iter.max = 100)$cluster)
  se$cluster.tsne <- as.factor(kmeans(SingleCellExperiment::reducedDim(se, "TSNE")[,1:2], k_clusters, iter.max = 100)$cluster)
  se$cluster.pca <- as.factor(kmeans(SingleCellExperiment::reducedDim(se, "PCA")[,1:2], k_clusters, iter.max = 100)$cluster)

  return(se)
}
#' Wrapper function for plotting dimensionality plots
#' plot Dimensionality Reduction with ggplot
#' @param se SummarizedExperiment Object with reducedDim data
#' @param redDim.method a single character parameter, either, "UMAP", "TSNE", "PCA"
#' @param colorColumns list of columns from colData adding colorings, every element with a different plot, default are clusters
#' @return a ggplot
#' @export
plotDimRed <- function(se, redDim.method, colorColumns = character()) {
  flattened.df <- as.data.frame(colData(se))
  redDF <- as.data.frame(SingleCellExperiment::reducedDim(se, redDim.method))

  clustername <- grep(tolower(redDim.method), colnames(flattened.df), value = TRUE)
  if (length(clustername) == 0) clustername <- colorColumns[1]

  p <- ggplot2::ggplot(redDF, aes(x = redDF[, 1], y = redDF[, 2], color = flattened.df[[clustername]])) +
    ggplot2::geom_point() +
    ggplot2::ggtitle(paste("Colored by", clustername))

  colorColumns <- colorColumns[colorColumns %in% colnames(flattened.df)]

  if (length(colorColumns) > 0) {
    plots <- lapply(colorColumns, function(cols) {
      if(is.numeric(flattened.df[[cols]])){
      ggplot2::ggplot(redDF, aes(x = redDF[, 1], y = redDF[, 2], color = scale(flattened.df[[cols]], center=TRUE))) +
        ggplot2::geom_point() +
        ggplot2::ggtitle(cols)+
       viridis::scale_colour_viridis()

      } else{
        ggplot2::ggplot(redDF, aes(x = redDF[, 1], y = redDF[, 2], color = flattened.df[[cols]])) +
          ggplot2::geom_point() +
          ggplot2::ggtitle(cols)
    }

    })
    return(ggpubr::ggarrange(plotlist = c(list(p), plots)))
  } else {
    return(p)
  }
}
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
  ggplot2::ggplot(melted.se, aes(x=melted.se[[rowCol]], y=value, color=Well)) +
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
#' \code{Well} + \code{Plate_ID}), the function extracts three metrics and stores
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
#' The assay matrix must have sweeps as rows and wells as columns. \code{V_Clamp}
#' must be present in \code{rowData(se)} and contain the holding potential for
#' each sweep, as produced by \code{prepareSE()}.
#'
#' @param se A \code{SingleCellExperiment} (or \code{SummarizedExperiment}) with
#'   step-wise voltage clamp data. Must contain \code{V_Clamp} in
#'   \code{rowData} and \code{Well} / \code{Plate_ID} in \code{colData}.
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
    rowDat.columns = "V_Clamp"
  )

  # Split once by Well + Plate_ID
  groups <- split(dat, interaction(dat$Well, dat$Plate_ID, drop = TRUE))

  res <- lapply(groups, function(subdat) {

    x_vals <- subdat$V_Clamp
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
  colnames(keys) <- c("Well", "Plate_ID")

  res$Well <- keys[, 1]
  res$Plate_ID <- keys[, 2]

  idx <- match(
    interaction(se$Well, se$Plate_ID),
    interaction(res$Well, res$Plate_ID)
  )

  colData(se)[[imax_col]] <- res$Imax[idx]
  colData(se)[[vhalf_col]] <- res$Vhalf[idx]
  colData(se)[[vmax_col]] <- res$Vmax[idx]

  se
}

fit_boltzmann_se <- function(
    se,
    assay = "Gmin",
    id_cols = c("Well", "Plate_ID"),
    vclamp_col = "V_Clamp",
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

  # Create the group key in the melted data (Well.Plate_ID)
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


