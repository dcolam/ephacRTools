#' @importFrom magrittr %>%
NULL
#' @import SingleCellExperiment
NULL
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
reducedDim.Cellwise <- function(se, assayList=c(), colNames=c(), scaling = "within", byRow=FALSE, method=c("pca", "tsne", "umap"), k_clusters=3){

  if(length(assayList) != 0){
  assayList <- assayList[assayList %in% assayNames(se)]
  pca_data <- lapply(assayList, function(x){
    if(scaling == "within"){
      #scale(t(assay(se, x)))
      temp <- sechm::safescale(t(assay(se, x)), byRow = byRow)
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

                       temp <- sechm::safescale(flattened.df[[x]])
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
#' Wrapper function to analyze IV-curves and extract Imax, Vmax and Vhalf.
#' Make sure that Vhalf is present in rowData
#' @param se SummarizedExperiment Object with reducedDim data
#' @param assay assay to be analyzed
#' @param inward boolean stating whether the IV-curve shows in- or outward current
#' @param getErev boolean to try to get Erev or value where y = 0
#' @return se with updated colData
#' @export
get_metric <- function(se, assay = "Minima", inward = TRUE) {
  suffix <- tolower(assay)

  imax_col <- paste0("Imax.", suffix)
  vhalf_col <- paste0("Vhalf.", suffix)
  vmax_col <- paste0("Vmax.", suffix)

  wells <- unique(se$Well)
  results <- data.frame(Well = wells,
                        Imax = NA_real_,
                        Vhalf = NA_real_,
                        Vmax = NA_real_)

  assay_data <- assay(se, assay)
  v_clamp <- rowData(se)$V_Clamp
  well_ids <- se$Well

  for (i in seq_along(wells)) {
    well <- wells[i]
    indices <- which(well_ids == well)

    x_vals <- v_clamp
    y_vals <- assay_data[, indices]

    # Skip if all NA
    if (all(is.na(y_vals))) next

    if (inward) {
      Imax <- min(y_vals, na.rm = TRUE)
      Vmax1 <- min(x_vals[y_vals == Imax])
    } else {
      Imax <- max(y_vals, na.rm = TRUE)
      Vmax1 <- max(x_vals[y_vals == Imax])
    }

    y_vals[!complete.cases(y_vals)] <- 1  # Avoid smooth.spline errors

    spl <- smooth.spline(x_vals, y = y_vals)
    fit <- predict(spl, seq(min(x_vals), max(x_vals), length.out = 100))

    Vhalf <- NA
    fit_sub <- fit$x[fit$x < Vmax1]
    pred_sub <- fit$y[fit$x < Vmax1]

    if (length(pred_sub) > 0) {
      idx <- which.min(abs(pred_sub - Imax / 2))
      Vhalf <- fit_sub[idx]
    }

    results$Imax[i] <- Imax
    results$Vhalf[i] <- Vhalf
    results$Vmax[i] <- Vmax1
  }

  coldata <- colData(se)
  coldata[[imax_col]] <- results$Imax[match(se$Well, results$Well)]
  coldata[[vhalf_col]] <- results$Vhalf[match(se$Well, results$Well)]
  coldata[[vmax_col]] <- results$Vmax[match(se$Well, results$Well)]

  colData(se) <- coldata
  return(se)
}


