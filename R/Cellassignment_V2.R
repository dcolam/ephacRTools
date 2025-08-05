#' Add your imaging columns to your SummarizedExperiment
#' @param SE Your SummarizedExperiment
#' @param imaging_df Your dataframe containing the imaging data
#' @param join_keys The columns you want to join by (match the wells by Well and Plate_ID)
#' @param columns_to_add Columns from your imaging dataset you're interested in adding to your SE
#' @return An SE with your added imaging columns
#' @export
add_imagingcols_toSE <- function(SE, imaging_df, join_keys = c("Well", "Plate_ID"), columns_to_add = NULL) {
  coldata <- as.data.frame(colData(SE))

  if (!is.null(columns_to_add)){
    missing_cols <- setdiff(columns_to_add, colnames(imaging_df))
    if (length(missing_cols) > 0) {
      stop("The following columns are missing in imaging_df", paste(missing_cols, collapse = ", "))
    }
    imaging_df <- imaging_df[, unique(c(join_keys, columns_to_add)), drop = FALSE]
  } else {
    imaging_df <- imaging_df[, unique(c(join_keys, setdiff(names(imaging_df), join_keys))), drop = FALSE]
  }

  merged_coldata <- left_join(coldata, imaging_df, by = join_keys)

  colData(SE) <- S4Vectors::DataFrame(merged_coldata)

  return(SE)

}
#' Assign cell colors to your wells based on the imaging data using just one value column
#' @param SE Your SummarizedExperiment
#' @param channels The different channels that your used in your experiment
#' @param chan2analyze The specific channels you're interested in using for the cell assignment (could also be all of them)
#' @param join_keys The columns you want to join by (match the wells by Well and Plate_ID)
#' @param value_col The feature in your data that you're using to base the cell assignment on (could be Mean, Area, IntDen..)
#' @param scale_suffix The suffix that you want to have in your log1p scaled channel column names
#' @param color_output The name you wan to give the column that contains your cell color assignments
#' @param classification_method Set whether you prefer the cell assignment to be done based on clustering or based on manually-set thresholds
#' @param threshold_mode Should you choose manual classification, here you can set whether you want to set the threshold for a given quantile or a hard-cut threshold
#' @param red_thresh The manual hard threshold you want to set for the red channel
#' @param green_thresh The manual hard threshold you want to set for the green channel
#' @param red_quant The quantile-based threshold you want to set for the red channel
#' @param green_quant The quantile-based threshold you want to set for the green channel
#' @param c3_quant The quantile-based threshold you want to set for your third channel
#' @param k_centers The number of centers you want for your k-means clustering
#' @param plot_umap Whether you want the function to return the umap plot from the clustering if you ran clustering classification
#' @return A SE with your scaled channel columns and cell assignments
#' @export
assign_cell_FINAL<- function(SE,
                               channels = c("C1", "C2", "C3"),
                               chan2analyze = c("C2", "C3"),
                               join_keys = c("Well", "Plate_ID"),
                               value_col = "Mean",
                               scale_suffix = "_scaled",
                               color_output = "Color",
                               classification_method = c("manual", "clustering"),
                               threshold_mode = c("manual", "quantile"),
                               red_thresh = 8.5,
                               green_thresh = 5.2,
                               c3_thresh = 5.6,
                               red_quant = 0.5,
                               green_quant = 0.75,
                               c3_quant = 0.7,
                               k_centers = 4,
                               plot_umap = FALSE) {

  classification_method <- match.arg(classification_method)

  meta_df <- as.data.frame(colData(SE))
  channel_data <- list()

  for (channel in channels) {
    if (!channel %in% colnames(colData(SE))) {
      stop("Channel ", channel, " not found in colData(SE).")
    }

    channel_df <- as.data.frame(colData(SE)[[channel]])
    full_df <- bind_cols(channel_df, meta_df[join_keys])

    summary_df <- full_df %>%
      group_by(across(all_of(join_keys))) %>%
      summarize(!!paste0("feature_", channel) := mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

    channel_data[[channel]] <- summary_df
  }

  # Step 2: Join summaries and scale
  well_summary <- Reduce(function(x, y) full_join(x, y, by = join_keys), channel_data) %>%
    mutate(across(starts_with("feature_"), ~replace_na(., 0))) %>%
    mutate(
      Row = substr(Well, 1, 1),
      Column = substr(Well, 2, 3))

  scaled_cols <- c()
  for (channel in channels) {

    if (!channel %in% colnames(colData(SE))) {
      stop("Channel ", channel, " not found in colData(SE).")
    }
    if (!"Mean" %in% colnames(colData(SE)[[channel]])) {
      stop("Channel ", channel, " does not contain 'Mean' column.")
    }

    feat_col <- paste0("feature_", channel)
    scaled_col <- paste0(feat_col, scale_suffix)

    log_vals <- log1p(well_summary[[feat_col]])

    cat(channel, ": raw min =", min(well_summary[[feat_col]], na.rm = TRUE),
        " | log1p max =", max(log_vals, na.rm = TRUE), "\n")

    scaled <- log_vals

    well_summary[[scaled_col]] <- scaled

    scaled_cols <- c(scaled_cols, scaled_col)
  }

  # Step 3: Classification
  if (classification_method == "manual") {
    red_feat <- paste0("feature_", chan2analyze[2], scale_suffix)
    green_feat <- paste0("feature_", chan2analyze[1], scale_suffix)

    red_vals <- well_summary[[red_feat]]
    green_vals <- well_summary[[green_feat]]

    if (length(chan2analyze > 2)) {
      c3_feat    <- paste0("feature_", chan2analyze[3], scale_suffix)
      c3_vals    <- well_summary[[c3_feat]]
    }

    #red_thresh <- quantile(red_vals[red_vals > 0], red_threshold_quantile, na.rm = TRUE)
    if (threshold_mode == "manual"){
      red_thresh <- red_thresh
      green_thresh <- green_thresh
      #c3_thresh <- c3_thresh
    } #else if (threshold_mode == "quantile"){
      #red_thresh   <- quantile(red_vals, red_quant, na.rm = TRUE)
      #green_thresh <- quantile(green_vals, green_quant, na.rm = TRUE)
      #c3_thresh <- quantile(c3_vals, c3_quant, na.rm = TRUE)

      #cat("Green threshold (quantile", green_quant, ") =", green_thresh, "\n")
      #cat("Red threshold (quantile", red_quant, ") =", red_thresh, "\n")
      #cat("C3 threshold (quantile", c3_quant, ") =", c3_thresh, "\n")
    #}

    well_summary[[color_output]] <- case_when(
      #red_vals < 0.05 & green_vals < 0.05 ~ "black",
      #red_vals >= red_thresh & green_vals >= green_thresh ~ "yellow",
      green_vals >= green_thresh ~ "green",  # Normal green case
      #red_vals >= red_thresh & green_vals < green_thresh & c3_vals >= c3_thresh ~ "green",  # override red to green
      red_vals >= red_thresh & green_vals >= green_thresh ~ "lightblue",
      red_vals >= red_thresh ~ "red",        # Normal red case
      TRUE ~ "black"                         # Default
    )

  } else if (classification_method == "clustering") {
    # UMAP + k-means clustering
    if (length(scaled_cols) < 2) {
      warning("Only one feature provided; skipping UMAP. Clustering will be done in 1D.")
      cluster_data <- well_summary %>%
        select(all_of(scaled_cols)) %>%
        as.matrix()

      set.seed(42)
      km <- kmeans(cluster_data, centers = k_centers, nstart = 50)
      well_summary$cluster <- factor(km$cluster)
      well_summary[[color_output]] <- as.factor(km$cluster)

      well_summary$UMAP1 <- NA
      well_summary$UMAP2 <- NA
    } else {
      cluster_data <- well_summary %>%
        select(all_of(scaled_cols)) %>%
        as.matrix()

      set.seed(42)
      umap_layout <- uwot::umap(cluster_data, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
      colnames(umap_layout) <- c("UMAP1", "UMAP2")
      umap_df <- as.data.frame(umap_layout)

      km <- kmeans(umap_df, centers = k_centers, nstart = 50)
      umap_df$cluster <- factor(km$cluster)

      well_summary$UMAP1 <- umap_df$UMAP1
      well_summary$UMAP2 <- umap_df$UMAP2
      well_summary$cluster <- umap_df$cluster
      well_summary[[color_output]] <- as.factor(km$cluster)

      if (plot_umap) {
        print(
          ggplot(well_summary, aes(x = UMAP1, y = UMAP2, color = .data[[color_output]])) +
            geom_point(alpha = 0.7, size = 2) +
            scale_color_brewer(palette = "Set1", name = color_output) +
            theme_minimal() +
            labs(
              title = paste("UMAP of Scaled Features"),
              subtitle = "Clusters assigned via k-means on UMAP",
              x = "UMAP1", y = "UMAP2"
            )
        )
      }
    }
  }

  # Step 4: Merge into colData
  col_df <- as.data.frame(colData(SE))
  new_columns <- c(scaled_cols, color_output, "UMAP1", "UMAP2", "cluster")
  new_columns <- intersect(new_columns, colnames(well_summary))
  join_cols <- unique(c(join_keys, new_columns))

  merge_df <- left_join(col_df, well_summary[, join_cols], by = join_keys)

  for (col in setdiff(colnames(merge_df), colnames(col_df))) {
    colData(SE)[[col]] <- merge_df[[col]]
  }

  return(SE)
}
#' Assign cell colors to your wells based on the imaging data using multiple value columns
#' @param SE Your SummarizedExperiment
#' @param channels The different channels that your used in your experiment
#' @param chan2analyze The specific channels you're interested in using for the cell assignment (could also be all of them)
#' @param join_keys The columns you want to join by (match the wells by Well and Plate_ID)
#' @param value_col The features in your data that you're using to base the cell assignment on (could be Mean and Area for instance)
#' @param scale_suffix The suffix that you want to have in your log1p scaled channel column names
#' @param color_output The name you wan to give the column that contains your cell color assignments
#' @param classification_method Set whether you prefer the cell assignment to be done based on clustering or based on manually-set thresholds
#' @param threshold_mode Should you choose manual classification, here you can set whether you want to set the threshold for a given quantile or a hard-cut threshold
#' @param red_thresh The manual hard threshold you want to set for the red channel
#' @param green_thresh The manual hard threshold you want to set for the green channel
#' @param c3_thresh The manual hard threshold you want to set for your third channel
#' @param red_quant The quantile-based threshold you want to set for the red channel
#' @param green_quant The quantile-based threshold you want to set for the green channel
#' @param c3_quant The quantile-based threshold you want to set for your third channel
#' @param k_centers The number of centers you want for your k-means clustering
#' @param plot_umap Whether you want the function to return the umap plot from the clustering if you ran clustering classification
#' @return A SE with your scaled channel columns and cell assignments
#' @export
assign_cell_FINAL_valcol <- function(SE,
                              channels = c("C1", "C2", "C3"),
                              chan2analyze = c("C2", "C3"),
                              join_keys = c("Well", "Plate_ID"),
                              value_cols = c("Mean", "Area"),
                              scale_suffix = "_scaled",
                              color_output = "Color",
                              classification_method = c("manual", "clustering"),
                              threshold_mode = c("manual", "quantile"),
                              red_thresh = 8.5,
                              green_thresh = 5.2,
                              c3_thresh = 5.6,
                              red_quant = 0.5,
                              green_quant = 0.75,
                              c3_quant = 0.7,
                              k_centers = 4,
                              plot_umap = FALSE) {

  classification_method <- match.arg(classification_method)

  meta_df <- as.data.frame(colData(SE))
  channel_data <- list()

  for (channel in channels) {
    if (!channel %in% colnames(colData(SE))) {
      stop("Channel ", channel, " not found in colData(SE).")
    }

    missing_vals <- setdiff(value_cols, colnames(colData(SE)[[channel]]))
    if (length(missing_vals) > 0) {
      stop("Channel ", channel, " is missing required value columns: ", paste(missing_vals, collapse = ", "))
    }

    channel_df <- as.data.frame(colData(SE)[[channel]])
    full_df <- bind_cols(channel_df, meta_df[join_keys])

    summary_df <- full_df %>%
      group_by(across(all_of(join_keys))) %>%
      summarise(
        across(all_of(value_cols), ~mean(.x, na.rm = TRUE),
               .names = paste0("feature_", channel, ".{.col}")),
        .groups = "drop"
      )

    channel_data[[channel]] <- summary_df
  }

  well_summary <- Reduce(function(x, y) full_join(x, y, by = join_keys), channel_data) %>%
    mutate(across(starts_with("feature_"), ~replace_na(., 0))) %>%
    mutate(
      Row = substr(Well, 1, 1),
      Column = substr(Well, 2, 3)
    )

  scaled_cols <- c()
  for (channel in channels) {
    for (val_col in value_cols) {
      feat_col <- paste0("feature_", channel, ".", val_col)
      scaled_col <- paste0(feat_col, scale_suffix)

      if (!feat_col %in% colnames(well_summary)) {
        next
      }

      log_vals <- log1p(well_summary[[feat_col]])

      cat(channel, ".", val_col, ": raw min =", min(well_summary[[feat_col]], na.rm = TRUE),
          " | log1p max =", max(log_vals, na.rm = TRUE), "\n")

      scaled <- log_vals
      well_summary[[scaled_col]] <- scaled
      scaled_cols <- c(scaled_cols, scaled_col)
    }
  }

  if (classification_method == "manual") {
    # Assuming user wants to classify based on first value_col (e.g., "Mean")
    val_metric <- value_cols[1]

    red_feat <- paste0("feature_", chan2analyze[2], ".", val_metric, scale_suffix)
    green_feat <- paste0("feature_", chan2analyze[1], ".", val_metric, scale_suffix)

    red_vals <- well_summary[[red_feat]]
    green_vals <- well_summary[[green_feat]]

    if (length(chan2analyze) > 2) {
      c3_feat <- paste0("feature_", chan2analyze[3], ".", val_metric, scale_suffix)
      c3_vals <- well_summary[[c3_feat]]
    }

    if (threshold_mode == "manual") {
      # thresholds already provided
    } else if (threshold_mode == "quantile") {
      red_thresh <- quantile(red_vals, red_quant, na.rm = TRUE)
      green_thresh <- quantile(green_vals, green_quant, na.rm = TRUE)
      if (exists("c3_vals")) c3_thresh <- quantile(c3_vals, c3_quant, na.rm = TRUE)
    }

    well_summary[[color_output]] <- case_when(
      green_vals >= green_thresh ~ "green",
      red_vals >= red_thresh & green_vals >= green_thresh ~ "lightblue",
      red_vals >= red_thresh ~ "red",
      TRUE ~ "black"
    )
  }

  if (classification_method == "clustering") {
    if (length(scaled_cols) < 2) {
      warning("Only one feature provided; skipping UMAP. Clustering will be done in 1D.")
      cluster_data <- well_summary %>%
        select(all_of(scaled_cols)) %>%
        as.matrix()

      set.seed(42)
      km <- kmeans(cluster_data, centers = k_centers, nstart = 50)
      well_summary$cluster <- factor(km$cluster)
      well_summary[[color_output]] <- as.factor(km$cluster)

      well_summary$UMAP1 <- NA
      well_summary$UMAP2 <- NA
    } else {
      cluster_data <- well_summary %>%
        select(all_of(scaled_cols)) %>%
        as.matrix()

      set.seed(42)
      umap_layout <- uwot::umap(cluster_data, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
      colnames(umap_layout) <- c("UMAP1", "UMAP2")
      umap_df <- as.data.frame(umap_layout)

      km <- kmeans(umap_df, centers = k_centers, nstart = 50)
      umap_df$cluster <- factor(km$cluster)

      well_summary$UMAP1 <- umap_df$UMAP1
      well_summary$UMAP2 <- umap_df$UMAP2
      well_summary$cluster <- umap_df$cluster
      well_summary[[color_output]] <- as.factor(km$cluster)

      if (plot_umap) {
        print(
          ggplot(well_summary, aes(x = UMAP1, y = UMAP2, color = .data[[color_output]])) +
            geom_point(alpha = 0.7, size = 2) +
            scale_color_brewer(palette = "Set1", name = color_output) +
            theme_minimal() +
            labs(
              title = paste("UMAP of Scaled Features"),
              subtitle = "Clusters assigned via k-means on UMAP",
              x = "UMAP1", y = "UMAP2"
            )
        )
      }
    }
  }

  # Step 4: Merge into colData
  col_df <- as.data.frame(colData(SE))
  new_columns <- c(scaled_cols, color_output, "UMAP1", "UMAP2", "cluster")
  new_columns <- intersect(new_columns, colnames(well_summary))
  join_cols <- unique(c(join_keys, new_columns))

  merge_df <- left_join(col_df, well_summary[, join_cols], by = join_keys)

  for (col in setdiff(colnames(merge_df), colnames(col_df))) {
    colData(SE)[[col]] <- merge_df[[col]]
  }

  return(SE)
}
#' Generate density or ECDF plots based on your SE's log1p scaled channels to look at the general distribution
#' @param SE Your SummarizedExperiment
#' @param channels The different channels that your used in your experiment
#' @param join_keys The columns you want to join by (match the wells by Well and Plate_ID)
#' @param value_col The feature in your data that you're using to base the cell assignment on (could be Mean, Area, IntDen..)
#' @param channel_colors The colors you want your channels to be in the plots
#' @return A density or ECDF plot
#' @export
plot_intensity_distribution <- function(SE,
                                        channels = c("C2", "C3"),
                                        value_col = "Mean",
                                        join_keys = c("Well", "Plate_ID"),
                                        plot_type = c("density", "ecdf"),
                                        channel_colors = c(C1 = "forestgreen", C2 = "lightcoral", C3 = "plum")
                                        ) {
  require(ggplot2)
  require(dplyr)

  plot_type <- match.arg(plot_type)

  meta_df <- as.data.frame(colData(SE))
  plot_data <- data.frame()

  channel_colors <- channel_colors[channels]

  for (channel in channels) {
    if (!channel %in% colnames(colData(SE))) {
      warning(paste("Channel", channel, "not found in colData. Skipping."))
      next
    }

    channel_df <- as.data.frame(colData(SE)[[channel]])
    full_df <- bind_cols(channel_df, meta_df[join_keys])

    summary_df <- full_df %>%
      group_by(across(all_of(join_keys))) %>%
      summarize(Intensity = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

    # log1p transform and scale
    log_vals <- log1p(summary_df$Intensity)
    scaled_vals <- log_vals

    df <- data.frame(
      Channel = channel,
      Scaled_Value = scaled_vals
    )
    plot_data <- bind_rows(plot_data, df)
  }

  p <- ggplot(plot_data, aes(x = Scaled_Value, color = Channel, fill = Channel))

  if (plot_type == "density") {
    p <- p +
      geom_density(alpha = 0.6) +
      labs(
        title = "Density Plot of log1p Scaled Feature Intensities",
        x = "Intensity", y = "Density"
      )
  } else if (plot_type == "ecdf") {
    p <- p +
      stat_ecdf(geom = "step", linewidth = 1) +
      labs(
        title = "ECDF Plot of log1p Scaled Feature Intensities",
        x = "Intensity", y = "Cumulative Probability"
      )
  }

  p <- p + theme_minimal() +
    scale_color_manual(values=channel_colors) +
    scale_fill_manual(values=channel_colors)

  print(p)
}
#' Generate a well plate plot that shows the color assignment for the given wells
#' @param SE Your SummarizedExperiment
#' @param color_col The column that contains your color assignment
#' @param title The title you want to give your plot
#' @param colors How you want the colors to be plotted
#' @return A well plate plot showing the color assignment
#' @export
plot_well_assignment <- function(SE,
                                 color_col = "Color_manual",
                                 title = NULL,
                                 colors = NULL) {
  se_df <- as.data.frame(colData(SE))
  se_df$Column <- as.numeric(se_df$Column)

  if (is.null(title)) {
    title <- paste("Well Plate Cell Assignment for", unique(se_df$Plate_ID))
  }

  # Get unique values in the color column
  color_values <- sort(unique(se_df[[color_col]]))

  # If colors aren't provided, generate distinct colors using RColorBrewer or default palette
  if (is.null(colors)) {
    if (requireNamespace("RColorBrewer", quietly = TRUE) && length(color_values) <= 9) {
      palette_colors <- RColorBrewer::brewer.pal(max(3, length(color_values)), "Set1")
    } else {
      palette_colors <- scales::hue_pal()(length(color_values))
    }
    colors <- setNames(palette_colors, color_values)
  } else {
    # If colors provided but missing some levels, assign defaults to missing ones
    missing_levels <- setdiff(color_values, names(colors))
    if (length(missing_levels) > 0) {
      palette_colors <- scales::hue_pal()(length(missing_levels))
      colors <- c(colors, setNames(palette_colors, missing_levels))
    }
  }

  ggplot(se_df, aes(x = Column, y = Row, fill = .data[[color_col]])) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = 1:24) +
    scale_y_discrete(limits = rev(sort(unique(se_df$Row)))) +
    geom_text(aes(label = Well), color = "black", size = 3) +
    labs(title = title, fill = color_col) +
    theme_minimal()
}
#' Generate a well plate plot that shows the whether the well is classified as induced, not-induced, ambiguous or unknown based on the well color assignment
#' @param SE Your SummarizedExperiment
#' @param color_col The column that contains your color assignment
#' @param values The colors you want your categories to be assigned
#' @param title The title you want to give your plot
#' @return A well plate plot showing the induction distribution
#' @export
plot_well_binary_assignment <- function(SE, color_col = "Color_binary", values= c(
                                          "red" = "red",
                                          "green" = "green",
                                          "navy" = "navy",
                                          "yellow" = "yellow",
                                          "orange" = "orange",
                                          "turquoise" = "turquoise",
                                          "white" = "white",
                                          "black" = "grey"
                                        ),
                                        title = NULL) {
  se_df <- as.data.frame(colData(SE))
  se_df$Column <- as.numeric(se_df$Column)

  if (is.null(title)) {
    title <- paste("Well Plate Cell Binary Assignment for", unique(se_df$Plate_ID))
  }

  ggplot(se_df, aes(x = Column, y = Row, fill = .data[[color_col]])) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_manual(
      values = values,
      na.value = "grey"
    ) +
    scale_x_continuous(breaks = 1:24) +
    scale_y_discrete(limits = rev(sort(unique(se_df$Row)))) +
    geom_text(aes(label = Well), color = "black", size = 3) +
    labs(title = title, x = "Column", y = "Row") +
    theme_minimal()
}



