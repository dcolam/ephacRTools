#' @importFrom magrittr %>%
NULL

# ============================================================================
# Particle-level cell classification pipeline
# ============================================================================

#' Filter particles by scaled intensity
#'
#' Scales the \code{Mean} intensity of each particle within groups (e.g.
#' per channel per plate), then sets numeric columns to \code{NA} for
#' particles whose intensity falls below the threshold.  This removes
#' background / out-of-focus particles before aggregation.
#'
#' @param df Data frame with unaggregated particle data.  Must contain at
#'   least a \code{Mean} column and the columns named in \code{group_vars}.
#' @param method Scaling strategy: \code{"zscore"} (default) standardises
#'   \code{Mean} (center=TRUE, scale=TRUE) per group and rejects particles
#'   with scaled value below \code{threshold} (default \code{0}, i.e. below
#'   the group mean); \code{"median_ratio"} rejects particles whose
#'   \code{Mean} is below \code{median(Mean) * threshold} per group
#'   (default threshold \code{1/3}).
#' @param threshold Cut-off value. \code{NULL} applies the method's default.
#' @param group_vars Character vector of column names used to define groups
#'   for scaling.  Default: \code{c("Channel_Name", "Plate_ID")}.
#' @param num_cols Numeric columns to set to \code{NA} for rejected
#'   particles.  Default: \code{c("Area", "Mean", "IntDen",
#'   "Number_of_Particles")}.
#'
#' @return The input data frame with numeric columns of rejected particles
#'   set to \code{NA}.  A \code{Mean_scaled} column is appended.
#'
#' @examples
#' \dontrun{
#' df_filtered <- filterParticles(df_raw, method = "zscore", threshold = 0)
#' }
#'
#' @importFrom dplyr group_by mutate ungroup across any_of
#' @export
filterParticles <- function(df,
                            method     = c("zscore", "median_ratio"),
                            threshold  = NULL,
                            group_vars = c("Channel_Name", "Plate_ID"),
                            num_cols   = c("Area", "Mean", "IntDen",
                                           "Number_of_Particles")) {

  method <- match.arg(method)

  if (method == "zscore") {
    if (is.null(threshold)) threshold <- 0
    df <- df %>%
      dplyr::group_by(dplyr::across(dplyr::any_of(group_vars))) %>%
      dplyr::mutate(
        Mean_scaled = as.numeric(scale(.data$Mean, center = TRUE, scale = TRUE))
      ) %>%
      dplyr::ungroup()
    reject <- !is.na(df$Mean_scaled) & df$Mean_scaled < threshold
  } else {
    # median_ratio: uncentered SD scaling — Mean / sd(Mean) per group.
    # Default threshold is median(Mean_scaled) / 3 computed per group.
    df <- df %>%
      dplyr::group_by(dplyr::across(dplyr::any_of(group_vars))) %>%
      dplyr::mutate(
        Mean_scaled = as.numeric(scale(.data$Mean, center = FALSE, scale = TRUE))
      ) %>%
      dplyr::ungroup()
    if (is.null(threshold)) {
      # Dynamic per-group threshold: median of scaled values / 3
      thr_df <- df %>%
        dplyr::group_by(dplyr::across(dplyr::any_of(group_vars))) %>%
        dplyr::mutate(.thr = median(.data$Mean_scaled, na.rm = TRUE) / 3) %>%
        dplyr::ungroup()
      reject <- !is.na(df$Mean_scaled) & df$Mean_scaled < thr_df$.thr
    } else {
      reject <- !is.na(df$Mean_scaled) & df$Mean_scaled < threshold
    }
  }

  # Blank out numeric columns for rejected particles
  cols_to_na <- intersect(num_cols, names(df))
  df[reject, cols_to_na] <- NA

  df
}


#' Aggregate filtered particle data per well and channel
#'
#' Groups particle-level rows and computes per-well summary statistics:
#' mean intensity, mean particle area, normalised area (occupancy), and
#' particle count.  Rows where \code{Mean} is \code{NA} (rejected by
#' \code{\link{filterParticles}}) are excluded from all summaries.
#'
#' @param df Data frame (optionally pre-filtered by \code{filterParticles}).
#'   Should contain \code{Mean}, \code{Area}, and \code{Selection_Area}
#'   columns, plus all columns listed in \code{group_vars}.  If a
#'   \code{CorrSel} column is present, only \code{"Hole_ROI"} rows are used.
#' @param group_vars Columns defining the aggregation groups.
#'   Default: \code{c("Channel_Name", "Plate_ID", "Well")}.
#'
#' @return A data frame with one row per group containing:
#'   \describe{
#'     \item{Mean_agg}{Mean of particle \code{Mean} intensities.}
#'     \item{Area_agg}{Mean of particle \code{Area} values.}
#'     \item{normArea}{Sum of particle areas divided by mean selection area
#'       (occupancy fraction, 0--1).}
#'     \item{n_particles}{Number of non-rejected particles.}
#'   }
#'
#' @examples
#' \dontrun{
#' agg <- aggregateParticles(df_filtered)
#' }
#'
#' @importFrom dplyr filter group_by summarise across any_of distinct left_join
#' @export
aggregateParticles <- function(df,
                               group_vars   = c("Channel_Name", "Plate_ID", "Well"),
                               agg_fun      = c("mean", "median", "sum"),
                               scale_within = "Channel_Name",
                               scale_center = FALSE) {

  agg_fun <- match.arg(agg_fun)
  .agg <- switch(agg_fun,
    mean   = function(x) mean(x,   na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE),
    sum    = function(x) sum(x,    na.rm = TRUE)
  )

  # Restrict to hole ROI if available
  if ("CorrSel" %in% names(df)) {
    df <- dplyr::filter(df, .data$CorrSel == "Hole_ROI")
  }

  # Snapshot all group combinations before dropping NAs — we need to keep
  # wells where every particle was filtered out so they score 0 -> "Negative"
  all_groups <- dplyr::distinct(df, dplyr::across(dplyr::any_of(group_vars)))

  agg <- df %>%
    dplyr::filter(!is.na(.data$Mean)) %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(group_vars))) %>%
    dplyr::summarise(
      Mean_agg    = .agg(.data$Mean),
      Area_agg    = .agg(.data$Area),
      normArea    = sum(.data$Area,             na.rm = TRUE) /
                    (mean(.data$Selection_Area, na.rm = TRUE) + .Machine$double.eps),
      n_particles = dplyr::n(),
      .groups = "drop"
    )

  # Re-join to restore groups that lost all particles
  out <- dplyr::left_join(all_groups, agg,
                          by = intersect(group_vars, names(all_groups)))

  # Wells where every particle was filtered: treat as genuinely absent (0),
  # not missing data.  This matches the original pipeline behaviour and
  # ensures those wells participate in per-channel SD computation as zeros,
  # so their final channel_score is 0 -> classified as Negative.
  stat_cols <- intersect(c("Mean_agg", "Area_agg", "normArea", "n_particles"),
                         names(out))
  for (col in stat_cols) out[[col]][is.na(out[[col]])] <- 0

  # Optional: scale aggregated metrics within scale_within groups
  # Produces Mean_z, Area_z, normArea_z alongside the raw aggregated values.
  if (!is.null(scale_within) && length(scale_within) > 0) {
    .scale_fn <- function(x) {
      x <- as.numeric(x)
      if (isTRUE(scale_center)) {
        m <- mean(x, na.rm = TRUE)
        s <- sd(x,   na.rm = TRUE)
        if (is.na(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
        (x - m) / (s + .Machine$double.eps)
      } else {
        s <- sd(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
        x / (s + .Machine$double.eps)
      }
    }
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::any_of(scale_within))) %>%
      dplyr::mutate(
        Mean_z     = .scale_fn(.data$Mean_agg),
        Area_z     = .scale_fn(.data$Area_agg),
        normArea_z = .scale_fn(.data$normArea)
      ) %>%
      dplyr::ungroup()
  }

  out
}


#' Compute soft classification scores for aggregated well-channel data
#'
#' Rescales each metric (uncentred) \strong{within each channel group}
#' separately, then combines them into a single \code{channel_score} using
#' the supplied weights.  Scoring each channel against its own distribution
#' ensures that a dim channel and a bright channel are both scored on a
#' 0–N scale relative to their own well populations, so cross-channel
#' comparisons inside \code{\link{classifyWells}} remain meaningful.
#' \code{NA}s (wells where all particles were filtered) are set to \code{0}
#' before combining, so those wells score 0 and land in \code{"Negative"}.
#'
#' @param agg_df Aggregated data frame returned by
#'   \code{\link{aggregateParticles}}.
#' @param weights Named numeric vector of relative weights for the three
#'   metrics \code{"Mean"}, \code{"Area"}, and \code{"normArea"}.
#'   Values are normalised internally so they sum to 1.
#'   Default: equal weights \code{c(Mean = 1, Area = 1, normArea = 1)}.
#' @param score_group_vars Columns that define independent scaling groups.
#'   Each channel (and plate, if present) is scaled against its own
#'   distribution of wells.
#'   Default: \code{c("Channel_Name", "Plate_ID")}.
#'
#' @return \code{agg_df} with additional columns \code{Mean_z},
#'   \code{Area_z}, \code{normArea_z}, and \code{channel_score}.
#'
#' @examples
#' \dontrun{
#' scored <- scoreParticles(agg_df, weights = c(Mean = 2, Area = 1, normArea = 1))
#' }
#'
#' @importFrom dplyr group_by mutate ungroup across any_of
#' @export
scoreParticles <- function(agg_df,
                           weights          = c(Mean = 1, Area = 1, normArea = 1),
                           score_group_vars = c("Channel_Name", "Plate_ID"),
                           center           = FALSE) {

  if (isTRUE(center)) {
    warning(paste(
      "center=TRUE uses z-score scaling which can assign near-average scores to",
      "empty wells (those set to 0 by aggregateParticles). Uncentered scoring",
      "(center=FALSE, the default) is recommended for most imaging workflows."
    ), call. = FALSE)
  }

  # Normalise weights to sum to 1
  w <- weights[c("Mean", "Area", "normArea")]
  w[is.na(w)] <- 0
  if (sum(w) == 0) stop("weights must have at least one positive value")
  w <- w / sum(w)

  # Scaling function: uncentred (÷ sd) or centred (z-score) within a group
  .scale_fn <- if (isTRUE(center)) {
    function(x) {
      x <- as.numeric(x)
      m <- mean(x, na.rm = TRUE)
      s <- sd(x,   na.rm = TRUE)
      if (is.na(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
      (x - m) / (s + .Machine$double.eps)
    }
  } else {
    function(x) {
      x <- as.numeric(x)
      s <- sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
      x / (s + .Machine$double.eps)
    }
  }

  # If _z columns are already present (computed by aggregateParticles),
  # only recompute when score_group_vars differs from what was used there.
  # As a simple heuristic: recompute when score_group_vars has >1 element
  # (i.e. Channel × Plate), skip when they are already available.
  has_z <- all(c("Mean_z", "Area_z", "normArea_z") %in% names(agg_df))
  if (!has_z) {
    agg_df <- agg_df %>%
      dplyr::group_by(dplyr::across(dplyr::any_of(score_group_vars))) %>%
      dplyr::mutate(
        Mean_z     = .scale_fn(.data$Mean_agg),
        Area_z     = .scale_fn(.data$Area_agg),
        normArea_z = .scale_fn(.data$normArea)
      ) %>%
      dplyr::ungroup()
  }

  # Replace NA with 0 before combining
  mz <- ifelse(is.na(agg_df$Mean_z),     0, agg_df$Mean_z)
  az <- ifelse(is.na(agg_df$Area_z),     0, agg_df$Area_z)
  nz <- ifelse(is.na(agg_df$normArea_z), 0, agg_df$normArea_z)

  agg_df$channel_score <- w["Mean"] * mz + w["Area"] * az + w["normArea"] * nz

  agg_df
}


#' Classify wells as positive or negative per imaging channel
#'
#' Pivots the scored data to wide format (one row per well-plate), computes
#' a per-well maximum score across channels, and marks each channel as
#' positive if its score is within \code{delta} of the maximum \emph{and}
#' its occupancy (\code{normArea}) exceeds \code{min_area}.  A combined
#' \code{Classification} string is assembled from the positive channel
#' labels (e.g. \code{"GFP+ mScarlett+"}, \code{"Negative"}).
#'
#' @param score_df Scored data frame from \code{\link{scoreParticles}}.
#'   Must contain \code{Well}, \code{Plate_ID}, \code{Channel_Name},
#'   \code{channel_score}, and \code{normArea} columns.
#' @param delta Dominance tolerance: a channel is called positive only if
#'   \code{channel_score >= max_score - delta}.  Default: \code{0.5}.
#' @param min_area Minimum occupancy (normArea) required for a positive
#'   call.  Default: \code{0.1} (10\%).
#' @param channel_labels Optional named character vector mapping
#'   \code{Channel_Name} values to display labels used in the
#'   \code{Classification} string, e.g.
#'   \code{c(C1 = "GFP", C2 = "mScarlett")}.  If \code{NULL} (default),
#'   the raw \code{Channel_Name} values are used.
#'
#' @return A data frame with one row per well-plate containing:
#'   \describe{
#'     \item{Well, Plate_ID}{Well and plate identifiers.}
#'     \item{<channel>_score}{Per-channel composite score.}
#'     \item{<channel>_normArea}{Per-channel occupancy fraction.}
#'     \item{<channel>_Mean_z, <channel>_Area_z, <channel>_normArea_z}{
#'       Per-channel z-scored metrics (included automatically when present in
#'       \code{score_df}, i.e. when called after \code{\link{scoreParticles}}).}
#'     \item{<channel>_positive}{Logical positivity flag per channel.}
#'     \item{max_score}{Maximum score across all channels for that well.}
#'     \item{Classification}{Combined label string.}
#'   }
#'
#' @examples
#' \dontrun{
#' cls <- classifyWells(scored,
#'                      delta          = 0.5,
#'                      min_area       = 0.1,
#'                      channel_labels = c(C1 = "GFP", C2 = "mScarlett"))
#' }
#'
#' @importFrom dplyr select left_join all_of across any_of
#' @importFrom tidyr pivot_wider
#' @export
classifyWells <- function(score_df,
                          delta          = 0.5,
                          min_area       = 0.1,
                          channel_labels = NULL) {

  # Apply optional channel label remapping
  if (!is.null(channel_labels)) {
    score_df$Channel_Name <- ifelse(
      score_df$Channel_Name %in% names(channel_labels),
      channel_labels[score_df$Channel_Name],
      score_df$Channel_Name
    )
  }

  channels <- unique(score_df$Channel_Name)
  id_vars  <- intersect(c("Well", "Plate_ID"), names(score_df))

  # Pivot scores and normAreas to wide
  wide_score <- tidyr::pivot_wider(
    dplyr::select(score_df,
                  dplyr::all_of(c(id_vars, "Channel_Name", "channel_score"))),
    names_from  = "Channel_Name",
    values_from = "channel_score",
    names_glue  = "{Channel_Name}_score"
  )

  wide_area <- tidyr::pivot_wider(
    dplyr::select(score_df,
                  dplyr::all_of(c(id_vars, "Channel_Name", "normArea"))),
    names_from  = "Channel_Name",
    values_from = "normArea",
    names_glue  = "{Channel_Name}_normArea"
  )

  wide <- dplyr::left_join(wide_score, wide_area, by = id_vars)

  # Pivot z-scored metrics if present (produced by scoreParticles / aggregateParticles)
  z_cols <- intersect(c("Mean_z", "Area_z", "normArea_z"), names(score_df))
  if (length(z_cols) > 0) {
    wide_z <- tidyr::pivot_wider(
      dplyr::select(score_df,
                    dplyr::all_of(c(id_vars, "Channel_Name", z_cols))),
      names_from  = "Channel_Name",
      values_from = dplyr::all_of(z_cols),
      names_glue  = "{Channel_Name}_{.value}"
    )
    wide <- dplyr::left_join(wide, wide_z, by = id_vars)
  }

  # Per-well maximum score across all channels
  score_cols <- intersect(paste0(channels, "_score"), names(wide))
  wide$max_score <- do.call(
    pmax,
    c(lapply(wide[, score_cols, drop = FALSE],
             function(x) ifelse(is.na(x), -Inf, x)),
      list(na.rm = TRUE))
  )

  # Positivity flags
  for (ch in channels) {
    sc_col  <- paste0(ch, "_score")
    na_col  <- paste0(ch, "_normArea")
    pos_col <- paste0(ch, "_positive")
    if (!sc_col %in% names(wide)) next
    sc <- ifelse(is.na(wide[[sc_col]]), -Inf, wide[[sc_col]])
    na <- ifelse(is.na(wide[[na_col]]),  0,   wide[[na_col]])
    wide[[pos_col]] <- sc >= (wide$max_score - delta) & na > min_area
  }

  # Assemble classification string
  pos_cols <- intersect(paste0(channels, "_positive"), names(wide))

  wide$Classification <- apply(
    wide[, pos_cols, drop = FALSE], 1,
    function(row) {
      pos_ch <- channels[which(as.logical(row))]
      if (length(pos_ch) == 0) return("Negative")
      if (length(pos_ch) == length(channels) && length(channels) > 2)
        return("Multiple+")
      paste(paste0(pos_ch, "+"), collapse = " ")
    }
  )

  wide
}


#' Full particle-level cell classification pipeline
#'
#' Convenience wrapper that runs
#' \code{\link{filterParticles}} ->
#' \code{\link{aggregateParticles}} ->
#' \code{\link{scoreParticles}} ->
#' \code{\link{classifyWells}}
#' in sequence.
#'
#' @param df Unaggregated particle data frame (one row per particle per
#'   ROI).  Typically the output of \code{\link{prepareImgDF}} with
#'   \code{aggregate = FALSE} and \code{cleanNames = FALSE}, filtered to
#'   \code{Image_Type == "fluor"}.
#' @param filter_method Passed to \code{\link{filterParticles}}.
#'   \code{"zscore"} (default) or \code{"median_ratio"}.
#' @param filter_threshold Intensity cut-off passed to
#'   \code{\link{filterParticles}}.  \code{NULL} uses the method default.
#' @param filter_group_vars Grouping columns for intensity scaling.
#'   Default: \code{c("Channel_Name", "Plate_ID")}.
#' @param agg_group_vars Grouping columns for aggregation.
#'   Default: \code{c("Channel_Name", "Plate_ID", "Well")}.
#' @param weights Named weight vector for scoring.
#'   Default: equal weights for \code{Mean}, \code{Area}, \code{normArea}.
#' @param score_group_vars Columns used as independent scaling groups in
#'   \code{\link{scoreParticles}}.
#'   Default: \code{c("Channel_Name", "Plate_ID")}.
#' @param delta Dominance tolerance for \code{\link{classifyWells}}.
#'   Default: \code{0.5}.
#' @param min_area Minimum occupancy threshold.  Default: \code{0.1}.
#' @param channel_labels Optional named vector mapping raw channel names to
#'   display labels, e.g. \code{c(C1 = "GFP", C2 = "mScarlett")}.
#' @param return_scores Controls what is returned.  \code{FALSE} (default)
#'   returns the wide-format classification table.  \code{TRUE} returns the
#'   long-format scored data frame (output of \code{\link{scoreParticles}}).
#'   \code{"both"} returns a named list with elements \code{classification}
#'   (wide) and \code{scores} (long), which can be passed directly to
#'   \code{\link{mergeClassificationToSE}} to also embed z-scored metrics.
#'
#' @return A wide-format classification data frame; the long-format scored
#'   data frame if \code{return_scores = TRUE}; or a named list
#'   \code{list(classification, scores)} if \code{return_scores = "both"}.
#'
#' @examples
#' \dontrun{
#' # Standard usage
#' cls <- classifyImgParticles(
#'   df_raw,
#'   filter_method  = "zscore",
#'   delta          = 0.5,
#'   min_area       = 0.1,
#'   channel_labels = c(C1 = "GFP", C2 = "mScarlett")
#' )
#' se <- mergeClassificationToSE(se, cls)
#'
#' # Also embed per-channel z-scored metrics into colData
#' res <- classifyImgParticles(df_raw, return_scores = "both")
#' se  <- mergeClassificationToSE(se, res$classification, scores = res$scores)
#' }
#'
#' @export
classifyImgParticles <- function(df,
                                 filter_method     = c("zscore", "median_ratio"),
                                 filter_threshold  = NULL,
                                 filter_group_vars = c("Channel_Name", "Plate_ID"),
                                 agg_group_vars    = c("Channel_Name", "Plate_ID", "Well"),
                                 agg_fun           = c("mean", "median", "sum"),
                                 scale_within      = "Channel_Name",
                                 scale_center      = FALSE,
                                 weights           = c(Mean = 1, Area = 1, normArea = 1),
                                 score_group_vars  = c("Channel_Name", "Plate_ID"),
                                 center            = FALSE,
                                 delta             = 0.5,
                                 min_area          = 0.1,
                                 channel_labels    = NULL,
                                 return_scores     = FALSE) {

  filter_method <- match.arg(filter_method)

  df_filt   <- filterParticles(df,
                                method     = filter_method,
                                threshold  = filter_threshold,
                                group_vars = filter_group_vars)

  df_agg    <- aggregateParticles(df_filt, group_vars   = agg_group_vars,
                                  agg_fun      = agg_fun,
                                  scale_within = scale_within,
                                  scale_center = scale_center)

  df_scored <- scoreParticles(df_agg, weights = weights,
                               score_group_vars = score_group_vars,
                               center = center)

  if (isTRUE(return_scores)) return(df_scored)

  cls <- classifyWells(df_scored,
                       delta          = delta,
                       min_area       = min_area,
                       channel_labels = channel_labels)

  if (identical(return_scores, "both"))
    return(list(classification = cls, scores = df_scored))

  cls
}


#' Merge image classification results into a SummarizedExperiment
#'
#' Joins the output of \code{\link{classifyImgParticles}} (or
#' \code{\link{classifyWells}}) to the \code{colData} of a
#' \code{SummarizedExperiment} by \code{Well} and \code{Plate_ID}.
#' Optionally also embeds per-channel z-scored metrics from the long-format
#' \code{\link{scoreParticles}} output.
#'
#' @param se A \code{SummarizedExperiment} whose \code{colData} contains
#'   \code{Well} and \code{Plate_ID} columns.
#' @param classification Data frame returned by
#'   \code{\link{classifyImgParticles}} or \code{\link{classifyWells}}.
#'   Must contain at least \code{Well}, \code{Plate_ID}, and
#'   \code{Classification}.
#' @param col_name Name of the colData column to write results into.  Use
#'   \code{NULL} (default) to add all new columns flat into \code{colData}.
#'   If a string is supplied, all new columns are stored as a nested
#'   \code{DataFrame} under that name.
#' @param scores Optional long-format data frame returned by
#'   \code{\link{scoreParticles}} (or the \code{$scores} element of
#'   \code{classifyImgParticles(..., return_scores = "both")}).  When
#'   supplied, the per-channel z-scored metrics (\code{Mean_z},
#'   \code{Area_z}, \code{normArea_z}) are pivoted to wide format and merged
#'   alongside the classification columns, producing columns such as
#'   \code{GFP_Mean_z}, \code{GFP_Area_z}, \code{GFP_normArea_z}.
#'
#' @return The updated \code{se} object.
#'
#' @examples
#' \dontrun{
#' # Classification only
#' se <- mergeClassificationToSE(se, cls)
#' colData(se)$Classification
#'
#' # Classification + per-channel z-scored metrics
#' res <- classifyImgParticles(df_raw, return_scores = "both")
#' se  <- mergeClassificationToSE(se, res$classification, scores = res$scores)
#' colData(se)$GFP_Mean_z
#'
#' # Store everything as a nested DataFrame
#' se <- mergeClassificationToSE(se, cls, col_name = "img_classification")
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join select all_of across any_of
#' @importFrom tidyr pivot_wider
#' @export
mergeClassificationToSE <- function(se, classification, col_name = NULL,
                                    scores = NULL) {

  cd     <- as.data.frame(SummarizedExperiment::colData(se))
  joined <- dplyr::left_join(cd, classification, by = c("Well", "Plate_ID"))

  # Optionally pivot z-scored metrics from the long-format scored data frame
  if (!is.null(scores)) {
    z_cols  <- intersect(c("Mean_z", "Area_z", "normArea_z"), names(scores))
    id_vars <- intersect(c("Well", "Plate_ID"), names(scores))
    if (length(z_cols) > 0 && "Channel_Name" %in% names(scores)) {
      scores_wide <- tidyr::pivot_wider(
        dplyr::select(scores,
                      dplyr::all_of(c(id_vars, "Channel_Name", z_cols))),
        names_from  = "Channel_Name",
        values_from = dplyr::all_of(z_cols),
        names_glue  = "{Channel_Name}_{.value}"
      )
      joined <- dplyr::left_join(joined, scores_wide, by = id_vars)
    }
  }

  new_cols <- setdiff(names(joined), names(cd))

  if (is.null(col_name)) {
    for (col in new_cols) {
      SummarizedExperiment::colData(se)[[col]] <- joined[[col]]
    }
  } else {
    SummarizedExperiment::colData(se)[[col_name]] <-
      S4Vectors::DataFrame(joined[, new_cols, drop = FALSE])
  }

  se
}
