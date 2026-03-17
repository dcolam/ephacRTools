#' @importFrom magrittr %>%
NULL

# ---------------------------------------------------------------------------
# previewImgDB — scan databases before import to discover column structure
# ---------------------------------------------------------------------------

#' Scan Cluster Analysis databases to discover available columns and values
#'
#' Queries one or more SQLite databases to return the metadata needed to
#' configure the import (column names, distinct Selection values with counts,
#' Image_ID rank distribution per well, available measurement columns, etc.).
#' This is used by the \code{tinySEV} Import Data UI to populate dynamic
#' configuration inputs before the actual import.
#'
#' @param paths Character vector of paths to \code{.db} SQLite files.
#' @param analysis \code{"pa"} (Particle Analysis, default) or
#'   \code{"coloc"} (Colocalization).
#'
#' @return A named list:
#'   \describe{
#'     \item{pa_cols}{Character vector — all columns in the analysis table
#'       (union across all DBs, excluding \code{PA_ID}).}
#'     \item{meas_cols}{Character vector — columns available in the
#'       measurement table (excluding \code{PA_ID} and \code{Label}).}
#'     \item{selections}{Named integer vector — Selection value → total row
#'       count, sorted descending.}
#'     \item{channels}{Character vector — distinct Channel_Name values.}
#'     \item{n_rows}{Integer — total particle rows across all DBs.}
#'     \item{n_wells}{Integer — approximate number of unique wells.}
#'     \item{image_ranks}{data.frame with columns \code{rank} and \code{n} —
#'       how many Image_IDs appear at each rank position within a well
#'       (rank 1 = lowest Image_ID per well, typically fluorescence).}
#'   }
#'
#' @export
previewImgDB <- function(paths, analysis = "pa") {
  pa_tbl   <- if (analysis == "pa") "Particle_Analysis_Table" else "Coloc_Analysis_Table"
  meas_tbl <- if (analysis == "pa") "PA_Measurement_Tables"   else "Coloc_Measurement_Tables"

  all_pa_cols   <- character(0)
  all_meas_cols <- character(0)
  sel_counts    <- list()
  channels      <- character(0)
  total_rows    <- 0L
  all_wells     <- character(0)
  rank_counts   <- list()

  per_db <- list()

  for (p in paths) {
    con <- DBI::dbConnect(RSQLite::SQLite(), p)
    db_rank_counts <- list()
    tryCatch({
      all_pa_cols   <- union(all_pa_cols,
                             setdiff(DBI::dbListFields(con, pa_tbl),   "PA_ID"))
      all_meas_cols <- union(all_meas_cols,
                             setdiff(DBI::dbListFields(con, meas_tbl), c("PA_ID", "Label")))

      # Selection value counts (global + per-DB)
      sel <- DBI::dbGetQuery(con,
        paste0("SELECT Selection, COUNT(*) AS n FROM ", pa_tbl,
               " GROUP BY Selection ORDER BY n DESC LIMIT 200"))
      db_sel <- setNames(as.integer(sel$n), sel$Selection)
      for (i in seq_len(nrow(sel))) {
        k <- sel$Selection[i]
        sel_counts[[k]] <- (sel_counts[[k]] %||% 0L) + as.integer(sel$n[i])
      }

      # Channels
      ch <- DBI::dbGetQuery(con,
        paste0("SELECT DISTINCT Channel_Name FROM ", pa_tbl))$Channel_Name
      channels <- union(channels, ch)

      # Row count
      n <- DBI::dbGetQuery(con,
        paste0("SELECT COUNT(*) AS n FROM ", pa_tbl))$n
      total_rows <- total_rows + as.integer(n)

      # Wells
      w <- DBI::dbGetQuery(con,
        paste0("SELECT DISTINCT Well FROM ", pa_tbl))$Well
      all_wells <- union(all_wells, w)

      # Image_ID rank per well (rank 1 = smallest Image_ID within a well)
      img <- DBI::dbGetQuery(con,
        paste0("SELECT Well, Image_ID FROM ", pa_tbl,
               " GROUP BY Well, Image_ID ORDER BY Well, CAST(Image_ID AS REAL)"))
      if (nrow(img) > 0) {
        img$rank <- ave(as.numeric(img$Image_ID), img$Well,
                        FUN = function(x) rank(x, ties.method = "first"))
        for (r in unique(img$rank)) {
          k <- as.character(as.integer(r))
          rank_counts[[k]]    <- (rank_counts[[k]]    %||% 0L) + sum(img$rank == r)
          db_rank_counts[[k]] <- (db_rank_counts[[k]] %||% 0L) + sum(img$rank == r)
        }
      }

      # Per-DB record
      db_rk <- as.integer(names(db_rank_counts))
      per_db[[length(per_db) + 1]] <- list(
        name        = basename(p),
        n_rows      = as.integer(n),
        n_wells     = length(w),
        image_ranks = if (length(db_rank_counts) > 0)
          data.frame(rank = db_rk[order(db_rk)],
                     n    = unlist(db_rank_counts)[order(db_rk)],
                     stringsAsFactors = FALSE)
        else data.frame(rank = integer(0), n = integer(0)),
        top_selections = head(db_sel, 10)
      )
    }, error = function(e) warning("previewImgDB: error reading ", p, ": ", e$message))
    DBI::dbDisconnect(con)
  }

  sel_sorted <- sort(unlist(sel_counts), decreasing = TRUE)

  rank_keys <- as.integer(names(rank_counts))
  rank_df   <- if (length(rank_counts) > 0)
    data.frame(rank = rank_keys[order(rank_keys)],
               n    = unlist(rank_counts)[order(rank_keys)],
               stringsAsFactors = FALSE)
  else data.frame(rank = integer(0), n = integer(0))

  list(
    pa_cols     = all_pa_cols,
    meas_cols   = all_meas_cols,
    selections  = sel_sorted,
    channels    = channels,
    n_rows      = total_rows,
    n_wells     = length(all_wells),
    image_ranks = rank_df,
    per_db      = per_db        # list: one entry per uploaded DB
  )
}

# ---------------------------------------------------------------------------
# prepareSingleImgDF — import a single .db file
# ---------------------------------------------------------------------------

#' Prepare imaging-results table from a single Cluster Analysis SQLite database
#'
#' Joins the analysis table with its measurement table, normalises well labels,
#' computes per-well Image_ID rank (rank 1 = lowest Image_ID per well, typically
#' the fluorescence image), and optionally filters to specified selection values.
#'
#' @param pathDB Path to a single SQLite \code{.db} file.
#' @param analysis \code{"pa"} (Particle Analysis, default) or \code{"coloc"}.
#' @param plate_col Name of the column containing the plate identifier.
#'   Default: \code{"Plate_ID"}.
#' @param well_col Name of the column containing the well identifier.
#'   Default: \code{"Well"}.
#' @param channel_col Name of the column containing the channel name.
#'   Default: \code{"Channel_Name"}.
#' @param selection_col Name of the column containing the selection/ROI name.
#'   Default: \code{"Selection"}.
#' @param extra_id_cols Additional columns from the analysis table to retain
#'   (e.g. \code{c("Species", "Date")}). Default: none.
#' @param meas_cols Measurement columns to import from the measurement table.
#'   Default: \code{c("Area","Mean","IntDen","Number_of_Particles")}.
#'   \code{"Number_of_Particles"} and \code{"Selection_Area"} are always kept
#'   from the analysis table regardless.
#' @param fluor_rank Integer. Image_ID rank within a well that corresponds to
#'   the fluorescence image (rank 1 = lowest Image_ID, typically fluorescence).
#'   Set to \code{NULL} to skip rank computation and keep all images.
#'   Default: \code{NULL}.
#' @param roi_selections Character vector of \code{Selection} values to keep.
#'   \code{NULL} retains all rows. Default: \code{NULL}.
#' @param scale_num Logical; add scaled versions of numeric columns. Default: \code{FALSE}.
#' @param scale_cols Columns to scale (defaults to \code{meas_cols} if \code{NULL}).
#' @param scale_fun Scaling function. Default: z-score (\code{scale(x,TRUE,TRUE)}).
#' @param aggregate Logical; if \code{TRUE}, average rows within
#'   Channel × Well × Plate × Selection groups. Default: \code{FALSE}.
#' @param cleanNames Logical; normalise well labels (strip site suffix, pad
#'   column digits). Default: \code{TRUE}.
#'
#' @return A data.frame with one row per particle (or per group if
#'   \code{aggregate = TRUE}).
#'
#' @seealso \code{\link{previewImgDB}}, \code{\link{prepareImgDF}}
#' @export
prepareSingleImgDF <- function(pathDB,
                               analysis      = c("pa", "coloc"),
                               plate_col     = "Plate_ID",
                               well_col      = "Well",
                               channel_col   = "Channel_Name",
                               selection_col = "Selection",
                               extra_id_cols = character(0),
                               meas_cols     = c("Area", "Mean", "IntDen",
                                                 "Number_of_Particles"),
                               fluor_rank     = NULL,
                               roi_selections = NULL,
                               apply_corr_sel = TRUE,
                               corr_sel_filter = NULL,
                               scale_num     = FALSE,
                               scale_cols    = NULL,
                               scale_fun     = function(x)
                                 as.numeric(scale(x, TRUE, TRUE)),
                               aggregate     = FALSE,
                               cleanNames    = TRUE) {

  analysis <- match.arg(analysis)
  pa_tbl   <- if (analysis == "pa") "Particle_Analysis_Table" else "Coloc_Analysis_Table"
  meas_tbl <- if (analysis == "pa") "PA_Measurement_Tables"   else "Coloc_Measurement_Tables"
  id_col   <- if (analysis == "pa") "PA_ID"                   else "COLOC_ID"

  con <- DBI::dbConnect(RSQLite::SQLite(), pathDB)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  tbl <- DBI::dbGetQuery(con, paste0(
    "SELECT * FROM ", pa_tbl, " AS pa",
    " JOIN ", meas_tbl, " AS meas ON meas.", id_col, " = pa.", id_col))

  DBI::dbDisconnect(con)

  # Remove duplicate join key column
  dup <- which(names(tbl) == id_col)
  if (length(dup) > 1) tbl <- tbl[, -dup[-1], drop = FALSE]

  # ── Rename user-specified columns to standard names ──────────────────────
  col_map <- setNames(c(plate_col, well_col, channel_col, selection_col),
                      c("Plate_ID",  "Well",   "Channel_Name", "Selection"))
  for (std in names(col_map)) {
    src <- col_map[[std]]
    if (src != std && src %in% names(tbl))
      names(tbl)[names(tbl) == src] <- std
  }

  # ── Select columns ────────────────────────────────────────────────────────
  always  <- c("Plate_ID", "Well", "Channel_Name", "Selection",
               "Image_ID", "Selection_Area", "Number_of_Particles")
  meas_extra <- setdiff(meas_cols, c("Number_of_Particles", "Selection_Area"))
  keep    <- unique(c(always, extra_id_cols, meas_extra))
  tbl     <- dplyr::select(tbl, tidyselect::any_of(keep))

  # ── Well normalisation ────────────────────────────────────────────────────
  if (isTRUE(cleanNames)) {
    # Strip site suffix: "B10-1" → "B10"
    tbl$Well <- vapply(tbl$Well,
      function(x) strsplit(x, "-")[[1]][1], character(1))
    # Normalise to "A01" format
    row_ltr  <- substr(tbl$Well, 1, 1)
    col_num  <- suppressWarnings(as.integer(substr(tbl$Well, 2, 3)))
    tbl$Well <- paste0(row_ltr, sprintf("%02d", col_num))
    # Plate_ID: strip \r / whitespace
    tbl$Plate_ID <- vapply(tbl$Plate_ID,
      function(x) trimws(strsplit(x, "\r")[[1]][1]), character(1))
  }

  # ── CorrSel labeling: Hole_ROI (min area) / background_ROI (max area) ─────
  # For each unique Selection value, the smallest Selection_Area = Hole_ROI,
  # the largest = background_ROI.  This resolves ambiguity when the same
  # Selection name (e.g. "14618.vsi - 009 BF") appears at two different areas.
  if (isTRUE(apply_corr_sel) &&
      "Selection" %in% names(tbl) &&
      "Selection_Area" %in% names(tbl)) {
    tbl$CorrSel <- NA_character_
    for (sel in unique(tbl$Selection)) {
      idx   <- tbl$Selection == sel
      areas <- tbl$Selection_Area[idx]
      sel_min <- min(areas, na.rm = TRUE)
      sel_max <- max(areas, na.rm = TRUE)
      tbl$CorrSel[idx & tbl$Selection_Area == sel_min] <- "Hole_ROI"
      tbl$CorrSel[idx & tbl$Selection_Area == sel_max] <- "background_ROI"
    }
  }

  # ── Filter by CorrSel ─────────────────────────────────────────────────────
  if (!is.null(corr_sel_filter) && length(corr_sel_filter) > 0 &&
      "CorrSel" %in% names(tbl)) {
    tbl <- tbl[tbl$CorrSel %in% corr_sel_filter, , drop = FALSE]
  }

  # ── Image_ID rank within well ─────────────────────────────────────────────
  if ("Image_ID" %in% names(tbl) && !is.null(fluor_rank)) {
    tbl$Image_ID_num <- suppressWarnings(as.numeric(tbl$Image_ID))
    tbl <- dplyr::group_by(tbl, Plate_ID, Well) %>%
      dplyr::mutate(Image_rank = dplyr::dense_rank(Image_ID_num)) %>%
      dplyr::ungroup() %>%
      as.data.frame()
    tbl$Image_Type <- ifelse(tbl$Image_rank == as.integer(fluor_rank),
                             "fluor", "bf")
  }

  # ── Filter by ROI selections ──────────────────────────────────────────────
  if (!is.null(roi_selections) && length(roi_selections) > 0 &&
      "Selection" %in% names(tbl)) {
    tbl <- tbl[tbl$Selection %in% roi_selections, , drop = FALSE]
  }

  # ── Aggregate (optional) ─────────────────────────────────────────────────
  if (isTRUE(aggregate)) {
    grp <- intersect(c("Plate_ID", "Well", "Channel_Name", "Selection",
                       "Image_Type", "Image_rank", extra_id_cols), names(tbl))
    tbl <- ag(tbl, cols = grp, fun = mean)
  }

  # ── Optional scaling ──────────────────────────────────────────────────────
  if (isTRUE(scale_num)) {
    if (is.null(scale_cols))
      scale_cols <- intersect(meas_extra, names(tbl))
    for (col in scale_cols) {
      if (!col %in% names(tbl)) next
      new_col       <- paste0(col, "_Scaled")
      tbl[[new_col]] <- scale_fun(tbl[[col]])
      tbl <- dplyr::relocate(tbl, dplyr::all_of(new_col),
                              .after = dplyr::all_of(col))
    }
  }

  tbl
}

# ---------------------------------------------------------------------------
# prepareImgDF — wrapper for one or multiple .db files
# ---------------------------------------------------------------------------

#' Prepare imaging-results from one or more Cluster Analysis SQLite databases
#'
#' Calls \code{\link{prepareSingleImgDF}} for each file and row-binds the
#' results.
#'
#' @inheritParams prepareSingleImgDF
#' @param pathDB Character vector of paths to one or more \code{.db} files.
#' @param coloc_cols Extra columns added when \code{analysis = "coloc"}.
#'   Default: \code{c("Second_Channel","Mask_Area")}.
#'
#' @return A combined data.frame. If \code{fluor_rank} is not \code{NULL},
#'   an \code{Image_Type} column (\code{"fluor"} / \code{"bf"}) is present.
#'
#' @seealso \code{\link{previewImgDB}}, \code{\link{prepareSingleImgDF}}
#' @export
prepareImgDF <- function(pathDB,
                         analysis        = "pa",
                         plate_col       = "Plate_ID",
                         well_col        = "Well",
                         channel_col     = "Channel_Name",
                         selection_col   = "Selection",
                         extra_id_cols   = character(0),
                         meas_cols       = c("Area", "Mean", "IntDen",
                                             "Number_of_Particles"),
                         coloc_cols      = c("Second_Channel", "Mask_Area"),
                         fluor_rank      = NULL,
                         fluor_ranks     = NULL,
                         roi_selections  = NULL,
                         apply_corr_sel  = TRUE,
                         corr_sel_filter = NULL,
                         scale_num       = FALSE,
                         scale_cols      = NULL,
                         scale_fun       = function(x) as.numeric(scale(x, TRUE, TRUE)),
                         aggregate       = FALSE,
                         cleanNames      = TRUE) {

  if ("coloc" %in% analysis)
    extra_id_cols <- union(extra_id_cols, coloc_cols)

  # fluor_ranks: named integer vector "1","2",... → per-DB rank override
  # fluor_rank:  global default used when a DB has no entry in fluor_ranks
  call_single <- function(p, idx) {
    this_rank <- if (!is.null(fluor_ranks) &&
                     as.character(idx) %in% names(fluor_ranks))
      fluor_ranks[[as.character(idx)]]
    else
      fluor_rank
    out <- tryCatch(
      prepareSingleImgDF(p,
        analysis        = analysis,
        plate_col       = plate_col,
        well_col        = well_col,
        channel_col     = channel_col,
        selection_col   = selection_col,
        extra_id_cols   = extra_id_cols,
        meas_cols       = meas_cols,
        fluor_rank      = this_rank,
        roi_selections  = roi_selections,
        apply_corr_sel  = apply_corr_sel,
        corr_sel_filter = corr_sel_filter,
        scale_num       = scale_num,
        scale_cols      = scale_cols,
        scale_fun       = scale_fun,
        aggregate       = aggregate,
        cleanNames      = cleanNames),
      error = function(e) {
        warning("prepareSingleImgDF failed for ", p, ": ", e$message); NULL }
    )
    out
  }

  if (length(pathDB) == 1) {
    df <- call_single(pathDB, 1L)
  } else {
    dfs  <- lapply(seq_along(pathDB), function(i) call_single(pathDB[i], i))
    dfs  <- Filter(Negate(is.null), dfs)
    if (length(dfs) == 0) stop("All databases failed to import.")
    names(dfs) <- paste0("DB", seq_along(dfs))
    df <- dplyr::bind_rows(dfs, .id = "column_label")
    if (!"Plate_ID" %in% names(df))
      df$Plate_ID <- df$column_label
  }

  df
}
#' Clean and normalize dataframe. Adds column and row identifiers and finds the
#' selections with the Hole-ROI
#' @param df dataframe of image results
#' @param channels The names of the channels you are using (in the right order)
#' @return A dataframe
#' @export
df_cleaned <- function(df, channels = c("Green", "Red", "ROMK")){

  df$Well_clean <- sapply(df$Well, function(x){

    unlist(stringr::str_split(x, "-"))[1]

  })

  df$Row <- sapply(df$Well_clean, function(x){

    stringr::str_sub(x, 1, 1)

  })

  df$Column <- sapply(df$Well_clean, function(x){

    stringr::str_sub(x, 2, 3)

  })

  df$CorrSel <- NA
  for(sel in unique(df$Selection)){
    selMin <- min(subset(df, Selection == sel)$Selection_Area)
    selMax <- max(subset(df, Selection == sel)$Selection_Area)

    df[df$Selection == sel & df$Selection_Area == selMin,"CorrSel"] <- "Hole_ROI"
    df[df$Selection == sel & df$Selection_Area == selMax,"CorrSel"] <- "background_ROI"
  }
  # Image_Type is no longer derived from Image_ID parity (the odd/even pattern
  # is not reliable across datasets).  Use prepareSingleImgDF(fluor_rank=) for
  # correct assignment; this column is left as-is here for backward compat.
  if (!"Image_Type" %in% names(df)) df$Image_Type <- NA_character_

  df$Channel <- ifelse(df$Channel_Name == "BFP", channels[1],
                        ifelse(df$Channel_Name == "mCherry", channels[2],
                               ifelse(df$Channel_Name == "GFP", channels[3], NA)))
  df$Well <- paste(df$Row, stringr::str_pad(df$Column, 2, pad = "0"), sep="")

  return(df)
}
#' Merge together the imaging-results into the Column Data of the SE
#' @param se SummarizedExperiment Object with the Ephys-Data
#' @param df_img DataFrame with imaging results returned by prepareImgDF()
#' @param tableType Indicate whether the merge should be using particle analysis data or colocalization
#' @param Selection Indicate whether the focus should be on the hole or the background
#' @param suffix Indicate how the new columns should be identified at the end of the name
#' @return A dataframe
#' @export
mergeSEandImg <- function(se, df_img, tableType = "pa", selType = c("Hole_ROI", "background_ROI"), suffix = "hole"){
  # if (Selection == "Hole_ROI"){
  #   df_img <- subset(df_img, Image_Type == "fluor" & CorrSel == "Hole_ROI")
  # }
  # if (Selection == "background_ROI"){
  #   df_img <- subset(df_img, Image_Type == "fluor" & CorrSel == "background_ROI")
  # }

  df_img <- subset(df_img, Image_Type == "fluor" & CorrSel %in% selType)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  # Loop through each channel
  if(tableType == "pa"){
    channels <- unique(df_img$Channel_Name)
    for (channel in channels) {
      # Subset df_img for current channel
      df_channel <- df_img %>%
        dplyr::filter(Channel_Name == channel) %>%
        dplyr::select(-Channel_Name)  # optional: remove the channel label
      # Perform join
      joined <- cd %>%
        dplyr::left_join(df_channel, by = c("Well", "Plate_ID"))
      # Extract just the new columns (everything except original colData)
      new_cols <- dplyr::setdiff(names(joined), names(cd))
      # Create a DataFrame object from just the new data
      channel_data <- S4Vectors::DataFrame(joined[, new_cols])
      SummarizedExperiment::colData(se)[[paste(channel, suffix, sep=".")]] <- channel_data
    }

  }else{
    channels <- unique(df_img$Channel_Name)

    for (channel in channels) {
      second_channels <- unique(subset(df_img, Channel_Name == channel)$Second_Channel)
      for (second_channel in second_channels){
      # Subset df_img for current channel
      df_channel <- df_img %>%
        dplyr::filter(Channel_Name == channel, Second_Channel == second_channel) %>%
        dplyr::select(-Channel_Name, -Second_Channel)  # optional: remove the channel label
      # Perform join
      joined <- cd %>%
        dplyr::left_join(df_channel, by = c("Well", "Plate_ID"))
      # Extract just the new columns (everything except original colData)
      new_cols <- dplyr::setdiff(names(joined), names(cd))
      # Create a DataFrame object from just the new data
      channel_data <- S4Vectors::DataFrame(joined[, new_cols])
      # Assign to colData(se), one column per channel, as a nested DataFrame
      SummarizedExperiment::colData(se)[[paste(channel, second_channel, sep=".")]] <- channel_data
      }
    }

  }
  return(se)
}

#' Add JPG thumbnail file paths to a SummarizedExperiment colData
#'
#' Scans \code{folder} recursively for JPG thumbnails whose filenames follow
#' the Cluster Analysis export convention:
#' \preformatted{..._<PlateID>_<Well>-<site>_<Channel>_<Class>_crop.jpg}
#' (e.g. \code{Exp_18T39265_A01-1_BF_class1_crop.jpg}).
#'
#' Files are matched to SE wells by Plate_ID (via sub-directory names that
#' contain the plate ID as a substring) and well label.  The resolved paths
#' are stored as a nested \code{\link[S4Vectors]{DataFrame}} inside a single
#' colData column (default: \code{"thumbnails"}).  The nested DataFrame has
#' one row per colData row (well \eqn{\times} plate) and one column per
#' unique \strong{Channel.Class} combination found in the folder
#' (e.g. \code{"BF.class1"}, \code{"GFP.class1"}, \code{"mCherry.class2"}).
#' Cells contain absolute file paths; \code{NA} where no matching thumbnail
#' exists.
#'
#' @param se A \code{SummarizedExperiment} (or subclass) object whose
#'   \code{colData} contains at least \code{Well} and \code{Plate_ID} columns.
#' @param folder Path to the parent folder containing the JPG thumbnails
#'   (searched recursively).
#' @param col_name Name of the colData column to write the nested DataFrame
#'   into. Default: \code{"thumbnails"}.
#' @param pattern Regex pattern used to list image files.
#'   Default: \code{"\\\\.(jpg|jpeg)$"}.
#'
#' @return The updated \code{se} object with a nested DataFrame stored in
#'   \code{colData(se)[[col_name]]}.
#'
#' @details
#' \strong{Plate matching}: Sub-directories inside \code{folder} whose
#' \code{basename} contains a \code{Plate_ID} from the SE as a fixed
#' substring are matched to that plate.  If no sub-directory matches but
#' only one plate is present in the SE, the folder itself is used as a
#' fallback.  Files not traceable to a known plate are silently dropped.
#'
#' \strong{Well normalisation}: Wells are standardised to the format
#' \code{"A01"} (upper-case letter + two-digit column, e.g. \code{"B03"}).
#'
#' \strong{Site de-duplication}: When multiple sites exist for the same
#' well/channel/class (e.g. \code{A01-1}, \code{A01-2}) the first
#' alphabetically is kept, consistent with the tinySEV image plate viewer.
#'
#' @examples
#' \dontrun{
#' se <- addThumbnailPaths(se, folder = "path/to/thumbnails")
#' # Inspect stored paths for the first well
#' colData(se)[["thumbnails"]][1, ]
#' }
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @export
addThumbnailPaths <- function(se,
                              folder,
                              col_name = "thumbnails",
                              pattern  = "\\.(jpg|jpeg)$") {

  folder <- normalizePath(folder, winslash = "/", mustWork = TRUE)

  # ------------------------------------------------------------------
  # 1. List all matching image files
  # ------------------------------------------------------------------
  all_files <- list.files(folder,
                          pattern     = pattern,
                          full.names  = TRUE,
                          recursive   = TRUE,
                          ignore.case = TRUE)
  if (length(all_files) == 0) {
    warning("No JPG/JPEG files found in: ", folder)
    return(se)
  }

  # ------------------------------------------------------------------
  # 2. Parse filename metadata
  #    Pattern: ..._<Well>-<site>_<Channel>_<Class>_crop.jpg
  # ------------------------------------------------------------------
  .parse_thumb <- function(fp) {
    bn <- basename(fp)
    m  <- regexpr("_([A-P][0-9]{1,2}-[0-9]+)_([^_]+)_([^_]+)_crop\\.jpe?g$",
                  bn, perl = TRUE, ignore.case = TRUE)
    if (m[1] == -1)
      return(data.frame(file = fp, well = NA_character_,
                        channel = NA_character_, img_class = NA_character_,
                        stringsAsFactors = FALSE))
    tok   <- sub("^_", "", regmatches(bn, m))
    parts <- strsplit(tok, "_")[[1]]     # well-site, channel, class, "crop.jpg"
    well_site <- parts[1]
    channel   <- if (length(parts) >= 2) parts[2] else NA_character_
    img_class <- if (length(parts) >= 3) parts[3] else NA_character_
    raw_well  <- sub("-.*$", "", well_site)   # strip site suffix
    row_ltr   <- sub("^([A-P]).*$",       "\\1", raw_well)
    col_num   <- suppressWarnings(
      as.integer(sub("^[A-P]([0-9]{1,2})$", "\\1", raw_well)))
    well <- if (!is.na(col_num) && col_num >= 1 && col_num <= 24)
      paste0(row_ltr, sprintf("%02d", col_num)) else NA_character_
    data.frame(file = fp, well = well, channel = channel,
               img_class = img_class, stringsAsFactors = FALSE)
  }

  meta <- do.call(rbind, lapply(all_files, .parse_thumb))
  meta <- meta[!is.na(meta$well), , drop = FALSE]

  if (nrow(meta) == 0) {
    warning("No thumbnails matched the expected filename convention ",
            "(*_<Well>-<site>_<Channel>_<Class>_crop.jpg).")
    return(se)
  }

  # ------------------------------------------------------------------
  # 3. Resolve Plate_ID for each file via sub-directory names
  # ------------------------------------------------------------------
  cd        <- as.data.frame(SummarizedExperiment::colData(se))
  plate_ids <- unique(cd$Plate_ID)

  subdirs    <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
  subnames   <- basename(subdirs)

  dir_plate_map <- list()   # dir_full_path -> plate_id
  for (pid in plate_ids) {
    idx <- which(grepl(pid, subnames, fixed = TRUE))
    if (length(idx) > 0)
      dir_plate_map[[subdirs[idx[1]]]] <- pid
  }
  # Fallback: single plate, use the folder itself
  if (length(dir_plate_map) == 0 && length(plate_ids) == 1)
    dir_plate_map[[folder]] <- plate_ids[1]

  .file_to_plate <- function(fp) {
    fp_norm <- normalizePath(fp, winslash = "/", mustWork = FALSE)
    for (d in names(dir_plate_map)) {
      d_norm <- normalizePath(d, winslash = "/", mustWork = FALSE)
      if (startsWith(fp_norm, paste0(d_norm, "/")))
        return(dir_plate_map[[d]])
    }
    if (length(plate_ids) == 1) return(plate_ids[1])
    NA_character_
  }

  meta$plate_id <- vapply(meta$file, .file_to_plate, character(1))
  meta          <- meta[!is.na(meta$plate_id), , drop = FALSE]

  if (nrow(meta) == 0) {
    warning("No thumbnails could be matched to a Plate_ID in the SE. ",
            "Ensure sub-folder names contain the Plate_ID as a substring.")
    return(se)
  }

  # ------------------------------------------------------------------
  # 4. De-duplicate: one file per well x plate x channel x class
  #    (keep first occurrence — alphabetically earliest site)
  # ------------------------------------------------------------------
  meta     <- meta[order(meta$file), ]
  meta_dd  <- meta[!duplicated(meta[, c("plate_id", "well", "channel", "img_class")]),
                   , drop = FALSE]

  # ------------------------------------------------------------------
  # 5. Build nested S4Vectors::DataFrame
  #    Rows  = one per colData row
  #    Cols  = one per Channel.Class combination (e.g. "BF.class1")
  # ------------------------------------------------------------------
  chan_cls <- sort(unique(paste(meta_dd$channel, meta_dd$img_class, sep = ".")))

  thumb_mat <- matrix(NA_character_,
                      nrow = nrow(cd),
                      ncol = length(chan_cls),
                      dimnames = list(NULL, chan_cls))

  for (i in seq_len(nrow(meta_dd))) {
    r   <- meta_dd[i, ]
    idx <- which(cd$Well == r$well & cd$Plate_ID == r$plate_id)
    if (length(idx) == 0) next
    key <- paste(r$channel, r$img_class, sep = ".")
    thumb_mat[idx, key] <- r$file
  }

  SummarizedExperiment::colData(se)[[col_name]] <-
    S4Vectors::DataFrame(as.data.frame(thumb_mat, stringsAsFactors = FALSE))

  n_matched <- sum(!is.na(thumb_mat[, 1]))
  message(sprintf(
    "Stored thumbnail paths in colData column '%s'.\n  %d / %d wells matched.\n  Channel \u00d7 Class keys: %s",
    col_name, n_matched, nrow(cd), paste(chan_cls, collapse = ", ")
  ))
  return(se)
}



#' Display thumbnail images for selected wells from a SummarizedExperiment
#'
#' Retrieves JPG thumbnail paths stored by \code{\link{addThumbnailPaths}} and
#' renders them as a grid.  For \code{show = "both"} (default), BF and fluoro
#' images are placed \strong{side by side} within each well panel.  Multiple
#' wells are tiled according to \code{ncol} / \code{nrow}.
#'
#' @param se A \code{SummarizedExperiment} with a nested thumbnail
#'   \code{DataFrame} in \code{colData} (added by \code{addThumbnailPaths}).
#' @param wells Character vector of well identifiers in
#'   \code{"<Well><sep><Plate_ID>"} format,
#'   e.g. \code{c("B10.18T39383", "C05.18T39383")}.
#' @param img_class Class label to display, e.g. \code{"class1"}.
#' @param show Which images to show: \code{"both"} (BF + fluoro side by side,
#'   default), \code{"BF"}, or \code{"fluoro"}.
#' @param col_name Name of the colData column holding the thumbnail paths.
#'   Default: \code{"thumbnails"}.
#' @param sep Separator between Well and Plate_ID in \code{wells}.
#'   Default: \code{"."}.
#' @param ncol Number of columns in the well grid. \code{NULL} (default) puts
#'   all wells in a single row.
#' @param nrow Number of rows in the well grid. \code{NULL} (default) is
#'   derived from \code{ncol}.
#' @param coldata_vars Optional character vector of \code{colData} columns
#'   (e.g. \code{c("Condition", "QC")}) whose values are appended to each
#'   well's title label.
#' @param label_size Font size (pts) for panel labels. Default: \code{8}.
#' @param bg Background colour for missing-image panels.
#'   Default: \code{"#1a1a1a"}.
#'
#' @return Invisibly returns the assembled \code{gtable}; called for its
#'   side-effect of rendering the plot.
#'
#' @examples
#' \dontrun{
#' # Default: BF | fluoro side by side for each well
#' plotThumbnails(se, wells = c("B10.18T39383", "C05.18T39383"),
#'                img_class = "class1")
#'
#' # 3 wells in a 2-column grid, show Condition in label
#' plotThumbnails(se,
#'                wells        = c("B10.18T39383", "B11.18T39383", "B12.18T39383"),
#'                img_class    = "class1",
#'                ncol         = 2,
#'                coldata_vars = "Condition")
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom jpeg readJPEG
#' @importFrom grid rasterGrob rectGrob textGrob unit gpar grid.newpage grid.draw
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
plotThumbnails <- function(se,
                           wells,
                           img_class,
                           show         = c("both", "BF", "fluoro"),
                           col_name     = "thumbnails",
                           sep          = ".",
                           ncol         = NULL,
                           nrow         = NULL,
                           coldata_vars = NULL,
                           label_size   = 8,
                           bg           = "#1a1a1a") {

  show <- match.arg(show)

  # ------------------------------------------------------------------
  # Retrieve thumbnail lookup table
  # ------------------------------------------------------------------
  thumb_df <- SummarizedExperiment::colData(se)[[col_name]]
  if (is.null(thumb_df))
    stop("colData column '", col_name, "' not found. Run addThumbnailPaths() first.")
  thumb_df <- as.data.frame(thumb_df, stringsAsFactors = FALSE)

  cd    <- as.data.frame(SummarizedExperiment::colData(se))
  avail <- colnames(thumb_df)

  # ------------------------------------------------------------------
  # Identify BF and fluoro column names for the requested class
  # ------------------------------------------------------------------
  bf_col     <- paste0("BF.", img_class)
  fluoro_col <- avail[grepl(paste0("\\.", img_class, "$"), avail) &
                        !grepl("^BF\\.", avail)]

  if (show %in% c("both", "BF") && !bf_col %in% avail)
    stop("Column '", bf_col, "' not found in '", col_name,
         "'. Available: ", paste(avail, collapse = ", "))
  if (show %in% c("both", "fluoro") && length(fluoro_col) == 0)
    stop("No non-BF column for class '", img_class,
         "' found. Available: ", paste(avail, collapse = ", "))
  if (length(fluoro_col) > 1) {
    message("Multiple fluoro columns for class '", img_class, "': ",
            paste(fluoro_col, collapse = ", "), " - using: ", fluoro_col[1])
    fluoro_col <- fluoro_col[1]
  }

  # ------------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------------

  # Load jpeg -> rasterGrob; dark placeholder on failure / missing file
  .load_grob <- function(fp) {
    if (length(fp) == 0 || is.na(fp) || !nzchar(fp) || !file.exists(fp))
      return(grid::rectGrob(gp = grid::gpar(fill = bg, col = NA)))
    img <- tryCatch(jpeg::readJPEG(fp), error = function(e) NULL)
    if (is.null(img))
      return(grid::rectGrob(gp = grid::gpar(fill = "#3a0000", col = NA)))
    grid::rasterGrob(img, interpolate = TRUE,
                     width  = grid::unit(1, "npc"),
                     height = grid::unit(1, "npc"))
  }

  # Small channel sub-label above an image grob
  .sublabel <- function(g, txt) {
    gridExtra::arrangeGrob(
      g,
      top     = grid::textGrob(txt, gp = grid::gpar(fontsize = label_size - 1,
                                                     col      = "#888888")),
      padding = grid::unit(2, "pt")
    )
  }

  # Full well title: well ID + optional colData values
  .well_title <- function(wid, cd_row) {
    lbl <- wid
    if (!is.null(coldata_vars) && length(coldata_vars) > 0 && !is.null(cd_row)) {
      extra <- vapply(coldata_vars, function(v) {
        if (v %in% colnames(cd_row))
          paste0(v, ": ", cd_row[[v]][1])
        else ""
      }, character(1))
      extra <- extra[nzchar(extra)]
      if (length(extra) > 0)
        lbl <- paste0(lbl, "\n", paste(extra, collapse = "  |  "))
    }
    lbl
  }

  # ------------------------------------------------------------------
  # Build one panel grob per well
  # ------------------------------------------------------------------
  panels      <- list()
  found_wells <- character(0)

  for (wid in wells) {
    idx <- regexpr(sep, wid, fixed = TRUE)
    if (idx == -1)
      stop("'", wid, "' does not contain separator '", sep, "'.")
    w   <- substr(wid, 1, idx - 1)
    pid <- substr(wid, idx + nchar(sep), nchar(wid))

    cd_idx <- which(cd$Well == w & cd$Plate_ID == pid)
    if (length(cd_idx) == 0) {
      warning("Well '", w, "' / Plate_ID '", pid,
              "' not found in SE colData - skipped.")
      next
    }

    row     <- thumb_df[cd_idx[1], , drop = FALSE]
    cd_row  <- cd[cd_idx[1], , drop = FALSE]
    title   <- .well_title(wid, cd_row)

    if (show == "both") {
      bf_g  <- .sublabel(.load_grob(row[[bf_col]]),     "BF")
      fl_g  <- .sublabel(.load_grob(row[[fluoro_col]]), "fluoro")
      panel <- gridExtra::arrangeGrob(
        bf_g, fl_g,
        ncol    = 2,
        top     = grid::textGrob(title,
                                 gp = grid::gpar(fontsize = label_size,
                                                 col      = "#333333",
                                                 fontface = "bold")),
        padding = grid::unit(3, "pt")
      )
    } else {
      fp    <- if (show == "BF") row[[bf_col]] else row[[fluoro_col]]
      panel <- gridExtra::arrangeGrob(
        .load_grob(fp),
        top     = grid::textGrob(title,
                                 gp = grid::gpar(fontsize = label_size,
                                                 col      = "#333333",
                                                 fontface = "bold")),
        padding = grid::unit(3, "pt")
      )
    }

    panels      <- c(panels, list(panel))
    found_wells <- c(found_wells, wid)
  }

  if (length(panels) == 0) {
    message("No valid wells to display.")
    return(invisible(NULL))
  }

  # ------------------------------------------------------------------
  # Resolve ncol / nrow defaults and render
  # ------------------------------------------------------------------
  n <- length(panels)

  if (is.null(ncol) && is.null(nrow)) {
    ncol <- n          # single row by default
    nrow <- 1L
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }

  p <- do.call(gridExtra::arrangeGrob, c(panels, list(ncol = ncol, nrow = nrow)))

  grid::grid.newpage()
  grid::grid.draw(p)
  invisible(p)
}
