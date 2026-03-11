#' @importFrom magrittr %>%
NULL
#' Prepare Imaging-results tables from Cluster-Analysis SQLite databases
#' @param pathDB Path to SQlite-DB
#' @param analysis pa or coloc, which table to extract
#' @param id_cols Columns which metadata about the image and measurement
#' @param num_cols Numeric Columns about the particle metrices
#' @param scale_num Boolean parameter to include scaled numeric metrices (Default FALSE)
#' @param scale_fun Scale function passed through the numeric columns
#' @param aggregate Boolean; if TRUE (default) rows are averaged across sites before cleaning
#' @param cleanNames Boolean; if TRUE (default) runs df_cleaned() and drops duplicated column suffixes
#' @return A dataframe
#' @export
prepareSingleImgDF <- function(pathDB,
                               analysis   = c("pa", "coloc"),
                               id_cols    = c("Date","Plate_ID","Well",
                                              "Image_ID","Channel_Name",
                                              "Selection","Selection_Area"),
                               num_cols   = c("Area","Mean","IntDen"),
                               scale_num  = FALSE,
                               scale_cols = NULL,
                               scale_fun  = function(x)
                                 as.numeric(scale(x, TRUE, TRUE)),
                               aggregate  = TRUE,
                               cleanNames = TRUE) {

  analysis <- match.arg(analysis)

  ## ---------- helper that does your existing pipeline ------------------
  process_tbl <- function(tbl) {
    tbl <- dplyr::select(tbl, tidyselect::any_of(c(id_cols, num_cols)))

    if (aggregate) {
      id_cols <- intersect(names(tbl), id_cols)
      tbl <- ag(tbl, cols = id_cols, fun = mean)
    }

    if (isTRUE(cleanNames)) {
      tbl <- df_cleaned(tbl)
      tbl <- tbl[ , !grepl("(\\.1|\\.\\.\\.[0-9]+)$", names(tbl)) ]
    }

    # optional scaling
    if (isTRUE(scale_num)) {
      if (is.null(scale_cols)){
        scale_cols <- intersect(num_cols, names(tbl))}
      for (col in scale_cols) {
        new_col <- paste0(col, "_Scaled")
        tbl[[new_col]] <- scale_fun(tbl[[col]])
        tbl <- dplyr::relocate(tbl, dplyr::all_of(new_col),
                               .after = dplyr::all_of(col))
      }
    }
      return(tbl)
  }
    con <- DBI::dbConnect(RSQLite::SQLite(), pathDB)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    if (analysis == "pa") {
      tbl <- DBI::dbGetQuery(con, "
        SELECT *
        FROM Particle_Analysis_Table  AS pa
        JOIN  PA_Measurement_Tables   AS meas
             ON meas.PA_ID = pa.PA_ID")
  } else {                          # analysis == "coloc"
    tbl <- DBI::dbGetQuery(con, "
        SELECT *
        FROM Coloc_Analysis_Table  AS ca
        JOIN  Coloc_Measurement_Tables  AS meas
             ON meas.COLOC_ID = ca.COLOC_ID")
  }
    process_tbl(tbl)
}
#' Prepare Imaging-results tables from Cluster-Analysis SQLite databases
#' @param pathDB Path to SQlite-DB
#' @param analysis pa or coloc, which table to extract, single option only
#' @param id_cols Columns which metadata about the image and measurement
#' @param num_cols Numeric Columns about the particle metrices
#' @param coloc_cols Colocalisation specific columns to include, default "Second_Channel","Mask_Area"
#' @param scale_num Boolean parameter to include scaled numeric metrices (Default FALSE)
#' @param scale_fun Scale function passed through the numeric columns
#' @param aggregate Boolean; if TRUE (default) rows are averaged across sites before cleaning
#' @param cleanNames Boolean; if TRUE (default) runs df_cleaned() and drops duplicated column suffixes
#' @return A dataframe
#' @export
prepareImgDF <- function(pathDB,
                               analysis   = "pa",
                               id_cols    = c("Date","Plate_ID","Well",
                                              "Image_ID","Channel_Name",
                                              "Selection","Selection_Area"),
                               num_cols   = c("Area","Mean","IntDen",
                                              "Number_of_Particles"),
                               coloc_cols = c("Second_Channel","Mask_Area"),
                              new_channel_names = c("blue", "green", "red", "farred"),
                               scale_num  = FALSE,
                               scale_cols = NULL,
                               scale_fun  = function(x)
                                 as.numeric(scale(x, TRUE, TRUE)),
                               aggregate  = TRUE,
                               cleanNames = TRUE){

  if("coloc" %in% analysis){
    id_cols <- c(id_cols,
                 coloc_cols)
  }

  if(length(pathDB) == 1){
    df <- prepareSingleImgDF(pathDB, analysis=analysis,
                             id_cols = id_cols,
                             num_cols = num_cols,
                             scale_num=scale_num,
                             scale_cols=scale_cols,
                             scale_fun=scale_fun,
                             aggregate=aggregate,
                             cleanNames=cleanNames)
  }else{
    dfs <- lapply(pathDB, function(x) {
      out <- prepareSingleImgDF(x, analysis=analysis,
                                id_cols = id_cols,
                                num_cols = num_cols,
                                scale_num=scale_num,
                                scale_cols=scale_cols,
                                scale_fun=scale_fun,
                                aggregate=aggregate,
                                cleanNames=cleanNames)
      if (!is.data.frame(out)) {
        warning(paste("prepareSingleImgDF failed for", x))
        return(NULL)
      }
      return(out)
    })
      #safe_names <- lapply(pathDB, function(x){basename(x)})
      #print(safe_names)
      #names(dfs) <- safe_names
      dfs <- Filter(Negate(is.null), dfs)
      names(dfs) <- paste0("DB", seq_along(pathDB))
      df <- dplyr::bind_rows(dfs, .id = "column_label")
    }

    if (!"Plate_ID" %in% colnames(df)) {
      df$Plate_ID <- df$column_label
    }else{
      df$Plate_ID <- sapply(df$Plate_ID, function(x){
        unlist(stringr::str_split(x, "\\r"))[1]
      })
    }
  df$Image_ID <- as.numeric(df$Image_ID)
  df$Image_Type <- ifelse(df$Image_ID %% 2 != 0, "fluor", "bf")

  return(df)
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
    #sel <- "4617.vsi - 283 BF"
    selMin <- min(subset(df, Selection == sel)$Selection_Area)
    selMax <- max(subset(df, Selection == sel)$Selection_Area)

    df[df$Selection == sel & df$Selection_Area == selMin,"CorrSel"] <- "Hole_ROI"
    df[df$Selection == sel & df$Selection_Area == selMax,"CorrSel"] <- "background_ROI"
  }
  df$Image_ID <- as.numeric(df$Image_ID)
  df$Image_Type <- ifelse(as.numeric(df$Image_ID) %% 2 != 0, "fluor", "bf")

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
