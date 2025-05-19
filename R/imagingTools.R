#' @importFrom magrittr %>%
NULL
#' Prepare Imaging-results tables from Cluster-Analysis SQLite databases
#' @param pathDB Path to SQlite-DB
#' @param analysis pa or coloc, which table to extract
#' @param id_cols Columns which metadata about the image and measurement
#' @param num_cols Numeric Columns about the particle metrices
#' @param scale_num Boolean parameter to include scaled numeric metrices (Default FALSE)
#' @param scale_fun Scale function passed through the numeric columns
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
                                    as.numeric(scale(x, TRUE, TRUE))) {

    analysis <- match.arg(analysis)

    ## ---------- helper that does your existing pipeline ------------------
    process_tbl <- function(tbl) {

      tbl <- dplyr::select(tbl, tidyselect::any_of(c(id_cols, num_cols)))

      tbl <- ag(tbl, cols = id_cols, fun = mean)   # your ag()
      tbl <- df_cleaned(tbl)                       # your ROI labelling
      tbl <- tbl[ , !grepl("(\\.1|\\.\\.\\.[0-9]+)$", names(tbl)) ]
      tbl <- dplyr::filter(tbl, CorrSel == "Hole_ROI")

      # optional scaling
      if (isTRUE(scale_num)) {
        if (is.null(scale_cols))
          scale_cols <- intersect(num_cols, names(tbl))
        for (col in scale_cols) {
          new_col <- paste0(col, "_Scaled")
          tbl[[new_col]] <- scale_fun(tbl[[col]])
          tbl <- dplyr::relocate(tbl, dplyr::all_of(new_col),
                                 .after = dplyr::all_of(col))
        }
      }

      tbl
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
#' @param analysis pa or coloc, which table to extract
#' @param id_cols Columns which metadata about the image and measurement
#' @param num_cols Numeric Columns about the particle metrices
#' @param scale_num Boolean parameter to include scaled numeric metrices (Default FALSE)
#' @param scale_fun Scale function passed through the numeric columns
#' @return A dataframe
#' @export
prepareImgDF <- function(pathDB,
                               analysis   = c("pa", "coloc"),
                               id_cols    = c("Date","Plate_ID","Well",
                                              "Image_ID","Channel_Name",
                                              "Selection","Selection_Area"),
                               num_cols   = c("Area","Mean","IntDen"),
                               scale_num  = FALSE,
                               scale_cols = NULL,
                               scale_fun  = function(x)
                                 as.numeric(scale(x, TRUE, TRUE))){


  if(length(pathDB) == 1){
    df <- prepareSingleImgDF(pathDB, analysis=analysis,
                             id_cols = id_cols,
                             num_cols = num_cols,
                             scale_num=scale_num,
                             scale_cols=scale_cols,
                             scale_fun=scale_fun)
  }else{
    dfs <- lapply(pathDB, function(x) {
      prepareSingleImgDF(x, analysis=analysis,
                         id_cols = id_cols,
                         num_cols = num_cols,
                         scale_num=scale_num,
                         scale_cols=scale_cols,
                         scale_fun=scale_fun)
    })

    names(dfs) <- pathDB
    df <- dplyr::bind_rows(dfs, .id = "column_label")
    df$Plate_ID <- sapply(df$Plate_ID, function(x){
      unlist(stringr::str_split(x, "\\r"))[1]
    })
  }
  return(df)
}
#' Clean and normalize dataframe. Adds column and row identifiers and finds the
#' selections with the Hole-ROI
#' @param df dataframe of image results
#' @return A dataframe
#' @export
df_cleaned <- function(df){

  df$Well_clean <- sapply(df$Well, function(x){

    unlist(str_split(x, "-"))[1]

  })

  df$Row <- sapply(df$Well_clean, function(x){

    str_sub(x, 1, 1)

  })

  df$Column <- sapply(df$Well_clean, function(x){

    str_sub(x, 2, 3)

  })

  df$CorrSel <- NA
  for(sel in unique(df$Selection)){
    #sel <- "4617.vsi - 283 BF"
    selMin <- min(subset(df, Selection == sel)$Selection_Area)
    selMax <- max(subset(df, Selection == sel)$Selection_Area)

    df[df$Selection == sel & df$Selection_Area == selMin,"CorrSel"] <- "Hole_ROI"
    df[df$Selection == sel & df$Selection_Area == selMax,"CorrSel"] <- "background_ROI"
  }

  df$Image_Type <- ifelse(df$Image_ID %% 2 != 0, "fluor", "bf")

  df$Channel <- ifelse(df$Channel_Name == "C1", "DAPI",
                       ifelse(df$Channel_Name == "C2", "Green",
                              ifelse(df$Channel_Name == "C3", "Red", NA)))
  df$Well <- paste(df$Row, stringr::str_pad(df$Column, 2, pad = "0"), sep="")

  df
}
#' Merge together the imgaging-results into the Column Data of the SE
#' @param se SummarizedExperiment Object with the Ephys-Data
#' @param df_img DataFrame with imaging results returned by prepareImgDF()
#' @return A dataframe
#' @export
mergeSEandImg <- function(se, df_img){
  df_img <- subset(df_img, Image_Type == "fluor")
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  # Loop through each channel
  channels <- unique(df_img$Channel_Name)

  for (channel in channels) {
    # Subset df_img for current channel
    df_channel <- df_img %>%
      filter(Channel_Name == channel) %>%
      select(-Channel_Name)  # optional: remove the channel label

    # Perform join
    joined <- cd %>%
      dplyr::left_join(df_channel, by = c("Well", "Plate_ID"))

    # Extract just the new columns (everything except original colData)
    new_cols <- setdiff(names(joined), names(cd))

    # Create a DataFrame object from just the new data
    channel_data <- DataFrame(joined[, new_cols])

    # Assign to colData(se), one column per channel, as a nested DataFrame
    SummarizedExperiment::colData(se)[[channel]] <- channel_data
  }
  return(se)
}

