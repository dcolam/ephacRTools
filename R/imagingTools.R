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
#' @param analysis pa or coloc, which table to extract, single option only
#' @param id_cols Columns which metadata about the image and measurement
#' @param num_cols Numeric Columns about the particle metrices
#' @param coloc_cols Colocalisation specific columns to include, default "Second_Channel","Mask_Area"
#' @param scale_num Boolean parameter to include scaled numeric metrices (Default FALSE)
#' @param scale_fun Scale function passed through the numeric columns
#' @return A dataframe
#' @export
prepareImgDF <- function(pathDB,
                               analysis   = "pa",
                               id_cols    = c("Date","Plate_ID","Well",
                                              "Image_ID","Channel_Name",
                                              "Selection","Selection_Area"),
                               num_cols   = c("Area","Mean","IntDen"),
                                coloc_cols = c("Second_Channel","Mask_Area"),
                               scale_num  = FALSE,
                               scale_cols = NULL,
                               scale_fun  = function(x)
                                 as.numeric(scale(x, TRUE, TRUE))){

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
#' Merge together the imaging-results into the Column Data of the SE
#' @param se SummarizedExperiment Object with the Ephys-Data
#' @param df_img DataFrame with imaging results returned by prepareImgDF()
#' @return A dataframe
#' @export
mergeSEandImg <- function(se, df_img, tableType = "pa"){
  df_img <- subset(df_img, Image_Type == "fluor")
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  # Loop through each channel
  if(tableType == "pa"){
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

  }else{
    channels <- unique(df_img$Channel_Name)

    for (channel in channels) {
      second_channels <- unique(subset(df_img, Channel_Name == channel)$Second_Channel)
      for (second_channel in second_channels){
      # Subset df_img for current channel
      df_channel <- df_img %>%
        filter(Channel_Name == channel, Second_Channel == second_channel) %>%
        select(-Channel_Name, -Second_Channel)  # optional: remove the channel label
      # Perform join
      joined <- cd %>%
        dplyr::left_join(df_channel, by = c("Well", "Plate_ID"))
      # Extract just the new columns (everything except original colData)
      new_cols <- setdiff(names(joined), names(cd))
      # Create a DataFrame object from just the new data
      channel_data <- DataFrame(joined[, new_cols])
      # Assign to colData(se), one column per channel, as a nested DataFrame
      SummarizedExperiment::colData(se)[[paste(channel, second_channel, sep=".")]] <- channel_data
      }
    }

  }
  return(se)
}
#' Function that generates all the files paths for the images
#' @param parent_folder Path to where your various experimental data is stored
#' @param idx The well(s) you want to look at
#' @param plate_ID The plate ID(s) in question
#' @param location The location of the plate_ID in the folder names (we assume here that it is always in the same spot)
#' @return A brightened image
#' @export
image_paths <- function(parent_folder, idx, plate_ID, location) {
  all_dirs <- list.dirs(parent_folder, full.names = TRUE, recursive = FALSE)
  matched_dir <- NULL
  short_idx <- gsub("^([A-Z])0*", "\\1", idx)
  for (d in all_dirs) {
    dir_name <- basename(d)
    parts <- strsplit(dir_name, "_")[[1]]
    if(length(parts) >= location && parts[location] == plate_ID) {
      particle_path <- file.path(d, "Particle_Analysis")
      if (!dir.exists(particle_path)) next

      subdirs <- list.dirs(particle_path, full.names = TRUE, recursive = FALSE)

      for (sub in subdirs) {
        sub_dir_name <- basename(sub)
        sub_parts <- strsplit(sub_dir_name, "_")[[1]]
        if (length(sub_parts) >= location && sub_parts[location] == plate_ID) {
          matched_dir <- sub
          break
        }
      }
      if (is.null(matched_dir)) break
    }
  }
  if (is.null(matched_dir)) {
    stop(paste("No folder with the plate_ID", plate_ID, "in the", location,"th position found"))
    return(NA)
  }
  img.list <- list.files(matched_dir, pattern = "\\.tif$", recursive = TRUE,
                         full.names = TRUE)
  pattern <- paste0(short_idx, "-\\d+")
  imgs <- img.list[grepl(pattern, img.list)]
  if (length(imgs) == 0) {
    warning(paste("No matching images found for well", idx, "in", matched_dir))
    return(NA)
  }
  return(imgs)
}
#' Helper function for imageval to brighten the images
#' @param img.path Path to where your particle analysis images are located
#' @param factor Factor by which the image should be brightened
#' @return A brightened image
#' @export
brighten_image <- function(img, factor = 1.5) {
  img_normalized <- (img - min(img)) / (max(img) - min(img))
  img_brightened <- img_normalized * factor
  img_brightened <- pmin(img_brightened, 1)
  return(img_brightened)
}
#' Plot the images with the different channels for a given well and plate
#' @param se_imagepath Path to your summarized experiment
#' @param idx Well number you want to look at, in the form of "H14" or "H09"
#' @param plate_ID Plate number you are interested in exploring
#' @return A plot of four images (brightfield, green and red channels and overlay of the two channels on top of the BF)
#' @export
imageval <- function(se, idx, plate_ID) {
  load_bright <- function(file, factor)
    brighten_image(readImage(file), factor = factor)

  short_idx <- gsub("^([A-Z])0*", "\\1", idx)

  # helper function: take one channel and put it in the desired colour slot
  make_grob <- function(img, src_slice = NULL, colour = c("red","green","blue")) {
    if (is.null(src_slice)) {                        # full RGB
      x <- normalize(img)
    } else {                                         # monochrome as chosen colour
      colour  <- match.arg(colour)
      chan    <- normalize(img[,,src_slice])
      rgb_arr <- array(0, dim = c(dim(chan), 3))
      rgb_arr[,, match(colour, c("red","green","blue")) ] <- chan
      x <- rgb_arr
    }
    rasterGrob(x, interpolate = TRUE)
  }

  ## --- Load images from the matched directory ---
  imgs <- se$Image_paths[idx]

  ## ---- bright‑field / BF image (first file) ---------------------------
  all_imgs <- imgs[[1]]
  bf_img <- all_imgs[grepl("BF\\.tif$", all_imgs)]
  img1 <- load_bright(bf_img, factor = 2)
  bf_channel <- normalize(img1)
  img1_grob <- make_grob(rotate(img1, -90))

  ## ---- fluorescence image (second file) -------------------------------
  fluorescent_img <- all_imgs[grepl("nm\\.tif$", all_imgs)]

  if (length(fluorescent_img) != 1) {
    stop(paste("Expected exactly 1 fluorescent image ending in 'nm.tif', but found", length(fluorescent_img)))
  }

  img2 <- load_bright(fluorescent_img, factor = 40)
  img2_grob_green <- make_grob(img2, src_slice = 2, colour = "green")  # plane 2 → green
  img2_grob_red   <- make_grob(img2, src_slice = 3, colour = "red")    # plane 3 → red

  ## composite: R=red, G=green, B=bright‑field
  comp_rgb <- array(0, dim = c(dim(img2[,,3]), 3))
  comp_rgb[,,1] <- normalize(img2[,,3])
  comp_rgb[,,2] <- normalize(img2[,,2])
  comp_rgb[,,3] <- bf_channel
  img2_grob_color <- rasterGrob(comp_rgb, interpolate = TRUE)

  ## ---- arrange --------------------------------------------------------
  grid.arrange(img1_grob, img2_grob_color, img2_grob_green, img2_grob_red, ncol = 2)

}


