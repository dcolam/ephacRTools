library(DBI)
library(RSQLite)
library(dplyr)
library(EBImage)
library(grid)
library(gridExtra)
library(glue)
library(ggplot2)
library(viridis)
library(RBioFormats)
library(RImageJROI)
library(tibble)
library(tidyr)
library(pROC)
library(umap)
library(cluster)
library(stringr)
library(SummarizedExperiment)


pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}


## Load Packages
packages <- c("RSQLite", "ggplot2", "readxl", "stringr", "reshape2", "lmerTest", 
              "emmeans", "data.table", "umap", "dplyr", "purrr")
for(package in packages){
  pkgTest(package)
}


# ────────────────────────────────────────────────────────────────────────
## Aggregation function based on columns
ag <- function(df, cols, fun=mean) {
  c <- colnames(df[,unlist(lapply(df, is.numeric))])
  print(c)
  df <- aggregate(df[,c], by=as.list(df[,cols]), FUN=fun, na.rm=TRUE)
  return(df[, colSums(is.na(df)) != nrow(df)])
}

# ────────────────────────────────────────────────────────────────────────
get_metric <- function(df, well, parameter = "Minima", plot=FALSE){

  x <- subset(df, Well == well)
  x <- x[order(x$V_Clamp),]
  Imax <- min(x[,parameter])
  Vmax1 <- min(x[x[,parameter]==Imax,]$V_Clamp)
  
  x_axis <- seq(-80, 20, length=100)
  x[!complete.cases(x[,parameter]),parameter] <- 1
  spl <- smooth.spline(x$V_Clamp, y=x[,parameter])
  fit.pred <- predict(spl, data.frame(V_Clamp=x_axis))
  fit.pred <- as.data.frame(fit.pred)
  names(fit.pred)[which(names(fit.pred) == "V_Clamp.1")] <- "fit.pred"
  Vhalf <- fit.pred[fit.pred$V_Clamp < Vmax1,]
  Vhalf <- Vhalf[which.min(abs(Vhalf$fit.pred- Imax/2)),]$V_Clamp

  if(plot){
    print("Plotting..")
    plot(x$V_Clamp, x[,parameter]*10^12)
    lines(x$V_Clamp, x[,parameter]*10^12)
    abline(v=Vmax1, col=2)
    abline(v=Vhalf, col=3)
    abline(h=Imax*10^12/2)
    lines(x=fit.pred$V_Clamp, y=fit.pred$fit.pred*10^12, col=2)
  }
  return(list(Imax=Imax, Vhalf=Vhalf, Vmax=Vmax1))
}

# ────────────────────────────────────────────────────────────────────────
prepareDF <- function(pathDF){
  
  df <- as.data.frame(read_excel(pathDF, sheet="OA Export", col_types = "text"))
  df$`\r` <- NULL
  names(df)[1:2] <- c("Well", "QC")
  df
  
  df <- df[-1,]
  sweeps <- grep("Sweep \\d", colnames(df), value=TRUE)
  no.sweeps <- unique(sapply(sweeps, FUN=function(s){
    unlist(str_split(s, " "))[2]
  }))
  
  volt <- df[which(df$Well == "Sweep Voltage"),]
  df <- df[-1,]
  volt <- volt[, grep("Compound", names(volt))]
  
  new.cols <- sapply(grep(no.sweeps[1], sweeps, value=T), function(x){
    unlist(str_split(x, " "))[3]
  })
  new.cols <- c("Well", "QC", new.cols, "Sweep", "V_Clamp")
  new.df <- data.frame(matrix(ncol=length(new.cols),nrow=0, dimnames=list(NULL, new.cols)))
  
  for(s in no.sweeps){
    cols <- c("Well", "QC", grep(s, sweeps, value=T))
    temp <- df[,cols]
    temp$Sweep <- s
    temp$V_Clamp <- volt[,grep(s, names(volt), value=T)]
    
    colnames(temp) <- colnames(new.df)
    new.df <- rbind(new.df,temp)
  }
  
  new.df$V_Clamp <- as.numeric(gsub("m", "", new.df$V_Clamp))
  for(cols in colnames(new.df)){
    tryCatch(expr = {
      recoverCol <- new.df[,cols]
      new.df[,cols] <- as.numeric(new.df[,cols])
    }, warning = function(w){
      new.df[,cols] <- new.df[,cols]
    })
  }
  
  
  return(new.df)
  
}

# ────────────────────────────────────────────────────────────────────────
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
  df
}

# ────────────────────────────────────────────────────────────────────────
extract_imaging_cols <- function(pathDB,
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

imageval_alt <- function(img.path, idx) {
  
  load_bright <- function(file, factor)
    brighten_image(readImage(file), factor = factor)
  
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
  
  ## --------------------------------------------------------------------
  img.list <- list.files(img.path, pattern = "\\.tif$", recursive = TRUE,
                         full.names = TRUE)
  imgs <- img.list[grepl(idx, img.list)]
  
  ## ---- bright‑field / BF image (first file) ---------------------------
  img1 <- load_bright(imgs[1], factor = 2)
  bf_channel <- normalize(img1)              
  img1_grob <- make_grob(rotate(img1, 180)) 
  
  ## ---- fluorescence image (second file) -------------------------------
  img2 <- load_bright(imgs[2], factor = 40)
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


assigncell <- function(PAoutput_path, analysis_type) {
  
  df <- extract_imaging_cols(
    pathDB    = PAoutput_path,
    analysis = analysis_type,
    id_cols = c("Date", "Plate_ID", "Well", "Image_ID", "Channel_Name","Selection", "Selection_Area"),
    num_cols = c("Area", "Mean", "IntDen"),
    scale_num = FALSE
    #scale_cols = c("Area", "Mean", "IntDen")
  )
  
  df$Image_Type <- ifelse(df$Image_ID %% 2 != 0, "fluor", "bf")
  
  df$Channel <- ifelse(df$Channel_Name == "C1", "DAPI", 
                       ifelse(df$Channel_Name == "C2", "Green",
                              ifelse(df$Channel_Name == "C3", "Red", NA)))
  
  df_filtered <- df %>% 
    filter(Image_Type == "fluor", Channel_Name %in% c("C2", "C3"), !is.na(CorrSel), CorrSel == "Hole_ROI")
  
  well_summary <- df_filtered %>%
    group_by(Date, Plate_ID, Image_ID, Well, Channel_Name) %>%
    summarize(Mean = mean(Mean, na.rm = TRUE),
              Area = mean(Area, na.rm = TRUE),
              IntDen = mean(IntDen, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Channel_Name, values_from = c(Mean, Area, IntDen), names_sep = "_")
  
  well_summary <- well_summary %>%
    dplyr::rename(
      Mean_Red = Mean_C3,
      Mean_Green = Mean_C2,
      Area_Red = Area_C3,
      Area_Green = Area_C2,
      IntDen_Red = IntDen_C3,
      IntDen_Green = IntDen_C2
    )
  
  well_summary$Mean_Red[is.na(well_summary$Mean_Red)] <- 0
  well_summary$Mean_Green[is.na(well_summary$Mean_Green)] <- 0
  
  well_summary$Row <- substring(well_summary$Well, 1, 1)
  well_summary$Column <- as.numeric(gsub("^[A-Z]+(\\d+)-.*", "\\1", well_summary$Well))
  
  well_summary$Well <- sub("-\\d+$", "", well_summary$Well)
  
  well_summary <- well_summary %>%
    group_by(Column) %>%
    mutate(
      red_range = max(Mean_Red, na.rm = TRUE) - min(Mean_Red, na.rm = TRUE),
      green_range = max(Mean_Green, na.rm = TRUE) - min(Mean_Green, na.rm = TRUE),
      Mean_Red_scaled = ifelse(red_range == 0, 0, (Mean_Red - min(Mean_Red, na.rm = TRUE)) / red_range),
      Mean_Green_scaled = ifelse(green_range == 0, 0, (Mean_Green - min(Mean_Green, na.rm = TRUE)) / green_range)
    ) %>%
    ungroup()
  
  min_active <- 0.05
  
  non_zero_red <- well_summary$Mean_Red_scaled[well_summary$Mean_Red_scaled > 0]
  non_zero_green <- well_summary$Mean_Green_scaled[well_summary$Mean_Green_scaled > 0]
  
  red_threshold <- if (length(non_zero_red) > 0) {
    quantile(non_zero_red, 0.4, na.rm = TRUE)
  } else {
    0 
  }
  
  green_threshold <- 0.1
  
  well_summary <- well_summary %>%
    mutate(Color = case_when(
      (Mean_Red_scaled < min_active & Mean_Green_scaled < min_active) ~ "black",
      (Mean_Red_scaled >= red_threshold & Mean_Green_scaled >= green_threshold) ~ "yellow",
      (Mean_Green_scaled >= green_threshold) ~ "green",
      (Mean_Red_scaled >= red_threshold) ~ "red",
      TRUE ~ "black"
    )) 
  
  return (well_summary)
  
}


