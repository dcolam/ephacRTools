#' @importFrom magrittr %>%
NULL
#' @import data.table
NULL
#' Prepare a DataControl Excel File into a DataFrame in long format
#' Make sure that it contains Online Analysis Results as well as QCs per sweep:
#' - Seal Resistance
#' - Series Resistance
#' - Capacitance
#' Make sure that the Excel file contains a "Nanion Chip Barcode" column
#' @param path Path to DataControl Excel
#'
#' @return A cleaned data.frame in long format
#' @export
prepareDF <- function(pathToDF) {

  cat("📦 Starting prepareDF\n")
  cat("📁 Path:", pathToDF, "\n")
  cat("📄 File exists:", file.exists(pathToDF), "\n")
  cat("🧠 Memory (start):", format(utils::object.size(ls(envir = environment())), units = "auto"), "\n")

  # Read file
  df <- readxl::read_excel(pathToDF, sheet = "OA Export")
  df <- as.data.frame(df)
  cat("✅ Full Excel loaded\n")
  cat("🧠 Memory (after read):", format(object.size(df), units = "auto"), "\n")

  # Clean up unwanted column
  if ("\r" %in% colnames(df)) df$`\r` <- NULL

  # Standardize columns
  names(df)[1:2] <- c("Well", "QC")
  df$Well <- sapply(df$Well, function(x) unlist(stringr::str_split(x, "\\r"))[1])

  # Ensure barcode column exists
  if (!("Nanion Chip Barcode" %in% colnames(df))) {
    df[,"Nanion Chip Barcode"] <- "NoPlateID"
  }

  # Handle voltage sweep info
  volt_steps <- FALSE
  if ("Sweep Voltage" %in% df$Well) {
    volt <- df[grepl("Sweep Voltage", df$Well), ]
    df <- df[!grepl("Sweep Voltage", df$Well), ]
    volt <- volt[, grep("Compound", names(volt))]
    volt_steps <- TRUE
  } else {
    volt <- df[1, ]
    volt[1, ] <- "NAm"
    volt <- volt[, grep("Compound", names(volt))]
  }

  # Filter to acceptable wells
  acceptable_wells <- as.vector(outer(LETTERS[1:16], sprintf("%02d", 1:24), paste0))
  df <- df[df$Well %in% acceptable_wells,]

  # Sweep parsing
  sweeps <- grep("Sweep \\d", colnames(df), value = TRUE)
  no.sweeps <- unique(sapply(sweeps, function(s) unlist(stringr::str_split(s, " "))[2]))

  new.cols <- sapply(grep(no.sweeps[1], sweeps, value = TRUE), function(x) {
    unlist(stringr::str_split(x, " "))[3]
  })
  new.cols <- c("Well", "QC", "Plate_ID", new.cols, "Sweep", "V_Clamp")
  new.df <- data.frame(matrix(ncol = length(new.cols), nrow = 0, dimnames = list(NULL, new.cols)))

  # Rebuild the long-format dataframe
  for (s in no.sweeps) {
    cols <- c("Well", "QC", "Nanion Chip Barcode", grep(s, sweeps, value = TRUE))
    tempdf <- df[, cols]
    tempdf$Sweep <- s
    tempdf$V_Clamp <- volt[, grep(s, names(volt), value = TRUE)]
    colnames(tempdf) <- colnames(new.df)
    new.df <- rbind(new.df, tempdf)
  }

  # Optional numeric conversion for V_Clamp
  if (volt_steps) {
    new.df$V_Clamp <- as.numeric(gsub("m", "", new.df$V_Clamp))
  }

  # Standardize Plate_ID column
  new.df$Plate_ID <- sapply(new.df$Plate_ID, function(x) unlist(stringr::str_split(x, "\\r"))[1])

  # Re-type using hablar
  new.df <- new.df %>% hablar::retype()

  cat("📈 Preview of final dataframe:\n")
  print(head(new.df, 3))
  cat("🧠 Memory (before cleanup):", format(object.size(new.df), units = "auto"), "\n")

  # Cleanup: Remove all but new.df
  rm(list = setdiff(ls(), "new.df"))
  gc()
  cat("✅ prepareDF complete\n")
  cat("🧠 Memory (after cleanup):", format(object.size(new.df), units = "auto"), "\n")

  return(new.df)
}

#' Prepare several DataControl Excel File into a DataFrame in long format
#' Make sure that it contains Online Analysis Results as well as QCs per sweep:
#' - Seal Resistance
#' - Series Resistance
#' - Capacitance
#' Make sure that the Excel file contains a "Nanion Chip Barcode" column
#' @param vector list of Paths to DataControl Excel files
#'
#' @return A cleaned data.frame in long format
#' @export
prepareMultipleDFs <- function(pathList){

  dfs <- lapply(pathList, function(x) {
    print(paste("File:", x, "exists:", file.exists(x)))
    df <- prepareDF(as.character(x))
    if (is.null(df)) {
      warning(paste("Skipping file due to read error:", x))
    }
    return(df)
  })
  dfs <- dfs[!sapply(dfs, is.null)]
  if (length(dfs) == 0) {
    stop("All uploaded files failed to read. Please check file format.")
  }
  print("Excels Loaded")
  safe_names <- lapply(pathList, function(x){basename(x)})
  print(safe_names)
  names(dfs) <- safe_names

  df <- dplyr::bind_rows(dfs, .id = "column_label")
  if (!"Plate_ID" %in% colnames(df)) {
    df$Plate_ID <- df$column_label
  }else{
    df$Plate_ID <- sapply(df$Plate_ID, function(x){
      unlist(stringr::str_split(x, "\\r"))[1]
    })
  }

  return(df)
}


#' Prepare SummarizedExperiment Object from DataControl Excel-file
#' @param pathDF Path to DataControl Excel files
#' @param conditionColumns array of columns that describe experimental conditions
#'
#' @return A SummarizedExperiment Object with OAs as assays
#' @export
prepareSE <- function(pathDF, conditionColumns= c("Compound")){

#pathDF <- "data-raw/iNeurons/IV neurons_14.28.48_18T39265_LC_new_LC.xlsx"
#pathDF <- l_files
#length(pathDF)

if(length(pathDF) > 1){
    df <- prepareMultipleDFs(pathDF)
    print("multi")
}else{
      df <- prepareDF(as.character(pathDF))
      print("single")
}

  df <- df %>% hablar::retype()

  numeric_cols <- names(df)[sapply(df, is.numeric)]
  numeric_cols <- numeric_cols[!(numeric_cols %in% c("Sweep", "V_Clamp"))]
  description_cols <- colnames (df)[!(colnames (df) %in% numeric_cols)]
  description_cols <- description_cols[description_cols != "Sweep"]

  assays <- lapply(numeric_cols, \(cols) {
    x <- reshape2::dcast(df, Well*Plate_ID ~Sweep, value.var = cols)
    m <- x[ , !(names(x) %in% c("Well", "Plate_ID"))] |> as.matrix()
    #names(m) <-interaction(x$Well, x$Plate_ID)
    t(m)
    #colnames(m) <- x$Well
    #m
  })

  cols <- reshape2::dcast(df, Well*Plate_ID ~Sweep, value.var = "Well")
  cols <- interaction(cols$Well, cols$Plate_ID)
  names(assays) <- numeric_cols

  df <- data.table::as.data.table(df)

  cd <- S4Vectors::DataFrame(unique(df[, .(Well, QC, Plate_ID)]))
  rownames(cd) <- interaction(cd$Well, cd$Plate_ID)
  cd <- cd[cols,]

  cd$Row <- sapply(cd$Well, function(x){

    stringr::str_sub(x, 1, 1)

  })

  cd$Column <- sapply(cd$Well, function(x){

    stringr::str_sub(x, 2, 3)

  })

  rd <- S4Vectors::DataFrame(unique(df[, .(Sweep)]))
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                   rowData = rd,
                                                   colData = cd)

  colnames(se) <- cd$Well


  description_cols <- description_cols[!(description_cols %in% names(cd))]

  descr <-lapply(description_cols, function(var) {
    x<- reshape2::dcast(df, Sweep ~ Well+Plate_ID, value.var = var)
    x[, !(names(x) %in% c("Sweep"))]
  })

  names(descr) <- description_cols


  for (colname in names(descr)){


    mat <- descr[[colname]]  # descr[[colname]] is a DataFrame or matrix-like

    # Check if all rows have the same values across columns
    same_across <- apply(mat, 1, function(x) length(unique(x)) == 1)

    if (all(same_across)) {
      # Collapse to a single column with unique values per row
      collapsed <- apply(mat, 1, function(x) unique(x))
      SummarizedExperiment::rowData(se)[[colname]] <- collapsed
    } else {
      # Retain the full DataFrame if values vary
      SummarizedExperiment::rowData(se)[[colname]] <- mat

    }
  }

  se <- SingleCellExperiment::SingleCellExperiment(
    assays=assays(se),
    colData=colData(se),
    rowData= rowData(se)
  )
  return(se)


}




