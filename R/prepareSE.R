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

  # Read file
  sheets <- readxl::excel_sheets(pathToDF)
  sheets <- sheets[grepl("Export", sheets)]
  sheets <- sheets[length(sheets)]  # prefer last Export sheet (main data, not header)
  df <- readxl::read_excel(pathToDF, sheet = sheets, col_types = "text")
  df <- as.data.frame(df)
  # Clean up unwanted column
  if ("\r" %in% colnames(df)) df$`\r` <- NULL

  # Standardize columns
  names(df)[1:2] <- c("well_id", "qc")
  df$well_id <- sapply(df$well_id, function(x) unlist(stringr::str_split(x, "\\r"))[1])

  # Ensure barcode column exists
  if (!("Nanion Chip Barcode" %in% colnames(df))) {
    df[,"Nanion Chip Barcode"] <- "NoPlateID"
  }

  # Handle voltage sweep info
  volt_steps <- FALSE
  if ("Sweep Voltage" %in% df$well_id | "Abs" %in% df$well_id) {
    volt <- df[grepl("Sweep Voltage", df$well_id), ]
    df <- df[!grepl("Sweep Voltage", df$well_id), ]
    volt <- volt[, grep("Compound", names(volt))]
    volt_steps <- TRUE
  } else {
    volt <- df[1, ]
    volt[1, ] <- "NAm"
    volt <- volt[, grep("Compound", names(volt))]
  }

  # Filter to acceptable wells
  acceptable_wells <- as.vector(outer(LETTERS[1:16], sprintf("%02d", 1:24), paste0))
  df <- df[df$well_id %in% acceptable_wells,]

  # Sweep parsing
  sweeps <- grep("Sweep \\d", colnames(df), value = TRUE)
  no.sweeps <- unique(sapply(sweeps, function(s) unlist(stringr::str_split(s, " "))[2]))

  new.cols <- sapply(grep(no.sweeps[1], sweeps, value = TRUE), function(x) {
    unlist(stringr::str_split(x, " "))[3]
  })
  new.cols <- c("well_id", "qc", "plate_id", new.cols, "sweep", "v_clamp_mV")
  new.df <- data.frame(matrix(ncol = length(new.cols), nrow = 0, dimnames = list(NULL, new.cols)))

  # Rebuild the long-format dataframe
  for (s in no.sweeps) {
    cols <- c("well_id", "qc", "Nanion Chip Barcode", grep(s, sweeps, value = TRUE))
    tempdf <- df[, cols]
    try({
    tempdf$sweep <- s
    tempdf$v_clamp_mV <- volt[, grep(s, names(volt), value = TRUE)]
    colnames(tempdf) <- colnames(new.df)
    new.df <- rbind(new.df, tempdf)
    })
  }

  # Optional numeric conversion for v_clamp_mV
  if (volt_steps) {
    new.df$v_clamp_mV <- as.numeric(gsub("m", "", new.df$v_clamp_mV))
  }

  # Standardize plate_id column
  new.df$plate_id <- sapply(new.df$plate_id, function(x) unlist(stringr::str_split(x, "\\r"))[1])

  # Re-type using hablar
  new.df <- new.df %>% hablar::retype()

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
prepareMultipleDFs <- function(pathList, progress_callback = NULL){

  n <- length(pathList)
  if (is.null(progress_callback) && n > 1) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    progress_callback <- function(i, total, name) {
      utils::setTxtProgressBar(pb, i)
      if (i == total) close(pb)
    }
  }

  dfs <- lapply(seq_along(pathList), function(i) {
    x <- pathList[[i]]
    df <- prepareDF(as.character(x))
    if (is.null(df)) {
      warning(paste("Skipping file due to read error:", x))
    }
    if (!is.null(progress_callback)) progress_callback(i, length(pathList), basename(x))
    return(df)
  })
  dfs <- dfs[!sapply(dfs, is.null)]
  if (length(dfs) == 0) {
    stop("All uploaded files failed to read. Please check file format.")
  }
  safe_names <- lapply(pathList, function(x){basename(x)})
  names(dfs) <- safe_names

  df <- dplyr::bind_rows(dfs, .id = "column_label")
  if (!"plate_id" %in% colnames(df)) {
    df$plate_id <- df$column_label
  } else {
    df$plate_id <- sapply(df$plate_id, function(x){
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
prepareSE <- function(pathDF, conditionColumns= c("Compound"), progress_callback = NULL){

if(length(pathDF) > 1){
    df <- prepareMultipleDFs(pathDF, progress_callback = progress_callback)
}else{
    if (!is.null(progress_callback)) progress_callback(1L, 1L, basename(pathDF))
    df <- prepareDF(as.character(pathDF))
}

  df <- df %>% hablar::retype()

  numeric_cols <- names(df)[sapply(df, is.numeric)]
  numeric_cols <- numeric_cols[!(numeric_cols %in% c("sweep", "v_clamp_mV"))]
  description_cols <- colnames(df)[!(colnames(df) %in% numeric_cols)]
  description_cols <- description_cols[description_cols != "sweep"]

  assays <- lapply(numeric_cols, \(cols) {
    x <- reshape2::dcast(df, well_id*plate_id ~sweep, value.var = cols)
    m <- x[, !(names(x) %in% c("well_id", "plate_id"))] |> as.matrix()
    t(m)
  })

  cols <- reshape2::dcast(df, well_id*plate_id ~sweep, value.var = "well_id")
  cols <- interaction(cols$well_id, cols$plate_id)
  names(assays) <- numeric_cols

  df <- data.table::as.data.table(df)

  cd <- S4Vectors::DataFrame(unique(df[, .(well_id, qc, plate_id)]))
  rownames(cd) <- interaction(cd$well_id, cd$plate_id)
  cd <- cd[cols,]

  cd$row <- sapply(cd$well_id, function(x) {
    stringr::str_sub(x, 1, 1)
  })

  cd$col <- sapply(cd$well_id, function(x) {
    stringr::str_sub(x, 2, 3)
  })

  rd <- S4Vectors::DataFrame(unique(df[, .(sweep)]))
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                   rowData = rd,
                                                   colData = cd)

  colnames(se) <- as.character(interaction(cd$well_id, cd$plate_id))

  description_cols <- description_cols[!(description_cols %in% names(cd))]

  descr <- lapply(description_cols, function(var) {
    x <- reshape2::dcast(df, sweep ~ well_id+plate_id, value.var = var)
    x[, !(names(x) %in% c("sweep"))]
  })

  names(descr) <- description_cols

  for (colname in names(descr)){

    mat <- descr[[colname]]

    # Check if all rows have the same values across columns
    same_across <- apply(mat, 1, function(x) length(unique(x)) == 1)

    if (all(same_across)) {
      collapsed <- apply(mat, 1, function(x) unique(x))
      SummarizedExperiment::rowData(se)[[colname]] <- collapsed
    } else {
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
