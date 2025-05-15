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
prepareDF <- function(pathDF){

  df <- as.data.frame(readxl::read_excel(pathDF, sheet="OA Export", col_types = "text"))
  df$`\r` <- NULL
  names(df)[1:2] <- c("Well", "QC")
  df

  df <- df[-1,]
  sweeps <- grep("Sweep \\d", colnames(df), value=TRUE)
  no.sweeps <- unique(sapply(sweeps, FUN=function(s){
    unlist(stringr::str_split(s, " "))[2]
  }))

  volt <- df[which(df$Well == "Sweep Voltage"),]

  if(nrow(volt) != 0){
    df <- df[-1,]
    volt <- volt[, grep("Compound", names(volt))]
    volt_steps <- TRUE
  }else{
    #volt <- NA
    volt[1,] <- "NAm"
    volt_steps <- FALSE

  }
  new.cols <- sapply(grep(no.sweeps[1], sweeps, value=T), function(x){
    unlist(stringr::str_split(x, " "))[3]
  })

  if(!volt_steps){
    new.cols <- c("Well", "QC","Plate_ID", new.cols, "Sweep")
    }else{
  new.cols <- c("Well", "QC","Plate_ID", new.cols, "Sweep", "V_Clamp")
  }
  new.df <- data.frame(matrix(ncol=length(new.cols),nrow=0, dimnames=list(NULL, new.cols)))

  for(s in no.sweeps){
    cols <- c("Well", "QC", "Nanion Chip Barcode", grep(s, sweeps, value=T))
    temp <- df[,cols]
    temp$Sweep <- s
    if(volt_steps){
    temp$V_Clamp <- volt[,grep(s, names(volt), value=T)]
    }
    colnames(temp) <- colnames(new.df)
    new.df <- rbind(new.df,temp)
  }

  if(volt_steps){
    new.df$V_Clamp <- as.numeric(gsub("m", "", new.df$V_Clamp))
    for(cols in colnames(new.df)){
      tryCatch(expr = {
        recoverCol <- new.df[,cols]
        new.df[,cols] <- as.numeric(new.df[,cols])
      }, warning = function(w){
        new.df[,cols] <- new.df[,cols]
      })
    }
  }

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
prepareMultipleDFs <- function(pathDF.list){

  dfs <- lapply(l_files, function(x) {
    df <- prepareDF(x)
    return(df)
  })

  names(dfs) <- l_files
  df <- dplyr::bind_rows(dfs, .id = "column_label")
  head(df)
  df$Plate_ID <- sapply(df$Plate_ID, function(x){
    unlist(stringr::str_split(x, "\\r"))[1]
  })

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
df <- prepareDF(pathDF)

df <- df %>% hablar::retype()
numeric_cols <- names(df)[sapply(df, is.numeric)]
numeric_cols <- numeric_cols[!(numeric_cols %in% c("Sweep", "V_Clamp"))]
description_cols <- colnames (df)[!(colnames (df) %in% numeric_cols)]
description_cols <- description_cols[description_cols != "Sweep"]

assays <- lapply(numeric_cols, \(numeric_cols) {
  x <- reshape2::dcast(df, Well ~ Sweep, value.var = numeric_cols)
  m <- x[ , !(names(x) %in% c("Well", "QC"))] |> as.matrix()
  t(m)
})

names(assays) <- numeric_cols

df <- data.table::as.data.table(df)

cd <- S4Vectors::DataFrame(unique(df[, .(Well, QC, Plate_ID)]))


rd <- S4Vectors::DataFrame(unique(df[, .(Sweep)]))
se <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                           rowData = rd,
                           colData = cd)
colnames(se) <- cd$Well

type(description_cols)

description_cols <- description_cols[!(description_cols %in% names(cd))]

descr <-lapply(description_cols, function(var) {
   x<- reshape2::dcast(df, Sweep ~ Well, value.var = var)
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

return(se)

}
