#' @importFrom magrittr %>%
NULL
#' Prepare the excel files from the Syncropatch Online Analysis
#' @param pathToExcel Path to the excel file in question
#' @return A data frame
#' @export
prep_excel <- function(pathToExcel){

  df <- as.data.frame(read_excel(pathToExcel, sheet="OA Export", col_types = "text"))
  df$`\r` <- NULL
  names(df)[1:2] <- c("Well", "QC")
  df

  df <- df[-1,]
  sweeps <- grep("Sweep \\d", colnames(df), value=TRUE)
  no.sweeps <- unique(sapply(sweeps, FUN=function(s){
    unlist(str_split(s, " "))[2]
  }))

  new.cols <- sapply(grep(no.sweeps[1], sweeps, value=T), function(x){
    unlist(str_split(x, " "))[3]
  })
  new.cols <- c("Well", "QC","Plate_ID", new.cols, "Sweep")
  new.df <- data.frame(matrix(ncol=length(new.cols),nrow=0, dimnames=list(NULL, new.cols)))

  for(s in no.sweeps){
    cols <- c("Well", "QC", "Nanion Chip Barcode", grep(s, sweeps, value=T))
    temp <- df[,cols]
    temp$Sweep <- s

    colnames(temp) <- colnames(new.df)
    new.df <- rbind(new.df,temp)
  }

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
#' Prepare the newly created data frame from the excel file: incorporation of the conditions, row and column columns etc
#' @param wd_path Path to the excel files location (including the meta data file)
#' @param df The data frame that you just created from the excel file you are looking at
#' @return An adjusted data frame
#' @export
prep_df <- function(wd_path, df) {
  l_files <- list.files(path = wd_path ,pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)

  l_files <- l_files[!grepl("~", l_files)]
  meta <- l_files[grepl("Experiment", l_files)]
  l_files <- l_files[!grepl("Experiment", l_files)]
  length(l_files)

  df$Plate_ID <- sapply(df$Plate_ID, function(x){
    unlist(str_split(x, "\\r"))[1]
  })

  meta <- as.data.frame(read_excel(meta))

  df <- merge(df, meta, by="Plate_ID")

  df$Row <- sapply(df$Well, function(x){
    str_sub(x, 1, 1)
  })

  df$Column <- sapply(df$Well, function(x){
    as.numeric(str_sub(x, 2, 3))
  })

  letter2number <- function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}

  df$RowNum <- sapply(df$Row, function(c){letter2number(c)})

  return(df)

}
ag <- function(df, cols, fun=mean) {

  c <- colnames(df[,unlist(lapply(df, is.numeric))])
  print(c)
  df <- aggregate(df[,c], by=as.list(df[,cols]), FUN=fun, na.rm=TRUE)

  return(df[, colSums(is.na(df)) != nrow(df)])
}
#' Filter your newly adjusted data frame so that it doesn't include the "empty" wells
#' @param prepared_df The data frame that you just adjusted in prep_df
#' @param plate_ID The plate ID that you are focusing on
#' @return A filtered data frame
#' @export
filtered_df <- function(prepared_df, plate_ID) {
  # Filter condition: rows B–O and cols 2–23
  valid_rows <- prepared_df$Row %in% LETTERS[2:15]   # B (2) to O (15)
  valid_cols <- prepared_df$Column >= 2 & prepared_df$Column <= 23
  valid_indices <- valid_rows & valid_cols

  filtered_wells <- prepared_df$Well[valid_indices]

  prepared_df <- prepared_df %>%
    filter(Plate_ID == plate_ID) %>%
    filter(Well %in% filtered_wells)

  prepared_df$Compound <- ifelse(prepared_df$RowNum %in% 1:4 | prepared_df$RowNum %in% 9:12, "Na Addition", "Na Removal")

  prepared_df <- ag(subset(prepared_df), cols= c("Well", "QC", "Compound", "Conditions", "Plate_ID"))

  return(prepared_df)

}
#' Assign your wells to their corresponding groups (wildtype, hyperplasia, donors etc)
#' @param prepared_df The data frame that you just filtered in filtered_df
#' @param pattern The pattern by which you separated your groups in the experiment
#' @param manual_map A named vector like c("A01" = "HP", "A02" = "WT", ...)
#' @param block_size Parameter that determines after how many rows the condition alternates from HP to WT
#' @param cycle_pattern A vector containing a custom, repeating sequence like the triplet c("HP", "WT", "HP")
#' @return A data frame with group assignment
#' @export
group_assignment <- function(prepared_df, pattern = c("Conditions", "Alternating", "Block", "Manual", "Cycle"),
                             manual_map = NULL, block_size = NULL, cycle_pattern = NULL) {
  pattern <- match.arg(pattern)
  prepared_df$Group <- NA   #initialize the column

  if (pattern == "Conditions") {
    prepared_df$Group <- ifelse(df$Conditions %in% c("Donor", "Donor1", "Donor2"), "WT", "HP")
  }
  else if (pattern == "Alternating") {
    prepared_df$Group <- ifelse(seq_len(nrow(prepared_df)) %% 2 == 0, "WT", "HP")
  }
  else if (pattern == "Block") {
    if (is.null(block_size)) stop("block_size must be provided for pattern = 'Block'")
    reps <- rep(c("HP", "WT"), length.out = ceiling(nrow(prepared_df) / block_size))
    prepared_df$Group <- rep(reps, each = block_size, length.out = nrow(prepared_df))
  }
  else if (pattern == "Manual") {
    if (is.null(manual_map)) stop("manual_map must be provided for pattern = 'manual'")
    prepared_df <- manual_map[as.character(prepared_df$Well)]
  }
  else if (pattern == "Cycle") {
    if (is.null(cycle_pattern)) stop("Please provide cycle_pattern for 'cycle' pattern")
    prepared_df$Group <- rep(cycle_pattern, length.out = nrow(prepared_df))
  }

return(prepared_df)

}
