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
<<<<<<< HEAD
#' @return A filtered data frame
#' @export
filtered_df <- function(prepared_df, plate_ID) {
=======
#' @param ion The compound that was added or removed during the measurements
#' @return A filtered data frame
#' @export
filtered_df <- function(prepared_df, plate_ID, ion = c("Na", "K"), columns = c("Well", "QC", "Compound", "Conditions", "Plate_ID")) {
>>>>>>> fceaa49 (Updated scripts)
  # Filter condition: rows B–O and cols 2–23
  valid_rows <- prepared_df$Row %in% LETTERS[2:15]   # B (2) to O (15)
  valid_cols <- prepared_df$Column >= 2 & prepared_df$Column <= 23
  valid_indices <- valid_rows & valid_cols

  filtered_wells <- prepared_df$Well[valid_indices]

  prepared_df <- prepared_df %>%
    filter(Plate_ID == plate_ID) %>%
    filter(Well %in% filtered_wells)

<<<<<<< HEAD
  prepared_df$Compound <- ifelse(prepared_df$RowNum %in% 1:4 | prepared_df$RowNum %in% 9:12, "Na Addition", "Na Removal")

  prepared_df <- ag(subset(prepared_df), cols= c("Well", "QC", "Compound", "Conditions", "Plate_ID"))
=======
  if (ion == "Na") {
    prepared_df$Compound <- ifelse(prepared_df$RowNum %in% 1:4 | prepared_df$RowNum %in% 9:12, "Na Addition", "Na Removal")
  } else if (ion == "K") {
      mapping <- prepared_df %>%
        select(Well, Plate_ID, Compound) %>%
        distinct()
      prepared_df <- prepared_df %>%
        left_join(mapping, by = c("Well", "Plate_ID"), suffix = c("", ".mapped"), relationship = "many-to-many") %>%
        mutate(Compound = Compound.mapped) %>%
        select(-Compound.mapped)
    }


  prepared_df <- ag(subset(prepared_df), cols= columns)
>>>>>>> fceaa49 (Updated scripts)

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
<<<<<<< HEAD
group_assignment <- function(prepared_df, pattern = c("Conditions", "Alternating", "Block", "Manual", "Cycle"),
                             manual_map = NULL, block_size = NULL, cycle_pattern = NULL) {
  pattern <- match.arg(pattern)
=======
group_assignment <- function(data, se = c("Yes", "No"), pattern = c("Conditions", "Alternating", "Block", "Manual", "Cycle"),
                             manual_map = NULL, block_size = NULL, cycle_pattern = NULL) {
  se <- match.arg(se)
  pattern <- match.arg(pattern)

  if (se == "Yes") {
    se_obj <- data
    prepared_df <- as.data.frame(colData(se_obj))
  } else {
    prepared_df <- data
  }

>>>>>>> fceaa49 (Updated scripts)
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
<<<<<<< HEAD
    prepared_df <- manual_map[as.character(prepared_df$Well)]
  }
  else if (pattern == "Cycle") {
    if (is.null(cycle_pattern)) stop("Please provide cycle_pattern for 'cycle' pattern")
    prepared_df$Group <- rep(cycle_pattern, length.out = nrow(prepared_df))
  }

return(prepared_df)

=======
    prepared_df$Group <- manual_map[as.character(prepared_df$Well)]
  }
  else if (pattern == "Cycle") {
    if (is.null(cycle_pattern)) stop("Please provide cycle_pattern for 'cycle' pattern")
    prepared_df <- prepared_df[order(prepared_df$Well), ]
    prepared_df$Group <- rep(cycle_pattern, length.out = nrow(prepared_df))
  }

  if (se == "Yes") {
    colData(se_obj)$Group  <- prepared_df$Group
    return(se_obj)
  } else {
    return(prepared_df)
  }
>>>>>>> fceaa49 (Updated scripts)
}
#' Create a data frame with the feature you want to add to your summarized experiment's colData
#' @param se The summarized experiment you have created previously
#' @param assay_name The metric / features you want to add
#' @param FUN The method you would like to aggregate by
#' @param na.rm Remove NAs
#' @param periods Define your liquid periods
#' @param mode Aggregation mode by period or by last n sweeps
#' @param last.n The value of n if aggregating by the last n sweeps
#' @return A data frame with the added feature in colData
#' @export
rowAG <- function(
    se,
    assay_name = "Erev",  # which assay to aggregate
    FUN = mean,           # aggregation function (mean, median, sd, …)
    na.rm = TRUE,
    # if you want period‐based: supply a named list of sweep ranges
    periods = list(
      Add1 = 1:24,
      Add2 = 25:42,
      Add3 = 43:60,
      Add4 = 61:84
    ),
    # OR if you prefer “last-n”-Sweeps mode:
    mode = c("by.period", "last.n"),
    last.n = 3
) {
  mode <- match.arg(mode)

  mat <- assay(se, assay_name) * 1000
  rd  <- as.data.frame(rowData(se))
  cd  <- as.data.frame(colData(se))

  df_long <- mat %>%
    as.data.frame() %>%
    mutate(Sweep = rd$Sweep) %>%
    pivot_longer(-Sweep, names_to = "Well", values_to = assay_name) %>%
    left_join(cd, by = "Well")

  if (mode == "last.n") {
    # ─── take last n sweeps per Well ──────────────────────────────
    out <- df_long %>%
      group_by(Well, Plate_ID, QC, Compound, Group) %>%
      arrange(Sweep) %>%
      slice_tail(n = last.n) %>%
      summarise(
        !!paste0(assay_name, "_last", last.n) := FUN(.data[[assay_name]], na.rm = na.rm),
        .groups = "drop"
      )
    return(as.data.frame(out))
  }

  # ─── otherwise: by.period mode ────────────────────────────────
  # build a little table mapping each Sweep → period label
  period_map <- bind_rows(
    lapply(names(periods), function(prd) {
      data.frame(Sweep = periods[[prd]], period = prd, stringsAsFactors = FALSE)
    })
  )

  out <- df_long %>%
    inner_join(period_map, by = "Sweep") %>%
    group_by(Well, Plate_ID, QC, Compound, Group, period) %>%
    summarise(
      agg = FUN(.data[[assay_name]], na.rm = na.rm),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = period, values_from = agg)

  return(as.data.frame(out))
}
#' Create a SE with assays for the additions (Add1 - Add4)
#' @param se The summarized experiment you have created previously
#' @param assay_name The metric / features you want to have as assays
#' @param FUN The method you would like to aggregate by
#' @param na.rm Remove NAs
#' @param periods Define your liquid periods
#' @param mode Aggregation mode by period or by last n sweeps
#' @param last.n The value of n if aggregating by the last n sweeps
#' @return A summarized experiment with the assays you set
#' @export
row_aggregate_SE <- function(
    se,
    metrics   = c("Erev"),
    FUN       = mean,
    na.rm     = TRUE,
    periods   = list(
      Add1 = 1:24,
      Add2 = 25:42,
      Add3 = 43:60,
      Add4 = 61:84
    ),
    mode      = c("by.period","last.n"),
    last.n    = 3
) {
  mode <- match.arg(mode)

  # 1) grab metadata
  rd_in <- as.data.frame(rowData(se))
  rd_in$Sweep <- as.integer(rd_in$Sweep)

  cd_in <- as.data.frame(colData(se))
  # drop any old metric/period cols
  drop <- intersect(colnames(cd_in), c(metrics, names(periods)))
  if (length(drop)) cd_in <- cd_in[, setdiff(colnames(cd_in), drop), drop=FALSE]
  # set rownames to well IDs
  if ("Well" %in% colnames(cd_in)) rownames(cd_in) <- cd_in$Well
  else                          rownames(cd_in) <- colnames(se)

  # 2) sweep→period lookup
  period_map <- bind_rows(lapply(names(periods), function(lbl){
    data.frame(Sweep=as.integer(periods[[lbl]]),
               period=lbl,
               stringsAsFactors=FALSE)
  }))

  # 3) build one wells×periods matrix per metric
  assays_list <- lapply(metrics, function(met){
    mat <- assay(se, met)                # sweeps × wells

    df_long <- as.data.frame(mat) %>%
      mutate(Sweep = rd_in$Sweep) %>%
      pivot_longer(-Sweep, names_to="Well", values_to=met) %>%
      left_join(cd_in, by="Well") %>%
      inner_join(period_map, by="Sweep")

    # choose grouping strategy
    if (mode=="last.n") {
      df_p <- df_long %>%
        group_by(Well, period) %>%
        arrange(Sweep) %>%
        slice_tail(n = last.n) %>%
        summarise(val = FUN(.data[[met]], na.rm=na.rm), .groups="drop")
    } else {
      df_p <- df_long %>%
        group_by(Well, period) %>%
        summarise(val = FUN(.data[[met]], na.rm=na.rm), .groups="drop")
    }

    mat_wide <- df_p %>%
      pivot_wider(names_from=period, values_from=val) %>%
      arrange(Well)

    out <- as.matrix(mat_wide[,-1,drop=FALSE])  # wells × periods
    rownames(out) <- mat_wide$Well
    colnames(out) <- names(periods)
    out
  })
  names(assays_list) <- metrics

  # 4) transpose so rows=periods, cols=wells
  assays_flipped <- lapply(assays_list, t)

  # 5) new rowData = the 4 periods
  rn <- names(periods)
  new_rd <- DataFrame(period = rn, row.names = rn)

  # 6) new colData = well metadata in correct order
  wells  <- cd_in$Well
  new_cd <- DataFrame(cd_in[wells, , drop=FALSE], row.names = wells)

  # 7) assemble
  SummarizedExperiment(
    assays  = assays_flipped,
    rowData = new_rd,
    colData = new_cd
  )
}
