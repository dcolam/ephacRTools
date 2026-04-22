#' @importFrom MultiAssayExperiment MultiAssayExperiment
NULL

#' Assemble a MultiAssayExperiment from SCE objects with automatic shared colData
#'
#' Builds a \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} from a
#' named list of \code{SingleCellExperiment} objects.  The shared
#' \code{colData} is assembled automatically by joining well-level metadata
#' across experiments.
#'
#' \strong{Column promotion logic:}
#' \enumerate{
#'   \item Columns listed in \code{shared_cols} are always pulled from
#'     whichever experiment has them (first experiment wins for conflicts).
#'   \item When \code{auto_shared = TRUE}, any additional column that (a)
#'     appears in \emph{all} experiments and (b) has consistent values for
#'     every well is also promoted to shared \code{colData}.
#'   \item All remaining columns stay in each experiment's own \code{colData}
#'     and are accessible via \code{mae[["exp_name"]]}.
#' }
#'
#' @param experiments Named list of \code{SummarizedExperiment} /
#'   \code{SingleCellExperiment} objects.  Each must have \code{well_id} and
#'   \code{plate_id} in its \code{colData}.
#' @param coldata Optional \code{data.frame} or \code{DataFrame} of shared
#'   well-level metadata with \code{well_id.plate_id} row names.  When
#'   supplied, \code{shared_cols} and \code{auto_shared} are ignored.
#' @param shared_cols Character vector of column names to always promote to
#'   shared \code{colData}.  Default:
#'   \code{c("well_id", "plate_id", "row", "col")}.
#' @param auto_shared Logical; if \code{TRUE} (default), columns present in
#'   all experiments with consistent per-well values are automatically
#'   promoted to shared \code{colData}.
#' @param sep Separator for the \code{well_id.plate_id} primary key.
#'   Default \code{"."}.
#'
#' @return A \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}.
#'
#' @examples
#' \dontrun{
#' se_vc <- prepareSE("vc.xlsx")
#' se_cc <- prepareCCSE("sweep.parquet")
#'
#' # Auto-join shared columns
#' mae <- buildMAE(list(voltage_clamp = se_vc, current_clamp = se_cc))
#'
#' # Inspect what ended up shared
#' colnames(colData(mae))
#'
#' # Access experiment-specific colData
#' colnames(colData(mae[["voltage_clamp"]]))
#' }
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @export
buildMAE <- function(experiments,
                     coldata     = NULL,
                     shared_cols = c("well_id", "plate_id", "row", "col"),
                     auto_shared = TRUE,
                     sep         = ".") {

  if (!is.list(experiments) || is.null(names(experiments)) ||
      any(!nzchar(names(experiments))))
    stop("`experiments` must be a named list.")

  # ── Validate each experiment has well_id + plate_id ─────────────────────────
  for (nm in names(experiments)) {
    cd <- as.data.frame(SummarizedExperiment::colData(experiments[[nm]]))
    if (!all(c("well_id", "plate_id") %in% colnames(cd)))
      stop("Experiment '", nm, "' colData must contain 'well_id' and 'plate_id'.")
  }

  .primary_key <- function(se) {
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    paste(cd$well_id, cd$plate_id, sep = sep)
  }

  all_keys <- unique(unlist(lapply(experiments, .primary_key)))

  # ── Build shared colData ─────────────────────────────────────────────────────
  if (!is.null(coldata)) {
    # User-supplied: align to all_keys, fill missing with NA
    coldata <- as.data.frame(coldata)
    missing <- setdiff(all_keys, rownames(coldata))
    if (length(missing) > 0) {
      message("buildMAE: ", length(missing),
              " well(s) not in supplied coldata — filling with NA.")
      extra <- as.data.frame(
        matrix(NA, nrow = length(missing), ncol = ncol(coldata),
               dimnames = list(missing, colnames(coldata)))
      )
      coldata <- rbind(coldata[intersect(rownames(coldata), all_keys), ,
                               drop = FALSE], extra)
    }
    shared_cd <- S4Vectors::DataFrame(coldata[all_keys, , drop = FALSE])

  } else {
    # ── Auto-build from experiments ────────────────────────────────────────────

    # Collect one colData data.frame per experiment, keyed by primary key
    cds <- lapply(experiments, function(se) {
      df       <- as.data.frame(SummarizedExperiment::colData(se))
      df$.key  <- paste(df$well_id, df$plate_id, sep = sep)
      df[!duplicated(df$.key), ]
    })

    # Start with the baseline: well_id + plate_id for every key
    base_df <- data.frame(
      well_id  = sub(paste0(sep, ".*$"),           "", all_keys, perl = TRUE),
      plate_id = sub(paste0("^[^", sep, "]*", sep), "", all_keys, perl = TRUE),
      row.names = all_keys, stringsAsFactors = FALSE
    )

    # Collect additional shared columns
    promote <- unique(c(shared_cols, "well_id", "plate_id"))

    for (nm in names(cds)) {
      df   <- cds[[nm]]
      cols <- intersect(setdiff(promote, colnames(base_df)), colnames(df))
      if (length(cols) == 0) next
      m <- match(all_keys, df$.key)
      for (col in cols) base_df[[col]] <- df[[col]][m]
    }

    # ── Auto-promote columns consistent across all experiments ─────────────────
    if (auto_shared && length(experiments) > 1) {
      # Columns in ALL experiments (excluding already-promoted and key helpers)
      already <- c(colnames(base_df), ".key")
      common  <- Reduce(intersect, lapply(cds, function(df)
        setdiff(colnames(df), already)))

      for (col in common) {
        # Build a key→value map from each experiment; check consistency
        vals_list <- lapply(cds, function(df) {
          idx <- match(all_keys, df$.key)
          setNames(df[[col]][idx], all_keys)
        })
        # A column is consistent if all experiments agree on its value per key
        # (NA counts as "no data" — use the first non-NA value)
        merged_vals <- vals_list[[1]]
        consistent  <- TRUE
        for (v in vals_list[-1]) {
          conflict <- !is.na(merged_vals) & !is.na(v) & merged_vals != v
          if (any(conflict, na.rm = TRUE)) { consistent <- FALSE; break }
          # Fill NAs from subsequent experiments
          na_idx <- is.na(merged_vals)
          merged_vals[na_idx] <- v[na_idx]
        }
        if (consistent) {
          base_df[[col]] <- merged_vals
          message("buildMAE: auto-promoted '", col, "' to shared colData.")
        }
      }
    }

    shared_cd <- S4Vectors::DataFrame(base_df)
  }

  rownames(shared_cd) <- all_keys

  # ── Build sampleMap ──────────────────────────────────────────────────────────
  sm_df <- do.call(rbind, lapply(names(experiments), function(nm) {
    se   <- experiments[[nm]]
    data.frame(assay   = nm,
               primary = .primary_key(se),
               colname = colnames(se),
               stringsAsFactors = FALSE)
  }))
  sm_df$assay <- factor(sm_df$assay, levels = names(experiments))

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = experiments,
    colData     = shared_cd,
    sampleMap   = S4Vectors::DataFrame(sm_df)
  )
}
