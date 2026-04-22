#' @importFrom SingleCellExperiment SingleCellExperiment
NULL

#' Plot raw current-clamp traces
#'
#' Reads a raw trace Parquet produced by \code{\link{csvToParquet}} and plots
#' sweeps overlaid, coloured by sweep number (viridis gradient).  A small
#' stimulus inset shows the injected current waveform.  Pass a vector of well
#' IDs to \code{well_id} to get a faceted multi-well overview.
#'
#' @param parquet_path  Path to a raw trace Parquet.
#' @param well_id  Character vector of wells to plot (e.g. \code{c("A01","B02")}).
#'   \code{NULL} plots all wells in the parquet.
#' @param plate_id  Optional plate filter. \code{NULL} (default) includes all
#'   plates in the parquet.
#' @param sweeps    Integer vector of sweep numbers to include. \code{NULL} = all.
#' @param v_is_volts  Multiply \code{recorded} by 1000 (V → mV). Default \code{TRUE}.
#' @param stim_scale  Factor applied to \code{stimulus} for display.
#'   Default \code{1e12} (A → pA).
#' @param stim_unit  Label for the stimulus scalebar. Default \code{"pA"}.
#' @param alpha      Line transparency. Default \code{0.65}.
#' @param linewidth  Line width for voltage traces. Default \code{0.45}.
#' @param palette    Viridis palette name. Default \code{"plasma"}.
#' @param color      Fixed line colour (e.g. \code{"black"}). Overrides
#'   \code{palette} when set. Default \code{NULL}.
#' @param scalebar   Replace main-panel axes with voltage + time scalebars.
#'   Default \code{FALSE}.
#' @param inset_x    Horizontal anchor of the stimulus inset as a 0–1 fraction
#'   of the x range. \code{0} = left, \code{1} = right. Default \code{1}.
#' @param inset_y    Vertical anchor of the stimulus inset as a 0–1 fraction
#'   of the y range. \code{0} = bottom, \code{1} = top. Default \code{0}.
#' @param inset_size Inset width as a fraction of the x range. Default \code{0.27}.
#' @param ncol       Number of columns when faceting multiple wells. Default \code{4}.
#'
#' @return A \code{ggplot} object (single well) or a faceted \code{ggplot}
#'   (multiple wells, no inset).
#'
#' @examples
#' \dontrun{
#' # Single well
#' plotCCTraces("18T37419.parquet", "G10")
#'
#' # Multiple wells — faceted grid
#' plotCCTraces("18T37419.parquet", c("G10", "G11", "H01", "H02"))
#'
#' # Minimal style
#' plotCCTraces("18T37419.parquet", "G10", scalebar = TRUE, color = "black")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs theme_classic
#'   theme_void theme element_blank element_line element_text element_rect
#'   margin annotation_custom ggplotGrob annotate scale_color_viridis_c
#' @export
plotCCTraces <- function(parquet_path,
                         well_id    = NULL,
                         plate_id   = NULL,
                         sweeps     = NULL,
                         v_is_volts = TRUE,
                         stim_scale = 1e12,
                         stim_unit  = "pA",
                         alpha      = 0.65,
                         linewidth  = 0.45,
                         palette    = "plasma",
                         color      = NULL,
                         scalebar   = FALSE,
                         inset_x    = 1,
                         inset_y    = 0,
                         inset_size    = 0.38,
                         inset_alpha    = 0,
                         scalebar_lwd   = 0.6,
                         scalebar_size  = 2.8,
                         scale_main     = c(NA, NA),
                         scale_inset    = NA,
                         ncol          = 4) {

  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' required. Install with: install.packages('arrow')")

  # ── Load data ────────────────────────────────────────────────────────────────
  ds <- arrow::open_dataset(parquet_path)
  if (!is.null(plate_id)) ds <- dplyr::filter(ds, plate_id == !!plate_id)
  if (!is.null(well_id))  ds <- dplyr::filter(ds, well_id  %in% !!well_id)
  df <- dplyr::collect(ds) |> as.data.frame()

  if (nrow(df) == 0) stop("No data found for the requested wells / plate.")
  if (!is.null(sweeps)) df <- df[df$sweep %in% sweeps, ]

  df$v_mV    <- if (v_is_volts) df$recorded * 1000 else df$recorded
  df$stim    <- df$stimulus * stim_scale
  df$sweep_n <- as.numeric(df$sweep)

  facet_mode <- length(unique(df$well_id)) > 1
  sweep_range <- range(df$sweep_n)

  # ── Color setup ───────────────────────────────────────────────────────────────
  if (is.null(color)) {
    color_aes   <- ggplot2::aes(x = t_s, y = v_mV,
                                color = sweep_n, group = factor(sweep_n))
    color_scale <- ggplot2::scale_color_viridis_c(option = palette,
                                                   guide  = "none",
                                                   limits = sweep_range)
    stim_aes    <- ggplot2::aes(x = t_s, y = stim,
                                color = sweep_n, group = factor(sweep_n))
    stim_cscale <- ggplot2::scale_color_viridis_c(option = palette,
                                                   guide  = "none",
                                                   limits = sweep_range)
  } else {
    color_aes   <- ggplot2::aes(x = t_s, y = v_mV, group = factor(sweep_n))
    color_scale <- NULL
    stim_aes    <- ggplot2::aes(x = t_s, y = stim, group = factor(sweep_n))
    stim_cscale <- NULL
  }

  # ── Build main plot ───────────────────────────────────────────────────────────
  p_main <- ggplot2::ggplot(df, color_aes) +
    (if (!is.null(color))
      ggplot2::geom_line(alpha = alpha, linewidth = linewidth, color = color)
    else
      ggplot2::geom_line(alpha = alpha, linewidth = linewidth)) +
    color_scale

  if (facet_mode) {
    p_main <- p_main +
      ggplot2::facet_wrap(~ well_id, ncol = ncol, scales = "free_y")
  }

  protocol <- if ("protocol" %in% names(df)) df$protocol[1] else ""
  title    <- if (facet_mode) {
    if (nzchar(protocol)) protocol else NULL
  } else {
    paste0(df$well_id[1], "  •  ", df$plate_id[1],
           if (nzchar(protocol)) paste0("  •  ", protocol) else "")
  }

  if (scalebar) {
    p_main <- p_main +
      ggplot2::labs(x = NULL, y = NULL, title = title) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 9, hjust = 0.5,
                                             colour = "grey40"),
        strip.text   = ggplot2::element_text(size = 8, colour = "grey30")
      )
  } else {
    p_main <- p_main +
      ggplot2::labs(x = "Time (s)", y = "Voltage (mV)", title = title) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        axis.line  = ggplot2::element_line(linewidth = 0.4),
        plot.title = ggplot2::element_text(size = 9, hjust = 0.5,
                                           colour = "grey40"),
        strip.text = ggplot2::element_text(size = 8, colour = "grey30")
      )
  }

  # Faceted mode: return here — inset doesn't work across facets
  if (facet_mode) return(p_main)

  # ── Stimulus inset (single-well only) ────────────────────────────────────────
  stim_range <- range(df$stim)
  if (!is.na(scale_inset[1])) {
    bar_val <- scale_inset[1]
  } else {
    raw_scale <- diff(stim_range) / 4
    magnitude <- 10 ^ floor(log10(max(raw_scale, 1e-9)))
    bar_val   <- magnitude * round(raw_scale / magnitude)
    if (bar_val <= 0) bar_val <- magnitude
  }

  bar_x  <- max(df$t_s) - diff(range(df$t_s)) * 0.04
  bar_y0 <- stim_range[1] + diff(stim_range) * 0.05
  bar_y1 <- bar_y0 + bar_val

  p_stim <- ggplot2::ggplot(df, stim_aes) +
    (if (!is.null(color))
      ggplot2::geom_line(linewidth = linewidth * 0.75, alpha = 0.85, color = color)
    else
      ggplot2::geom_line(linewidth = linewidth * 0.75, alpha = 0.85)) +
    stim_cscale +
    ggplot2::annotate("segment",
                      x = bar_x, xend = bar_x, y = bar_y0, yend = bar_y1,
                      colour = "black", linewidth = scalebar_lwd) +
    ggplot2::annotate("text",
                      x = bar_x, y = (bar_y0 + bar_y1) / 2,
                      label = paste0(bar_val, " ", stim_unit),
                      hjust = -0.15, size = scalebar_size * 0.65,
                      colour = "black") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.25))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(
        fill   = scales::alpha("white", inset_alpha),
        colour = NA)
    )

  xr       <- range(df$t_s)
  yr       <- range(df$v_mV)
  xw       <- diff(xr) * inset_size
  yh       <- diff(yr) * inset_size * 0.55
  x_anchor <- xr[1] + inset_x * diff(xr)
  y_anchor <- yr[1] + inset_y * diff(yr)
  xmin_i   <- x_anchor - (if (inset_x >= 0.5) xw else 0)
  xmax_i   <- xmin_i + xw
  ymin_i   <- y_anchor - (if (inset_y >= 0.5) yh else 0)
  ymax_i   <- ymin_i + yh

  out <- p_main +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::annotation_custom(
      grob = ggplot2::ggplotGrob(p_stim),
      xmin = xmin_i, xmax = xmax_i,
      ymin = ymin_i, ymax = ymax_i
    )

  # ── Scalebars on main panel ───────────────────────────────────────────────────
  if (scalebar) {
    if (!is.na(scale_main[2])) {
      v_bar <- scale_main[2]
    } else {
      v_raw <- diff(yr) / 5
      v_mag <- 10 ^ floor(log10(max(v_raw, 1e-9)))
      v_bar <- v_mag * round(v_raw / v_mag)
      if (v_bar <= 0) v_bar <- v_mag
    }

    if (!is.na(scale_main[1])) {
      t_bar <- scale_main[1]
    } else {
      t_raw <- diff(xr) / 6
      t_mag <- 10 ^ floor(log10(max(t_raw, 1e-9)))
      t_bar <- t_mag * round(t_raw / t_mag)
      if (t_bar <= 0) t_bar <- t_mag
    }
    t_label <- if (t_bar < 1) paste0(round(t_bar * 1000), " ms") else paste0(t_bar, " s")

    sb_x1 <- xr[2] - diff(xr) * 0.04
    sb_x0 <- sb_x1 - t_bar
    sb_y0 <- yr[1] + diff(yr) * 0.04
    sb_y1 <- sb_y0 + v_bar

    out <- out +
      ggplot2::annotate("segment", x = sb_x0, xend = sb_x1,
                        y = sb_y0, yend = sb_y0,
                        colour = "black", linewidth = scalebar_lwd) +
      ggplot2::annotate("segment", x = sb_x1, xend = sb_x1,
                        y = sb_y0, yend = sb_y1,
                        colour = "black", linewidth = scalebar_lwd) +
      ggplot2::annotate("text", x = (sb_x0 + sb_x1) / 2, y = sb_y0,
                        label = t_label, vjust = 1.6, size = scalebar_size,
                        colour = "black") +
      ggplot2::annotate("text", x = sb_x1, y = (sb_y0 + sb_y1) / 2,
                        label = paste0(v_bar, " mV"), hjust = -0.15,
                        size = scalebar_size, colour = "black")
  }

  out
}

#' Build a current-clamp analysis configuration
#'
#' Constructs a named list of parameters that mirrors the Python
#' \code{EfelConfig} dataclass used by \code{ap_analysis_v2.py}.  Pass the
#' result to \code{\link{runAPAnalysis}}.
#'
#' @param v_is_volts Logical; \code{TRUE} (default) if the recorded signal is
#'   in Volts and should be multiplied by 1000 to obtain mV.
#' @param smooth_ms Gaussian smoothing window in ms before dV/dt computation.
#'   Default \code{0.2}.
#' @param dvdt_thr_mV_per_ms dV/dt threshold for spike detection (mV/ms).
#'   Default \code{25}.
#' @param refractory_ms Minimum inter-spike interval in ms. Default \code{2}.
#' @param peak_search_ms Window after dV/dt crossing to find voltage peak (ms).
#'   Default \code{2}.
#' @param min_peak_voltage_mV Minimum voltage for a detected peak to count as a
#'   spike (mV). Default \code{-20}.
#' @param min_height_mV Minimum spike height (mV). \code{NULL} (default) lets
#'   the noise-adaptive floor take over.
#' @param prominence_mV Minimum spike prominence (mV). \code{NULL} (default).
#' @param baseline_win_s Baseline window duration before stim onset (s).
#'   Default \code{0.05}.
#' @param steady_win_s Steady-state window duration before stim end (s).
#'   Default \code{0.05}.
#' @param after_win_s Post-stimulus window duration (s). Default \code{0.05}.
#' @param after_delay_s Delay after stim end before measuring post-stim window
#'   (s). Default \code{0}.
#' @param stim_frac Fraction of stimulus amplitude used to detect stim window.
#'   Default \code{0.2}.
#' @param stim_min_dur_s Minimum detected stim duration (s). Default \code{0.02}.
#' @param stim_start_ignore_ms Ignore spikes within this window at stim onset
#'   (ms). Default \code{2}.
#' @param noise_win_s Window used to estimate baseline noise (s).
#'   Default \code{0.05}.
#' @param noise_mult_height Noise multiplier for adaptive height threshold.
#'   Default \code{6}.
#' @param noise_mult_prom Noise multiplier for adaptive prominence threshold.
#'   Default \code{5}.
#' @param min_height_floor_mV Hard minimum for noise-adaptive height (mV).
#'   Default \code{6}.
#' @param prominence_floor_mV Hard minimum for noise-adaptive prominence (mV).
#'   Default \code{4}.
#' @param phase_window_after_peak_ms Phase-plane analysis window after spike
#'   peak (ms). Default \code{3}.
#'
#' @return A named list of class \code{"cc_config"}.
#'
#' @examples
#' cfg <- cc_config(dvdt_thr_mV_per_ms = 30, min_height_mV = 30)
#'
#' @export
cc_config <- function(
    v_is_volts              = TRUE,
    smooth_ms               = 0.2,
    dvdt_thr_mV_per_ms      = 25.0,
    refractory_ms           = 2.0,
    peak_search_ms          = 2.0,
    min_peak_voltage_mV     = -20.0,
    min_height_mV           = NULL,
    prominence_mV           = NULL,
    baseline_win_s          = 0.05,
    steady_win_s            = 0.05,
    after_win_s             = 0.05,
    after_delay_s           = 0.0,
    stim_frac               = 0.2,
    stim_min_dur_s          = 0.02,
    stim_start_ignore_ms    = 2.0,
    noise_win_s             = 0.05,
    noise_mult_height       = 6.0,
    noise_mult_prom         = 5.0,
    min_height_floor_mV     = 6.0,
    prominence_floor_mV     = 4.0,
    phase_window_after_peak_ms = 3.0
) {
  cfg <- list(
    v_is_volts              = v_is_volts,
    smooth_ms               = smooth_ms,
    dvdt_thr_mV_per_ms      = dvdt_thr_mV_per_ms,
    refractory_ms           = refractory_ms,
    peak_search_ms          = peak_search_ms,
    min_peak_voltage_mV     = min_peak_voltage_mV,
    min_height_mV           = min_height_mV,
    prominence_mV           = prominence_mV,
    baseline_win_s          = baseline_win_s,
    steady_win_s            = steady_win_s,
    after_win_s             = after_win_s,
    after_delay_s           = after_delay_s,
    stim_frac               = stim_frac,
    stim_min_dur_s          = stim_min_dur_s,
    stim_start_ignore_ms    = stim_start_ignore_ms,
    noise_win_s             = noise_win_s,
    noise_mult_height       = noise_mult_height,
    noise_mult_prom         = noise_mult_prom,
    min_height_floor_mV     = min_height_floor_mV,
    prominence_floor_mV     = prominence_floor_mV,
    phase_window_after_peak_ms = phase_window_after_peak_ms
  )
  class(cfg) <- "cc_config"
  cfg
}


# ── Internal: load ap_analysis_v2 Python module via reticulate ────────────────

.load_ap_module <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required. Install with: install.packages('reticulate')")

  python_dir <- system.file("python", package = "ephacRTools")
  if (!nzchar(python_dir) || !dir.exists(python_dir)) {
    dev_path <- file.path("inst", "python")
    if (dir.exists(dev_path)) python_dir <- normalizePath(dev_path)
  }
  if (!nzchar(python_dir) || !dir.exists(python_dir))
    stop("Could not locate inst/python/. ",
         "Install the package or run from the package root directory.")

  reticulate::py_run_string(
    paste0("import sys; sys.path.insert(0, '",
           gsub("\\\\", "/", python_dir), "')")
  )
  tryCatch(
    reticulate::import("ap_analysis_v2"),
    error = function(e) stop(
      "Failed to import ap_analysis_v2.py: ", e$message,
      "\nEnsure Python has: numpy, pandas, polars, pyarrow, scipy, matplotlib"
    )
  )
}


#' Convert a folder of per-well DataControl CSV files to a single Parquet file
#'
#' Step 1 of the current-clamp pipeline.  Each CSV in \code{input_dir} must
#' follow the DataControl naming convention
#' \code{<date>_<protocol>_<time>_<well>[_<plate>].csv}, e.g.
#' \code{241119_VI 5pA_16.10.11_A01__.csv}.  All CSVs are concatenated into
#' one long-format Parquet with columns \code{date}, \code{protocol},
#' \code{plate_id}, \code{well_id}, \code{sweep}, \code{sweep_time_s},
#' \code{t_s}, \code{recorded}, \code{stimulus}.
#'
#' @param input_dir  Path to a folder containing \code{*.csv} files (one per
#'   well, one folder per plate).
#' @param output_parquet  Destination \code{.parquet} path.  Parent directories
#'   are created automatically.
#' @param plate_id  Optional character string to override the plate ID parsed
#'   from filenames (useful when the plate barcode is not embedded in the
#'   filename).  Defaults to the folder name.
#' @param verbose_every Print progress every N files. \code{0} = silent.
#'   Default \code{20}.
#'
#' @return Path to the written Parquet file (invisibly).
#'
#' @examples
#' \dontrun{
#' csvToParquet(
#'   input_dir      = "data-raw/CC/241119/22W14827/260304",
#'   output_parquet = "data-raw/CC/241119/22W14827/260304.parquet"
#' )
#' }
#'
#' @export
csvToParquet <- function(input_dir,
                         output_parquet,
                         plate_id     = NULL,
                         verbose_every = 20L) {

  ap_mod <- .load_ap_module()

  input_dir      <- normalizePath(input_dir,      mustWork = TRUE)
  output_parquet <- normalizePath(output_parquet, mustWork = FALSE)

  ap_mod$folder_csvs_to_single_parquet(
    input_dir      = input_dir,
    output_parquet = output_parquet,
    plate_id       = if (is.null(plate_id)) reticulate::py_none() else plate_id,
    verbose_every  = as.integer(verbose_every)
  )

  message("Parquet written to: ", output_parquet)
  invisible(output_parquet)
}


#' Run the AP analysis pipeline on a raw trace Parquet file
#'
#' Step 2 of the current-clamp pipeline.  Calls the Python
#' \code{analyze_parquet_iterative()} function on \code{raw_parquet_path} and
#' saves the sweep-level and per-AP Parquet tables to \code{out_dir}.
#'
#' @param raw_parquet_path  Path to a raw trace Parquet file produced by
#'   \code{\link{csvToParquet}} (columns: \code{plate_id}, \code{well_id},
#'   \code{sweep}, \code{t_s}, \code{recorded}, \code{stimulus}).
#' @param out_dir  Directory where output Parquets will be written.  Created
#'   if it does not exist.  Two files are written:
#'   \code{<stem>_sweep.parquet} and \code{<stem>_ap.parquet}.
#' @param cfg  A \code{\link{cc_config}} object.  \code{NULL} uses Python
#'   defaults.
#' @param plates  Optional character vector of plate IDs to process
#'   (subset of plates in the Parquet).
#' @param wells   Optional character vector of well IDs to process.
#' @param verbose_every  Print progress every N wells. Default \code{25}.
#'
#' @return Named character vector with elements \code{sweep} and \code{ap}
#'   giving the paths of the written Parquet files (invisibly).
#'
#' @examples
#' \dontrun{
#' paths <- runAPAnalysis(
#'   raw_parquet_path = "data-raw/CC/260304.parquet",
#'   out_dir          = "data-raw/CC/260304_results",
#'   cfg              = cc_config(dvdt_thr_mV_per_ms = 30)
#' )
#' se_cc <- prepareCCSE(paths["sweep"], ap_parquet_path = paths["ap"])
#' }
#'
#' @export
runAPAnalysis <- function(raw_parquet_path,
                          out_dir,
                          cfg           = NULL,
                          plates        = NULL,
                          wells         = NULL,
                          verbose_every = 25L) {

  ap_mod <- .load_ap_module()

  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' is required. Install with: install.packages('arrow')")

  raw_parquet_path <- normalizePath(raw_parquet_path, mustWork = TRUE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  stem    <- tools::file_path_sans_ext(basename(raw_parquet_path))
  sw_path <- file.path(out_dir, paste0(stem, "_sweep.parquet"))
  ap_path <- file.path(out_dir, paste0(stem, "_ap.parquet"))

  # Build Python EfelConfig
  py_cfg <- if (is.null(cfg)) {
    ap_mod$EfelConfig()
  } else {
    cfg_clean <- cfg
    class(cfg_clean) <- NULL
    cfg_clean <- lapply(cfg_clean, function(x) if (is.null(x)) reticulate::py_none() else x)
    do.call(ap_mod$EfelConfig, cfg_clean)
  }

  message("Running AP analysis on: ", basename(raw_parquet_path))
  ap_mod$analyze_parquet_iterative(
    parquet_path      = raw_parquet_path,
    out_ap_parquet    = ap_path,
    out_sweep_parquet = sw_path,
    cfg               = py_cfg,
    plates            = if (!is.null(plates)) reticulate::r_to_py(plates) else reticulate::py_none(),
    wells             = if (!is.null(wells))  reticulate::r_to_py(wells)  else reticulate::py_none(),
    verbose_every     = as.integer(verbose_every)
  )

  message("Sweep parquet: ", sw_path)
  message("AP parquet:    ", ap_path)

  invisible(c(sweep = sw_path, ap = ap_path))
}


# Column roles in the sweep parquet
.CC_WELL_COLS  <- c("well_id", "plate_id", "date", "protocol", "cell_uid")
.CC_SWEEP_META <- c("sweep", "stim_amp_cmd", "stim_start_s", "stim_end_s",
                    "baseline_t0_s", "baseline_t1_s",
                    "steady_t0_s",   "steady_t1_s",
                    "after_t0_s",    "after_t1_s")
.CC_DROP       <- c("sweep_uid", "ap_uid")


#' Prepare a SingleCellExperiment from pre-computed current-clamp sweep parquets
#'
#' Step 3 of the current-clamp pipeline.  Reads one or more sweep-level Parquet
#' files produced by \code{\link{runAPAnalysis}} and assembles them into a
#' \code{SingleCellExperiment} compatible with \code{\link{buildMAE}}.
#'
#' The output \code{SingleCellExperiment} has:
#' \itemize{
#'   \item \strong{rows} = sweeps (stimulus metadata in \code{rowData}).
#'   \item \strong{cols} = \code{well_id.plate_id} keys, matching
#'     \code{prepareSE()} output.
#'   \item \strong{assays} = \code{baseline_v_mV}, \code{steady_v_mV},
#'     \code{after_v_mV}, \code{ap_count}, \code{ap_freq_hz},
#'     \code{mean_ap_amplitude_mV}, \code{mean_isi_s},
#'     \code{mean_dvdt_max_mV_per_ms}, etc.
#' }
#'
#' @param sweep_parquet_path  Character vector of paths to sweep-level Parquet
#'   files.  Multiple files are row-bound before processing.
#' @param ap_parquet_path  Optional character vector of paths to per-spike
#'   Parquet files.  Attached as \code{metadata(se)$ap_data}.
#' @param progress_callback  Optional function \code{f(i, total, name)} called
#'   after each file is read.
#'
#' @return A \code{SingleCellExperiment}.
#'
#' @examples
#' \dontrun{
#' paths <- runAPAnalysis("260304.parquet", out_dir = "260304_results")
#' se_cc <- prepareCCSE(paths["sweep"], ap_parquet_path = paths["ap"])
#' }
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
#' Extract per-well current-clamp features from a CC SingleCellExperiment
#'
#' Derives a summary feature table from a \code{SingleCellExperiment} produced
#' by \code{\link{prepareCCSE}}.  All features come from data already computed
#' by the Python AP-analysis pipeline — no raw traces are required.
#'
#' \strong{Sweep-parquet features (always computed):}
#' \itemize{
#'   \item \code{resting_vm_mV} — mean baseline voltage at the most
#'     hyperpolarised (zero or smallest) current step.
#'   \item \code{input_resistance_MOhm} — slope of the linear fit
#'     ΔV_steady / ΔI across all subthreshold sweeps (ap_count == 0);
#'     requires ≥ 3 subthreshold sweeps.
#'   \item \code{rheobase_pA} — commanded current at the first sweep with at
#'     least one AP.
#'   \item \code{max_firing_rate_Hz} — peak \code{ap_freq_hz} across all sweeps.
#'   \item \code{fi_slope_Hz_per_pA} — slope of the linear fit
#'     ap_freq_hz ~ stim_amp_cmd across all firing sweeps (ap_count > 0);
#'     requires ≥ 2 firing sweeps.
#'   \item \code{depolarization_block} — logical; \code{TRUE} when ap_count
#'     drops to < 50 \% of its maximum in a later sweep.
#'   \item \code{depolarization_block_stim_pA} — commanded current (pA) at the
#'     first sweep where depolarisation block is detected.
#' }
#'
#' \strong{AP-parquet features (requires \code{ap_parquet_path} or AP data in
#' \code{metadata(se_cc)$ap_data}):}
#' \itemize{
#'   \item \code{first_ap_amplitude_mV} — amplitude of the first AP at
#'     rheobase.
#'   \item \code{time_to_first_spike_ms} — latency from stimulus onset to the
#'     first AP at rheobase (ms).
#'   \item \code{ap_amplitude_tapering} — ratio of the last to first AP
#'     amplitude in the sweep with the most APs (< 1 = adapting, > 1 =
#'     facilitation).
#'   \item \code{isi_adaptation_ratio} — ratio of the last to first ISI in the
#'     best (most APs) sweep.
#' }
#'
#' @param se_cc A \code{SingleCellExperiment} from \code{\link{prepareCCSE}}.
#' @param ap_parquet_path Optional path(s) to per-spike Parquet file(s)
#'   produced by \code{\link{runAPAnalysis}}.  If \code{NULL} and AP data are
#'   stored in \code{metadata(se_cc)$ap_data}, those are used automatically.
#'
#' @return A \code{data.frame} with one row per well and columns
#'   \code{well_id}, \code{plate_id}, plus the features listed above.
#'
#' @examples
#' \dontrun{
#' paths  <- runAPAnalysis("260304.parquet", out_dir = "260304_results")
#' se_cc  <- prepareCCSE(paths["sweep"], ap_parquet_path = paths["ap"])
#' feat   <- extractCCFeatures(se_cc)
#' }
#'
#' @importFrom SummarizedExperiment assay assayNames rowData colData
#' @importFrom S4Vectors metadata
#' @export
extractCCFeatures <- function(se_cc, ap_parquet_path = NULL) {

  # ── 1. Pull sweep metadata and assays ────────────────────────────────────────
  rd <- as.data.frame(SummarizedExperiment::rowData(se_cc))
  cd <- as.data.frame(SummarizedExperiment::colData(se_cc))

  .get <- function(nm) {
    if (nm %in% SummarizedExperiment::assayNames(se_cc))
      SummarizedExperiment::assay(se_cc, nm)
    else
      NULL
  }

  a_baseline <- .get("baseline_v_mV")
  a_steady   <- .get("steady_v_mV")
  a_ap_count <- .get("ap_count")
  a_ap_freq  <- .get("ap_freq_hz")

  stim_A  <- if ("stim_amp_cmd" %in% names(rd)) rd$stim_amp_cmd else NULL
  stim_pA <- if (!is.null(stim_A)) stim_A * 1e12 else NULL

  # ── 2. Per-well sweep features ───────────────────────────────────────────────
  well_keys <- colnames(se_cc)

  sweep_feats <- lapply(well_keys, function(wk) {
    row <- list(well_key = wk)

    # Resting Vm: mean baseline_v_mV at zero/most-negative current
    if (!is.null(a_baseline) && !is.null(stim_pA)) {
      min_stim <- min(stim_pA, na.rm = TRUE)
      zero_idx <- which(stim_pA <= max(0, min_stim + abs(min_stim) * 0.05))
      if (length(zero_idx) == 0) zero_idx <- which.min(abs(stim_pA))
      row$resting_vm_mV <- mean(a_baseline[zero_idx, wk], na.rm = TRUE)
    }

    # Input resistance: ΔV_steady / ΔI across subthreshold sweeps
    # R [MΩ] = (ΔV [mV] / ΔI [pA]) × 1000
    if (!is.null(a_steady) && !is.null(a_ap_count) && !is.null(stim_pA)) {
      sub_idx <- which(a_ap_count[, wk] == 0 & !is.na(stim_pA))
      if (length(sub_idx) >= 3) {
        I_v <- stim_pA[sub_idx]
        V_v <- a_steady[sub_idx, wk]
        ok  <- is.finite(I_v) & is.finite(V_v)
        if (sum(ok) >= 3 && var(I_v[ok]) > 0) {
          fit <- tryCatch(lm(V_v[ok] ~ I_v[ok]), error = function(e) NULL)
          if (!is.null(fit)) {
            row$input_resistance_MOhm <- unname(coef(fit)[2]) * 1000
          }
        }
      }
    }

    # Rheobase: stim_amp_cmd at first firing sweep
    if (!is.null(a_ap_count) && !is.null(stim_pA)) {
      fire_idx <- which(a_ap_count[, wk] > 0 & is.finite(stim_pA))
      if (length(fire_idx) > 0) {
        row$rheobase_pA <- stim_pA[min(fire_idx)]
      }
    }

    # Max firing rate
    if (!is.null(a_ap_freq)) {
      mx <- suppressWarnings(max(a_ap_freq[, wk], na.rm = TRUE))
      row$max_firing_rate_Hz <- if (is.finite(mx)) mx else NA_real_
    }

    # FI slope: linear fit ap_freq_hz ~ stim_amp_pA over firing sweeps
    if (!is.null(a_ap_freq) && !is.null(a_ap_count) && !is.null(stim_pA)) {
      fire_idx <- which(a_ap_count[, wk] > 0)
      I_v <- stim_pA[fire_idx]
      F_v <- a_ap_freq[fire_idx, wk]
      ok  <- is.finite(I_v) & is.finite(F_v)
      if (sum(ok) >= 2 && var(I_v[ok]) > 0) {
        fit <- tryCatch(lm(F_v[ok] ~ I_v[ok]), error = function(e) NULL)
        if (!is.null(fit)) row$fi_slope_Hz_per_pA <- unname(coef(fit)[2])
      }
    }

    # Depolarisation block: ap_count drops to < 50 % of peak in a later sweep
    if (!is.null(a_ap_count) && !is.null(stim_pA)) {
      counts   <- a_ap_count[, wk]
      max_cnt  <- suppressWarnings(max(counts, na.rm = TRUE))
      max_idx  <- which.max(counts)
      row$depolarization_block <- FALSE
      if (is.finite(max_cnt) && max_cnt > 0 && max_idx < length(counts)) {
        later <- counts[(max_idx + 1):length(counts)]
        block_sub <- which(later < max_cnt * 0.5)
        if (length(block_sub) > 0) {
          row$depolarization_block <- TRUE
          block_i <- max_idx + block_sub[1]
          if (!is.null(stim_pA) && is.finite(stim_pA[block_i]))
            row$depolarization_block_stim_pA <- stim_pA[block_i]
        }
      }
    }

    as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
  })

  result <- do.call(function(...) {
    all_nms <- unique(unlist(lapply(list(...), names)))
    rows <- lapply(list(...), function(x) {
      missing_nms <- setdiff(all_nms, names(x))
      if (length(missing_nms))
        x[missing_nms] <- NA
      x[all_nms]
    })
    do.call(rbind, rows)
  }, sweep_feats)

  # Split well_key → well_id + plate_id
  parts           <- strsplit(result$well_key, ".", fixed = TRUE)
  result$well_id  <- vapply(parts, `[`, character(1), 1)
  result$plate_id <- vapply(parts, `[`, character(1), 2)
  result$well_key <- NULL

  # ── 3. AP-parquet features ────────────────────────────────────────────────────
  ap_data <- NULL
  if (!is.null(ap_parquet_path)) {
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("Package 'arrow' required for ap_parquet_path.")
    ap_dfs  <- lapply(ap_parquet_path, function(p)
      tryCatch(arrow::read_parquet(p),
               error = function(e) { warning("Could not read: ", p); NULL }))
    ap_dfs  <- Filter(Negate(is.null), ap_dfs)
    if (length(ap_dfs) > 0) ap_data <- dplyr::bind_rows(ap_dfs)
  } else {
    ap_data <- S4Vectors::metadata(se_cc)$ap_data
  }

  if (!is.null(ap_data) && nrow(ap_data) > 0) {

    rheobase_sweep <- if (!is.null(stim_A) && !is.null(a_ap_count)) {
      # row index (1-based) of rheobase sweep per well
      vapply(well_keys, function(wk) {
        fire_idx <- which(a_ap_count[, wk] > 0 & is.finite(stim_A))
        if (length(fire_idx) > 0 && "sweep" %in% names(rd))
          rd$sweep[min(fire_idx)]
        else
          NA_real_
      }, numeric(1))
    } else {
      setNames(rep(NA_real_, length(well_keys)), well_keys)
    }

    ap_feats <- lapply(well_keys, function(wk) {
      parts <- strsplit(wk, ".", fixed = TRUE)[[1]]
      wid <- parts[1]; pid <- parts[2]
      row <- list(well_id = wid, plate_id = pid)

      well_ap <- ap_data[ap_data$well_id == wid &
                           ap_data$plate_id == pid, , drop = FALSE]
      if (nrow(well_ap) == 0) return(as.data.frame(row, stringsAsFactors = FALSE))

      # -- Rheobase-sweep AP features --
      rheo_sw <- rheobase_sweep[wk]
      if (!is.na(rheo_sw) && "sweep" %in% names(well_ap)) {
        rheo_ap <- well_ap[well_ap$sweep == rheo_sw, , drop = FALSE]
        if (nrow(rheo_ap) >= 1) {
          if ("peak_time_s" %in% names(rheo_ap))
            rheo_ap <- rheo_ap[order(rheo_ap$peak_time_s), , drop = FALSE]
          if ("ap_amplitude_mV" %in% names(rheo_ap))
            row$first_ap_amplitude_mV <- rheo_ap$ap_amplitude_mV[1]
          # time_to_first_spike is a per-sweep value duplicated on every AP row
          if ("time_to_first_spike_ms" %in% names(rheo_ap))
            row$time_to_first_spike_ms <- rheo_ap$time_to_first_spike_ms[1]
        }
      }

      # -- Train features from the sweep with the most APs --
      if ("sweep" %in% names(well_ap)) {
        sweep_cnt <- table(well_ap$sweep)
        best_sw   <- as.integer(names(which.max(sweep_cnt)))
        best_ap   <- well_ap[well_ap$sweep == best_sw, , drop = FALSE]
        if ("peak_time_s" %in% names(best_ap))
          best_ap <- best_ap[order(best_ap$peak_time_s), , drop = FALSE]
        n_ap <- nrow(best_ap)

        if (n_ap >= 2) {
          if ("ap_amplitude_mV" %in% names(best_ap)) {
            amps <- best_ap$ap_amplitude_mV
            row$ap_amplitude_tapering <- amps[n_ap] / amps[1]
          }
          if ("peak_time_s" %in% names(best_ap)) {
            isis <- diff(best_ap$peak_time_s)
            if (length(isis) >= 1)
              row$isi_adaptation_ratio <- isis[length(isis)] / isis[1]
          }
        }
      }

      as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
    })

    # Bind AP features — columns may differ across wells
    all_ap_nms <- unique(unlist(lapply(ap_feats, names)))
    ap_feat_df <- do.call(rbind, lapply(ap_feats, function(x) {
      miss <- setdiff(all_ap_nms, names(x))
      if (length(miss)) x[miss] <- NA
      x[all_ap_nms]
    }))

    result <- merge(result, ap_feat_df, by = c("well_id", "plate_id"), all.x = TRUE)
  }

  # Reorder: identifiers first
  lead <- intersect(c("well_id", "plate_id"), names(result))
  result[, c(lead, setdiff(names(result), lead)), drop = FALSE]
}


prepareCCSE <- function(sweep_parquet_path,
                        ap_parquet_path   = NULL,
                        progress_callback = NULL) {

  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' is required. Install with: install.packages('arrow')")

  # ── 1. Read and combine sweep parquets ──────────────────────────────────────
  n <- length(sweep_parquet_path)
  dfs <- lapply(seq_along(sweep_parquet_path), function(i) {
    p  <- sweep_parquet_path[[i]]
    df <- tryCatch(
      arrow::read_parquet(p),
      error = function(e) {
        warning("Failed to read '", basename(p), "': ", e$message)
        NULL
      }
    )
    if (!is.null(progress_callback)) progress_callback(i, n, basename(p))
    df
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) stop("No sweep parquet files could be read.")
  df <- dplyr::bind_rows(dfs)

  # ── 2. Identify column roles ─────────────────────────────────────────────────
  non_assay <- c(.CC_WELL_COLS, .CC_SWEEP_META, .CC_DROP)
  assay_cols <- setdiff(
    names(df)[vapply(df, is.numeric, logical(1))],
    non_assay
  )
  if (length(assay_cols) == 0)
    stop("No numeric assay columns found in sweep parquet.")

  # ── 3. Primary column key ────────────────────────────────────────────────────
  df$`.col_key` <- paste(df$well_id, df$plate_id, sep = ".")

  # ── 4. Build assay matrices [sweeps × wells] ─────────────────────────────────
  assay_list <- lapply(assay_cols, function(col) {
    wide <- reshape2::dcast(df, sweep ~ `.col_key`, value.var = col,
                            fun.aggregate = mean)
    sweep_order <- wide$sweep
    mat <- as.matrix(wide[, -1, drop = FALSE])
    rownames(mat) <- paste0("Sweep", sweep_order)
    mat
  })
  names(assay_list) <- assay_cols

  col_names <- colnames(assay_list[[1]])

  # ── 5. rowData ────────────────────────────────────────────────────────────────
  rd_cols <- intersect(.CC_SWEEP_META, colnames(df))
  rd_raw  <- unique(as.data.frame(df[, rd_cols, drop = FALSE]))
  rd_raw  <- rd_raw[order(rd_raw$sweep), , drop = FALSE]
  rownames(rd_raw) <- paste0("Sweep", rd_raw$sweep)
  rd <- S4Vectors::DataFrame(rd_raw)

  # ── 6. colData ────────────────────────────────────────────────────────────────
  cd_cols <- intersect(.CC_WELL_COLS, colnames(df))
  cd_raw  <- as.data.frame(df[!duplicated(df$`.col_key`), c(".col_key", cd_cols), drop = FALSE])
  cd_raw  <- cd_raw[match(col_names, cd_raw$`.col_key`), , drop = FALSE]
  rownames(cd_raw) <- col_names
  cd_raw$`.col_key` <- NULL
  cd <- S4Vectors::DataFrame(cd_raw)

  # ── 7. Assemble SCE ───────────────────────────────────────────────────────────
  se <- SingleCellExperiment::SingleCellExperiment(
    assays  = assay_list,
    rowData = rd,
    colData = cd
  )
  colnames(se) <- col_names

  # ── 8. Attach per-spike data as metadata (optional) ──────────────────────────
  if (!is.null(ap_parquet_path)) {
    ap_dfs <- lapply(seq_along(ap_parquet_path), function(i) {
      p  <- ap_parquet_path[[i]]
      tryCatch(
        arrow::read_parquet(p),
        error = function(e) {
          warning("Failed to read AP parquet '", basename(p), "': ", e$message)
          NULL
        }
      )
    })
    ap_dfs <- Filter(Negate(is.null), ap_dfs)
    if (length(ap_dfs) > 0) {
      ap_combined <- dplyr::bind_rows(ap_dfs)
      S4Vectors::metadata(se)$ap_data <- ap_combined
    }
  }

  se
}


#' Run AP detection on a single well for interactive preview
#'
#' Loads one well from a raw trace Parquet produced by
#' \code{\link{csvToParquet}} and runs the same AP-detection pipeline as
#' \code{\link{runAPAnalysis}} — but returns the results directly instead of
#' writing Parquet files.  Intended for the \code{ccPreviewApp} Shiny app and
#' other interactive use.
#'
#' @param parquet_path  Path to a raw trace Parquet file.
#' @param well_id  Well identifier (e.g. \code{"G10"}).
#' @param plate_id  Plate identifier (e.g. \code{"18T37419"}).
#' @param sweeps  Integer vector of sweep numbers to analyse.  \code{NULL}
#'   (default) analyses all sweeps for the well.
#' @param cfg  A \code{\link{cc_config}} object.  \code{NULL} uses Python
#'   defaults.
#'
#' @return A named list with two \code{data.frame}s:
#' \describe{
#'   \item{\code{sweeps}}{One row per sweep: \code{sweep}, \code{stim_amp_cmd},
#'     \code{ap_count}, \code{ap_freq_hz}, \code{baseline_v_mV},
#'     \code{steady_v_mV}, etc.}
#'   \item{\code{aps}}{One row per detected AP: \code{sweep},
#'     \code{peak_time_s}, \code{peak_v_mV}, \code{ap_amplitude_mV}, etc.
#'     Empty \code{data.frame} when no APs are detected.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- previewAPDetection("18T37419.parquet", "G10", "18T37419")
#' result$sweeps  # per-sweep summary
#' result$aps     # individual APs
#' }
#'
#' @export
previewAPDetection <- function(parquet_path, well_id, plate_id,
                               sweeps = NULL, cfg = NULL) {
  ap_mod       <- .load_ap_module()
  parquet_path <- normalizePath(parquet_path, mustWork = TRUE)

  py_cfg <- if (is.null(cfg)) {
    ap_mod$EfelConfig()
  } else {
    cfg_clean        <- cfg
    class(cfg_clean) <- NULL
    cfg_clean        <- lapply(cfg_clean, function(x)
      if (is.null(x)) reticulate::py_none() else x)
    do.call(ap_mod$EfelConfig, cfg_clean)
  }

  result <- ap_mod$detect_aps_preview(
    parquet_path = parquet_path,
    well_id      = well_id,
    plate_id     = plate_id,
    sweeps       = if (!is.null(sweeps))
      reticulate::r_to_py(as.integer(sweeps))
    else
      reticulate::py_none(),
    cfg          = py_cfg
  )

  list(
    sweeps = as.data.frame(result[["sweeps"]]),
    aps    = as.data.frame(result[["aps"]])
  )
}


#' Quick scan of a raw trace Parquet for AP-positive wells
#'
#' Computes the maximum recorded voltage and sweep count per well using Arrow
#' lazy evaluation (pure R — no Python required).  Wells where the maximum
#' voltage exceeds \code{threshold_mV} are flagged as \code{ap_likely = TRUE}.
#'
#' @param parquet_path  Path (or character vector of paths) to raw trace
#'   Parquet file(s).
#' @param v_is_volts  Logical; multiply \code{recorded} by 1000 to convert
#'   V → mV.  Default \code{TRUE}.
#' @param threshold_mV  Voltage threshold (mV) above which a well is considered
#'   AP-likely.  Default \code{0}.
#'
#' @return A \code{data.frame} with columns \code{plate_id}, \code{well_id},
#'   \code{parquet_path}, \code{max_v_mV}, \code{n_sweeps}, \code{ap_likely};
#'   sorted AP-likely wells first (descending max_v_mV).
#'
#' @importFrom dplyr group_by summarise collect n_distinct bind_rows filter
#' @export
scanParquetForAPs <- function(parquet_path,
                              v_is_volts   = TRUE,
                              threshold_mV = 0,
                              sweeps       = NULL) {
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' required. Install with: install.packages('arrow')")

  scale      <- if (isTRUE(v_is_volts)) 1000 else 1
  sweep_ints <- if (!is.null(sweeps) && length(sweeps) > 0)
    as.integer(sweeps[is.finite(as.numeric(sweeps))]) else NULL

  results <- lapply(parquet_path, function(p) {
    p <- normalizePath(p, mustWork = TRUE)
    tryCatch({
      ds <- arrow::open_dataset(p)
      if (!is.null(sweep_ints))
        ds <- dplyr::filter(ds, sweep %in% sweep_ints)
      df <- ds |>
        dplyr::mutate(v_mV = recorded * scale) |>
        dplyr::group_by(plate_id, well_id, sweep) |>
        dplyr::summarise(
          sweep_max_v_mV = max(v_mV, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::group_by(plate_id, well_id) |>
        dplyr::summarise(
          max_v_mV    = max(sweep_max_v_mV, na.rm = TRUE),
          n_sweeps    = dplyr::n(),
          n_ap_sweeps = sum(sweep_max_v_mV > threshold_mV, na.rm = TRUE),
          .groups     = "drop"
        ) |>
        dplyr::collect() |>
        as.data.frame()
      df$ap_likely    <- df$n_ap_sweeps > 0L
      df$parquet_path <- p
      df
    }, error = function(e) {
      warning("Could not scan '", basename(p), "': ", e$message)
      NULL
    })
  })

  out <- dplyr::bind_rows(Filter(Negate(is.null), results))
  if (nrow(out) == 0) return(out)
  out[order(-out$n_ap_sweeps, -out$max_v_mV), ]
}
