# ccPreviewApp.R — Interactive AP detection preview
#
# Input:  raw trace Parquets produced by csvToParquet()
#         (columns: plate_id, well_id, sweep, t_s, recorded, stimulus)
# NOT the sweep/AP parquets from runAPAnalysis() — those are the output.
#
# Launch:
#   shiny::runApp("app/ccPreviewApp.R")    # from package root
#   shiny::runApp("ccPreviewApp.R")        # from app/ directory

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinyFiles)
  library(ggplot2)
  library(dplyr)
  library(arrow)
  library(jsonlite)
})

# Load package (dev vs installed)
if (file.exists("DESCRIPTION")) {
  suppressMessages(devtools::load_all("."))
} else if (file.exists("../DESCRIPTION")) {
  suppressMessages(devtools::load_all(".."))
} else {
  library(ephacRTools)
}

# ── Helpers ───────────────────────────────────────────────────────────────────

.roots <- c(Home = path.expand("~"), `C:` = "C:/", `D:` = "D:/")

.make_well_key  <- function(well_id, plate_id, path)
  paste(well_id, plate_id, path, sep = "|||")

.parse_well_key <- function(v) {
  parts <- strsplit(v, "|||", fixed = TRUE)[[1]]
  list(well_id = parts[1], plate_id = parts[2], path = parts[3])
}

.opt_val <- function(x)
  if (is.numeric(x) && length(x) == 1 && !is.na(x)) x else NULL

.well_rc <- function(well_id) {
  list(row = match(substr(well_id, 1, 1), LETTERS),
       col = as.integer(substr(well_id, 2, 3)))
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = "AP Detection Preview",
  theme = bs_theme(bootswatch = "flatly"),

  sidebar = sidebar(
    width = 340,

    # ── Step 1: CSV → Parquet ─────────────────────────────────────────────────
    accordion(
      id = "acc_csv",
      open = FALSE,
      accordion_panel(
        "Step 1 — CSV → Parquet",
        helpText("Convert DataControl CSVs (one file per well) to a single Parquet."),
        shinyDirButton("csv_in_dir",  "Input folder (CSVs)",
                       title = "Folder containing per-well CSV files",
                       buttonType = "default btn-sm w-100"),
        verbatimTextOutput("csv_in_label", placeholder = TRUE),
        shinyDirButton("csv_out_dir", "Output folder",
                       title = "Where to save the output Parquet",
                       buttonType = "default btn-sm w-100"),
        verbatimTextOutput("csv_out_label", placeholder = TRUE),
        textInput("csv_plate_id", "Plate ID (optional)", placeholder = "e.g. Plate1"),
        actionButton("csv_convert_btn", "Convert & add to queue",
                     icon = icon("exchange-alt"), class = "btn-primary btn-sm w-100")
      )
    ),
    hr(),

    # ── Step 2: File browser ──────────────────────────────────────────────────
    h6("Step 2 — Load Parquets"),
    layout_columns(
      shinyFilesButton("add_files",  "Add files",
                       title = "Select raw trace Parquets", multiple = TRUE,
                       buttonType = "default btn-sm"),
      shinyDirButton("add_folder", "Add folder",
                     title = "Add all .parquets in folder",
                     buttonType = "default btn-sm"),
      col_widths = c(6, 6)
    ),
    uiOutput("file_list_ui"),
    actionButton("clear_files", "Clear", class = "btn-outline-secondary btn-sm"),
    hr(),

    # ── Step 3: Pre-scan ──────────────────────────────────────────────────────
    h6("Step 3 — Scan & Preview"),
    layout_columns(
      numericInput("scan_threshold", "AP threshold (mV)", 0, step = 5),
      checkboxInput("v_is_volts", "V → mV", TRUE),
      col_widths = c(7, 5)
    ),
    textInput("scan_sweeps", "Sweeps to scan",
              placeholder = "e.g.  8, 9, 10, 11  — blank = all"),
    helpText("Scans voltage max AND pre-detects APs on likely wells with default config."),
    actionButton("scan_btn", "Scan & pre-detect", icon = icon("search"),
                 class = "btn-primary btn-sm w-100"),
    hr(),

    # ── Plate grid ────────────────────────────────────────────────────────────
    uiOutput("plate_selector_ui"),
    div(style = "position:relative;",
      plotOutput("plate_grid", height = "180px", click = "plate_click",
                 hover = hoverOpts("plate_hover", delay = 50, delayType = "debounce")),
      uiOutput("plate_tooltip")
    ),
    uiOutput("selected_well_label_ui"),
    hr(),

    # ── Sweep navigator ───────────────────────────────────────────────────────
    uiOutput("sweep_ui"),
    layout_columns(
      checkboxInput("quick_mode",    "Detect on this sweep only", TRUE),
      checkboxInput("show_all_sweeps", "Overlay all",             FALSE),
      col_widths = c(6, 6)
    ),
    hr(),

    # ── Run full analysis ─────────────────────────────────────────────────────
    h6("Step 4 — Run Full Analysis"),
    shinyDirButton("out_dir_btn", "Choose output folder",
                   title = "Output for sweep + AP parquets",
                   buttonType = "default btn-sm w-100"),
    verbatimTextOutput("out_dir_label", placeholder = TRUE),
    actionButton("run_btn", "Run on all loaded parquets",
                 icon = icon("play"), class = "btn-success btn-sm w-100")
  ),

  # ── Main panel ────────────────────────────────────────────────────────────────
  layout_columns(
    card(
      card_header(textOutput("plot_title", inline = TRUE)),
      plotOutput("trace_plot", height = "480px"),
      card_footer(textOutput("ap_count_msg"))
    ),
    card(
      navset_tab(
        nav_panel(
          "Sweep Summary",
          tableOutput("sweep_table"),
          helpText("stim_pA = commanded current  •  ap_count = detected spikes")
        ),
        nav_panel(
          "Phase Plane",
          plotOutput("phase_plot", height = "420px"),
          helpText("dV/dt vs V for selected sweep.  ▲ = dV/dt max,  ▼ = dV/dt min per AP.")
        )
      )
    ),
    col_widths = c(8, 4)
  ),

  # ── Detection config (below plots) ────────────────────────────────────────────
  accordion(
    id = "acc_cfg",
    open = FALSE,
    accordion_panel(
      "Detection Config",
      layout_columns(
        sliderInput("dvdt_thr",      "dV/dt threshold (mV/ms)", 5,  100, 25,  step = 1),
        sliderInput("smooth_ms",     "Smoothing (ms)",           0,    2,  0.2, step = 0.05),
        sliderInput("refractory_ms", "Refractory (ms)",        0.5,   10,  2,   step = 0.5),
        sliderInput("min_peak_v",    "Min peak voltage (mV)",  -60,   20, -20,  step = 5),
        col_widths = c(3, 3, 3, 3)
      ),
      layout_columns(
        numericInput("min_height", "Min height (mV)", NA, step = 1),
        numericInput("prominence", "Prominence (mV)", NA, step = 1),
        helpText(style = "margin-top:1.8em;",
                 "Leave height / prominence blank for adaptive threshold."),
        layout_columns(
          downloadButton("export_cfg", "Export config",
                         icon = icon("download"), class = "btn-sm btn-outline-secondary w-100"),
          actionButton("reset_cfg", "Defaults",
                       icon = icon("undo"), class = "btn-sm btn-outline-secondary w-100"),
          col_widths = c(6, 6)
        ),
        col_widths = c(2, 2, 4, 4)
      ),
      fileInput("import_cfg", NULL, accept = ".json",
                buttonLabel = "Load config", placeholder = "No file selected",
                width = "300px")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  parquet_paths     <- reactiveVal(character(0))
  scan_results      <- reactiveVal(NULL)
  scan_ap_cache     <- reactiveVal(NULL)   # well_key → previewAPDetection() result
  cfg_dirty         <- reactiveVal(FALSE)  # TRUE once user edits any config slider
  out_dir_path      <- reactiveVal(NULL)
  selected_well_key <- reactiveVal(NULL)
  csv_in_path       <- reactiveVal(NULL)
  csv_out_path      <- reactiveVal(NULL)

  shinyFileChoose(input, "add_files",   roots = .roots, filetypes = "parquet")
  shinyDirChoose(input,  "add_folder",  roots = .roots)
  shinyDirChoose(input,  "out_dir_btn", roots = .roots)
  shinyDirChoose(input,  "csv_in_dir",  roots = .roots)
  shinyDirChoose(input,  "csv_out_dir", roots = .roots)

  # ── File management ────────────────────────────────────────────────────────
  observeEvent(input$add_files, {
    req(is.list(input$add_files))
    new <- parseFilePaths(.roots, input$add_files)$datapath
    parquet_paths(unique(c(parquet_paths(), new[file.exists(new)])))
  })

  observeEvent(input$add_folder, {
    req(is.list(input$add_folder))
    folder <- parseDirPath(.roots, input$add_folder)
    req(nzchar(folder))
    found <- list.files(folder, pattern = "\\.parquet$", full.names = TRUE)
    if (!length(found)) { showNotification("No .parquet files found.", type = "warning"); return() }
    parquet_paths(unique(c(parquet_paths(), found)))
    showNotification(sprintf("Added %d parquet(s).", length(found)), type = "message")
  })

  observeEvent(input$clear_files, { parquet_paths(character(0)); scan_results(NULL) })

  observeEvent(input$out_dir_btn, {
    req(is.list(input$out_dir_btn))
    d <- parseDirPath(.roots, input$out_dir_btn)
    if (nzchar(d)) out_dir_path(d)
  })

  output$file_list_ui <- renderUI({
    paths <- parquet_paths()
    if (!length(paths)) return(helpText("No files loaded."))
    do.call(tags$ul, c(
      list(style = "font-size:0.8em; max-height:80px; overflow-y:auto; padding-left:1em;"),
      lapply(paths, function(p) tags$li(basename(p)))
    ))
  })

  output$out_dir_label <- renderText({
    d <- out_dir_path(); if (is.null(d)) "(none selected)" else d
  })

  # Mark config dirty when user touches any slider
  observeEvent(
    list(input$dvdt_thr, input$smooth_ms, input$refractory_ms,
         input$min_peak_v, input$min_height, input$prominence),
    { cfg_dirty(TRUE) },
    ignoreInit = TRUE
  )

  # ── Scan + pre-detect ────────────────────────────────────────────────────────
  observeEvent(input$scan_btn, {
    paths <- parquet_paths()
    if (!length(paths)) {
      showNotification("Add at least one raw trace Parquet first.", type = "error"); return()
    }

    scan_sw <- if (nzchar(trimws(input$scan_sweeps %||% ""))) {
      as.integer(na.omit(suppressWarnings(
        as.integer(trimws(strsplit(input$scan_sweeps, "[,;\\s]+", perl = TRUE)[[1]])))))
    } else NULL

    withProgress(message = "Scanning wells…", value = 0.1, {
      result <- tryCatch(
        scanParquetForAPs(paths,
                          v_is_volts   = isTRUE(input$v_is_volts),
                          threshold_mV = input$scan_threshold %||% 0,
                          sweeps       = scan_sw),
        error = function(e) {
          showNotification(paste("Scan failed:", e$message), type = "error"); NULL
        }
      )
      if (is.null(result)) return()

      n_likely <- sum(result$ap_likely)
      showNotification(
        sprintf("%d wells — %d AP-likely. Pre-detecting…", nrow(result), n_likely),
        type = "message", duration = 3
      )
      scan_results(result)

      # Auto-select best well immediately so plate renders
      best <- result[which.max(result$n_ap_sweeps), ]
      best_key <- .make_well_key(best$well_id, best$plate_id, best$parquet_path)
      selected_well_key(best_key)

      # Pre-run AP detection on all AP-likely wells with default config
      if (n_likely > 0) {
        likely_rows <- result[result$ap_likely, , drop = FALSE]
        cache <- list()
        default_cfg <- cc_config(v_is_volts = isTRUE(input$v_is_volts))

        incProgress(0.1)
        for (i in seq_len(nrow(likely_rows))) {
          incProgress(0.8 / nrow(likely_rows),
                      detail = likely_rows$well_id[i],
                      message = sprintf("Pre-detecting AP (%d/%d)…", i, nrow(likely_rows)))
          w   <- likely_rows[i, ]
          key <- .make_well_key(w$well_id, w$plate_id, w$parquet_path)
          cache[[key]] <- tryCatch(
            previewAPDetection(w$parquet_path, w$well_id, w$plate_id,
                               sweeps = scan_sw, cfg = default_cfg),
            error = function(e) NULL
          )
        }
        scan_ap_cache(cache)
        cfg_dirty(FALSE)   # cache matches current (default) config
        showNotification(
          sprintf("Done. %d AP-likely wells pre-detected.", n_likely),
          type = "message", duration = 4
        )
      }
    })
  })

  # ── Plate selector (when multiple plates loaded) ───────────────────────────
  output$plate_selector_ui <- renderUI({
    scan_df <- scan_results(); req(scan_df)
    plates <- unique(scan_df$plate_id)
    if (length(plates) <= 1) return(NULL)
    selectInput("current_plate", "Plate", choices = plates)
  })

  current_plate <- reactive({
    scan_df <- scan_results(); req(scan_df)
    plates <- unique(scan_df$plate_id)
    req(length(plates) > 0)
    if (length(plates) == 1 || is.null(input$current_plate)) plates[1] else input$current_plate
  })

  # ── Plate grid ──────────────────────────────────────────────────────────────
  output$plate_grid <- renderPlot({
    scan_df <- scan_results(); req(scan_df)
    cp      <- current_plate(); req(cp)
    df <- scan_df[scan_df$plate_id == cp, , drop = FALSE]
    req(nrow(df) > 0)

    df$row_idx <- match(substr(df$well_id, 1, 1), LETTERS)
    df$col_num <- as.integer(substr(df$well_id, 2, 3))

    sel <- selected_well_key()
    df$is_sel <- if (!is.null(sel)) {
      si <- .parse_well_key(sel)
      df$well_id == si$well_id & df$plate_id == si$plate_id
    } else FALSE

    n_rows <- max(df$row_idx, na.rm = TRUE)
    n_cols <- max(df$col_num, na.rm = TRUE)

    p <- ggplot2::ggplot(df, ggplot2::aes(col_num, row_idx, fill = n_ap_sweeps)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
      ggplot2::scale_fill_viridis_c(option = "plasma", name = "AP\nsweeps") +
      ggplot2::scale_y_reverse(breaks = seq_len(n_rows),
                                labels = LETTERS[seq_len(n_rows)]) +
      ggplot2::scale_x_continuous(breaks = seq(2, n_cols, by = 2),
                                   expand = ggplot2::expansion(0)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 7) +
      ggplot2::theme(
        legend.position      = "right",
        legend.key.height    = ggplot2::unit(0.5, "cm"),
        legend.key.width     = ggplot2::unit(0.25, "cm"),
        legend.title         = ggplot2::element_text(size = 6),
        axis.text            = ggplot2::element_text(size = 6),
        panel.grid           = ggplot2::element_blank()
      )

    # Highlight selected well
    sel_df <- df[df$is_sel, , drop = FALSE]
    if (nrow(sel_df) > 0) {
      p <- p + ggplot2::geom_tile(
        data = sel_df,
        ggplot2::aes(col_num, row_idx),
        fill = NA, colour = "white", linewidth = 1.2, inherit.aes = FALSE
      )
    }
    p
  }, bg = "transparent")

  # Hover tooltip on plate grid
  output$plate_tooltip <- renderUI({
    h <- input$plate_hover
    if (is.null(h)) return(NULL)
    scan_df <- scan_results(); req(scan_df)
    cp <- current_plate(); req(cp)
    df <- scan_df[scan_df$plate_id == cp, , drop = FALSE]
    df$row_idx <- match(substr(df$well_id, 1, 1), LETTERS)
    df$col_num <- as.integer(substr(df$well_id, 2, 3))
    col_near <- round(h$x)
    row_near <- round(h$y)
    hit <- df[df$col_num == col_near & df$row_idx == row_near, , drop = FALSE]
    if (nrow(hit) == 0) return(NULL)
    lbl <- sprintf("%s | %d AP sweeps | %.0f mV max",
                   hit$well_id[1], as.integer(hit$n_ap_sweeps[1]), hit$max_v_mV[1])
    left_px <- as.integer(round(h$coords_css$x)) + 8L
    top_px  <- as.integer(round(h$coords_css$y)) - 28L
    div(
      style = sprintf(paste0(
        "position:absolute;left:%dpx;top:%dpx;",
        "background:rgba(0,0,0,0.75);color:white;",
        "padding:2px 6px;border-radius:3px;",
        "font-size:11px;pointer-events:none;white-space:nowrap;"
      ), left_px, top_px),
      lbl
    )
  })

  # Click on plate grid → select well
  observeEvent(input$plate_click, {
    scan_df <- scan_results(); req(scan_df)
    cp      <- current_plate(); req(cp)

    col_click <- round(input$plate_click$x)
    row_click  <- round(input$plate_click$y)
    if (anyNA(c(col_click, row_click))) return()
    if (row_click < 1 || row_click > 16 || col_click < 1 || col_click > 24) return()

    clicked_well <- sprintf("%s%02d", LETTERS[row_click], col_click)
    df_cp  <- scan_df[scan_df$plate_id == cp, , drop = FALSE]
    hit    <- df_cp[df_cp$well_id == clicked_well, , drop = FALSE]
    if (nrow(hit) > 0)
      selected_well_key(.make_well_key(hit$well_id[1], hit$plate_id[1], hit$parquet_path[1]))
  })

  output$selected_well_label_ui <- renderUI({
    sel <- selected_well_key()
    if (is.null(sel)) return(helpText("Click a well in the grid to select it."))
    si <- .parse_well_key(sel)
    tags$p(style = "font-size:0.82em; margin:2px 0 0;",
           tags$b(si$well_id), " | ", si$plate_id)
  })

  # ── Well info from reactive val ────────────────────────────────────────────
  well_info <- reactive({
    sel <- selected_well_key(); req(sel)
    .parse_well_key(sel)
  })

  # ── Load raw traces for selected well ─────────────────────────────────────
  well_traces <- reactive({
    info <- well_info(); req(info$path, file.exists(info$path))
    arrow::open_dataset(info$path) |>
      dplyr::filter(well_id  == !!info$well_id,
                    plate_id == !!info$plate_id) |>
      dplyr::collect() |>
      dplyr::mutate(v_mV = recorded * (if (isTRUE(input$v_is_volts)) 1000 else 1))
  })

  # ── Sweep slider ──────────────────────────────────────────────────────────
  output$sweep_ui <- renderUI({
    df <- well_traces(); req(df)
    sweeps <- sort(unique(df$sweep))
    tagList(
      h6("Sweep"),
      sliderInput("selected_sweep", NULL,
                  min = min(sweeps), max = max(sweeps),
                  value = max(sweeps), step = 1, ticks = FALSE)
    )
  })

  # ── AP detection (Python, debounced) ──────────────────────────────────────
  cfg_inputs <- reactive({
    list(dvdt_thr      = input$dvdt_thr,
         smooth_ms     = input$smooth_ms,
         refractory_ms = input$refractory_ms,
         min_peak_v    = input$min_peak_v,
         min_height    = input$min_height,
         prominence    = input$prominence,
         quick         = input$quick_mode,
         sweep         = input$selected_sweep)
  }) |> debounce(500)

  ap_result <- reactive({
    info <- well_info()
    cfg  <- cfg_inputs()
    req(info, cfg)

    key   <- .make_well_key(info$well_id, info$plate_id, info$path)
    cache <- scan_ap_cache()

    # Return cached result instantly when config is still at scan defaults
    if (!isTRUE(cfg_dirty()) && !is.null(cache) && !is.null(cache[[key]])) {
      return(cache[[key]])
    }

    scan_sw_str <- trimws(input$scan_sweeps %||% "")
    scan_sw <- if (nzchar(scan_sw_str))
      as.integer(na.omit(suppressWarnings(
        as.integer(trimws(strsplit(scan_sw_str, "[,;\\s]+", perl = TRUE)[[1]])))))
    else NULL

    sweeps_arg <- if (isTRUE(cfg$quick) && !is.null(cfg$sweep))
      as.integer(cfg$sweep) else scan_sw

    tryCatch(
      previewAPDetection(
        parquet_path = info$path,
        well_id      = info$well_id,
        plate_id     = info$plate_id,
        sweeps       = sweeps_arg,
        cfg = cc_config(
          v_is_volts          = isTRUE(input$v_is_volts),
          dvdt_thr_mV_per_ms  = cfg$dvdt_thr,
          smooth_ms           = cfg$smooth_ms,
          refractory_ms       = cfg$refractory_ms,
          min_peak_voltage_mV = cfg$min_peak_v,
          min_height_mV       = .opt_val(cfg$min_height),
          prominence_mV       = .opt_val(cfg$prominence)
        )
      ),
      error = function(e) {
        showNotification(paste("Detection failed:", conditionMessage(e)),
                         type = "warning", duration = 10)
        NULL
      }
    )
  })

  # ── Trace plot ─────────────────────────────────────────────────────────────
  output$plot_title <- renderText({
    info <- well_info(); req(info)
    paste0(info$well_id, "  •  ", info$plate_id, "  •  ", basename(info$path))
  })

  output$trace_plot <- renderPlot({
    df       <- well_traces(); req(df)
    ap       <- ap_result()
    show_all <- isTRUE(input$show_all_sweeps)
    sw       <- if (show_all || is.null(input$selected_sweep)) NULL else input$selected_sweep
    plot_df  <- if (is.null(sw)) df else df[df$sweep == sw, ]

    # Base traces
    if (show_all) {
      p <- ggplot2::ggplot(plot_df,
                           ggplot2::aes(t_s, v_mV,
                                        colour = as.numeric(sweep),
                                        group  = factor(sweep))) +
        ggplot2::geom_line(alpha = 0.55, linewidth = 0.4) +
        ggplot2::scale_color_viridis_c(option = "plasma", guide = "none")
    } else {
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(t_s, v_mV)) +
        ggplot2::geom_line(colour = "steelblue", linewidth = 0.5)
    }

    # Stimulus window shading
    if (!is.null(ap) && nrow(ap$sweeps) > 0) {
      sw_row <- if (!show_all && !is.null(sw))
        ap$sweeps[ap$sweeps$sweep == sw, , drop = FALSE]
      else
        ap$sweeps[1, , drop = FALSE]
      if (all(c("stim_start_s", "stim_end_s") %in% names(sw_row))) {
        st <- sw_row[is.finite(sw_row$stim_start_s) & is.finite(sw_row$stim_end_s), ]
        if (nrow(st) > 0)
          p <- p + ggplot2::annotate("rect",
                                     xmin = st$stim_start_s[1], xmax = st$stim_end_s[1],
                                     ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "gold")
      }
    }

    # AP markers
    if (!is.null(ap) && nrow(ap$aps) > 0) {
      ap_plot <- ap$aps[!is.na(ap$aps$peak_time_s), ]
      if (!show_all && !is.null(sw))
        ap_plot <- ap_plot[ap_plot$sweep == sw, , drop = FALSE]
      if (nrow(ap_plot) > 0)
        p <- p + ggplot2::geom_point(
          data = ap_plot,
          ggplot2::aes(x = peak_time_s, y = peak_v_mV),
          colour = "firebrick", fill = "firebrick",
          shape = 25, size = 2.2, inherit.aes = FALSE
        )
    }

    p + ggplot2::labs(x = "Time (s)", y = "Voltage (mV)") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35))
  })

  # ── Phase plane plot ──────────────────────────────────────────────────────────
  output$phase_plot <- renderPlot({
    df    <- well_traces(); req(df)
    sw    <- input$selected_sweep; req(sw)
    df_sw <- df[df$sweep == sw & !is.na(df$v_mV), ]
    req(nrow(df_sw) > 3)
    df_sw <- df_sw[order(df_sw$t_s), ]

    dt_ms  <- diff(df_sw$t_s) * 1000          # seconds → ms
    dvdt   <- diff(df_sw$v_mV) / dt_ms         # mV/ms
    v_mid  <- (df_sw$v_mV[-1] + df_sw$v_mV[-nrow(df_sw)]) / 2

    # Drop any large time-gap steps (shouldn't exist within one sweep, but be safe)
    keep <- dt_ms > 0 & dt_ms < 5
    phase_df <- data.frame(v_mV = v_mid[keep], dvdt = dvdt[keep])

    p <- ggplot2::ggplot(phase_df, ggplot2::aes(v_mV, dvdt)) +
      ggplot2::geom_path(colour = "steelblue", alpha = 0.7, linewidth = 0.45) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "dashed") +
      ggplot2::labs(x = "Voltage (mV)", y = "dV/dt (mV/ms)") +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35))

    ap <- ap_result()
    if (!is.null(ap) && nrow(ap$aps) > 0) {
      ap_sw <- ap$aps[ap$aps$sweep == sw, , drop = FALSE]
      ap_sw <- ap_sw[!is.na(ap_sw$dvdt_max_v_mV) & !is.na(ap_sw$dvdt_max_mV_per_ms), , drop = FALSE]
      if (nrow(ap_sw) > 0) {
        p <- p +
          ggplot2::geom_point(
            data = ap_sw,
            ggplot2::aes(x = dvdt_max_v_mV, y = dvdt_max_mV_per_ms),
            colour = "firebrick", shape = 24, size = 2.8, fill = "firebrick",
            inherit.aes = FALSE
          ) +
          ggplot2::geom_point(
            data = ap_sw,
            ggplot2::aes(x = dvdt_min_v_mV, y = dvdt_min_mV_per_ms),
            colour = "navy", shape = 25, size = 2.8, fill = "navy",
            inherit.aes = FALSE
          )
      }
    }
    p
  })

  output$ap_count_msg <- renderText({
    ap <- ap_result(); req(ap)
    sw    <- input$selected_sweep
    total <- nrow(ap$aps)
    if (!isTRUE(input$show_all_sweeps) && !is.null(sw)) {
      n <- sum(ap$aps$sweep == sw, na.rm = TRUE)
      sprintf("Sweep %s: %d AP(s)   |   Total across analysed sweeps: %d", sw, n, total)
    } else {
      sprintf("Total APs across analysed sweeps: %d", total)
    }
  })

  # ── Sweep summary table ────────────────────────────────────────────────────
  output$sweep_table <- renderTable({
    ap <- ap_result(); req(ap, nrow(ap$sweeps) > 0)
    cols <- intersect(c("sweep", "stim_amp_cmd", "ap_count", "ap_freq_hz",
                        "baseline_v_mV", "steady_v_mV"), names(ap$sweeps))
    out <- ap$sweeps[, cols, drop = FALSE]
    if ("stim_amp_cmd" %in% names(out)) {
      out$stim_pA <- round(out$stim_amp_cmd * 1e12, 1); out$stim_amp_cmd <- NULL
    }
    for (col in c("ap_freq_hz", "baseline_v_mV", "steady_v_mV")) {
      if (col %in% names(out)) out[[col]] <- round(out[[col]], 1)
    }
    out
  }, digits = 1, striped = TRUE, hover = TRUE, spacing = "s",
  width = "100%", align = "r")

  # ── Run full analysis ──────────────────────────────────────────────────────
  observeEvent(input$run_btn, {
    out_dir <- out_dir_path()
    if (is.null(out_dir) || !nzchar(out_dir)) {
      showNotification("Choose an output folder first.", type = "error"); return()
    }
    paths <- parquet_paths()
    if (!length(paths)) {
      showNotification("No parquet files loaded.", type = "error"); return()
    }
    cfg <- cc_config(
      v_is_volts          = isTRUE(input$v_is_volts),
      dvdt_thr_mV_per_ms  = input$dvdt_thr,
      smooth_ms           = input$smooth_ms,
      refractory_ms       = input$refractory_ms,
      min_peak_voltage_mV = input$min_peak_v,
      min_height_mV       = .opt_val(input$min_height),
      prominence_mV       = .opt_val(input$prominence)
    )
    withProgress(message = "Running AP analysis…", value = 0, {
      for (i in seq_along(paths)) {
        incProgress(1 / length(paths), detail = basename(paths[i]))
        tryCatch(
          runAPAnalysis(paths[i], out_dir = out_dir, cfg = cfg),
          error = function(e)
            showNotification(paste0(basename(paths[i]), ": ", e$message),
                             type = "error", duration = 15)
        )
      }
    })
    showNotification(paste0("Done! Results in: ", out_dir),
                     type = "message", duration = 10)
  })

  # ── Config export / import / reset ────────────────────────────────────────
  .cfg_defaults <- list(
    dvdt_thr = 25, smooth_ms = 0.2, refractory_ms = 2, min_peak_v = -20,
    min_height = NA_real_, prominence = NA_real_,
    v_is_volts = TRUE, scan_threshold = 0
  )

  output$export_cfg <- downloadHandler(
    filename = function()
      paste0("cc_config_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"),
    content = function(file) {
      cfg <- list(
        dvdt_thr      = input$dvdt_thr,
        smooth_ms     = input$smooth_ms,
        refractory_ms = input$refractory_ms,
        min_peak_v    = input$min_peak_v,
        min_height    = if (is.na(input$min_height))    NULL else input$min_height,
        prominence    = if (is.na(input$prominence))    NULL else input$prominence,
        v_is_volts    = isTRUE(input$v_is_volts),
        scan_threshold = input$scan_threshold %||% 0
      )
      jsonlite::write_json(cfg, file, auto_unbox = TRUE, null = "null", pretty = TRUE)
    }
  )

  observeEvent(input$import_cfg, {
    req(input$import_cfg)
    cfg <- tryCatch(jsonlite::read_json(input$import_cfg$datapath, simplifyVector = TRUE),
                    error = function(e) NULL)
    if (is.null(cfg)) {
      showNotification("Invalid config JSON.", type = "error"); return()
    }
    if (!is.null(cfg$dvdt_thr))       updateSliderInput(session,  "dvdt_thr",       value = cfg$dvdt_thr)
    if (!is.null(cfg$smooth_ms))      updateSliderInput(session,  "smooth_ms",      value = cfg$smooth_ms)
    if (!is.null(cfg$refractory_ms))  updateSliderInput(session,  "refractory_ms",  value = cfg$refractory_ms)
    if (!is.null(cfg$min_peak_v))     updateSliderInput(session,  "min_peak_v",     value = cfg$min_peak_v)
    updateNumericInput(session, "min_height",     value = cfg$min_height    %||% NA)
    updateNumericInput(session, "prominence",     value = cfg$prominence    %||% NA)
    updateNumericInput(session, "scan_threshold", value = cfg$scan_threshold %||% 0)
    if (!is.null(cfg$v_is_volts))
      updateCheckboxInput(session, "v_is_volts", value = isTRUE(cfg$v_is_volts))
    showNotification("Config loaded.", type = "message")
  })

  observeEvent(input$reset_cfg, {
    d <- .cfg_defaults
    updateSliderInput(session,   "dvdt_thr",       value = d$dvdt_thr)
    updateSliderInput(session,   "smooth_ms",       value = d$smooth_ms)
    updateSliderInput(session,   "refractory_ms",   value = d$refractory_ms)
    updateSliderInput(session,   "min_peak_v",      value = d$min_peak_v)
    updateNumericInput(session,  "min_height",      value = d$min_height)
    updateNumericInput(session,  "prominence",      value = d$prominence)
    updateCheckboxInput(session, "v_is_volts",      value = d$v_is_volts)
    updateNumericInput(session,  "scan_threshold",  value = d$scan_threshold)
    showNotification("Config reset to defaults.", type = "message")
  })

  # ── CSV → Parquet ─────────────────────────────────────────────────────────
  observeEvent(input$csv_in_dir, {
    req(is.list(input$csv_in_dir))
    d <- parseDirPath(.roots, input$csv_in_dir)
    if (nzchar(d)) csv_in_path(d)
  })

  observeEvent(input$csv_out_dir, {
    req(is.list(input$csv_out_dir))
    d <- parseDirPath(.roots, input$csv_out_dir)
    if (nzchar(d)) csv_out_path(d)
  })

  output$csv_in_label  <- renderText(csv_in_path()  %||% "(none selected)")
  output$csv_out_label <- renderText(csv_out_path() %||% "(none selected)")

  observeEvent(input$csv_convert_btn, {
    inp <- csv_in_path();  out <- csv_out_path()
    if (is.null(inp) || !nzchar(inp)) {
      showNotification("Select CSV input folder first.", type = "error"); return()
    }
    if (is.null(out) || !nzchar(out)) {
      showNotification("Select output folder first.", type = "error"); return()
    }
    plate_id_arg <- trimws(input$csv_plate_id %||% "")
    out_parquet  <- file.path(out, paste0(basename(inp), ".parquet"))
    withProgress(message = "Converting CSVs to Parquet…", {
      tryCatch({
        csvToParquet(inp, out_parquet,
                     plate_id = if (nzchar(plate_id_arg)) plate_id_arg else NULL)
        parquet_paths(unique(c(parquet_paths(), out_parquet)))
        showNotification(
          paste0("Done: ", basename(out_parquet), "  — added to queue."),
          type = "message", duration = 8
        )
      }, error = function(e)
        showNotification(paste("Conversion failed:", e$message),
                         type = "error", duration = 15))
    })
  })
}

shinyApp(ui, server)
