#' CC Preview Shiny Module
#'
#' Interactive AP detection preview — plate grid, trace/phase-plane tabs,
#' pre-scan caching, and SE generation.  Embedded in tinySEV as the
#' "Current Clamp" tab.
#'
#' @param id Module namespace ID.
#' @name ccPreview
NULL

# ── Internal helpers ──────────────────────────────────────────────────────────

.cc_roots <- c(Home = path.expand("~"), `C:` = "C:/", `D:` = "D:/")

.cc_make_well_key  <- function(well_id, plate_id, path)
  paste(well_id, plate_id, path, sep = "|||")

.cc_parse_well_key <- function(v) {
  p <- strsplit(v, "|||", fixed = TRUE)[[1]]
  list(well_id = p[1], plate_id = p[2], path = p[3])
}

.cc_opt_val <- function(x)
  if (is.numeric(x) && length(x) == 1 && !is.na(x)) x else NULL

.cc_cfg_defaults <- list(
  dvdt_thr = 25, smooth_ms = 0.2, refractory_ms = 2, min_peak_v = -20,
  min_height = NA_real_, prominence = NA_real_,
  v_is_volts = TRUE, scan_threshold = 0
)

.cc_color_vars <- c(
  "AP sweeps (scan)"   = "n_ap_sweeps",
  "Total APs detected" = "total_aps",
  "Max freq. (Hz)"     = "max_freq_hz",
  "Baseline Vm (mV)"   = "baseline_vm_mV",
  "Mean AP amp. (mV)"  = "mean_ap_amp_mV"
)

# Bootstrap-3-safe collapsible panel (works in shinydashboard)
.cc_collapse_panel <- function(title, ..., btn_id, collapsed = TRUE) {
  panel_id <- paste0(btn_id, "_body")
  tagList(
    tags$a(
      href          = paste0("#", panel_id),
      `data-toggle` = "collapse",
      style         = paste0(
        "display:block; font-weight:600; font-size:0.88em; color:#337ab7;",
        " padding:4px 2px; margin-bottom:4px; text-decoration:none;"),
      title
    ),
    tags$div(
      id    = panel_id,
      class = if (collapsed) "collapse" else "collapse in",
      tags$div(style = "padding:4px 2px 6px 2px;", ...)
    )
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────

#' @rdname ccPreview
#' @export
ccPreviewUI_csv <- function(id) {
  ns <- NS(id)
  tagList(
    h4("CSV → Parquet Converter"),
    wellPanel(
      style = "padding:12px;",
      helpText("Convert DataControl per-well CSV files into a single Parquet file."),
      shinyFiles::shinyDirButton(ns("csv_in_dir"), "Input folder (CSVs)",
        title = "Folder with per-well CSV files",
        buttonType = "default btn-sm btn-block"),
      verbatimTextOutput(ns("csv_in_label"), placeholder = TRUE),
      shinyFiles::shinyDirButton(ns("csv_out_dir"), "Output folder",
        title = "Where to save the Parquet",
        buttonType = "default btn-sm btn-block"),
      verbatimTextOutput(ns("csv_out_label"), placeholder = TRUE),
      textInput(ns("csv_plate_id"), "Plate ID (optional)",
                placeholder = "e.g. Plate1"),
      actionButton(ns("csv_convert_btn"), "Convert & add to queue",
                   icon  = icon("exchange-alt"),
                   class = "btn-primary btn-sm btn-block")
    )
  )
}

#' @rdname ccPreview
#' @export
ccPreviewUI_detect <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(3,
      wellPanel(
        style = "overflow-y:auto; max-height:calc(100vh - 80px); padding:10px;",

        h6("Step 2 — Load Parquets"),
        tags$div(style = "display:flex; gap:4px; margin-bottom:4px;",
          shinyFiles::shinyFilesButton(ns("add_files"), "Add files",
            title = "Select raw trace Parquets", multiple = TRUE,
            buttonType = "default btn-sm"),
          shinyFiles::shinyDirButton(ns("add_folder"), "Add folder",
            title = "Add all .parquets in folder",
            buttonType = "default btn-sm")
        ),
        uiOutput(ns("file_list_ui")),
        actionButton(ns("clear_files"), "Clear",
                     class = "btn-outline-secondary btn-sm"),
        tags$hr(),

        h6("Step 3 — Scan & Pre-detect"),
        fluidRow(
          column(7, numericInput(ns("scan_threshold"), "AP threshold (mV)", 0, step = 5)),
          column(5, checkboxInput(ns("v_is_volts"), "V → mV", TRUE))
        ),
        textInput(ns("scan_sweeps"), "Sweeps to scan",
                  placeholder = "e.g. 8,9,10 — blank = all"),
        helpText("Scans voltage max then pre-detects APs on likely wells."),
        actionButton(ns("scan_btn"), "Scan & pre-detect", icon = icon("search"),
                     class = "btn-primary btn-sm btn-block"),
        tags$hr(),

        uiOutput(ns("sweep_ui")),
        checkboxInput(ns("show_all_sweeps"), "Overlay all sweeps", FALSE),
        tags$hr(),

        h6("Step 4 — Run Full Analysis"),
        shinyFiles::shinyDirButton(ns("out_dir_btn"), "Choose output folder",
          title = "Output folder for sweep + AP parquets",
          buttonType = "default btn-sm btn-block"),
        verbatimTextOutput(ns("out_dir_label"), placeholder = TRUE),
        tags$div(style = "display:flex; gap:4px;",
          actionButton(ns("run_btn"),      "Run",  icon = icon("play"),  class = "btn-success btn-sm"),
          actionButton(ns("stop_run_btn"), "Stop", icon = icon("stop"),  class = "btn-danger btn-sm")
        ),
        uiOutput(ns("run_status_ui")),
        tags$div(
          style = paste0("margin-top:4px; max-height:160px; overflow-y:auto;",
                         " background:#1e1e1e; border-radius:3px; padding:4px 6px;"),
          verbatimTextOutput(ns("run_log_out"), placeholder = FALSE)
        )
      )
    ),

    column(9,
      fluidRow(
        column(8,
          tags$div(style = "position:relative;",
            plotOutput(ns("plate_grid"), height = "360px",
                       click = ns("plate_click"),
                       hover = hoverOpts(ns("plate_hover"), delay = 50,
                                         delayType = "debounce")),
            uiOutput(ns("plate_tooltip"))
          )
        ),
        column(4,
          uiOutput(ns("plate_selector_ui")),
          uiOutput(ns("plate_color_ui")),
          uiOutput(ns("selected_well_label_ui"))
        )
      ),

      tags$hr(style = "margin:6px 0;"),

      tabsetPanel(type = "tabs",
        tabPanel("Trace",
          plotOutput(ns("trace_plot"), height = "380px"),
          tags$div(style = "font-size:0.82em; color:#666; padding:3px 6px;",
                   textOutput(ns("ap_count_msg"), inline = TRUE))
        ),
        tabPanel("Phase Plane",
          plotOutput(ns("phase_plot"), height = "380px"),
          helpText("dV/dt vs V.  ▲ = dV/dt max per AP,  ▼ = dV/dt min per AP.")
        ),
        tabPanel("Sweep Summary",
          tableOutput(ns("sweep_table")),
          helpText("stim_pA = commanded current  •  ap_count = detected spikes")
        )
      ),

      tags$br(),

      .cc_collapse_panel(
        "▶ Detection Config", btn_id = paste0(ns("acc_cfg")), collapsed = TRUE,
        fluidRow(
          column(3, sliderInput(ns("dvdt_thr"),      "dV/dt thr. (mV/ms)", 5, 100, 25,  step = 1)),
          column(3, sliderInput(ns("smooth_ms"),     "Smoothing (ms)",      0,   2,  0.2, step = 0.05)),
          column(3, sliderInput(ns("refractory_ms"), "Refractory (ms)",   0.5,  10,  2,   step = 0.5)),
          column(3, sliderInput(ns("min_peak_v"),    "Min peak V (mV)",   -60,  20, -20,  step = 5))
        ),
        fluidRow(
          column(2, numericInput(ns("min_height"), "Min height (mV)", NA, step = 1)),
          column(2, numericInput(ns("prominence"), "Prominence (mV)", NA, step = 1)),
          column(4, helpText(style = "margin-top:1.8em;",
                   "Leave height/prominence blank for adaptive threshold.")),
          column(2,
            tags$div(style = "margin-top:1.5em;",
              tags$div(downloadButton(ns("export_cfg"), "Export",
                                      class = "btn-sm btn-outline-secondary btn-block")),
              tags$div(style = "margin-top:4px;",
                actionButton(ns("reset_cfg"), "Defaults", icon = icon("undo"),
                             class = "btn-sm btn-outline-secondary btn-block"))
            )
          ),
          column(2,
            fileInput(ns("import_cfg"), NULL, accept = ".json",
                      buttonLabel = "Load", placeholder = "config .json")
          )
        )
      )
    )
  )
}

#' @rdname ccPreview
#' @export
ccPreviewUI_build <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Build SingleCellExperiment"),
    wellPanel(
      style = "max-width:600px; padding:14px;",
      helpText("Assemble a SingleCellExperiment from AP analysis output parquets."),
      tags$div(style = "display:flex; gap:4px; margin-bottom:6px;",
        actionButton(ns("se_scan_out"), "From output folder",
                     icon = icon("folder-open"), class = "btn-default btn-sm"),
        shinyFiles::shinyFilesButton(ns("se_add_sweep"), "+ Sweep",
          title = "Add sweep parquet(s)", multiple = TRUE,
          buttonType = "default btn-sm", filetypes = "parquet"),
        shinyFiles::shinyFilesButton(ns("se_add_ap"), "+ AP",
          title = "Add AP parquet(s)", multiple = TRUE,
          buttonType = "default btn-sm", filetypes = "parquet")
      ),
      tags$label("Sweep parquet(s)",
                 class = "control-label", style = "font-size:0.85em;"),
      selectizeInput(ns("se_sweep_files"), NULL, choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "None — add above",
                       plugins = list("remove_button"))),
      tags$label("AP parquet(s) (optional)",
                 class = "control-label", style = "font-size:0.85em;"),
      selectizeInput(ns("se_ap_files"), NULL, choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "Optional",
                       plugins = list("remove_button"))),
      textInput(ns("se_name"), "Object name",
                placeholder = "e.g. CC_experiment_1"),
      actionButton(ns("build_se_btn"), "Build SE",
                   icon = icon("cube"), class = "btn-primary btn-block"),
      uiOutput(ns("se_status_ui")),
      tags$hr(style = "margin:8px 0;"),
      downloadButton(ns("download_se"), "Download SE (.rds)",
                     class = "btn-sm btn-outline-secondary btn-block")
    )
  )
}

#' @rdname ccPreview
#' @export
ccPreviewUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(

      # ── Left control column ─────────────────────────────────────────────────
      column(3,
        wellPanel(
          style = "overflow-y:auto; max-height:calc(100vh - 80px); padding:10px;",

          # Step 1 ── CSV → Parquet (collapsible)
          .cc_collapse_panel(
            "▶ Step 1 — CSV → Parquet", btn_id = "acc_csv", collapsed = TRUE,
            helpText("Convert DataControl CSVs (one per well) to a single Parquet."),
            shinyFiles::shinyDirButton(ns("csv_in_dir"), "Input folder (CSVs)",
              title = "Folder with per-well CSV files",
              buttonType = "default btn-sm btn-block"),
            verbatimTextOutput(ns("csv_in_label"), placeholder = TRUE),
            shinyFiles::shinyDirButton(ns("csv_out_dir"), "Output folder",
              title = "Where to save the Parquet",
              buttonType = "default btn-sm btn-block"),
            verbatimTextOutput(ns("csv_out_label"), placeholder = TRUE),
            textInput(ns("csv_plate_id"), "Plate ID (optional)",
                      placeholder = "e.g. Plate1"),
            actionButton(ns("csv_convert_btn"), "Convert & add to queue",
                         icon = icon("exchange-alt"),
                         class = "btn-primary btn-sm btn-block")
          ),
          tags$hr(),

          # Step 2 ── Load parquets
          h6("Step 2 — Load Parquets"),
          tags$div(style = "display:flex; gap:4px; margin-bottom:4px;",
            shinyFiles::shinyFilesButton(ns("add_files"), "Add files",
              title = "Select raw trace Parquets", multiple = TRUE,
              buttonType = "default btn-sm"),
            shinyFiles::shinyDirButton(ns("add_folder"), "Add folder",
              title = "Add all .parquets in folder",
              buttonType = "default btn-sm")
          ),
          uiOutput(ns("file_list_ui")),
          actionButton(ns("clear_files"), "Clear",
                       class = "btn-outline-secondary btn-sm"),
          tags$hr(),

          # Step 3 ── Scan & pre-detect
          h6("Step 3 — Scan & Pre-detect"),
          fluidRow(
            column(7, numericInput(ns("scan_threshold"), "AP threshold (mV)", 0, step = 5)),
            column(5, checkboxInput(ns("v_is_volts"), "V → mV", TRUE))
          ),
          textInput(ns("scan_sweeps"), "Sweeps to scan",
                    placeholder = "e.g. 8,9,10 — blank = all"),
          helpText("Scans voltage max then pre-detects APs on likely wells (default config)."),
          actionButton(ns("scan_btn"), "Scan & pre-detect", icon = icon("search"),
                       class = "btn-primary btn-sm btn-block"),
          tags$hr(),

          # Sweep controls (rendered after scan)
          uiOutput(ns("sweep_ui")),
          checkboxInput(ns("show_all_sweeps"), "Overlay all sweeps", FALSE),
          tags$hr(),

          # Step 4 ── Run full analysis
          h6("Step 4 — Run Full Analysis"),
          shinyFiles::shinyDirButton(ns("out_dir_btn"), "Choose output folder",
            title = "Output folder for sweep + AP parquets",
            buttonType = "default btn-sm btn-block"),
          verbatimTextOutput(ns("out_dir_label"), placeholder = TRUE),
          tags$div(style = "display:flex; gap:4px;",
            actionButton(ns("run_btn"), "Run",
                         icon = icon("play"), class = "btn-success btn-sm"),
            actionButton(ns("stop_run_btn"), "Stop",
                         icon = icon("stop"), class = "btn-danger btn-sm")
          ),
          uiOutput(ns("run_status_ui")),
          tags$div(
            style = paste0("margin-top:4px; max-height:160px; overflow-y:auto;",
                           " background:#1e1e1e; border-radius:3px; padding:4px 6px;"),
            verbatimTextOutput(ns("run_log_out"), placeholder = FALSE)
          ),
          tags$hr(),

          # Step 5 ── Build SE (collapsible)
          .cc_collapse_panel(
            "▶ Step 5 — Build SE", btn_id = "acc_se", collapsed = TRUE,
            helpText("Assemble a SingleCellExperiment from run outputs."),
            tags$div(style = "display:flex; gap:4px; margin-bottom:6px;",
              actionButton(ns("se_scan_out"), "From output folder",
                           icon = icon("folder-open"),
                           class = "btn-default btn-sm"),
              shinyFiles::shinyFilesButton(ns("se_add_sweep"), "+ Sweep",
                title = "Add sweep parquet(s)", multiple = TRUE,
                buttonType = "default btn-sm", filetypes = "parquet"),
              shinyFiles::shinyFilesButton(ns("se_add_ap"), "+ AP",
                title = "Add AP parquet(s)", multiple = TRUE,
                buttonType = "default btn-sm", filetypes = "parquet")
            ),
            tags$label("Sweep parquet(s)",
                       class = "control-label", style = "font-size:0.85em;"),
            selectizeInput(ns("se_sweep_files"), NULL, choices = NULL,
              multiple = TRUE,
              options = list(placeholder = "None — add above",
                             plugins = list("remove_button"))),
            tags$label("AP parquet(s) (optional)",
                       class = "control-label", style = "font-size:0.85em;"),
            selectizeInput(ns("se_ap_files"), NULL, choices = NULL,
              multiple = TRUE,
              options = list(placeholder = "Optional",
                             plugins = list("remove_button"))),
            textInput(ns("se_name"), "Object name",
                      placeholder = "e.g. CC_experiment_1"),
            actionButton(ns("build_se_btn"), "Build SE",
                         icon = icon("cube"),
                         class = "btn-primary btn-sm btn-block"),
            uiOutput(ns("se_status_ui")),
            tags$hr(style = "margin:6px 0;"),
            tags$div(downloadButton(ns("download_se"), "Download SE (.rds)",
                                    class = "btn-sm btn-outline-secondary btn-block"))
          )
        ) # /wellPanel
      ), # /column 3

      # ── Main content ─────────────────────────────────────────────────────────
      column(9,

        # Plate grid + controls
        fluidRow(
          column(8,
            tags$div(style = "position:relative;",
              plotOutput(ns("plate_grid"), height = "360px",
                         click = ns("plate_click"),
                         hover = hoverOpts(ns("plate_hover"), delay = 50,
                                           delayType = "debounce")),
              uiOutput(ns("plate_tooltip"))
            )
          ),
          column(4,
            uiOutput(ns("plate_selector_ui")),
            uiOutput(ns("plate_color_ui")),
            uiOutput(ns("selected_well_label_ui"))
          )
        ),

        tags$hr(style = "margin:6px 0;"),

        # Plots (tabset)
        tabsetPanel(type = "tabs",
          tabPanel("Trace",
            plotOutput(ns("trace_plot"), height = "380px"),
            tags$div(style = "font-size:0.82em; color:#666; padding:3px 6px;",
                     textOutput(ns("ap_count_msg"), inline = TRUE))
          ),
          tabPanel("Phase Plane",
            plotOutput(ns("phase_plot"), height = "380px"),
            helpText("dV/dt vs V.  ▲ = dV/dt max per AP,  ▼ = dV/dt min per AP.")
          ),
          tabPanel("Sweep Summary",
            tableOutput(ns("sweep_table")),
            helpText("stim_pA = commanded current  •  ap_count = detected spikes")
          )
        ),

        tags$br(),

        # Detection Config (collapsible, below plots)
        .cc_collapse_panel(
          "▶ Detection Config", btn_id = "acc_cfg", collapsed = TRUE,
          fluidRow(
            column(3, sliderInput(ns("dvdt_thr"),      "dV/dt thr. (mV/ms)", 5, 100, 25,  step = 1)),
            column(3, sliderInput(ns("smooth_ms"),     "Smoothing (ms)",      0,   2,  0.2, step = 0.05)),
            column(3, sliderInput(ns("refractory_ms"), "Refractory (ms)",   0.5,  10,  2,   step = 0.5)),
            column(3, sliderInput(ns("min_peak_v"),    "Min peak V (mV)",   -60,  20, -20,  step = 5))
          ),
          fluidRow(
            column(2, numericInput(ns("min_height"), "Min height (mV)", NA, step = 1)),
            column(2, numericInput(ns("prominence"), "Prominence (mV)", NA, step = 1)),
            column(4, helpText(style = "margin-top:1.8em;",
                     "Leave height/prominence blank for adaptive threshold.")),
            column(2,
              tags$div(style = "margin-top:1.5em;",
                tags$div(downloadButton(ns("export_cfg"), "Export",
                                        class = "btn-sm btn-outline-secondary btn-block")),
                tags$div(style = "margin-top:4px;",
                  actionButton(ns("reset_cfg"), "Defaults", icon = icon("undo"),
                               class = "btn-sm btn-outline-secondary btn-block"))
              )
            ),
            column(2,
              fileInput(ns("import_cfg"), NULL, accept = ".json",
                        buttonLabel = "Load", placeholder = "config .json")
            )
          )
        )

      ) # /column 9
    ) # /fluidRow
  )
}


# ── Server ────────────────────────────────────────────────────────────────────

#' @rdname ccPreview
#' @export
ccPreviewServer <- function(id, ses = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    parquet_paths     <- reactiveVal(character(0))
    scan_results      <- reactiveVal(NULL)
    scan_ap_cache     <- reactiveVal(NULL)  # used only for plate-grid metrics
    suggested_sweep   <- reactiveVal(NULL)
    out_dir_path      <- reactiveVal(NULL)
    selected_well_key <- reactiveVal(NULL)
    csv_in_path       <- reactiveVal(NULL)
    csv_out_path      <- reactiveVal(NULL)
    se_sweep_opts     <- reactiveVal(character(0))  # named: basename → full path
    se_ap_opts        <- reactiveVal(character(0))
    se_object         <- reactiveVal(NULL)
    # Each element: list(job, log_file, path, done, ok)
    run_jobs          <- reactiveVal(list())
    run_log_text      <- reactiveVal("")
    run_active        <- reactiveVal(FALSE)

    shinyFiles::shinyFileChoose(input, "add_files",   roots = .cc_roots, filetypes = "parquet")
    shinyFiles::shinyDirChoose(input,  "add_folder",  roots = .cc_roots)
    shinyFiles::shinyDirChoose(input,  "out_dir_btn", roots = .cc_roots)
    shinyFiles::shinyDirChoose(input,  "csv_in_dir",  roots = .cc_roots)
    shinyFiles::shinyDirChoose(input,  "csv_out_dir", roots = .cc_roots)
    shinyFiles::shinyFileChoose(input, "se_add_sweep",roots = .cc_roots, filetypes = "parquet")
    shinyFiles::shinyFileChoose(input, "se_add_ap",   roots = .cc_roots, filetypes = "parquet")

    # ── File management ────────────────────────────────────────────────────────
    observeEvent(input$add_files, {
      req(is.list(input$add_files))
      new <- shinyFiles::parseFilePaths(.cc_roots, input$add_files)$datapath
      parquet_paths(unique(c(parquet_paths(), new[file.exists(new)])))
    })
    observeEvent(input$add_folder, {
      req(is.list(input$add_folder))
      folder <- shinyFiles::parseDirPath(.cc_roots, input$add_folder)
      req(nzchar(folder))
      found <- list.files(folder, pattern = "\\.parquet$", full.names = TRUE)
      if (!length(found)) { showNotification("No .parquet files found.", type = "warning"); return() }
      parquet_paths(unique(c(parquet_paths(), found)))
      showNotification(sprintf("Added %d parquet(s).", length(found)), type = "message")
    })
    observeEvent(input$clear_files, { parquet_paths(character(0)); scan_results(NULL); scan_ap_cache(NULL) })
    observeEvent(input$out_dir_btn, {
      req(is.list(input$out_dir_btn))
      d <- shinyFiles::parseDirPath(.cc_roots, input$out_dir_btn)
      if (nzchar(d)) out_dir_path(d)
    })
    output$file_list_ui <- renderUI({
      paths <- parquet_paths()
      if (!length(paths)) return(helpText("No files loaded."))
      do.call(tags$ul, c(
        list(style = "font-size:0.8em; max-height:70px; overflow-y:auto; padding-left:1em;"),
        lapply(paths, function(p) tags$li(basename(p)))
      ))
    })
    output$out_dir_label <- renderText(out_dir_path() %||% "(none selected)")

    # ── CSV → Parquet ──────────────────────────────────────────────────────────
    observeEvent(input$csv_in_dir, {
      req(is.list(input$csv_in_dir))
      d <- shinyFiles::parseDirPath(.cc_roots, input$csv_in_dir)
      if (nzchar(d)) csv_in_path(d)
    })
    observeEvent(input$csv_out_dir, {
      req(is.list(input$csv_out_dir))
      d <- shinyFiles::parseDirPath(.cc_roots, input$csv_out_dir)
      if (nzchar(d)) csv_out_path(d)
    })
    output$csv_in_label  <- renderText(csv_in_path()  %||% "(none selected)")
    output$csv_out_label <- renderText(csv_out_path() %||% "(none selected)")
    observeEvent(input$csv_convert_btn, {
      inp <- csv_in_path(); out <- csv_out_path()
      if (is.null(inp) || !nzchar(inp)) { showNotification("Select CSV input folder.", type="error"); return() }
      if (is.null(out) || !nzchar(out)) { showNotification("Select output folder.", type="error"); return() }
      plate_id_arg <- trimws(input$csv_plate_id %||% "")
      out_parquet  <- file.path(out, paste0(basename(inp), ".parquet"))
      withProgress(message = "Converting CSVs to Parquet…", {
        tryCatch({
          csvToParquet(inp, out_parquet,
                       plate_id = if (nzchar(plate_id_arg)) plate_id_arg else NULL)
          parquet_paths(unique(c(parquet_paths(), out_parquet)))
          showNotification(paste0("Done: ", basename(out_parquet), " — added to queue."),
                           type = "message", duration = 8)
        }, error = function(e)
          showNotification(paste("Conversion failed:", e$message), type = "error", duration = 15))
      })
    })

    # ── Scan + pre-detect ──────────────────────────────────────────────────────
    observeEvent(input$scan_btn, {
      paths <- parquet_paths()
      if (!length(paths)) {
        showNotification("Load at least one raw trace Parquet first.", type = "error"); return()
      }
      scan_sw <- if (nzchar(trimws(input$scan_sweeps %||% "")))
        as.integer(na.omit(suppressWarnings(
          as.integer(trimws(strsplit(input$scan_sweeps, "[,;\\s]+", perl = TRUE)[[1]])))))
      else NULL

      withProgress(message = "Scanning wells…", value = 0.05, {
        result <- tryCatch(
          scanParquetForAPs(paths,
                            v_is_volts   = isTRUE(input$v_is_volts),
                            threshold_mV = input$scan_threshold %||% 0,
                            sweeps       = scan_sw),
          error = function(e) { showNotification(paste("Scan failed:", e$message), type="error"); NULL }
        )
        if (is.null(result)) return()

        n_likely <- sum(result$ap_likely)
        scan_results(result)

        # Auto-select well with most AP sweeps
        best      <- result[which.max(result$n_ap_sweeps), ]
        best_key  <- .cc_make_well_key(best$well_id, best$plate_id, best$parquet_path)
        selected_well_key(best_key)

        if (n_likely == 0) {
          showNotification(sprintf("%d wells scanned — none above threshold.", nrow(result)),
                           type = "warning")
          return()
        }

        # Pre-detect APs on all likely wells with default config
        likely_rows  <- result[result$ap_likely, , drop = FALSE]
        default_cfg  <- cc_config(v_is_volts = isTRUE(input$v_is_volts))
        cache        <- list()
        incProgress(0.05)

        for (i in seq_len(nrow(likely_rows))) {
          incProgress(0.85 / nrow(likely_rows),
                      message = sprintf("Pre-detecting APs (%d/%d)…", i, nrow(likely_rows)),
                      detail  = likely_rows$well_id[i])
          w   <- likely_rows[i, ]
          key <- .cc_make_well_key(w$well_id, w$plate_id, w$parquet_path)
          cache[[key]] <- tryCatch(
            previewAPDetection(w$parquet_path, w$well_id, w$plate_id,
                               sweeps = scan_sw, cfg = default_cfg),
            error = function(e) NULL
          )
        }
        scan_ap_cache(cache)

        # Suggest the sweep with most APs for the best well
        bc <- cache[[best_key]]
        if (!is.null(bc) && nrow(bc$sweeps) > 0 && "ap_count" %in% names(bc$sweeps)) {
          suggested_sweep(as.integer(bc$sweeps$sweep[which.max(bc$sweeps$ap_count)]))
        }

        showNotification(
          sprintf("%d wells, %d AP-likely — pre-detection done.", nrow(result), n_likely),
          type = "message", duration = 4
        )
      })
    })

    # Update suggested sweep when a new well is selected
    observeEvent(selected_well_key(), {
      key   <- selected_well_key(); req(key)
      cache <- scan_ap_cache()
      if (!is.null(cache) && !is.null(cache[[key]])) {
        sw_data <- cache[[key]]$sweeps
        if (nrow(sw_data) > 0 && "ap_count" %in% names(sw_data))
          suggested_sweep(as.integer(sw_data$sweep[which.max(sw_data$ap_count)]))
        else suggested_sweep(NULL)
      } else {
        suggested_sweep(NULL)
      }
    })

    # ── scan_summary: scan_results + per-well metrics extracted from cache ─────
    scan_summary <- reactive({
      result <- scan_results(); req(result)
      cache  <- scan_ap_cache()
      if (is.null(cache) || length(cache) == 0) return(result)

      extras <- lapply(names(cache), function(key) {
        res <- cache[[key]]
        if (is.null(res) || nrow(res$sweeps) == 0) return(NULL)
        info <- .cc_parse_well_key(key)
        sw   <- res$sweeps
        freq_vals <- sw$ap_freq_hz[is.finite(sw$ap_freq_hz)]
        amp_vals  <- if ("mean_ap_amplitude_mV" %in% names(sw))
          sw$mean_ap_amplitude_mV[is.finite(sw$mean_ap_amplitude_mV)] else numeric(0)
        base_vals <- if ("baseline_v_mV" %in% names(sw))
          sw$baseline_v_mV[sw$ap_count == 0 & is.finite(sw$baseline_v_mV)] else numeric(0)
        data.frame(
          well_id        = info$well_id,
          plate_id       = info$plate_id,
          total_aps      = sum(sw$ap_count, na.rm = TRUE),
          max_freq_hz    = if (length(freq_vals)) max(freq_vals) else NA_real_,
          baseline_vm_mV = if (length(base_vals)) mean(base_vals) else NA_real_,
          mean_ap_amp_mV = if (length(amp_vals))  mean(amp_vals)  else NA_real_,
          stringsAsFactors = FALSE
        )
      })
      extras_df <- dplyr::bind_rows(Filter(Negate(is.null), extras))
      if (nrow(extras_df) == 0) return(result)
      merge(result, extras_df, by = c("well_id", "plate_id"), all.x = TRUE)
    })

    # ── Plate selector & color chooser ─────────────────────────────────────────
    output$plate_selector_ui <- renderUI({
      scan_df <- scan_results(); req(scan_df)
      plates  <- unique(scan_df$plate_id)
      if (length(plates) <= 1) return(NULL)
      selectInput(ns("current_plate"), "Plate", choices = plates)
    })

    output$plate_color_ui <- renderUI({
      df    <- scan_summary(); req(df)
      avail <- names(.cc_color_vars)[.cc_color_vars %in% names(df)]
      if (length(avail) == 0) return(NULL)
      selectInput(ns("plate_color_var"), "Color by",
                  choices  = .cc_color_vars[avail],
                  selected = .cc_color_vars[avail][1])
    })

    current_plate <- reactive({
      scan_df <- scan_results(); req(scan_df)
      plates  <- unique(scan_df$plate_id)
      req(length(plates) > 0)
      if (length(plates) == 1 || is.null(input$current_plate)) plates[1]
      else input$current_plate
    })

    # ── Plate grid ─────────────────────────────────────────────────────────────
    output$plate_grid <- renderPlot({
      df     <- scan_summary(); req(df)
      cp     <- current_plate(); req(cp)
      df     <- df[df$plate_id == cp, , drop = FALSE]
      req(nrow(df) > 0)

      color_var <- input$plate_color_var %||% "n_ap_sweeps"
      if (!color_var %in% names(df)) color_var <- "n_ap_sweeps"

      df$row_idx  <- match(substr(df$well_id, 1, 1), LETTERS)
      df$col_num  <- as.integer(substr(df$well_id, 2, 3))
      df$.fill    <- df[[color_var]]

      sel <- selected_well_key()
      df$is_sel <- if (!is.null(sel)) {
        si <- .cc_parse_well_key(sel)
        df$well_id == si$well_id & df$plate_id == si$plate_id
      } else FALSE

      n_rows <- max(df$row_idx, na.rm = TRUE)
      n_cols <- max(df$col_num, na.rm = TRUE)
      lbl    <- names(.cc_color_vars)[.cc_color_vars == color_var]
      if (!length(lbl)) lbl <- color_var

      p <- ggplot2::ggplot(df, ggplot2::aes(col_num, row_idx, fill = .fill)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
        ggplot2::scale_fill_viridis_c(option = "plasma", name = lbl,
                                      na.value = "grey80") +
        ggplot2::scale_y_reverse(breaks = seq_len(n_rows),
                                  labels = LETTERS[seq_len(n_rows)]) +
        ggplot2::scale_x_continuous(breaks = seq(2, n_cols, by = 2),
                                    expand = ggplot2::expansion(0)) +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 7) +
        ggplot2::theme(
          legend.position   = "right",
          legend.key.height = ggplot2::unit(0.4, "cm"),
          legend.key.width  = ggplot2::unit(0.2, "cm"),
          legend.title      = ggplot2::element_text(size = 6),
          axis.text         = ggplot2::element_text(size = 6),
          panel.grid        = ggplot2::element_blank()
        )

      sel_df <- df[df$is_sel, , drop = FALSE]
      if (nrow(sel_df) > 0)
        p <- p + ggplot2::geom_tile(
          data = sel_df, ggplot2::aes(col_num, row_idx),
          fill = NA, colour = "white", linewidth = 1.2, inherit.aes = FALSE
        )
      p
    }, bg = "transparent")

    output$plate_tooltip <- renderUI({
      h <- input$plate_hover; if (is.null(h)) return(NULL)
      df  <- scan_summary(); req(df)
      cp  <- current_plate(); req(cp)
      df  <- df[df$plate_id == cp, , drop = FALSE]
      df$row_idx <- match(substr(df$well_id, 1, 1), LETTERS)
      df$col_num <- as.integer(substr(df$well_id, 2, 3))
      hit <- df[df$col_num == round(h$x) & df$row_idx == round(h$y), , drop = FALSE]
      if (nrow(hit) == 0) return(NULL)

      color_var <- input$plate_color_var %||% "n_ap_sweeps"
      val_str <- if (color_var %in% names(hit) && !is.na(hit[[color_var]][1]))
        sprintf("%.1f", hit[[color_var]][1]) else "—"
      lbl_name <- names(.cc_color_vars)[.cc_color_vars == color_var]
      if (!length(lbl_name)) lbl_name <- color_var

      lbl <- sprintf("%s  |  %s: %s", hit$well_id[1], lbl_name, val_str)
      div(style = sprintf(paste0(
            "position:absolute;left:%dpx;top:%dpx;",
            "background:rgba(0,0,0,0.78);color:white;",
            "padding:2px 7px;border-radius:3px;",
            "font-size:11px;pointer-events:none;white-space:nowrap;"),
          as.integer(round(h$coords_css$x)) + 8L,
          as.integer(round(h$coords_css$y)) - 28L),
        lbl)
    })

    observeEvent(input$plate_click, {
      scan_df <- scan_results(); req(scan_df)
      cp      <- current_plate(); req(cp)
      col_c   <- round(input$plate_click$x)
      row_c   <- round(input$plate_click$y)
      if (anyNA(c(col_c, row_c))) return()
      if (row_c < 1 || row_c > 16 || col_c < 1 || col_c > 24) return()
      clicked <- sprintf("%s%02d", LETTERS[row_c], col_c)
      df_cp   <- scan_df[scan_df$plate_id == cp, , drop = FALSE]
      hit     <- df_cp[df_cp$well_id == clicked, , drop = FALSE]
      if (nrow(hit) > 0)
        selected_well_key(.cc_make_well_key(hit$well_id[1], hit$plate_id[1], hit$parquet_path[1]))
    })

    output$selected_well_label_ui <- renderUI({
      sel <- selected_well_key()
      if (is.null(sel)) return(helpText("Click a well to select."))
      si  <- .cc_parse_well_key(sel)
      tags$p(style = "font-size:0.82em; margin:2px 0 0;",
             tags$b(si$well_id), " | ", si$plate_id)
    })

    # ── Traces + sweep slider ──────────────────────────────────────────────────
    well_info <- reactive({
      sel <- selected_well_key(); req(sel)
      .cc_parse_well_key(sel)
    })

    well_traces <- reactive({
      info <- well_info(); req(info$path, file.exists(info$path))
      arrow::open_dataset(info$path) |>
        dplyr::filter(well_id  == !!info$well_id,
                      plate_id == !!info$plate_id) |>
        dplyr::collect() |>
        dplyr::mutate(v_mV = recorded * (if (isTRUE(input$v_is_volts)) 1000 else 1))
    })

    output$sweep_ui <- renderUI({
      df    <- well_traces(); req(df)
      sweeps <- sort(unique(df$sweep))
      sugg   <- suggested_sweep()
      init   <- if (!is.null(sugg) && sugg %in% sweeps) sugg else max(sweeps)
      tagList(
        h6("Sweep"),
        sliderInput(ns("selected_sweep"), NULL,
                    min = min(sweeps), max = max(sweeps),
                    value = init, step = 1, ticks = FALSE)
      )
    })

    # ── AP detection (always live on selected sweep) ───────────────────────────
    # Debounce config inputs so sliders don't trigger on every pixel drag.
    # sweep + config together so any change re-runs detection once settled.
    cfg_inputs <- reactive({
      list(dvdt_thr      = input$dvdt_thr,
           smooth_ms     = input$smooth_ms,
           refractory_ms = input$refractory_ms,
           min_peak_v    = input$min_peak_v,
           min_height    = input$min_height,
           prominence    = input$prominence,
           show_all      = input$show_all_sweeps,
           sweep         = input$selected_sweep)
    }) |> debounce(400)

    ap_result <- reactive({
      info <- well_info(); cfg <- cfg_inputs(); req(info, cfg)

      # When "overlay all" is checked run on all scan sweeps; otherwise
      # run only on the selected sweep for speed.
      sweeps_arg <- if (isTRUE(cfg$show_all)) {
        scan_sw_str <- trimws(input$scan_sweeps %||% "")
        if (nzchar(scan_sw_str))
          as.integer(na.omit(suppressWarnings(
            as.integer(trimws(strsplit(scan_sw_str, "[,;\\s]+", perl = TRUE)[[1]])))))
        else NULL
      } else {
        if (!is.null(cfg$sweep)) as.integer(cfg$sweep) else NULL
      }

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
            min_height_mV       = .cc_opt_val(cfg$min_height),
            prominence_mV       = .cc_opt_val(cfg$prominence)
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
    output$trace_plot <- renderPlot({
      df       <- well_traces(); req(df)
      ap       <- ap_result()
      show_all <- isTRUE(input$show_all_sweeps)
      sw       <- if (show_all || is.null(input$selected_sweep)) NULL
                  else input$selected_sweep
      plot_df  <- if (is.null(sw)) df else df[df$sweep == sw, ]

      if (show_all) {
        p <- ggplot2::ggplot(plot_df,
               ggplot2::aes(t_s, v_mV, colour = as.numeric(sweep),
                            group = factor(sweep))) +
          ggplot2::geom_line(alpha = 0.55, linewidth = 0.4) +
          ggplot2::scale_color_viridis_c(option = "plasma", guide = "none")
      } else {
        info <- well_info()
        title <- paste0(info$well_id, "  •  ", info$plate_id,
                        "  •  ", basename(info$path))
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(t_s, v_mV)) +
          ggplot2::geom_line(colour = "steelblue", linewidth = 0.5) +
          ggplot2::labs(title = title) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 9, hjust = 0))
      }

      # Stimulus shading
      if (!is.null(ap) && nrow(ap$sweeps) > 0) {
        sw_row <- if (!show_all && !is.null(sw))
          ap$sweeps[ap$sweeps$sweep == sw, , drop = FALSE]
        else ap$sweeps[1, , drop = FALSE]
        if (all(c("stim_start_s", "stim_end_s") %in% names(sw_row))) {
          st <- sw_row[is.finite(sw_row$stim_start_s) & is.finite(sw_row$stim_end_s), ]
          if (nrow(st) > 0)
            p <- p + ggplot2::annotate("rect",
                       xmin = st$stim_start_s[1], xmax = st$stim_end_s[1],
                       ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "gold")
        }
      }

      # AP markers
      if (!is.null(ap) && nrow(ap$aps) > 0) {
        ap_plot <- ap$aps[!is.na(ap$aps$peak_time_s), , drop = FALSE]
        if (!show_all && !is.null(sw))
          ap_plot <- ap_plot[ap_plot$sweep == sw, , drop = FALSE]
        if (nrow(ap_plot) > 0)
          p <- p + ggplot2::geom_point(
            data = ap_plot, ggplot2::aes(x = peak_time_s, y = peak_v_mV),
            colour = "firebrick", fill = "firebrick",
            shape = 25, size = 2.2, inherit.aes = FALSE
          )
      }

      p + ggplot2::labs(x = "Time (s)", y = "Voltage (mV)") +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35))
    })

    # ── Phase plane ────────────────────────────────────────────────────────────
    output$phase_plot <- renderPlot({
      df    <- well_traces(); req(df)
      sw    <- input$selected_sweep; req(sw)
      df_sw <- df[df$sweep == sw & !is.na(df$v_mV), ]
      req(nrow(df_sw) > 3)
      df_sw <- df_sw[order(df_sw$t_s), ]
      dt_ms <- diff(df_sw$t_s) * 1000
      dvdt  <- diff(df_sw$v_mV) / dt_ms
      v_mid <- (df_sw$v_mV[-1] + df_sw$v_mV[-nrow(df_sw)]) / 2
      keep  <- dt_ms > 0 & dt_ms < 5
      phase_df <- data.frame(v_mV = v_mid[keep], dvdt = dvdt[keep])

      p <- ggplot2::ggplot(phase_df, ggplot2::aes(v_mV, dvdt)) +
        ggplot2::geom_path(colour = "steelblue", alpha = 0.7, linewidth = 0.45) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                            linewidth = 0.3, linetype = "dashed") +
        ggplot2::labs(x = "Voltage (mV)", y = "dV/dt (mV/ms)") +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35))

      ap <- ap_result()
      if (!is.null(ap) && nrow(ap$aps) > 0) {
        ap_sw <- ap$aps[ap$aps$sweep == sw &
                          !is.na(ap$aps$dvdt_max_v_mV) &
                          !is.na(ap$aps$dvdt_max_mV_per_ms), , drop = FALSE]
        if (nrow(ap_sw) > 0)
          p <- p +
            ggplot2::geom_point(data = ap_sw,
              ggplot2::aes(dvdt_max_v_mV, dvdt_max_mV_per_ms),
              colour = "firebrick", shape = 24, size = 2.8, fill = "firebrick",
              inherit.aes = FALSE) +
            ggplot2::geom_point(data = ap_sw,
              ggplot2::aes(dvdt_min_v_mV, dvdt_min_mV_per_ms),
              colour = "navy", shape = 25, size = 2.8, fill = "navy",
              inherit.aes = FALSE)
      }
      p
    })

    # ── AP count footer + sweep table ──────────────────────────────────────────
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

    output$sweep_table <- renderTable({
      ap <- ap_result(); req(ap, nrow(ap$sweeps) > 0)
      cols <- intersect(c("sweep", "stim_amp_cmd", "ap_count", "ap_freq_hz",
                          "baseline_v_mV", "steady_v_mV",
                          "mean_ap_amplitude_mV", "mean_dvdt_max_mV_per_ms"),
                        names(ap$sweeps))
      out <- ap$sweeps[, cols, drop = FALSE]
      if ("stim_amp_cmd" %in% names(out)) {
        out$stim_pA <- round(out$stim_amp_cmd * 1e12, 1); out$stim_amp_cmd <- NULL
      }
      for (col in setdiff(names(out), "sweep"))
        if (is.numeric(out[[col]])) out[[col]] <- round(out[[col]], 2)
      out
    }, digits = 2, striped = TRUE, hover = TRUE, spacing = "s",
    width = "100%", align = "r")

    # ── Run full analysis — one callr::r_bg per parquet (parallel) ───────────────
    output$run_log_out <- renderText(run_log_text())

    output$run_status_ui <- renderUI({
      jobs <- run_jobs()
      if (!isTRUE(run_active()) || !length(jobs)) return(NULL)
      n_done <- sum(sapply(jobs, function(j) isTRUE(j$done)))
      tags$p(style = "font-size:0.8em; color:#e67e22; margin:3px 0 0;",
             icon("circle-notch"),
             sprintf(" Running… %d / %d plate(s) done", n_done, length(jobs)))
    })

    # Poll all jobs while active
    observe({
      if (!isTRUE(run_active())) return()
      invalidateLater(800, session)
      jobs <- run_jobs()
      if (!length(jobs)) return()

      # Update done/ok flags for finished jobs
      jobs <- lapply(jobs, function(j) {
        if (isTRUE(j$done)) return(j)
        if (!j$job$is_alive()) {
          j$done <- TRUE
          j$ok   <- identical(tryCatch(j$job$get_exit_status(), error = function(e) -1L), 0L)
        }
        j
      })
      run_jobs(jobs)

      # Merge log tails — one section per parquet
      log_sections <- lapply(jobs, function(j) {
        lines <- if (file.exists(j$log_file))
          readLines(j$log_file, warn = FALSE) else character(0)
        status <- if (isTRUE(j$done)) (if (isTRUE(j$ok)) "[done]" else "[FAILED]")
                  else "[running]"
        hdr <- sprintf("=== %s %s ===", basename(j$path), status)
        paste(c(hdr, tail(lines, 60)), collapse = "\n")
      })
      run_log_text(paste(log_sections, collapse = "\n\n"))

      # All finished?
      if (all(sapply(jobs, function(j) isTRUE(j$done)))) {
        run_active(FALSE)
        n_ok  <- sum(sapply(jobs, function(j) isTRUE(j$ok)))
        n_err <- length(jobs) - n_ok
        out_dir <- out_dir_path()
        if (!is.null(out_dir)) {
          .cc_add_se_sweep(list.files(out_dir, pattern = "_sweep\\.parquet$", full.names = TRUE))
          .cc_add_se_ap(   list.files(out_dir, pattern = "_ap\\.parquet$",    full.names = TRUE))
        }
        if (n_err == 0)
          showNotification(sprintf("All %d plate(s) complete!", n_ok),
                           type = "message", duration = 6)
        else
          showNotification(sprintf("%d plate(s) done, %d failed — see log.", n_ok, n_err),
                           type = "error", duration = 12)
      }
    })

    observeEvent(input$stop_run_btn, {
      jobs <- run_jobs()
      alive <- Filter(function(j) !isTRUE(j$done) && j$job$is_alive(), jobs)
      for (j in alive) j$job$kill()
      run_active(FALSE)
      run_jobs(list())
      if (length(alive)) showNotification(
        sprintf("Stopped %d running job(s).", length(alive)), type = "warning")
    })

    observeEvent(input$run_btn, {
      if (isTRUE(run_active())) {
        showNotification("Analysis already running.", type = "warning"); return()
      }
      out_dir <- out_dir_path()
      if (is.null(out_dir) || !nzchar(out_dir)) {
        showNotification("Choose an output folder first.", type = "error"); return()
      }
      paths <- parquet_paths()
      if (!length(paths)) {
        showNotification("No parquet files loaded.", type = "error"); return()
      }

      cfg_args <- list(
        v_is_volts          = isTRUE(input$v_is_volts),
        dvdt_thr_mV_per_ms  = input$dvdt_thr  %||% 25,
        smooth_ms           = input$smooth_ms  %||% 0.2,
        refractory_ms       = input$refractory_ms %||% 2,
        min_peak_voltage_mV = input$min_peak_v %||% -20,
        min_height_mV       = .cc_opt_val(input$min_height),
        prominence_mV       = .cc_opt_val(input$prominence)
      )
      py_exe   <- tryCatch(reticulate::py_config()$python, error = function(e) "")
      pkg_root <- tryCatch(normalizePath(find.package("ephacRTools")), error = function(e) "")
      is_dev   <- nzchar(pkg_root) && dir.exists(file.path(pkg_root, "R"))

      .bg_func <- function(path, out_dir, cfg_args, py_exe, pkg_root, is_dev) {
        if (nzchar(py_exe)) Sys.setenv(RETICULATE_PYTHON = py_exe)
        Sys.setenv(PYTHONUNBUFFERED = "1")
        if (is_dev && nzchar(pkg_root)) pkgload::load_all(pkg_root, quiet = TRUE)
        else library(ephacRTools)
        cfg <- do.call(cc_config, cfg_args)
        runAPAnalysis(path, out_dir = out_dir, cfg = cfg)
      }

      jobs <- lapply(paths, function(p) {
        lf  <- tempfile(fileext = ".log")
        job <- tryCatch(
          callr::r_bg(
            func  = .bg_func,
            args  = list(path = p, out_dir = out_dir, cfg_args = cfg_args,
                         py_exe = py_exe, pkg_root = pkg_root, is_dev = is_dev),
            stdout = lf, stderr = lf, supervise = TRUE
          ),
          error = function(e) {
            writeLines(paste("Failed to start:", e$message), lf)
            NULL
          }
        )
        list(job = job, log_file = lf, path = p,
             done = is.null(job), ok = FALSE)
      })

      started <- sum(!sapply(jobs, function(j) isTRUE(j$done)))
      run_jobs(jobs)
      run_log_text("")
      run_active(TRUE)
      showNotification(
        sprintf("Started %d parallel job(s).", started),
        type = "message", duration = 4)
    })

    # ── Config export / import / reset ─────────────────────────────────────────
    output$export_cfg <- downloadHandler(
      filename = function()
        paste0("cc_config_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"),
      content = function(file) {
        jsonlite::write_json(
          list(dvdt_thr = input$dvdt_thr, smooth_ms = input$smooth_ms,
               refractory_ms = input$refractory_ms, min_peak_v = input$min_peak_v,
               min_height = if (is.na(input$min_height)) NULL else input$min_height,
               prominence = if (is.na(input$prominence)) NULL else input$prominence,
               v_is_volts = isTRUE(input$v_is_volts),
               scan_threshold = input$scan_threshold %||% 0),
          file, auto_unbox = TRUE, null = "null", pretty = TRUE
        )
      }
    )
    observeEvent(input$import_cfg, {
      req(input$import_cfg)
      cfg <- tryCatch(jsonlite::read_json(input$import_cfg$datapath,
                                          simplifyVector = TRUE),
                      error = function(e) NULL)
      if (is.null(cfg)) { showNotification("Invalid config JSON.", type = "error"); return() }
      if (!is.null(cfg$dvdt_thr))      updateSliderInput(session, "dvdt_thr",      value = cfg$dvdt_thr)
      if (!is.null(cfg$smooth_ms))     updateSliderInput(session, "smooth_ms",     value = cfg$smooth_ms)
      if (!is.null(cfg$refractory_ms)) updateSliderInput(session, "refractory_ms", value = cfg$refractory_ms)
      if (!is.null(cfg$min_peak_v))    updateSliderInput(session, "min_peak_v",    value = cfg$min_peak_v)
      updateNumericInput(session, "min_height",     value = cfg$min_height     %||% NA)
      updateNumericInput(session, "prominence",     value = cfg$prominence     %||% NA)
      updateNumericInput(session, "scan_threshold", value = cfg$scan_threshold %||% 0)
      if (!is.null(cfg$v_is_volts))
        updateCheckboxInput(session, "v_is_volts", value = isTRUE(cfg$v_is_volts))
      showNotification("Config loaded.", type = "message")
    })
    observeEvent(input$reset_cfg, {
      d <- .cc_cfg_defaults
      updateSliderInput(session,   "dvdt_thr",      value = d$dvdt_thr)
      updateSliderInput(session,   "smooth_ms",      value = d$smooth_ms)
      updateSliderInput(session,   "refractory_ms",  value = d$refractory_ms)
      updateSliderInput(session,   "min_peak_v",     value = d$min_peak_v)
      updateNumericInput(session,  "min_height",     value = d$min_height)
      updateNumericInput(session,  "prominence",     value = d$prominence)
      updateCheckboxInput(session, "v_is_volts",     value = d$v_is_volts)
      updateNumericInput(session,  "scan_threshold", value = d$scan_threshold)
      showNotification("Config reset to defaults.", type = "message")
    })

    # ── Step 5: Build SE ───────────────────────────────────────────────────────

    # Helpers that merge new paths into the named choice vectors and update UI
    .cc_add_se_sweep <- function(new_paths) {
      new_paths <- new_paths[file.exists(new_paths)]
      if (!length(new_paths)) return()
      cur   <- se_sweep_opts()
      added <- setdiff(new_paths, cur)
      if (!length(added)) return()
      merged <- c(cur, added)
      se_sweep_opts(merged)
      updateSelectizeInput(session, "se_sweep_files",
        choices  = setNames(merged, basename(merged)),
        selected = merged, server = FALSE)
    }

    .cc_add_se_ap <- function(new_paths) {
      new_paths <- new_paths[file.exists(new_paths)]
      if (!length(new_paths)) return()
      cur   <- se_ap_opts()
      added <- setdiff(new_paths, cur)
      if (!length(added)) return()
      merged <- c(cur, added)
      se_ap_opts(merged)
      updateSelectizeInput(session, "se_ap_files",
        choices  = setNames(merged, basename(merged)),
        selected = merged, server = FALSE)
    }

    # Manual file-picker additions
    observeEvent(input$se_add_sweep, {
      req(is.list(input$se_add_sweep))
      .cc_add_se_sweep(shinyFiles::parseFilePaths(.cc_roots, input$se_add_sweep)$datapath)
    })
    observeEvent(input$se_add_ap, {
      req(is.list(input$se_add_ap))
      .cc_add_se_ap(shinyFiles::parseFilePaths(.cc_roots, input$se_add_ap)$datapath)
    })

    # Scan output folder button
    observeEvent(input$se_scan_out, {
      out_dir <- out_dir_path()
      if (is.null(out_dir) || !nzchar(out_dir)) {
        showNotification("Set an output folder in Step 4 first.", type = "warning"); return()
      }
      sw <- list.files(out_dir, pattern = "_sweep\\.parquet$", full.names = TRUE)
      ap <- list.files(out_dir, pattern = "_ap\\.parquet$",    full.names = TRUE)
      if (!length(sw) && !length(ap)) {
        showNotification("No _sweep/_ap parquets found in output folder.", type = "warning"); return()
      }
      .cc_add_se_sweep(sw); .cc_add_se_ap(ap)
      showNotification(sprintf("Found %d sweep, %d AP parquet(s).",
                               length(sw), length(ap)), type = "message")
    })

    observeEvent(input$build_se_btn, {
      sw_paths <- input$se_sweep_files
      if (is.null(sw_paths) || !length(sw_paths)) {
        showNotification("Add at least one sweep parquet.", type = "error"); return()
      }
      ap_paths <- if (length(input$se_ap_files)) input$se_ap_files else NULL

      withProgress(message = "Building SingleCellExperiment…", {
        se <- tryCatch(
          prepareCCSE(sw_paths, ap_parquet_path = ap_paths),
          error = function(e) {
            showNotification(paste("Failed:", e$message), type = "error"); NULL
          }
        )
        if (!is.null(se)) {
          se_object(se)
          nm <- trimws(input$se_name %||% "CC_experiment")
          if (!nzchar(nm)) nm <- "CC_experiment"
          if (!is.null(ses)) {
            existing <- names(shiny::reactiveValuesToList(ses))
            if (nm %in% existing)
              nm <- make.unique(c(existing, nm), sep = " ")[length(existing) + 1L]
            ses[[nm]] <- se
          }
          showNotification(
            sprintf("SE built: %d wells × %d sweeps%s", ncol(se), nrow(se),
                    if (!is.null(ses)) paste0(" — added as '", nm, "'") else ""),
            type = "message", duration = 8
          )
        }
      })
    })

    output$se_status_ui <- renderUI({
      se <- se_object()
      if (is.null(se)) return(NULL)
      tags$p(style = "font-size:0.82em; color:#2a7; margin:4px 0 0;",
             icon("check-circle"),
             sprintf(" SE ready: %d wells, %d sweeps", ncol(se), nrow(se)))
    })

    output$download_se <- downloadHandler(
      filename = function() {
        nm <- trimws(input$se_name %||% "cc_se")
        if (!nzchar(nm)) nm <- "cc_se"
        paste0(nm, ".rds")
      },
      content = function(file) {
        se <- se_object()
        if (is.null(se)) {
          showNotification("Build the SE first.", type = "error"); return()
        }
        saveRDS(se, file)
      }
    )
  })
}
