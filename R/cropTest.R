library(shiny)
library(shinyFiles)

rows_384 <- LETTERS[1:16]   # A-P
cols_384 <- 1:24            # 1-24

# Row-major order: A01..A24, B01..B24, ..., P24
well_ids <- function() {
  as.vector(t(outer(rows_384, sprintf("%02d", cols_384), paste0)))
}

# Extract well like O9 from a string containing _O9-1_ etc.
# Returns "O09" (zero-padded) or NA
extract_well_384 <- function(x) {
  m <- regexpr("(?<=_)([A-P])([0-9]{1,2})(?:-[0-9]+)?(?=_)", x, perl = TRUE)
  if (m[1] == -1) return(NA_character_)
  token <- regmatches(x, m)

  row <- sub("^([A-P]).*$", "\\1", token)
  col <- as.integer(sub("^[A-P]([0-9]{1,2}).*$", "\\1", token))
  if (is.na(col) || col < 1 || col > 24) return(NA_character_)

  paste0(row, sprintf("%02d", col))
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* ---- Page: remove dead space ---- */
      body { overflow: hidden; }
      .container-fluid { padding-left: 10px; padding-right: 10px; }
      h2 { margin-top: 10px; margin-bottom: 10px; }

      /* ---- App layout ---- */
      .app-row {
        display: grid;
        grid-template-columns: 320px 1fr; /* compact sidebar + full main */
        gap: 12px;
        height: calc(100vh - 70px);      /* fill screen under title */
      }
      .sidebar {
        border: 1px solid #e6e6e6;
        border-radius: 12px;
        padding: 12px;
        background: #fff;
        overflow: auto;
      }
      .main {
        border: 1px solid #e6e6e6;
        border-radius: 12px;
        background: #fafafa;
        overflow: hidden;
        display: flex;
        flex-direction: column;
      }
      .main-header {
        padding: 10px 12px;
        border-bottom: 1px solid #e6e6e6;
        background: #fff;
        border-top-left-radius: 12px;
        border-top-right-radius: 12px;
      }
      .plate-wrap {
        flex: 1;
        overflow: auto;                 /* scrollbars live HERE */
        padding: 10px;
      }

      /* ---- Make the folder button big and full-width ---- */
      #imgdir {
        width: 100% !important;
        padding: 14px 12px !important;
        font-size: 15px !important;
        border-radius: 10px !important;
      }

      /* Make inputs full width */
      .shiny-input-container { width: 100% !important; }

      /* ---- Plate density ---- */
      :root{
        --gap: 2px;
        --cell: var(--cellUser, var(--cellAuto, 56px)); /* default looks good */
      }

      .plate-grid{
        width: max-content;             /* only as wide as needed */
        display: grid;
        grid-template-columns: repeat(24, var(--cell));
        grid-auto-rows: var(--cell);
        gap: var(--gap);
        align-content: start;
      }

      .well{
        width: var(--cell);
        height: var(--cell);
        border-radius: 7px;
        border: 1px solid #eee;
        background: #fff;
        cursor: pointer;
        position: relative;
        overflow: hidden;
        box-shadow: none;
      }
      .well img{
        width: 100%;
        height: 100%;
        object-fit: cover;
        display: block;
      }
      .missing{
        display:flex;
        align-items:center;
        justify-content:center;
        height:100%;
        color:#aaa;
        font-size:10px;
      }
      .well-label{
        position: absolute;
        bottom: 2px;
        right: 2px;
        background: rgba(0,0,0,0.55);
        color: #fff;
        padding: 1px 4px;
        font-size: 9px;
        border-radius: 6px;
        pointer-events: none;
      }

      .meta { color:#666; font-size:12px; }
      .small-mono { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; font-size: 12px; }
    ")),

    tags$script(HTML("
      function clickWell(wellId) {
        Shiny.setInputValue('well_click', wellId, {priority: 'event'});
      }

      // Auto-compute cell size to fit 24 cols into the available MAIN width,
      // but clamp to avoid microscopic tiles.
      function setCellAuto() {
        const plateWrap = document.querySelector('.plate-wrap');
        if (!plateWrap) return;

        const w = plateWrap.clientWidth;
        const gap = 2; // must match --gap
        let cell = Math.floor((w - (23 * gap) - 20) / 24); // -20 for padding/scrollbar

        // clamp
        cell = Math.max(34, Math.min(80, cell));
        document.documentElement.style.setProperty('--cellAuto', cell + 'px');
      }

      window.addEventListener('load', setCellAuto);
      window.addEventListener('resize', setCellAuto);
      $(document).on('shiny:value', function(){ setCellAuto(); });

      Shiny.addCustomMessageHandler('setCellSize', function(px){
        if (px === null) {
          document.documentElement.style.removeProperty('--cellUser');
          setCellAuto();
        } else {
          document.documentElement.style.setProperty('--cellUser', px + 'px');
        }
      });
    "))
  ),

  titlePanel("384-well plate viewer (folder picker)"),

  div(class = "app-row",
      div(class = "sidebar",
          shinyDirButton("imgdir", "Pick image folder…", "Select a folder"),
          div(style="height:8px;"),
          div(class="small-mono", verbatimTextOutput("dirpath")),

          checkboxInput("recursive", "Search recursively", TRUE),

          # visible and right next to the folder picker
          sliderInput("cell_px", "Zoom (well size)", min = 0, max = 140, value = 0, step = 1),
          div(class="meta", "Tip: 0 = auto-fit. Move slider to force a fixed size."),

          hr(),
          strong("Stats"),
          verbatimTextOutput("stats"),
          div(class="meta", "Well parsing: filenames containing e.g. _O9-1_ → O09")
      ),

      div(class = "main",
          div(class="main-header",
              tags$h4(style="margin:0;", "Plate (click a well to enlarge)")
          ),
          div(class="plate-wrap",
              uiOutput("plate_ui")
          )
      )
  )
)

server <- function(input, output, session) {

  roots <- c(
    Home = normalizePath("~", winslash = "/", mustWork = TRUE),
    Z = "Z:"
  )
  shinyDirChoose(input, "imgdir", roots = roots)

  resource_token <- "imgfolder"

  chosen_dir <- reactive({
    req(input$imgdir)
    p <- parseDirPath(roots, input$imgdir)
    if (length(p) != 1 || is.na(p) || !nzchar(p)) return(NULL)
    p <- normalizePath(p, winslash = "/", mustWork = FALSE)
    if (!dir.exists(p)) return(NULL)
    p
  })

  observeEvent(chosen_dir(), {
    d <- chosen_dir()
    req(!is.null(d))
    addResourcePath(resource_token, d)
  })

  output$dirpath <- renderText({
    d <- chosen_dir()
    if (is.null(d)) "(no folder selected)"
    else d
  })

  # Slider:
  # 0 => auto-fit
  observe({
    req(input$cell_px)
    if (input$cell_px == 0) session$sendCustomMessage("setCellSize", NULL)
    else session$sendCustomMessage("setCellSize", input$cell_px)
  })

  well_map <- reactive({
    d <- chosen_dir()
    ids <- well_ids()

    if (is.null(d)) {
      return(list(map = setNames(rep(NA_character_, length(ids)), ids),
                  files = character(0), parsed_ok = 0L))
    }

    files <- list.files(
      d,
      pattern = "\\.(jpg|jpeg)$",
      full.names = TRUE,
      recursive = isTRUE(input$recursive),
      ignore.case = TRUE
    )

    if (length(files) == 0) {
      return(list(map = setNames(rep(NA_character_, length(ids)), ids),
                  files = character(0), parsed_ok = 0L))
    }

    wells <- vapply(files, extract_well_384, character(1))
    ok <- !is.na(wells)
    parsed_ok <- sum(ok)

    m <- setNames(rep(NA_character_, length(ids)), ids)
    if (parsed_ok > 0) {
      files_ok <- files[ok]
      wells_ok <- wells[ok]
      keep <- !duplicated(wells_ok)
      m[wells_ok[keep]] <- files_ok[keep]
    }

    list(map = m, files = files, parsed_ok = parsed_ok)
  })

  output$stats <- renderText({
    wm <- well_map()
    paste0(
      "JPGs found: ", length(wm$files), "\n",
      "Parsed to valid wells: ", wm$parsed_ok, "\n",
      "Wells filled: ", sum(!is.na(wm$map)), " / 384"
    )
  })

  file_to_url <- function(abs_file) {
    d <- chosen_dir()
    req(!is.null(d))

    abs_file <- normalizePath(abs_file, winslash = "/", mustWork = FALSE)
    d <- normalizePath(d, winslash = "/", mustWork = FALSE)

    if (!grepl("/$", d)) d <- paste0(d, "/")
    if (!startsWith(abs_file, d)) return(NA_character_)

    rel <- substring(abs_file, nchar(d) + 1)
    rel <- gsub("^/+", "", rel)

    file.path(resource_token, rel)
  }

  output$plate_ui <- renderUI({
    wm <- well_map()
    m <- wm$map
    ids <- well_ids()

    tags$div(
      class = "plate-grid",
      lapply(ids, function(well) {
        fp <- m[[well]]
        if (is.na(fp)) {
          tags$div(
            class = "well",
            onclick = sprintf("clickWell('%s')", well),
            title = paste("Well", well, "(no image)"),
            tags$div(class = "missing", "No image"),
            tags$span(class = "well-label", well)
          )
        } else {
          url <- file_to_url(fp)
          if (is.na(url)) {
            tags$div(
              class = "well",
              onclick = sprintf("clickWell('%s')", well),
              title = paste("Well", well, "(path error)"),
              tags$div(class = "missing", "Path error"),
              tags$span(class = "well-label", well)
            )
          } else {
            tags$div(
              class = "well",
              onclick = sprintf("clickWell('%s')", well),
              title = paste("Well", well),
              tags$img(src = url, loading = "lazy"),
              tags$span(class = "well-label", well)
            )
          }
        }
      })
    )
  })

  observeEvent(input$well_click, {
    well <- input$well_click
    wm <- well_map()
    fp <- wm$map[[well]]

    if (is.na(fp)) {
      showModal(modalDialog(
        title = paste("Well", well),
        tags$p("No image mapped to this well."),
        easyClose = TRUE,
        size = "l"
      ))
      return()
    }

    url <- file_to_url(fp)
    if (is.na(url)) {
      showModal(modalDialog(
        title = paste("Well", well),
        tags$p("Could not map file path to a served URL."),
        tags$div(class = "meta", fp),
        easyClose = TRUE,
        size = "l"
      ))
      return()
    }

    showModal(modalDialog(
      title = paste("Well", well),
      tags$div(
        style = "text-align:center;",
        tags$img(
          src = url,
          style = "width: 750px; max-width: 95vw; height: auto; border-radius: 10px; border: 1px solid #ddd;"
        ),
        tags$div(class="meta", fp)
      ),
      easyClose = TRUE,
      size = "l",
      footer = modalButton("Close")
    ))
  })
}

shinyApp(ui, server)
