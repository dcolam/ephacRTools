library(shiny)
library(shinyFiles)
library(jsonlite)
library(jpeg)

rows_384 <- LETTERS[1:16]   # A-P
cols_384 <- 1:24            # 1-24

well_ids <- function() {
  as.vector(t(outer(rows_384, sprintf("%02d", cols_384), paste0)))
}

extract_well_384 <- function(x) {
  m <- regexpr("(?<=_)([A-P])([0-9]{1,2})(?:-[0-9]+)?(?=_)", x, perl = TRUE)
  if (m[1] == -1) return(NA_character_)
  token <- regmatches(x, m)
  row <- sub("^([A-P]).*$", "\\1", token)
  col <- as.integer(sub("^[A-P]([0-9]{1,2}).*$", "\\1", token))
  if (is.na(col) || col < 1 || col > 24) return(NA_character_)
  paste0(row, sprintf("%02d", col))
}

coords_to_well <- function(x, y) {
  col_num <- floor(x) + 1L
  row_idx <- 16L - floor(y)
  if (col_num < 1 || col_num > 24 || row_idx < 1 || row_idx > 16) return(NA_character_)
  paste0(rows_384[row_idx], sprintf("%02d", col_num))
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { overflow: hidden; margin: 0; padding: 0; }
      .container-fluid { padding: 6px 10px; }
      h2 { margin: 4px 0; font-size: 18px; }

      .app-row {
        display: grid;
        grid-template-columns: 260px 1fr;
        gap: 10px;
        height: calc(100vh - 54px);
      }

      .sidebar {
        border: 1px solid #e0e0e0;
        border-radius: 10px;
        padding: 10px;
        background: #fff;
        overflow-y: auto;
      }

      .well-info {
        font-family: ui-monospace, Menlo, Consolas, monospace;
        font-size: 13px;
        color: #333;
        min-height: 20px;
        padding: 2px 0;
      }

      /* Classifier tab */
      .classifier-wrap {
        padding: 10px;
        overflow-y: auto;
        height: calc(100vh - 110px);
      }
      .classifier-grid {
        display: grid;
        grid-template-columns: 1fr 300px;
        gap: 12px;
        align-items: start;
      }
      .card {
        background: white;
        border: 1px solid #e0e0e0;
        border-radius: 10px;
        padding: 12px;
      }
      .class-buttons { display: grid; gap: 8px; }
      .class-buttons button {
        width: 100%;
        padding: 10px;
        border-radius: 8px;
        font-size: 14px;
      }
      .img-stage {
        display: flex;
        justify-content: center;
        align-items: center;
        min-height: 380px;
        background: #f5f5f5;
        border: 1px dashed #ccc;
        border-radius: 10px;
      }
      .img-stage img {
        max-width: 100%;
        max-height: 65vh;
        border-radius: 8px;
      }
      .meta { color: #888; font-size: 12px; }
    "))
  ),

  titlePanel("384-well plate viewer + manual classifier"),

  div(class = "app-row",

    # --- Sidebar ---
    div(class = "sidebar",
      shinyDirButton("imgdir", "Pick image folder…", "Select a folder"),
      div(style = "margin-top: 6px; font-family: monospace; font-size: 11px; color: #666; word-break: break-all;",
          textOutput("dirpath")),
      checkboxInput("recursive", "Search recursively", TRUE),
      hr(),
      strong("Hover"),
      div(class = "well-info", textOutput("hover_info")),
      hr(),
      strong("Stats"),
      verbatimTextOutput("stats")
    ),

    # --- Main tabs ---
    tabsetPanel(
      id = "main_tabs",

      tabPanel("Plate viewer",
        plotOutput(
          "plate_plot",
          click = "plot_click",
          hover = hoverOpts("plot_hover", delay = 60, nullOutside = TRUE),
          width  = "100%",
          height = "calc(100vh - 108px)"
        )
      ),

      tabPanel("Manual classifier",
        div(class = "classifier-wrap",
          div(class = "classifier-grid",

            div(class = "card",
              div(class = "img-stage", uiOutput("classifier_image_ui")),
              div(style = "margin-top: 10px;",
                actionButton("skip_btn",  "Skip / Next",  class = "btn-default"),
                actionButton("undo_btn",  "Undo last",    class = "btn-warning"),
                span(class = "meta", style = "margin-left: 10px;",
                     textOutput("progress_txt", inline = TRUE))
              )
            ),

            div(class = "card",
              tags$h4(style = "margin-top: 0;", "Classes"),
              textInput("new_class", NULL, placeholder = "New class label"),
              actionButton("add_class",     "Add class",    class = "btn-primary"),
              actionButton("clear_classes", "Clear classes", class = "btn-danger"),
              hr(),
              uiOutput("class_buttons_ui"),
              hr(),
              tags$h4(style = "margin-top: 0;", "Results"),
              downloadButton("download_csv", "Download CSV"),
              actionButton("clear_results", "Clear results", class = "btn-danger"),
              div(style = "margin-top: 10px;", tableOutput("results_preview"))
            )
          )
        )
      )
    )
  )
)

# -------------------------------------------------------------------------
server <- function(input, output, session) {

  roots <- c(
    Home = normalizePath("~", winslash = "/", mustWork = TRUE),
    Z    = "Z:"
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
    req(chosen_dir())
    addResourcePath(resource_token, chosen_dir())
  })

  output$dirpath <- renderText({
    d <- chosen_dir()
    if (is.null(d)) "(no folder selected)" else d
  })

  # All JPGs in chosen folder
  all_files <- reactive({
    d <- chosen_dir()
    if (is.null(d)) return(character(0))
    list.files(d,
               pattern   = "\\.(jpg|jpeg)$",
               full.names = TRUE,
               recursive  = isTRUE(input$recursive),
               ignore.case = TRUE)
  })

  # well -> file path map
  well_map <- reactive({
    ids   <- well_ids()
    files <- all_files()
    m     <- setNames(rep(NA_character_, length(ids)), ids)
    if (length(files) == 0)
      return(list(map = m, parsed_ok = 0L, files = files))

    wells     <- vapply(files, extract_well_384, character(1))
    ok        <- !is.na(wells)
    parsed_ok <- sum(ok)
    if (parsed_ok > 0) {
      fw   <- files[ok]; ww <- wells[ok]; keep <- !duplicated(ww)
      m[ww[keep]] <- fw[keep]
    }
    list(map = m, parsed_ok = parsed_ok, files = files)
  })

  output$stats <- renderText({
    wm <- well_map()
    paste0(
      "JPGs found: ",          length(wm$files),      "\n",
      "Parsed to wells: ",     wm$parsed_ok,           "\n",
      "Wells with image: ",    sum(!is.na(wm$map)), " / 384"
    )
  })

  # Load all images into memory (cached per file list)
  loaded_images <- reactive({
    wm <- well_map()
    fps <- na.omit(unname(wm$map))
    imgs <- list()
    for (fp in fps) {
      imgs[[fp]] <- tryCatch(jpeg::readJPEG(fp), error = function(e) NULL)
    }
    imgs
  })

  # ---- Plate plot ----
  output$plate_plot <- renderPlot({
    wm   <- well_map()
    imgs <- loaded_images()
    m    <- wm$map

    par(mar = c(0, 0, 0, 0), bg = "#111111")
    plot(NULL,
         xlim = c(0, 24), ylim = c(0, 16),
         xlab = "", ylab = "", axes = FALSE,
         xaxs = "i", yaxs = "i")

    for (well in well_ids()) {
      row_idx <- which(rows_384 == substr(well, 1, 1))
      col_num <- as.integer(substr(well, 2, 3))
      x0 <- col_num - 1;  x1 <- col_num
      y0 <- 16 - row_idx; y1 <- y0 + 1

      fp  <- m[[well]]
      img <- if (!is.na(fp)) imgs[[fp]] else NULL
      if (!is.null(img)) {
        rasterImage(img, x0, y0, x1, y1, interpolate = TRUE)
      } else {
        rect(x0, y0, x1, y1, col = "#2a2a2a", border = NA)
      }
    }
  },
  width  = function() session$clientData$output_plate_plot_width,
  height = function() session$clientData$output_plate_plot_height,
  bg     = "#111111")

  # Hover -> well label in sidebar
  output$hover_info <- renderText({
    h <- input$plot_hover
    if (is.null(h)) return("")
    well <- coords_to_well(h$x, h$y)
    if (is.na(well)) "" else well
  })

  # Click -> enlarge modal (served via addResourcePath)
  file_to_url <- function(abs_file) {
    d <- chosen_dir(); req(!is.null(d))
    abs_file <- normalizePath(abs_file, winslash = "/", mustWork = FALSE)
    d        <- normalizePath(d,        winslash = "/", mustWork = FALSE)
    if (!endsWith(d, "/")) d <- paste0(d, "/")
    if (!startsWith(abs_file, d)) return(NA_character_)
    file.path(resource_token, substring(abs_file, nchar(d) + 1))
  }

  observeEvent(input$plot_click, {
    cl   <- input$plot_click
    well <- coords_to_well(cl$x, cl$y)
    if (is.na(well)) return()
    fp   <- well_map()$map[[well]]
    if (is.na(fp)) return()
    url  <- file_to_url(fp)
    if (is.na(url)) return()
    showModal(modalDialog(
      title = paste("Well", well),
      tags$div(style = "text-align: center;",
        tags$img(src = url,
                 style = "width: 600px; max-width: 90vw; height: auto; border-radius: 8px;")),
      easyClose = TRUE, size = "l", footer = modalButton("Close")
    ))
  })

  # -------------------------
  # Manual classifier tab
  # -------------------------

  rv <- reactiveValues(
    classes = c("Class A", "Class B"),
    shuffled = character(0),
    idx      = 0L,
    results  = data.frame(
      timestamp = character(0), file = character(0),
      well = character(0), label = character(0),
      stringsAsFactors = FALSE
    )
  )

  csv_path <- reactive({
    d <- chosen_dir()
    if (is.null(d)) return(NULL)
    file.path(d, "manual_classifications.csv")
  })

  # Load existing CSV when folder changes
  observeEvent(chosen_dir(), {
    cp <- csv_path()
    if (!is.null(cp) && file.exists(cp)) {
      existing <- tryCatch(read.csv(cp, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(existing) && nrow(existing) > 0) {
        for (col in c("timestamp", "file", "well", "label"))
          if (!col %in% colnames(existing)) existing[[col]] <- NA_character_
        rv$results <- existing[, c("timestamp", "file", "well", "label")]
      }
    } else {
      rv$results <- rv$results[0, ]
    }
  })

  # Build shuffle queue excluding already-annotated files
  observeEvent(all_files(), {
    files <- all_files()
    if (length(files) == 0) {
      rv$shuffled <- character(0); rv$idx <- 0L; return()
    }
    pool <- setdiff(files, rv$results$file)
    if (length(pool) == 0) pool <- files
    rv$shuffled <- sample(pool)
    rv$idx      <- 1L
  })

  current_file <- reactive({
    if (length(rv$shuffled) == 0 || rv$idx < 1) return(NA_character_)
    rv$shuffled[[rv$idx]]
  })

  advance <- function() {
    n <- length(rv$shuffled)
    if (n == 0) return()
    nxt <- rv$idx %% n + 1L
    if (nxt == 1L) rv$shuffled <- sample(rv$shuffled)
    rv$idx <- nxt
  }

  output$classifier_image_ui <- renderUI({
    if (length(rv$shuffled) == 0)
      return(div(class = "meta", "Pick a folder with JPGs first."))
    fp  <- current_file()
    if (is.na(fp)) return(div(class = "meta", "No image available."))
    url <- file_to_url(fp)
    if (is.na(url)) return(div(class = "meta", "Could not serve image URL."))
    tags$img(src = url)
  })

  output$class_buttons_ui <- renderUI({
    req(length(rv$classes) > 0)
    div(class = "class-buttons",
      lapply(rv$classes, function(lab) {
        tags$button(lab, class = "btn btn-success",
          onclick = sprintf(
            "Shiny.setInputValue('class_selected', %s, {priority: 'event'});",
            jsonlite::toJSON(lab, auto_unbox = TRUE)))
      })
    )
  })

  # Single observer — no stacking bug
  observeEvent(input$class_selected, {
    fp <- current_file(); req(!is.na(fp))
    row <- data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      file      = fp,
      well      = extract_well_384(basename(fp)),
      label     = input$class_selected,
      stringsAsFactors = FALSE
    )
    rv$results <- rbind(rv$results, row)
    cp <- csv_path()
    if (!is.null(cp))
      write.table(row, cp, sep = ",", row.names = FALSE,
                  col.names = !file.exists(cp), append = file.exists(cp))
    advance()
  }, ignoreInit = TRUE)

  observeEvent(input$add_class, {
    x <- trimws(input$new_class)
    if (nzchar(x) && !(x %in% rv$classes)) {
      rv$classes <- c(rv$classes, x)
      updateTextInput(session, "new_class", value = "")
    }
  })

  observeEvent(input$clear_classes, { rv$classes <- character(0) })
  observeEvent(input$skip_btn,      { advance() })

  observeEvent(input$undo_btn, {
    if (nrow(rv$results) == 0) return()
    last_file  <- rv$results$file[nrow(rv$results)]
    rv$results <- rv$results[-nrow(rv$results), , drop = FALSE]
    rv$shuffled <- c(last_file, rv$shuffled)
    rv$idx      <- 1L
    cp <- csv_path()
    if (!is.null(cp) && file.exists(cp)) {
      lines <- readLines(cp)
      if (length(lines) > 1) writeLines(lines[-length(lines)], cp)
    }
  })

  observeEvent(input$clear_results, {
    rv$results <- rv$results[0, ]
    cp <- csv_path()
    if (!is.null(cp) && file.exists(cp)) file.remove(cp)
  })

  output$progress_txt <- renderText({
    total     <- length(all_files())
    labeled   <- nrow(rv$results)
    remaining <- length(rv$shuffled) - rv$idx + 1L
    paste0("Labeled: ", labeled, " / ", total,
           if (remaining > 0) paste0("  |  Queue: ", remaining) else "")
  })

  output$results_preview <- renderTable({
    tail(rv$results[, c("timestamp", "well", "label")], 10)
  })

  output$download_csv <- downloadHandler(
    filename = function() paste0("manual_classification_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content  = function(file) write.csv(rv$results, file, row.names = FALSE)
  )
}

shinyApp(ui, server)
