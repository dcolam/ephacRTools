#' tinySEV.server
#'
#' @param objects A named list of (paths to)
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} objects.
#' @param uploadMaxSize The maximum upload size. Set to zero to disable upload.
#' @param maxPlot The maximum number of features to allow for plotting heatmaps.
#' @param feature.lists An optional named list of genes/features which will be
#'   flagged in the gene tab.
#' @param filelist A named list of downloadable files optionally associated with
#'   elements of `objects`, or a folder where to find these files
#' @param feature.listsTab Logical, whether to show the feature list tab
#' @param logins An optional dataframe containing possible logins. Must contain
#'   the columns "user" and "password_hash" (sodium-encoded). Not providing the
#'   argument disables login.
#'
#' @return A shiny server function.
#' @export
#' @import shiny ggplot2 SummarizedExperiment sechm waiter plotly shinyjs
#' @importFrom shinydashboard updateTabItems
#' @importFrom shinyjs showElement hideElement
#' @importFrom DT datatable renderDT
#' @importFrom ComplexHeatmap draw
#' @importFrom S4Vectors metadata
#' @importFrom shinyauthr loginServer
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath
tinySEV.server <- function(objects=NULL, uploadMaxSize=1000*1024^2, maxPlot=500,
                           feature.lists=list(), filelist=list(), logins=NULL,
                           feature.listsTab=TRUE){
  options(shiny.maxRequestSize=uploadMaxSize)

  if(!is.null(objects) && is.null(names(objects))){
    names(objects) <- paste("Object ",seq_along(objects))
    x <- sapply(objects, FUN=function(x){
      if(is.character(x))
        return(gsub("\\.SE\\.rds$|\\.rds$", "", basename(x), ignore.case=TRUE))
      return(NULL)
    })
    x[which(is.null(x))] <- names(objects)[which(is.null(x))]
    names(objects) <- make.unique(x, sep=" ")
  }

  if(!is.null(filelist) && length(filelist)>0){
    if(is.character(filelist)){
      filelist <- list.dirs(filelist, recursive=FALSE)
      filelist <- lapply(setNames(filelist,basename(filelist)), FUN=function(x){
        paste0(basename(x),"/",list.files(x))
      })
    }
  }

  grepGene <- function(x,g){
    if(!is.character(x)){
      g <- grepGene(row.names(x), g)
      return(x[g,drop=FALSE])
    }
    if(all(g %in% x)) return(g)
    g <- paste0("^",g,"\\.|^",g,"$|\\.",g,"$")
    g <- lapply(g,FUN=function(i) grep(i, x, value=TRUE, ignore.case=TRUE))
    return(unique(unlist(g)))
  }

  getDef <- function(se,var, choices=NULL){
    if(length(var)>1){
      y <- unlist(lapply(var, FUN=function(x) getDef(se,x)))
      y <- y[!sapply(y,is.null)]
      if(length(y)==0){
        if(is.null(choices)) return(NULL)
        return(choices[[1]])
      }
      return(y)
    }
    if(is.null(se@metadata$default_view[[var]])){
      if(is.null(choices)) return(NULL)
      return(choices[[1]])
    }
    se@metadata$default_view[[var]]
  }

  DEAselRender <- "function(item, escape){
    if (item.value === item.label) return '<div><b>' + item.label + '</b><div>';
    var oOut = '<div><b>' + item.value + ': </b><br/>';
    oOut = oOut + '<div style=\"padding: 0 10px 10px 20px;\">';
    return oOut + item.label + '</div></div>';
  }"

  function(input, output, session) {

    previous_sel <- reactiveVal(value=NULL)

    if(!is.null(logins)){
      credentials <- shinyauthr::loginServer(
        id = "login",
        data = logins,
        user_col = user,
        pwd_col = password_hash,
        sodium_hashed = TRUE,
      )

      observe({
        if (credentials()$user_auth) {
          shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
        } else {
          shinyjs::addClass(selector = "body", class = "sidebar-collapse")
        }
      })
    }

    mem_used <- reactiveVal()

    # Update on button click
    observeEvent(input$refreshMem, {
      used_mb <- round(sum(gc()[, 2]), 1)
      mem_used(paste(used_mb, "MB RAM used"))
    })

    # Initialize once at startup
    observe({
      gc()
      used_mb <- round(sum(gc()[, 2]), 1)
      mem_used(paste(used_mb, "MB RAM used"))
    })

    # Display the memory usage
    output$memoryUsage <- renderText({
      mem_used()
    })

    output$uploadMenu <- renderUI({
      if(!(uploadMaxSize>0)) return(NULL)
      menuSubItem("Upload object", tabName="tab_fileinput")
    })

    mergeFlists <- function(se){
      fl <- tryCatch(metadata(se)$feature.lists, error=function(e) NULL)
      if(is.null(fl)) fl <- list()
      fl <- c(fl, feature.lists[setdiff(names(feature.lists), names(fl))])
      fl <- fl[lengths(fl)>0]
      metadata(se)$feature.lists <- fl
      se
    }

    output$menu_genelist <- renderUI({
      if(!feature.listsTab) return(NULL)
      menuItem("Feature lists", tabName="tab_genelists")
    })
    output$maxGenes <- renderText(maxPlot)

    slider_initialized <- reactiveVal(FALSE)
    initialized <- reactiveVal(FALSE)


    SEinit <- function(x){
      if(!is.null(logins)) req(credentials()$user_auth)
      if(is.null(assayNames(x)))
        assayNames(x) <- paste0("assay",1:length(assays(x)))
      if(ncol(rowData(x))==0) rowData(x)$name <- row.names(x)

      updateSelectInput(session, "plate_id", choices=unique(colData(x)$Plate_ID),
                        selected=unique(colData(x)$Plate_ID)[1])
      updateSelectInput(session, "assay_id", choices=c(assayNames(x), colnames(as.data.frame(colData(x)))),
                        selected=assayNames(x)[1])
      updateSelectInput(session, "sweep_id", choices=unique(rowData(x)$Sweep),
                        selected=unique(rowData(x)$Sweep)[1])
      updateSelectizeInput(session, "sweep_group",
                               choices=unique(rowData(x)$Sweep))
      updateSelectizeInput(session, "group_by_meta",
                           choices=colnames(rowData(x)))
      updateSelectInput(session, "plate_id1", choices=unique(colData(x)$Plate_ID),
                        selected=unique(colData(x)$Plate_ID)[1])
      updateSelectInput(session, "plate_id3", choices=unique(colData(x)$Plate_ID),
                        selected=unique(colData(x)$Plate_ID)[1])
      updateSelectInput(session, "plate_id4", choices=unique(colData(x)$Plate_ID),
                        selected=unique(colData(x)$Plate_ID)[1])
      updateSelectInput(session, "assay_id1", choices=assayNames(x),
                        selected=assayNames(x)[1])
      updateSelectInput(session, "color_group1", choices=colnames(as.data.frame(colData(x))),
                        selected=FALSE)
      updateSelectInput(session, "facet_group1", choices=colnames(as.data.frame(colData(x))),
                        selected=FALSE)
      updateSelectizeInput(session, "group_by_meta1",
                           choices=colnames(rowData(x)))
      coldat <- colnames(as.data.frame(colData(x))[ , !purrr::map_lgl(as.data.frame(colData(x)), is.numeric)])

      updateSelectInput(session, "condition", choices=coldat)



        slider_vals <- round(assay(x, assayNames(x)[1]), 2)

        updateSliderInput(session, "selected_slider",
                          min = min(slider_vals, na.rm = TRUE),
                          max = max(slider_vals, na.rm = TRUE),
                          value = c(min(slider_vals, na.rm = TRUE), max(slider_vals, na.rm = TRUE)),
                          step = 0.01)

        slider_initialized <- reactiveVal(FALSE)

        coldat <- colnames(as.data.frame(colData(x))[ , purrr::map_lgl(as.data.frame(colData(x)), is.numeric)])
        coldata <- as.data.frame(colData(x))
        updateSelectInput(session, "clusterAssay", choices = assayNames(x))
        updateSelectInput(session, "clusterColData", choices = coldat)
        updateSelectInput(session, "clustercolor1", choices = colnames(coldata), selected = "cluster.tsne")
        updateSelectInput(session, "clustercolor2", choices = colnames(coldata), selected = "cluster.umap")
        updateSelectInput(session, "clustercolor3", choices = colnames(coldata), selected = "cluster.pca")

        allWells <- list()
        for(plate in unique(coldata$Plate_ID)){

          allWells[[plate]] <- subset(coldata, Plate_ID == plate)$Well
        }

        if (!initialized()) {
        if(is.null(unlist(selected_wells$data))){
          selwel <- character(0)
        }else{
          selwel <- selected_wells$data[[unique(colData(x)$Plate_ID)[1]]]
          #selwel <- character(0)
        }
        updateSelectizeInput(session, "selected_well",
                             choices=allWells,
                             selected=selwel,
                             server=T)
        updateSelectizeInput(session, "selected_well1",
                             choices=allWells,
                             selected=selwel,
                             server=T)
        updateSelectizeInput(session, "clusteredwells",
                             choices=allWells,
                             selected=selwel,
                             server=T)
        initialized(TRUE)
        }
        updateSelectizeInput(session, "seDataset", choices = names(SEs), selected = input$object)
        updateSelectInput(session, "rdsObject", choices = names(SEs), selected = input$object)
        updateSelectInput(session, "img_coldata_var",
                          choices  = c("None", colnames(as.data.frame(colData(x)))),
                          selected = "None")
        updateSelectizeInput(session, "img_hover_vars",
                             choices  = colnames(as.data.frame(colData(x))),
                             selected = character(0))
        updateSelectizeInput(session, "img_seDataset",  choices=names(SEs), selected=input$object)
        updateSelectizeInput(session, "cls_seDataset",  choices=names(SEs), selected=input$object)

        # Customize Object Groups selectors
        all_cd_cols <- colnames(as.data.frame(colData(x)))
        num_cd_cols <- colnames(as.data.frame(colData(x))[ , purrr::map_lgl(as.data.frame(colData(x)), is.numeric), drop=FALSE])
        rd_cols     <- colnames(as.data.frame(rowData(x)))
        updateSelectInput(session, "cond_preview_col", choices = all_cd_cols, selected = all_cd_cols[1])
        updateSelectInput(session, "row_col_select",   choices = rd_cols,     selected = rd_cols[1])
        updateSelectInput(session, "lp_col",           choices = rd_cols,     selected = rd_cols[1])
        updateSelectInput(session, "ag_assays",        choices = assayNames(x))
        updateSelectInput(session, "transform_col",    choices = num_cd_cols, selected = if (length(num_cd_cols) > 0) num_cd_cols[length(num_cd_cols)] else NULL)
        updateSelectInput(session, "fw_col",           choices = all_cd_cols, selected = all_cd_cols[1])
        updateSelectInput(session, "dr_assays",        choices = assayNames(x))
        updateSelectInput(session, "dr_colnames",      choices = num_cd_cols)
        updateSelectInput(session, "dr_color_col",     choices = all_cd_cols, selected = all_cd_cols[1])

        # Auto-load image folder stored in SE metadata
        folder <- tryCatch(S4Vectors::metadata(x)$image_path_jpgs, error = function(e) NULL)
        if (!is.null(folder) && nzchar(folder) && dir.exists(folder)) {
          addResourcePath("imgplate", folder)
          img_parent(folder)
        }

      # Warn user before closing browser tab once data is loaded
      session$sendCustomMessage("setWarnBeforeClose", TRUE)

      x <- mergeFlists(x)
      return(x)
    }

    SEs <- reactiveValues()
    for(nn in names(objects)) SEs[[nn]] <- objects[[nn]]
    updateSelectInput(session, "object", choices=names(objects))


    SE <- reactive({
      if(is.null(input$object) || input$object=="" ||
         is.null(SEs[[input$object]])) return(NULL)
      if(is.character(fp <- SEs[[input$object]])){
        base <- gsub("\\.se\\.rds", "", fp, ignore.case=TRUE)
        if(file.exists(paste0(base, "assays.h5"))){
          SEs[[input$object]] <- loadHDF5SummarizedExperiment(dirname(fp),
                                                              prefix=basename(base))
        }else{
          SEs[[input$object]] <- readRDS(SEs[[input$object]])
        }
      }

    SEinit(SEs[[input$object]])
    })

    flists <- reactive({
      if(is.null(SE())) return(NULL)
      metadata(SE())$feature.lists
    })

    observeEvent(input$file, {
      tryCatch({
        if(!is.null(input$file)){
          x <- readRDS(input$file$datapath)
          if(is(x,"SummarizedExperiment")){
            SEname <- gsub("\\.SE\\.rds$|\\.rds$","",
                           basename(input$file$name), ignore.case=TRUE)
            SEs[[SEname]] <- x
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))
          }else{
            stop("The object is not a SummarizedExperiment!")
          }
        }}, error=function(e){
          print(conditionMessage(e))
          print(traceback())
          showModal(modalDialog(easyClose=TRUE, title="Error with upload",
                                "The file was not recognized. Are you sure that it is a R .rds file?",
                                tags$pre(e)))
        })
    })

    output$importXl <- renderTable({
      req(input$fileEphys)

      df <- input$fileEphys
      df$size_MB <- round(df$size / (1024^2), 2)  # Convert bytes to MB and round to 2 decimals
      df <- df[, c("name", "size_MB")]   # Select relevant columns

      df
    }, striped = TRUE, spacing = "s", bordered = TRUE)

    output$importDB <- renderTable({
      req(input$fileDB)

      df <- input$fileDB
      df$size_MB <- round(df$size / (1024^2), 2)  # Convert bytes to MB and round to 2 decimals
      df <- df[, c("name", "size_MB")]   # Select relevant columns

      df
    }, striped = TRUE, spacing = "s", bordered = TRUE)



    observeEvent(input$fileEphys, {
      req(input$fileEphys)
      tryCatch({
        if (!is.null(input$fileEphys)) {

          print(input$fileEphys)
          print(length(input$fileEphys$datapath))

          # Fix ifelse() misuse: use plain if()
          #l_files <- if (length(input$fileEphys$datapath) > 1) {
          #  as.list(input$fileEphys$datapath)
          #} else {
          #  input$fileEphys$datapath
          #}

          print(input$loadEphys)

          PathFiles <- input$fileEphys$datapath

          print(PathFiles)
          withProgress(message = 'Loading Excel-Files into SE', value = 0, {
            incProgress(0.5, detail = "This may take a while...")
            gc()
            tryCatch({
              x <- tryCatch({
                prepareSE(PathFiles)
              }, error = function(e) {
                stop(paste("prepareSE failed:", conditionMessage(e)))
              })
              gc()
              # rest of your logic
            }, error = function(e) {
              showModal(modalDialog("Error", tags$pre(conditionMessage(e))))
            })  # this is probably where it fails

            print(x)
            if (is(x, "SingleCellExperiment")) {
              SEname <- input$se_id
              SEs[[SEname]] <- x
              SEinit(SEs[[SEname]])
              incProgress(0.75, detail = "Updating UI")
              updateSelectInput(session, "object", selected = SEname,
                                choices = union(names(objects), names(SEs)))
            }
            incProgress(1, detail = "SE loaded")
          })

        } else {
          stop("The object is not a SummarizedExperiment!")
        }

      }, error = function(e) {
        # Show detailed error message
        print(conditionMessage(e))
        print(traceback())
        showModal(modalDialog(
          title = "Error with upload",
          easyClose = TRUE,
          "An error occurred while processing the uploaded file. Details:",
          tags$pre(conditionMessage(e)),   # show actual message
          tags$pre(paste(capture.output(traceback()), collapse = "\n"))  # full traceback
        ))
      })
    })

    observeEvent(input$mergeSE, {
      tryCatch({
        if(!is.null(input$seDataset) & !is.null(input$fileDB)){
          #x <- readRDS(input$fileEphys$datapath)
          #l_files <- ifelse(length(input$fileDB$datapath) > 1,
          #                  as.list(input$fileDB$datapath),
          #                  input$fileDB$datapath)
          l_files <- input$fileDB$datapath

          withProgress(message = 'Loading Imaging Results', value = 0, {
            incProgress(0.5, detail = "This may take a while..")

            req(input$tabletype)
            for(tabletype in input$tabletype){



                #df_img <- tryCatch({
                #  prepareImgDF(l_files, analysis = tabletype)
                #}, error = function(e) {
                #  stop(paste("prepareSE failed:", conditionMessage(e)))
                #})

                # rest of your logic
              #}, error = function(e) {
              #  showModal(modalDialog("Error", tags$pre(conditionMessage(e))))
              #})

              df_img <- prepareImgDF(l_files, analysis = tabletype)
              SEname <- input$seDataset


              SEs[[SEname]] <-  mergeSEandImg(SEs[[SEname]], df_img,
                                              tableType = tabletype)

              print("imported")

            #se <- SEs[[SEname]]
            #se <- mergeSEandImg(se, df_img)
            #SEs[[SEname]] <- se
            SEinit(SEs[[SEname]])
            incProgress(0.75, detail = "Updating UI")
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))

            incProgress(1, detail = "SE updated")
            }
          })
        }

      }, error = function(e){
        # Print error details to console/logs
        message("Caught an error:")
        message(conditionMessage(e))
        message("Traceback:")
        print(traceback())

        # Show detailed error to user in a modal dialog
        showModal(modalDialog(
          title = "Error with upload",
          easyClose = TRUE,
          "An error occurred while processing the uploaded file. Details:",
          tags$pre(conditionMessage(e)),  # show error message only
          tags$pre(paste(capture.output(traceback()), collapse = "\n"))  # show traceback
        ))
      })
    })

    # ── Scan reactive ────────────────────────────────────────────────────────
    img_scan_rv <- reactiveValues(preview = NULL)

    observeEvent(input$img_scan_btn, {
      req(input$img_fileDB)
      tryCatch({
        withProgress(message = "Scanning databases\u2026", value = 0.5, {
          prev <- previewImgDB(input$img_fileDB$datapath,
                               analysis = input$img_tabletype %||% "pa")
          img_scan_rv$preview <- prev
          incProgress(1)
        })
      }, error = function(e) {
        showModal(modalDialog(title = "Scan Error", easyClose = TRUE,
          tags$pre(conditionMessage(e))))
      })
    })

    # Reset scan when new files are uploaded
    observeEvent(input$img_fileDB, { img_scan_rv$preview <- NULL })

    # ── Dynamic configuration UI ─────────────────────────────────────────────
    output$img_config_ui <- renderUI({
      prev <- img_scan_rv$preview
      if (is.null(prev))
        return(helpText(style="color:#aaa; font-size:11px; padding:4px 0;",
                        "Upload .db files and click Scan to configure the import."))

      pa_cols   <- prev$pa_cols
      meas_cols <- prev$meas_cols
      sel_vals  <- prev$selections
      rank_df   <- prev$image_ranks

      # Auto-detect column defaults
      auto <- function(candidates, fallback = pa_cols[1])
        Filter(function(x) x %in% pa_cols, candidates)[1] %||% fallback
      auto_plate   <- auto(c("Plate_ID", "plate_id", "plate"))
      auto_well    <- auto(c("Well", "well", "Well_ID"))
      auto_channel <- auto(c("Channel_Name", "channel", "Channel"))
      auto_sel     <- auto(c("Selection", "selection", "ROI"))

      # Auto-detect ROI: prefer "allSelected"; else first non-BF value
      non_bf_names <- names(sel_vals)[!grepl("\\sBF$", names(sel_vals))]
      auto_roi <- if ("allSelected" %in% names(sel_vals)) "allSelected" else
                  if (length(non_bf_names) > 0) non_bf_names[1] else character(0)

      # Build labelled choices for selection checkboxes (top 20)
      top_sel      <- head(sel_vals, 20)
      sel_choices  <- setNames(names(top_sel),
                               paste0(names(top_sel), "  (n=", top_sel, ")"))

      # Extra optional columns (everything not already mandatory)
      mandatory    <- unique(c(auto_plate, auto_well, auto_channel, auto_sel,
                               "Image_ID", "Selection_Area", "Number_of_Particles"))
      extra_avail  <- setdiff(pa_cols, mandatory)

      # Rank hint text
      rank_hint <- if (nrow(rank_df) > 0)
        paste(paste0("rank ", rank_df$rank, " \u2192 ",
                     rank_df$n, " Image_IDs"), collapse = " | ")
      else "no rank info"

      tagList(
        # ── Summary ──────────────────────────────────────────────────────────
        tags$p(style="font-size:11px; color:#555; margin-bottom:6px;",
               icon("circle-info"),
               sprintf(" %d rows | %d wells | channels: %s",
                       prev$n_rows, prev$n_wells,
                       paste(prev$channels, collapse=", "))),
        tags$h5(style="margin-top:4px;", "1. Column Mapping"),
        fluidRow(
          column(6,
            selectInput("img_col_plate",   "Plate ID column:",
                        choices=pa_cols, selected=auto_plate),
            selectInput("img_col_channel", "Channel column:",
                        choices=pa_cols, selected=auto_channel)
          ),
          column(6,
            selectInput("img_col_well",      "Well column:",
                        choices=pa_cols, selected=auto_well),
            selectInput("img_col_selection", "Selection column:",
                        choices=pa_cols, selected=auto_sel)
          )
        ),
        if (length(extra_avail) > 0)
          checkboxGroupInput("img_extra_cols",
                             "Extra metadata columns to keep (optional):",
                             choices=extra_avail, inline=TRUE),
        hr(),
        # ── ROI labeling & filter ─────────────────────────────────────────────
        tags$h5("2. ROI Labeling & Filter"),
        checkboxInput("img_apply_corr_sel",
          "Label ROIs by area (adds CorrSel column: Hole_ROI / background_ROI)",
          value = TRUE),
        helpText(style="font-size:10px; margin-top:-4px;",
          "For each unique Selection value, smallest area \u2192 Hole_ROI, ",
          "largest \u2192 background_ROI. Resolves cases where the same ",
          "Selection name appears at two areas (e.g. BF hole vs inverted ROI)."),
        conditionalPanel("input.img_apply_corr_sel == true",
          checkboxGroupInput("img_corr_sel_filter",
            "Keep only CorrSel:",
            choices  = c("Hole_ROI", "background_ROI"),
            selected = "Hole_ROI",
            inline   = TRUE)
        ),
        tags$p(style="font-size:11px; color:#777; margin-top:6px; margin-bottom:2px;",
               "Additional filter by raw Selection name (optional):"),
        checkboxGroupInput("img_roi_selections", NULL,
                           choices  = sel_choices,
                           selected = auto_roi),
        hr(),
        # ── Fluorescence image detection ──────────────────────────────────────
        tags$h5("3. Fluorescence Image Detection"),
        helpText(style="font-size:10px;",
          "Images within each well are ranked by Image_ID (rank 1 = lowest). ",
          "Set the fluorescence rank per database — useful when different ",
          "databases in the same upload have different Image_ID patterns."),
        radioButtons("img_fluor_method", NULL,
          choices  = c("Keep only fluorescence (by rank)" = "rank",
                       "Keep all images"                  = "all"),
          selected = "rank", inline = TRUE),
        conditionalPanel("input.img_fluor_method == 'rank'",
          if (!is.null(prev$per_db) && length(prev$per_db) > 0) {
            tagList(
              tags$p(style="font-size:11px; color:#555; margin-bottom:4px;",
                     "Fluorescence rank per database:"),
              do.call(tagList, lapply(seq_along(prev$per_db), function(i) {
                db   <- prev$per_db[[i]]
                rdf  <- db$image_ranks
                hint <- if (nrow(rdf) > 0)
                  paste(paste0("rank ", rdf$rank, ": ", rdf$n, " imgs"),
                        collapse = " | ")
                else "no rank info"
                fluidRow(
                  column(5,
                    tags$p(style="font-size:11px; margin-top:8px; white-space:nowrap;
                                  overflow:hidden; text-overflow:ellipsis;",
                           tags$abbr(title=db$name, db$name))),
                  column(4,
                    tags$p(style="font-size:10px; color:#888; margin-top:8px;",
                           hint)),
                  column(3,
                    numericInput(paste0("img_fluor_rank_", i), NULL,
                                 value=2L, min=1L, max=10L, step=1L))
                )
              }))
            )
          } else {
            fluidRow(column(4,
              numericInput("img_fluor_rank_1", "Fluorescence = rank:",
                           value=2L, min=1L, max=10L, step=1L)))
          }
        ),
        hr(),
        # ── Measurement columns ───────────────────────────────────────────────
        tags$h5("4. Particle Measurements to Import"),
        checkboxGroupInput("img_meas_cols", NULL,
          choices  = meas_cols,
          selected = intersect(c("Area","Mean","StdDev","IntDen"), meas_cols),
          inline   = TRUE),
        helpText(style="font-size:10px;",
                 "Number_of_Particles and Selection_Area are always included."),
        hr(),
        # ── Import ────────────────────────────────────────────────────────────
        actionButton("img_mergeSE", "Import", class="btn-primary btn-block"),
        verbatimTextOutput("img_import_status"),
        helpText(style="font-size:11px; color:#aaa; margin-top:4px;",
                 "Imported particles are available in the Auto Classification pipeline."),
        hr(),
        tags$p(style="font-size:12px; color:#888; margin-bottom:4px;",
               "Optional \u2014 merge aggregated data into an SE object:"),
        fluidRow(
          column(8,
            selectizeInput("img_seDataset", "SE:", choices=c(), multiple=FALSE)
          ),
          column(4,
            br(),
            actionButton("img_connectSE", "Connect to SE",
                         class="btn-default btn-sm btn-block")
          )
        )
      )
    })

    # ── Import handler ───────────────────────────────────────────────────────
    observeEvent(input$img_mergeSE, {
      req(input$img_fileDB, input$img_tabletype)
      tryCatch({
        l_files         <- input$img_fileDB$datapath
        ttype           <- input$img_tabletype[[1]] %||% "pa"
        plate_col       <- input$img_col_plate     %||% "Plate_ID"
        well_col        <- input$img_col_well      %||% "Well"
        channel_col     <- input$img_col_channel   %||% "Channel_Name"
        selection_col   <- input$img_col_selection %||% "Selection"
        extra_cols      <- input$img_extra_cols    %||% character(0)
        meas_sel        <- input$img_meas_cols     %||% c("Area", "Mean", "IntDen")
        roi_sel         <- if (length(input$img_roi_selections) > 0)
                             input$img_roi_selections else NULL
        apply_cs        <- isTRUE(input$img_apply_corr_sel)
        corr_sel_filter <- if (apply_cs && length(input$img_corr_sel_filter) > 0)
                             input$img_corr_sel_filter else NULL

        # Per-DB fluor ranks — one numericInput per DB (img_fluor_rank_1, _2, ...)
        fluor_ranks <- if (identical(input$img_fluor_method, "rank")) {
          n_dbs <- length(l_files)
          ranks <- vapply(seq_len(n_dbs), function(i) {
            as.integer(input[[paste0("img_fluor_rank_", i)]] %||% 2L)
          }, integer(1))
          setNames(ranks, as.character(seq_len(n_dbs)))
        } else NULL

        withProgress(message = "Importing particle data", value = 0.3, {
          df <- prepareImgDF(
            l_files,
            analysis        = ttype,
            plate_col       = plate_col,
            well_col        = well_col,
            channel_col     = channel_col,
            selection_col   = selection_col,
            extra_id_cols   = extra_cols,
            meas_cols       = meas_sel,
            fluor_rank      = NULL,
            fluor_ranks     = fluor_ranks,
            roi_selections  = roi_sel,
            apply_corr_sel  = apply_cs,
            corr_sel_filter = corr_sel_filter,
            aggregate       = FALSE,
            cleanNames      = TRUE
          )
          cls_rv$particles  <- df
          cls_rv$filtered   <- NULL
          cls_rv$aggregated <- NULL
          cls_rv$scored     <- NULL
          cls_rv$classified <- NULL
          incProgress(1)
        })
      }, error = function(e) {
        showModal(modalDialog(title = "Import Error", easyClose = TRUE,
          tags$pre(conditionMessage(e))))
      })
    })

    output$img_import_status <- renderText({
      df <- cls_rv$particles
      if (is.null(df)) return("No data imported yet.")
      n_w   <- length(unique(paste(df$Plate_ID, df$Well)))
      n_ch  <- length(unique(df$Channel_Name))
      n_pl  <- length(unique(df$Plate_ID))
      sels   <- if ("Selection"  %in% names(df)) paste(unique(df$Selection),  collapse=", ") else "\u2014"
      itype  <- if ("Image_Type" %in% names(df)) paste(unique(df$Image_Type), collapse=", ") else "not set"
      corrsl <- if ("CorrSel"    %in% names(df)) paste(unique(df$CorrSel),    collapse=", ") else "not applied"
      paste0(nrow(df), " particles | ", n_w, " wells | ",
             n_ch, " channel(s) | ", n_pl, " plate(s)\n",
             "Selections: ", sels, "\n",
             "CorrSel:    ", corrsl, "\n",
             "Image types: ", itype, "\n",
             "Ready for Auto Classification.")
    })

    # ── RDS particle table loader ────────────────────────────────────────────
    observeEvent(input$img_load_rds, {
      req(input$img_rds_file)
      tryCatch({
        df <- readRDS(input$img_rds_file$datapath)
        if (!is.data.frame(df))
          stop("The .rds file must contain a data.frame.")
        cls_rv$particles  <- df
        cls_rv$filtered   <- NULL
        cls_rv$aggregated <- NULL
        cls_rv$scored     <- NULL
        cls_rv$classified <- NULL
      }, error = function(e) {
        showModal(modalDialog(title = "Load Error", easyClose = TRUE,
          tags$pre(conditionMessage(e))))
      })
    })

    output$img_rds_status <- renderText({
      req(input$img_rds_file)
      df <- cls_rv$particles
      if (is.null(df)) return("Not loaded.")
      n_w  <- length(unique(paste(df$Plate_ID %||% "?", df$Well %||% "?")))
      n_ch <- if ("Channel_Name" %in% names(df)) length(unique(df$Channel_Name)) else "?"
      n_pl <- if ("Plate_ID"    %in% names(df)) length(unique(df$Plate_ID))    else "?"
      paste0(nrow(df), " rows | ", n_w, " wells | ",
             n_ch, " channel(s) | ", n_pl, " plate(s)")
    })

    observeEvent(input$img_connectSE, {
      req(input$img_seDataset, input$img_fileDB, input$img_tabletype)
      tryCatch({
        l_files       <- input$img_fileDB$datapath
        plate_col     <- input$img_col_plate     %||% "Plate_ID"
        well_col      <- input$img_col_well      %||% "Well"
        channel_col   <- input$img_col_channel   %||% "Channel_Name"
        selection_col <- input$img_col_selection %||% "Selection"
        withProgress(message = "Connecting to SE", value = 0.3, {
          for (tabletype in input$img_tabletype) {
            df_img <- prepareImgDF(l_files,
              analysis      = tabletype,
              plate_col     = plate_col,
              well_col      = well_col,
              channel_col   = channel_col,
              selection_col = selection_col,
              aggregate     = TRUE,
              cleanNames    = TRUE)
            SEname <- input$img_seDataset
            SEs[[SEname]] <- mergeSEandImg(SEs[[SEname]], df_img, tableType=tabletype)
            SEinit(SEs[[SEname]])
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))
          }
          incProgress(1)
        })
      }, error = function(e) {
        showModal(modalDialog(title="Error", easyClose=TRUE,
          tags$pre(conditionMessage(e))))
      })
    })

    output$img_importDB_preview <- renderTable({
      req(input$img_fileDB)
      df <- input$img_fileDB
      df$size_MB <- round(df$size/(1024^2), 2)
      df[, c("name","size_MB")]
    }, striped=TRUE, spacing="s", bordered=TRUE)

    observeEvent(input$dataset_button, {
      tryCatch({
        if(!is.null(input$datasets)){

          choices <-list("Human Adrenal Glands" = "se_hAG",
                         "Primary Neurons"  = "se_pn",
                         "iPSC-Tricultures" = "se_iN",
                         "ROMK" = "se_romk")
          selected_label <- names(choices)[choices %in% input$datasets]
          toimport <- setdiff(selected_label, names(SEs))

          for(i in seq_along(toimport)){

            #dataset <- input$datasets[i]
            #SEname <- selected_label[i]
            SEname <- toimport[i]
            dataset <- choices[[SEname]]

            if(!(SEname %in% names(SEs))){
              initialized(FALSE)

            data(list = dataset, package = "ephacRTools", envir = environment())

            SEs[[SEname]] <- get(dataset, envir = environment())
            SEinit(SEs[[SEname]])
            print(SEname)
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))

            }
          }

        }}, error=function(e){
          print(conditionMessage(e))
          print(traceback())
          showModal(modalDialog(easyClose=TRUE, title="Error",
                                "Choose at least one dataset.",
                                tags$pre(e)))
        })
    })

    ############
    ### Overview tabs

    output$objOverview <- renderUI({
      if(!is.null(logins)) req(credentials()$user_auth)
      if(is.null(SE())){
        return(box(width=12, tags$p("No object loaded.")))
      }
      desImg <- ""
      ff <- NULL
      if(!is.null(filelist) && length(filelist[[input$object]])>0){
        ff <- filelist[[input$object]]
        if(length(wDes <- which(basename(ff)=="design.png"))>0){
          desImg <- tags$img(src=gsub("^www/","",ff[[head(wDes,1)]]))
          ff <- ff[-head(wDes,1)]
        }
      }
      md <- metadata(SE())
      md <- md[intersect(c("title","name","source"),names(md))]
      tagList(
        box(width=12, title="Object overview",
            tags$p(metadata(SE())$description),
            tagList(
              tags$p(
                paste("A SummarizedExperiment with", ncol(SE()), "samples and",
                      length(unique(SE()$Plate_ID)), "electrophysiology plates:")
              ),
              tags$ul(
                lapply(unique(SE()$Plate_ID), function(pid) {
                  tags$li(pid)
                })
              )
            ),
            tags$ul(lapply(names(md), FUN=function(x){
              tags$li(tags$b(x), tags$span(md[[x]]))
            })),
            tags$p(desImg)
        ),
        box(width=12, title="Associated files",
            tags$ul(lapply(ff, FUN=function(x){
              tags$li(tags$a(href=x, basename(x), download=gsub(" ","_",basename(x))))
            })))
      )
    })

    output$features <- renderDT({
      if(is.null(SE())) return(NULL)
      RD <- rowData(SE())
      RD <- RD[,unlist(sapply(RD, is.vector)),drop=FALSE]
      datatable( as.data.frame(RD), filter="top", class="compact",
                 options=list( pageLength=30, dom = "fltBip" ),
                 extensions=c("ColReorder") )
    }, server = TRUE)

    output$samples <- renderDT({
      if(is.null(SE())) return(NULL)
      datatable( as.data.frame(colData(SE())), filter="top", class="compact",
                 options=list( pageLength=30, dom = "fltBip" ),
                 extensions=c("ColReorder") )
    })


    ############
    ### BEGIN Plotting TABS

    observeEvent(input$group_by_meta, {
      req(SE(), input$group_by_meta)

      # Extract metadata column values from SE() — assuming SE() has colData with metadata
      meta_values <- NULL

      # Check if the metadata column exists in SE colData
      if (input$group_by_meta %in% colnames(SummarizedExperiment::rowData(SE()))) {
        meta_values <- unique(SummarizedExperiment::rowData(SE())[[input$group_by_meta]])
        meta_values <- sort(meta_values)  # optional: sort choices
      }

      # Update the sweep_group selectInput choices accordingly
      updateSelectInput(session, "sweep_group",
                        choices = meta_values,
                        selected = NULL)
      updateSelectInput(session, "sweep_id",
                        choices = meta_values,
                        selected = 1)
    })



        output$plate_view <- renderPlotly({

          req(input$assay_id)
          if(!is.null(SE())){

            if(!input$assay_id %in% assayNames(SE())){
              assayNameSham <- assayNames(SE())[1]
            }else{assayNameSham <-input$assay_id}
            assayName <- input$assay_id

            melted.dat <- sechm::meltSE(SE()[,SE()$Plate_ID == input$plate_id],
                                          features = row.names(rowData(SE())),
                                          assayName = assayNameSham,
                                          rowDat.columns = c(input$sweep_id, input$group_by_meta))

            melted.dat <- melted.dat[melted.dat[[input$group_by_meta]] %in% input$sweep_id, ]
            letter2number <- function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}
            melted.dat$RowNum <- sapply(melted.dat$Row, function(c){letter2number(c)})

            if(is.numeric(melted.dat[[assayName]])){
              numID <- TRUE
              if(input$assay_option == "raw"){
                melted.dat[[assayName]] <- melted.dat[[assayName]]
                legend <- assayName
              }
              if(input$assay_option == "log10"){
                melted.dat[[assayName]] <- log10(abs(melted.dat[[assayName]]))
                legend <- paste("log10(", assayName, ")", sep="")
              }
              if(input$assay_option == "scale"){
                melted.dat[[assayName]] <- scale(melted.dat[[assayName]], center = T)
                legend <- paste("Z-scaled(", assayName, ")", sep="")
              }

            } else{
              numID <- FALSE
              melted.dat[[assayName]] <- as.factor((melted.dat[[assayName]]))
              legend <- assayName
            }



            if(length(get_wells(input$plate_id)) > 0){
              melted.dat$is_selected <- ifelse(melted.dat$Well %in% get_wells(input$plate_id), TRUE, FALSE)

            }else{
              melted.dat$is_selected <- TRUE

            }
            #melted.dat$is_selected <- ifelse(melted.dat$Well %in% input$selected_well, TRUE, FALSE)

            melted.dat$key_combined <- paste(melted.dat$Well, melted.dat$Plate_ID, sep = ", ")
            p <- ggplot(melted.dat, aes(x = as.numeric(Column), y = Row, key = key_combined, fill = .data[[assayName]])) +
              geom_tile() +
              scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2)) +
              scale_x_continuous(breaks = 1:24) +
              scale_y_discrete(limits = rev) +
              geom_text(aes(label = paste(Row, Column, sep=""), alpha = is_selected), color = "white") +
              theme_minimal()


            if(numID){p <- p + scale_fill_viridis_c(option = "magma")}

            # Register click events and return the plot
            p_plotly <- ggplotly(p, source = "plate_plot") %>%
              layout(dragmode = "select") %>%  # Makes box selection the default
              config(
                displayModeBar = TRUE,
                modeBarButtonsToAdd = list('select2d'),  # Ensures the select tool is in the toolbar
                displaylogo = FALSE
              )
            #plotly::event_register(p_plotly, "plotly_click")
            #slider_initialized(TRUE)

            return(p_plotly)
          }
        })
        jqui_resizable(ui="#plate_view")



        ############
        ### Definition of Well Selection

        # --- Reactive storage ---
        selected_wells <- reactiveValues(
          data = data.frame(
            well = character(),
            plate_id = character(),
            stringsAsFactors = FALSE
          )
        )

        get_wells <- function(plate_id) {
          unique(selected_wells$data[selected_wells$data$plate_id %in% plate_id,]$well)
        }

        set_wells <- function(plate_id, wells) {
          current <- get_wells(plate_id)
          if (!identical(sort(current), sort(wells))) {
            selected_wells$data <- rbind(
              selected_wells$data[selected_wells$data$plate_id != plate_id, ],
              data.frame(well = wells, plate_id = plate_id, stringsAsFactors = FALSE)
            )
          }
        }

        update_all_select_inputs <- function() {
          if (!is.null(input$plate_id) && nzchar(input$plate_id)) {
            updateSelectizeInput(session, "selected_well", selected = get_wells(input$plate_id))
          }
          if (!is.null(input$plate_id1) && nzchar(input$plate_id1)) {
            updateSelectizeInput(session, "selected_well1", selected = get_wells(input$plate_id1))
          }
          if (!is.null(input$plate_id3) && nzchar(input$plate_id3)) {
            updateSelectizeInput(session, "clusteredwells", selected = get_wells(input$plate_id3))
          }
        }

        handle_plotly_interaction <- function(event_data, type = c("click", "selected"), source_id) {
          type <- match.arg(type)
          print(paste0("plotly_", type))
          event <- plotly::event_data(paste0("plotly_", type), source = source_id)
          req(event)
          print(event)
          df <- do.call(rbind, lapply(event$key, function(k) {
            parts <- strsplit(k, ",\\s*")[[1]]
            data.frame(well = parts[1], plate_id = parts[2], stringsAsFactors = FALSE)
          }))

          plate_ids <- unique(df$plate_id)
          #if (length(plate_ids) != 1) return(NULL)  # Only allow one plate per event

          for(plate in plate_ids){
            print(plate)

          current_wells <- get_wells(plate)

          if (type == "click") {
            for (well in df$well) {
              if (well %in% current_wells) {
                current_wells <- setdiff(current_wells, well)
              } else {
                current_wells <- union(current_wells, well)
              }
            }
          }

          if (type == "selected") {
            current_wells <- unique(subset(df, plate_id == plate)$well)
          }

          set_wells(plate, current_wells)

          #if (source_id == "cluster_plot") {
          #  updateSelectizeInput(session, "plate_id3", selected = plate_id)
          #}
          update_all_select_inputs()
          }

        }

        # --- Observe plot events ---
        observeEvent(plotly::event_data("plotly_click", source = "plate_plot"), {
          handle_plotly_interaction("plotly_click", "click", source_id = "plate_plot")
        })

        observeEvent(plotly::event_data("plotly_selected", source = "plate_plot"), {
          handle_plotly_interaction("plotly_selected", "selected", source_id = "plate_plot")
        })


        # --- Reset selection (only for current plate) ---
        observeEvent(input$reset_well, {
          if (!is.null(input$plate_id) && nzchar(input$plate_id)) {
            selected_wells$data <- subset(selected_wells$data, plate_id != input$plate_id)
            update_all_select_inputs()
          }
        })

        # --- Abstract observer for plate/well syncing ---
        observe_plate_select_sync <- function(plate_input, well_input) {
          observeEvent(input[[plate_input]], {
            plate_id <- input[[plate_input]]
            if (!is.null(plate_id) && nzchar(plate_id)) {
              updateSelectizeInput(session, well_input, selected = get_wells(plate_id))
            }
          })
          observeEvent(input[[well_input]], {
            plate_id <- input[[plate_input]]
            if (!is.null(plate_id) && nzchar(plate_id)) {
              set_wells(plate_id, input[[well_input]])
            }
          })
        }

        observe_plate_select_sync("plate_id", "selected_well")
        observe_plate_select_sync("plate_id1", "selected_well1")
        observe_plate_select_sync("plate_id3", "clusteredwells")

      ############
      ### BEGIN Plot Sweeps

      output$sweep_view <- renderPlotly({
        req(input$facets)
        if(!is.null(SE())){

          if(!input$all_plates){
          se <- SE()[,SE()$Plate_ID == input$plate_id1]
          }else{se <- SE()}


          #if(input$assay_id1 %in% assayNames(se)){
            assayNames <- input$assay_id1
          #}
          print(input$color_group1)
          print(input$group_by_meta1)
          if(input$color_group1 == ""){
            color_group <- input$group_by_meta1
          }else{
            color_group <- input$color_group1
          }




        if(length(get_wells(input$plate_id1)) > 0){
            se <- se[,se$Well %in% get_wells(input$plate_id1)]

        }



         # p <-  plotAssayVSSweeps(se, assayList = assayNames,
         #                    rowCol = input$group_by_meta1, colorGroup = color_group,
         #                    wrapFormula = variable~.)

            rowCol <-input$group_by_meta1

            melted.se <- sechm::meltSE(se, features=row.names(se),
                                       assayName=assayNames,
                                       rowDat.columns = rowCol)

            melted.se <-reshape2::melt(melted.se,
                                       measure.vars = assayNames,
                                       value.name = "assays", variable.name="vars")


            # if(input$assay_option1 == "raw"){
            #   melted.dat[[assayName]] <- melted.dat[[assayName]]
            #   legend <- assayName
            # }
            # if(input$assay_option1 == "log10"){
            #   melted.dat[[assayName]] <- log10(abs(melted.dat[[assayName]]))
            #   legend <- paste("log10(", assayName, ")", sep="")
            # }
            # if(input$assay_option1 == "scale"){
            #   melted.dat[[assayName]] <- scale(melted.dat[[assayName]], center = T)
            #   legend <- paste("Z-scaled(", assayName, ")", sep="")
            # }

              p <- ggplot2::ggplot(melted.se, aes(x=.data[[rowCol]], y=assays, color=.data[[color_group]])) +
                ggplot2::stat_summary(geom='errorbar',fun.data=mean_se, size=1, alpha=0.6) +
                ggplot2::stat_summary(geom='line', fun = "mean", size=1, alpha=1) +
                ggplot2::theme_minimal(base_size = 16) +
                ggplot2::ylab(input$y.labels) +
                ggplot2::xlab(input$x.labels) +
                ggplot2::geom_hline(yintercept=0, linetype="dashed")

              print(input$facet_group1)
                if(input$facet_group1 == ""){

                      facet_group1 <- "."
                }else{
                  facet_group1 <- as.character(input$facet_group1)
    }
                print(melted.se)
              if(input$invertFacet){
                facetForm <- paste(facet_group1, "vars", sep=" ~ ")
              }else{
                facetForm <- paste("vars", facet_group1, sep=" ~ ")
              }

              print(as.formula(facetForm))

                if(input$facets == "grid"){
                    p <- p + facet_grid(as.formula(facetForm, env = parent.frame()), scales="free")
                }
              if(input$facets == "wrap"){
                if(input$invertFacet){
                  p <- p + facet_wrap(as.formula(facetForm, env = parent.frame()),
                                      scales="free", ncol = length(unique(melted.se[["vars"]])))
                }else{
                  p <- p + facet_wrap(as.formula(facetForm, env = parent.frame()),
                                      scales="free", nrow = length(unique(melted.se[["vars"]])))
                }
              }
          ggplotly(p)
        }
      })
      jqui_resizable(ui="#sweep_view")


    ############
    ### BEGIN Clustering





  observeEvent(input$clustering, {

    if(!is.null(input$clusterAssay) | !is.null(input$clusterColData)){

      withProgress(message = 'Perform Clustering...', value = 0.2, {
    req(SE())

    se <- SE()
    print(unique(se$Plate_ID))
    incProgress(0.5)
    se <- reducedDim.Cellwise(se, assayList = input$clusterAssay,
                              colNames = input$clusterColData,
                              k_clusters = input$k_cluster)
    incProgress(0.5, message="Updating Results")

    SEs[[input$object]] <- se
    initialized(FALSE)
    SEinit(SEs[[input$object]])
    incProgress(0.75, message="Finished")

      })
    }
  })

  generate_cluster_plot <- function(SE, reduction_name, color_var, xlab, ylab, selected_wells, plot_source) {

    df <- cbind(as.data.frame(reducedDim(SE, reduction_name)),
                as.data.frame(colData(SE)))
    colnames(df)[1:2] <- c("X1", "X2")

    df$color_value <- df[[color_var]]
    df$hover_text <- paste("Well: ", df$Well,
                           "<br>Plate ID: ", df$Plate_ID,
                           "<br>", color_var, ": ", df$color_value)

    df$key_combined <- paste(df$Well, df$Plate_ID, sep=", ")
    p <- ggplot(df, aes(x = X1, y = X2, text = hover_text, key = key_combined)) +
      geom_point(aes(color = color_value)) +
      labs(x = xlab, y = ylab, color = color_var) +
      theme_minimal() +
      theme(plot.margin = margin(10, 10, 10, 10))

    if(is.numeric(df$color_value)){
      p <- p + viridis::scale_colour_viridis()
    }


    # if(length(get_wells(unique(df$Plate_ID))) > 0){
    #
    # highlighted_df <- df %>%
    #   dplyr::filter(Plate_ID %in% selected_wells$data$plate_id) %>%  # only plates with selected wells
    #   dplyr::rowwise() %>%
    #   dplyr::filter(Well %in% get_wells(Plate_ID)) %>%
    #   dplyr::ungroup()
    # print(highlighted_df)

    selected <- selected_wells$data

    if (!is.null(selected) && nrow(selected) > 0 && !is.null(df) && nrow(df) > 0) {
    highlighted_df <- dplyr::semi_join(df, selected, by = c("Plate_ID" = "plate_id", "Well" = "well"))
    print(highlighted_df)

      if (nrow(highlighted_df) > 0) {
        p <- p +
          geom_point(data = highlighted_df, aes(x = X1, y = X2),
                     shape = 21, fill = NA, color = "red", size = 3, stroke = 0.5,
                     inherit.aes = FALSE)
      }

}

    # Add independent title and legend
    ggplotly(p, tooltip = "text", source = "plate_plot") %>%
      layout(title = list(text = color_var, x = 0.5, xanchor = "center"),
             showlegend = TRUE)%>%
      layout(dragmode = "lasso") %>%  # Makes box selection the default
      config(
        displayModeBar = TRUE,
        modeBarButtonsToAdd = list('select2d'),  # Ensures the select tool is in the toolbar
        displaylogo = FALSE
      )

  }

  output$cluster_tsne_ui <- renderUI({

    if (!is.null(SE())) {
      if (length(reducedDims(SE())) > 1) {
        req(input$clustercolor1)
        n_cols <- 2
        plots <- lapply(input$clustercolor1, function(var) {
          plotlyOutput(outputId = paste0("tsne_plot_", var), height = "400px")
        })

        # Arrange into a grid
        rows <- split(plots, ceiling(seq_along(plots) / n_cols))

        tagList(
          lapply(rows, function(row) {
            fluidRow(
              lapply(row, function(plot) {
                column(width = 6, plot)
              })
            )
          })
        )
      } else {
        tagList(
          tags$p("No Clustering available yet! Generate above first"),
          withSpinner(plotlyOutput(outputId = "placeholder", height = "400px"))
        )
      }
    } else {
      tagList(
        tags$p("No Clustering available yet! Generate above first"),
        withSpinner(plotlyOutput(outputId = "placeholder", height = "400px"))
      )
    }

  })


  observe({
    req(input$clustercolor1)

    lapply(input$clustercolor1, function(var) {
      local({
        v <- var
        output[[paste0("tsne_plot_", v)]] <- renderPlotly({
          generate_cluster_plot(
            SE = SE(),
            reduction_name = "TSNE",
            color_var = v,
            xlab = "TSNE1",
            ylab = "TSNE2",
            selected_wells = selected_wells,
            plot_source = paste0("tsne_", v)
          )
        })
      })
    })


  })

  output$cluster_umap_ui <- renderUI({

    if (!is.null(SE())) {
      if (length(reducedDims(SE())) > 1) {
        req(input$clustercolor2)
        n_cols <- 2
        plots <- lapply(input$clustercolor2, function(var) {
          plotlyOutput(outputId = paste0("umap_plot_", var), height = "400px")
        })

        # Arrange into a grid
        rows <- split(plots, ceiling(seq_along(plots) / n_cols))

        tagList(
          lapply(rows, function(row) {
            fluidRow(
              lapply(row, function(plot) {
                column(width = 6, plot)
              })
            )
          })
        )
      } else {
        tagList(
          tags$p("No Clustering available yet! Generate above first"),
          withSpinner(plotlyOutput(outputId = "placeholder", height = "400px"))
        )
      }
    } else {
      tagList(
        tags$p("No Clustering available yet! Generate above first"),
        withSpinner(plotlyOutput(outputId = "placeholder", height = "400px"))
      )
    }

  })

  observe({
    req(input$clustercolor2)

    lapply(input$clustercolor2, function(var) {
      local({
        v <- var
        output[[paste0("umap_plot_", v)]] <- renderPlotly({
          generate_cluster_plot(
            SE = SE(),
            reduction_name = "UMAP",
            color_var = v,
            xlab = "UMAP1",
            ylab = "UMAP2",
            selected_wells = selected_wells,
            plot_source = paste0("umap_", v)
          )
        })
      })
    })
  })

  output$cluster_pca_ui <- renderUI({

    if (!is.null(SE())) {
      if (length(reducedDims(SE())) > 1) {
        req(input$clustercolor2)
        n_cols <- 2
        plots <- lapply(input$clustercolor3, function(var) {
          plotlyOutput(outputId = paste0("pca_plot_", var), height = "400px")
        })

        # Arrange into a grid
        rows <- split(plots, ceiling(seq_along(plots) / n_cols))

        tagList(
          lapply(rows, function(row) {
            fluidRow(
              lapply(row, function(plot) {
                column(width = 6, plot)
              })
            )
          })
        )
      } else {
        tagList(
          tags$p("No Clustering available yet! Generate above first"),
          withSpinner(plotlyOutput(outputId = "placeholder", height = "400px"))
        )
      }
    } else {
      tagList(
        tags$p("No Clustering available yet! Generate above first"),
        withSpinner(plotlyOutput(outputId = "placeholder", height = "400px"))
      )
    }

  })

  observe({
    req(input$clustercolor3)

    lapply(input$clustercolor3, function(var) {
      local({
        v <- var
        output[[paste0("pca_plot_", v)]] <- renderPlotly({
          generate_cluster_plot(
            SE = SE(),
            reduction_name = "PCA",
            color_var = v,
            xlab = "PC1",
            ylab = "PC2",
            selected_wells = selected_wells,
            plot_source = paste0("pca_", v)
          )
        })
      })
    })
  })

    #### End of Clustering
    ############
    ### BEGIN Data Manipulation Tab

  observeEvent(input$selected_slider, { req(slider_initialized())
    req(SE(), input$assay_id, input$plate_id)

    assayName <- input$assay_id

    # Melt and subset the data
    melted.dat <- sechm::meltSE(
      SE()[, SE()$Plate_ID == input$plate_id],
      features = row.names(rowData(SE())),
      assayName = assayName,
      rowDat.columns = c(input$sweep_id, input$group_by_meta)
    )

    # Clean and filter
    slider_vals <- melted.dat[[assayName]]
    slider_range <- input$selected_slider

    # Filter wells based on slider range
    matching_wells <- unique(melted.dat$Well[slider_vals >= slider_range[1] & slider_vals <= slider_range[2]])

    # Update selected_well input with matching wells
    #updateSelectInput(session, "selected_well", selected = matching_wells)
    #updateSelectInput(session, "selected_well1", selected = matching_wells)
  })

  # ============================================================
  # CUSTOMIZE OBJECT GROUPS — Define Conditions (tab_coldata)
  # ============================================================

  # Helper: make a plate grid plotly from colData
  .make_plate_grid <- function(coldata_df, fill_col = "Well") {
    rows <- LETTERS[1:16]
    cols <- sprintf("%02d", 1:24)
    grid <- expand.grid(Row = rows, Col = cols, stringsAsFactors = FALSE)
    grid$Well <- paste0(grid$Row, grid$Col)
    grid$Col_num <- as.integer(grid$Col)
    if (!is.null(coldata_df) && fill_col %in% colnames(coldata_df)) {
      grid <- merge(grid, coldata_df[, unique(c("Well", fill_col))], by = "Well", all.x = TRUE)
      fill_vals <- grid[[fill_col]]
    } else {
      fill_vals <- NA_character_
    }
    grid$fill_val <- as.character(fill_vals)
    p <- ggplot(grid, aes(x = Col_num, y = Row, fill = fill_val,
                          text = paste0(Well, ": ", fill_val))) +
      geom_tile(color = "grey70") +
      scale_x_continuous(breaks = 1:24, expand = c(0, 0)) +
      scale_y_discrete(limits = rev, expand = c(0, 0)) +
      geom_text(aes(label = paste0(Row, sub("^0", "", Col))),
                size = 2.2, color = "white") +
      theme_minimal(base_size = 11) +
      labs(fill = fill_col, x = "Column", y = NULL)
    ggplotly(p, tooltip = "text")
  }

  # Reactive: rule list stored as a list of lists
  cond_rules <- reactiveVal(list())
  cond_rule_counter <- reactiveVal(0L)

  observeEvent(input$cond_add_rule, {
    n <- cond_rule_counter() + 1L
    cond_rule_counter(n)
    rules <- cond_rules()
    rules[[as.character(n)]] <- list(id = n, plates = NULL, col_min = 1L, col_max = 24L, label = "")
    cond_rules(rules)
  })

  observeEvent(input$cond_clear_rules, {
    cond_rules(list())
    cond_rule_counter(0L)
  })

  output$cond_rules_ui <- renderUI({
    se <- SE()
    rules <- cond_rules()
    if (length(rules) == 0) return(helpText("Click 'Add Rule' to define a mapping."))
    plate_ids <- if (!is.null(se)) unique(as.character(colData(se)$Plate_ID)) else character(0)
    tagList(lapply(seq_along(rules), function(i) {
      r <- rules[[i]]
      rid <- as.character(r$id)
      wellPanel(
        style = "padding: 8px; margin-bottom: 6px;",
        fluidRow(
          column(12,
            tags$strong(paste("Rule", i)),
            actionButton(paste0("cond_rm_rule_", rid), "×",
                         class = "btn-danger btn-xs pull-right",
                         style = "margin-top:-2px;")
          )
        ),
        fluidRow(
          column(6,
            selectizeInput(paste0("cond_rule_plates_", rid), "Plate_ID(s):",
                           choices = plate_ids, multiple = TRUE,
                           selected = r$plates,
                           options = list(placeholder = "All plates"))
          ),
          column(3,
            numericInput(paste0("cond_rule_col_min_", rid), "Col min:", value = r$col_min, min = 1L, max = 24L, step = 1L)
          ),
          column(3,
            numericInput(paste0("cond_rule_col_max_", rid), "Col max:", value = r$col_max, min = 1L, max = 24L, step = 1L)
          )
        ),
        fluidRow(
          column(12,
            textInput(paste0("cond_rule_label_", rid), "Label value:", value = r$label, placeholder = "e.g. WT")
          )
        )
      )
    }))
  })

  # Remove individual rules
  observe({
    rules <- cond_rules()
    lapply(names(rules), function(rid) {
      btn_id <- paste0("cond_rm_rule_", rid)
      observeEvent(input[[btn_id]], {
        cur <- cond_rules()
        cur[[rid]] <- NULL
        cond_rules(cur)
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # Plate grid for condition preview
  output$plate_view_col <- renderPlotly({
    se <- SE()
    req(se, input$plate_id4)
    cd <- as.data.frame(colData(se))
    cd <- cd[cd$Plate_ID %in% input$plate_id4, , drop = FALSE]
    fill_col <- if (!is.null(input$cond_preview_col) && input$cond_preview_col %in% colnames(cd))
      input$cond_preview_col else "Well"
    .make_plate_grid(cd, fill_col)
  })

  output$cond_coldata_preview <- renderDT({
    se <- SE()
    req(se)
    cd <- as.data.frame(colData(se))
    if (!is.null(input$plate_id4) && nchar(input$plate_id4[1]) > 0)
      cd <- cd[cd$Plate_ID %in% input$plate_id4, , drop = FALSE]
    datatable(cd, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  # Collect current rule values from inputs
  .collect_rules <- function(rules) {
    lapply(names(rules), function(rid) {
      list(
        plates  = input[[paste0("cond_rule_plates_", rid)]],
        col_min = input[[paste0("cond_rule_col_min_", rid)]],
        col_max = input[[paste0("cond_rule_col_max_", rid)]],
        label   = input[[paste0("cond_rule_label_", rid)]]
      )
    })
  }

  observeEvent(input$cond_apply, {
    se <- SE()
    req(se, input$cond_new_col, nchar(trimws(input$cond_new_col)) > 0)
    rules <- .collect_rules(cond_rules())
    if (length(rules) == 0) {
      showNotification("No rules defined.", type = "warning")
      return()
    }
    tryCatch({
      cd <- as.data.frame(colData(se))
      new_col <- trimws(input$cond_new_col)
      cd[[new_col]] <- NA_character_
      for (r in rules) {
        plates  <- if (is.null(r$plates) || length(r$plates) == 0) unique(cd$Plate_ID) else r$plates
        col_min <- as.integer(r$col_min %||% 1L)
        col_max <- as.integer(r$col_max %||% 24L)
        label   <- r$label
        well_col <- as.integer(sub("^[A-Za-z]+", "", cd$Well))
        idx <- cd$Plate_ID %in% plates & !is.na(well_col) & well_col >= col_min & well_col <= col_max
        cd[[new_col]][idx] <- label
      }
      colData(se)[[new_col]] <- cd[[new_col]]
      SEs[[input$object]] <- se
      SEinit(SEs[[input$object]])
      updateSelectInput(session, "cond_preview_col",
                        choices = colnames(as.data.frame(colData(se))),
                        selected = new_col)
      showNotification(paste0("Column '", new_col, "' added to colData."), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)), type = "error")
    })
  })


  # ============================================================
  # CUSTOMIZE OBJECT GROUPS — Define Sweeps (tab_rowdata)
  # ============================================================

  row_recode_counter <- reactiveVal(0L)

  observeEvent(input$row_col_select, {
    se <- SE()
    req(se, input$row_col_select)
    rd <- as.data.frame(rowData(se))
    vals <- unique(as.character(rd[[input$row_col_select]]))
    updateSelectInput(session, "lp_val", choices = vals, selected = vals[1])
  })

  observeEvent(input$row_add_recode, {
    n <- row_recode_counter() + 1L
    row_recode_counter(n)
  })

  output$row_recode_ui <- renderUI({
    n <- row_recode_counter()
    se <- SE()
    if (n == 0 || is.null(se)) return(helpText("Click 'Add mapping' to add old → new value pairs."))
    rd <- as.data.frame(rowData(se))
    col <- input$row_col_select
    old_vals <- if (!is.null(col) && col %in% colnames(rd))
      unique(as.character(rd[[col]])) else character(0)

    tagList(lapply(seq_len(n), function(i) {
      fluidRow(
        column(5, selectInput(paste0("recode_old_", i), if (i == 1) "Old value" else NULL,
                              choices = old_vals, width = "100%")),
        column(5, textInput(paste0("recode_new_", i), if (i == 1) "New value" else NULL,
                            placeholder = "New label", width = "100%")),
        column(2, if (i == 1) br() else NULL,
               actionButton(paste0("recode_rm_", i), "×", class = "btn-danger btn-xs"))
      )
    }))
  })

  output$row_data_table <- renderDT({
    se <- SE()
    req(se)
    rd <- as.data.frame(rowData(se))
    datatable(rd, options = list(pageLength = 10, scrollX = TRUE), rownames = TRUE)
  })

  observeEvent(input$row_apply_recode, {
    se <- SE()
    req(se, input$row_col_select)
    n <- row_recode_counter()
    if (n == 0) { showNotification("No mappings defined.", type = "warning"); return() }
    tryCatch({
      col <- input$row_col_select
      rd_col <- as.character(rowData(se)[[col]])
      for (i in seq_len(n)) {
        old_v <- input[[paste0("recode_old_", i)]]
        new_v <- input[[paste0("recode_new_", i)]]
        if (!is.null(old_v) && !is.null(new_v) && nchar(trimws(new_v)) > 0)
          rd_col[rd_col == old_v] <- trimws(new_v)
      }
      rowData(se)[[col]] <- rd_col
      SEs[[input$object]] <- se
      SEinit(SEs[[input$object]])
      showNotification(paste0("Column '", col, "' recoded."), type = "message")
      row_recode_counter(0L)
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)), type = "error")
    })
  })

  observeEvent(input$lp_apply, {
    se <- SE()
    req(se, input$lp_col, input$lp_val, input$lp_new_col)
    tryCatch({
      lp_vals <- rowData(se)[[input$lp_col]]
      rowData(se)[[trimws(input$lp_new_col)]] <- lp_vals == input$lp_val
      SEs[[input$object]] <- se
      SEinit(SEs[[input$object]])
      showNotification(paste0("Logical column '", input$lp_new_col, "' created."), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)), type = "error")
    })
  })


  # ============================================================
  # CUSTOMIZE OBJECT GROUPS — Change Assays (tab_assays)
  # ============================================================

  # -- Sweep selector UI for colAG --
  output$ag_sweep_ui <- renderUI({
    se <- SE()
    req(se)
    mode <- input$ag_sweep_mode
    if (mode == "all") return(NULL)
    if (mode == "range") {
      nr <- nrow(se)
      return(fluidRow(
        column(6, numericInput("ag_sweep_min", "First sweep:", value = 1L, min = 1L, max = nr, step = 1L)),
        column(6, numericInput("ag_sweep_max", "Last sweep:", value = nr, min = 1L, max = nr, step = 1L))
      ))
    }
    if (mode == "logical") {
      rd <- as.data.frame(rowData(se))
      log_cols <- colnames(rd)[purrr::map_lgl(rd, is.logical)]
      if (length(log_cols) == 0) return(helpText("No logical columns in rowData. Create one in 'Define Sweeps'."))
      return(selectInput("ag_sweep_logical_col", "Logical column (TRUE = use sweep):", choices = log_cols))
    }
  })

  ag_status_msg <- reactiveVal("")
  output$ag_status <- renderText(ag_status_msg())

  observeEvent(input$ag_run, {
    se <- SE()
    req(se, input$ag_assays)
    tryCatch({
      sweeps <- switch(input$ag_sweep_mode,
        "all"     = row.names(se),
        "range"   = {
          idx <- seq(input$ag_sweep_min, input$ag_sweep_max)
          row.names(se)[idx[idx >= 1 & idx <= nrow(se)]]
        },
        "logical" = {
          req(input$ag_sweep_logical_col)
          row.names(se)[isTRUE(rowData(se)[[input$ag_sweep_logical_col]])]
        }
      )
      se <- colAG(se, assayList = input$ag_assays, sweeps = sweeps)
      SEs[[input$object]] <- se
      SEinit(SEs[[input$object]])
      # Refresh transform col selector
      cd <- as.data.frame(colData(se))
      num_cols <- colnames(cd)[purrr::map_lgl(cd, is.numeric)]
      updateSelectInput(session, "transform_col", choices = num_cols, selected = num_cols[length(num_cols)])
      added <- paste(paste0(input$ag_assays, "_mean"), collapse = ", ")
      ag_status_msg(paste0("Done. Added: ", added))
    }, error = function(e) {
      ag_status_msg(paste("Error:", conditionMessage(e)))
    })
  })

  output$ag_result_preview <- renderDT({
    se <- SE()
    req(se)
    cd <- as.data.frame(colData(se))
    mean_cols <- grep("_mean$", colnames(cd), value = TRUE)
    if (length(mean_cols) == 0) return(datatable(data.frame(message = "No _mean columns yet.")))
    datatable(cd[, c("Well", "Plate_ID", mean_cols), drop = FALSE],
              options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  # -- Transform Columns --

  transform_status_msg <- reactiveVal("")
  output$transform_status <- renderText(transform_status_msg())

  output$transform_preview <- renderPlotly({
    se <- SE()
    req(se, input$transform_col, input$transform_fn)
    col <- input$transform_col
    cd <- as.data.frame(colData(se))
    req(col %in% colnames(cd))
    x_raw <- cd[[col]]
    x_new <- tryCatch({
      switch(input$transform_fn,
        "*1000"  = x_raw * 1000,
        "/1000"  = x_raw / 1000,
        "*1e9"   = x_raw * 1e9,
        "/1e9"   = x_raw / 1e9,
        "log1p"  = log1p(x_raw),
        "exp"    = exp(x_raw),
        "abs"    = abs(x_raw),
        "negate" = -x_raw,
        "custom" = {
          req(input$transform_custom_expr)
          .x <- x_raw
          eval(parse(text = input$transform_custom_expr))
        }
      )
    }, error = function(e) NA_real_)
    df_plot <- data.frame(
      Before = x_raw, After = x_new,
      Well = cd$Well, Plate_ID = cd$Plate_ID
    )
    plot_ly(df_plot, x = ~Before, y = ~After, type = "scatter", mode = "markers",
            text = ~paste0(Well, " (", Plate_ID, ")"),
            hoverinfo = "text+x+y", marker = list(size = 6, opacity = 0.7)) %>%
      layout(xaxis = list(title = paste("Before:", col)),
             yaxis = list(title = "After transform"))
  })

  observeEvent(input$transform_apply, {
    se <- SE()
    req(se, input$transform_col, input$transform_fn)
    col <- input$transform_col
    cd <- as.data.frame(colData(se))
    req(col %in% colnames(cd))
    tryCatch({
      x_raw <- cd[[col]]
      x_new <- switch(input$transform_fn,
        "*1000"  = x_raw * 1000,
        "/1000"  = x_raw / 1000,
        "*1e9"   = x_raw * 1e9,
        "/1e9"   = x_raw / 1e9,
        "log1p"  = log1p(x_raw),
        "exp"    = exp(x_raw),
        "abs"    = abs(x_raw),
        "negate" = -x_raw,
        "custom" = {
          req(input$transform_custom_expr)
          .x <- x_raw
          eval(parse(text = input$transform_custom_expr))
        }
      )
      dest_col <- trimws(input$transform_new_col %||% "")
      if (nchar(dest_col) == 0) dest_col <- col
      colData(se)[[dest_col]] <- x_new
      SEs[[input$object]] <- se
      SEinit(SEs[[input$object]])
      transform_status_msg(paste0("Saved to '", dest_col, "'."))
    }, error = function(e) {
      transform_status_msg(paste("Error:", conditionMessage(e)))
    })
  })

  # -- Dimensionality Reduction --
  dr_status_msg <- reactiveVal("")
  output$dr_status <- renderText(dr_status_msg())

  observeEvent(input$dr_run, {
    se <- SE()
    req(se)
    an_sel <- if (length(input$dr_assays) > 0) input$dr_assays else c()
    cn_sel <- if (length(input$dr_colnames) > 0) input$dr_colnames else c()
    if (length(an_sel) == 0 && length(cn_sel) == 0) {
      dr_status_msg("Select at least one assay or column."); return()
    }
    scaling <- if (input$dr_scaling == "none") "global" else input$dr_scaling
    tryCatch({
      dr_status_msg("Running… (may take a moment)")
      se_new <- reducedDim.Cellwise(
        se, assayList = an_sel, colNames = cn_sel,
        scaling = scaling, center = isTRUE(input$dr_center),
        k_clusters = as.integer(input$dr_k)
      )
      SEs[[input$object]] <- se_new
      SEinit(SEs[[input$object]])
      dr_status_msg("Done. PCA, tSNE, and UMAP added to reducedDims.")
      cd <- as.data.frame(colData(se_new))
      updateSelectInput(session, "dr_color_col", choices = colnames(cd))
    }, error = function(e) {
      dr_status_msg(paste("Error:", conditionMessage(e)))
    })
  })

  output$dr_preview <- renderPlotly({
    se <- SE()
    req(se)
    rds <- tryCatch(SingleCellExperiment::reducedDims(se), error = function(e) NULL)
    req(!is.null(rds), "PCA" %in% names(rds))
    pca <- as.data.frame(rds[["PCA"]])
    req(ncol(pca) >= 2)
    cd <- as.data.frame(colData(se))
    df <- cbind(pca[, 1:2], cd[, intersect(colnames(cd), c("Well","Plate_ID", input$dr_color_col)), drop=FALSE])
    colnames(df)[1:2] <- c("PC1", "PC2")
    color_col <- input$dr_color_col
    if (!is.null(color_col) && color_col %in% colnames(df)) {
      plot_ly(df, x = ~PC1, y = ~PC2, color = ~get(color_col),
              type = "scatter", mode = "markers",
              text = if ("Well" %in% colnames(df)) ~Well else NULL,
              marker = list(size = 7)) %>%
        layout(title = "PCA (PC1 vs PC2)", coloraxis = list(colorbar = list(title = color_col)))
    } else {
      plot_ly(df, x = ~PC1, y = ~PC2, type = "scatter", mode = "markers",
              marker = list(size = 7)) %>%
        layout(title = "PCA (PC1 vs PC2)")
    }
  })


  # ============================================================
  # CUSTOMIZE OBJECT GROUPS — Filter Wells (tab_filter_wells)
  # ============================================================

  fw_filters <- reactiveVal(list())

  output$fw_val_ui <- renderUI({
    se <- SE()
    req(se, input$fw_col)
    cd <- as.data.frame(colData(se))
    col <- input$fw_col
    req(col %in% colnames(cd))
    vals <- unique(as.character(cd[[col]]))
    selectizeInput("fw_val", "Value(s):", choices = vals, multiple = TRUE, selected = vals[1])
  })

  observeEvent(input$fw_add, {
    req(input$fw_col, input$fw_val)
    filters <- fw_filters()
    key <- paste0(input$fw_col, "_", length(filters) + 1L)
    filters[[key]] <- list(col = input$fw_col, vals = input$fw_val, mode = input$fw_mode)
    fw_filters(filters)
  })

  output$fw_active_filters_ui <- renderUI({
    filters <- fw_filters()
    if (length(filters) == 0) return(helpText("No active filters."))
    tagList(lapply(names(filters), function(k) {
      f <- filters[[k]]
      verb <- if (f$mode == "keep") "IN" else "NOT IN"
      wellPanel(style = "padding:6px; margin-bottom:4px;",
        fluidRow(
          column(10, tags$small(paste0(f$col, " ", verb, " {", paste(f$vals, collapse=", "), "}"))),
          column(2,  actionButton(paste0("fw_rm_", k), "×", class = "btn-danger btn-xs"))
        )
      )
    }))
  })

  observe({
    filters <- fw_filters()
    lapply(names(filters), function(k) {
      btn_id <- paste0("fw_rm_", k)
      observeEvent(input[[btn_id]], {
        cur <- fw_filters()
        cur[[k]] <- NULL
        fw_filters(cur)
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # Reactive: wells to keep after applying all filters
  fw_wells_keep <- reactive({
    se <- SE()
    req(se)
    cd <- as.data.frame(colData(se))
    keep <- rep(TRUE, nrow(cd))
    filters <- fw_filters()
    for (f in filters) {
      col_vals <- as.character(cd[[f$col]])
      if (f$mode == "keep")
        keep <- keep & (col_vals %in% f$vals)
      else
        keep <- keep & !(col_vals %in% f$vals)
    }
    keep
  })

  output$fw_preview_n <- renderText({
    se <- SE()
    req(se)
    keep <- fw_wells_keep()
    paste0("Will keep ", sum(keep), " / ", length(keep), " wells.")
  })

  output$fw_plate_preview <- renderPlotly({
    se <- SE()
    req(se)
    cd <- as.data.frame(colData(se))
    keep <- fw_wells_keep()
    cd$Keep <- ifelse(keep, "Keep", "Remove")
    plates <- unique(cd$Plate_ID)
    cd <- cd[cd$Plate_ID == plates[1], , drop = FALSE]
    .make_plate_grid(cd, "Keep")
  })

  output$fw_coldata_preview <- renderDT({
    se <- SE()
    req(se)
    cd <- as.data.frame(colData(se))
    keep <- fw_wells_keep()
    cd$Keep <- ifelse(keep, "Keep", "Remove")
    datatable(cd[, c("Well","Plate_ID","Keep"), drop=FALSE],
              options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  fw_status_msg <- reactiveVal("")
  output$fw_status <- renderText(fw_status_msg())

  observeEvent(input$fw_apply, {
    se <- SE()
    req(se)
    keep <- fw_wells_keep()
    if (sum(keep) == 0) { showNotification("No wells would remain!", type = "error"); return() }
    tryCatch({
      se_filtered <- se[, keep]
      SEs[[input$object]] <- se_filtered
      SEinit(SEs[[input$object]])
      fw_filters(list())
      showNotification(paste0("SE filtered to ", sum(keep), " wells."), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)), type = "error")
    })
  })


  # Toggle visibility
  observeEvent(input$toggle_plate, {
    toggle("mini_plate_plot")
  })


  ## Download-Tab

  output$downloadRDS <- downloadHandler(
    filename = function() {
      paste(input$rdsObject, ".rds", sep="")
    },
    content = function(file) {
      saveRDS(SEs[[input$rdsObject]], file)
      # Data has been exported — safe to close without warning
      session$sendCustomMessage("setWarnBeforeClose", FALSE)
    }
  )



    observeEvent(input$quickStart, showModal(.getHelp("general")))
    observeEvent(input$help_SE, showModal(.getHelp("SE")))
    observeEvent(input$help_gassay, showModal(.getHelp("assay")))
    observeEvent(input$help_ggroup, showModal(.getHelp("group")))
    observeEvent(input$help_ggrid, showModal(.getHelp("grid")))
    observeEvent(input$help_gfreeaxes, showModal(.getHelp("grid")))
    observeEvent(input$help_hmassay, showModal(.getHelp("assay")))
    observeEvent(input$help_hmscale, showModal(.getHelp("scale")))
    observeEvent(input$help_hmtrim, showModal(.getHelp("scaletrim")))
    observeEvent(input$help_feature.lists, showModal(.getHelp("feature.lists")))

    # ---------------------------------------------------------------
    # Image Plate Viewer tab
    # ---------------------------------------------------------------

    .rows384 <- LETTERS[1:16]

    .extract_well <- function(bn) {
      m <- regexpr("(?<=_)([A-P])([0-9]{1,2})(?:-[0-9]+)?(?=_)", bn, perl = TRUE)
      if (m[1] == -1) return(NA_character_)
      tok <- regmatches(bn, m)
      row <- sub("^([A-P]).*$", "\\1", tok)
      col <- as.integer(sub("^[A-P]([0-9]{1,2}).*$", "\\1", tok))
      if (is.na(col) || col < 1 || col > 24) return(NA_character_)
      paste0(row, sprintf("%02d", col))
    }

    # Parse well, channel, and class from filenames like:
    #   ..._PLATEID_A01-1_BF_class1_crop.jpg
    .parse_img_meta <- function(bn) {
      m <- regexpr("_([A-P][0-9]{1,2}-[0-9]+)_([^_]+)_([^_]+)_crop\\.jpe?g$",
                   bn, perl = TRUE, ignore.case = TRUE)
      if (m[1] == -1)
        return(data.frame(well = NA_character_, channel = NA_character_,
                          img_class = NA_character_, stringsAsFactors = FALSE))
      tok   <- sub("^_", "", regmatches(bn, m))
      parts <- strsplit(tok, "_")[[1]]  # well-site, channel, class, "crop.jpg(e?)"
      well_site <- parts[1]
      channel   <- if (length(parts) >= 2) parts[2] else NA_character_
      img_class <- if (length(parts) >= 3) parts[3] else NA_character_
      raw <- sub("-.*$", "", well_site)
      row <- sub("^([A-P]).*$", "\\1", raw)
      col <- as.integer(sub("^[A-P]([0-9]{1,2})$", "\\1", raw))
      well <- if (!is.na(col) && col >= 1 && col <= 24)
        paste0(row, sprintf("%02d", col)) else NA_character_
      data.frame(well = well, channel = channel, img_class = img_class,
                 stringsAsFactors = FALSE)
    }

    # ---- client-side folder module (works locally + deployed) ----
    imgbrowser_data <- localImgBrowserServer("imgbrowser")

    # TRUE when the client-side module has supplied images
    img_client_mode <- reactive({
      d <- tryCatch(imgbrowser_data(), error = function(e) NULL)
      !is.null(d) && isTRUE(d$n > 0)
    })

    # ---- server-side folder picker (local / annotation tab) --------
    .img_roots <- c(
      Home = normalizePath("~", winslash = "/", mustWork = TRUE),
      C = "C:/", D = "D:/", Y = "Y:/", Z = "Z:/"
    )
    .img_roots <- .img_roots[c(TRUE, dir.exists(.img_roots[-1]))]

    shinyFiles::shinyDirChoose(input, "img_dir_btn", roots = .img_roots,
                               allowDirCreate = FALSE)

    img_parent <- reactiveVal(NULL)

    observeEvent(input$img_dir_btn, {
      p <- shinyFiles::parseDirPath(.img_roots, input$img_dir_btn)
      if (length(p) != 1 || !nzchar(p) || !dir.exists(p)) return()
      p <- normalizePath(p, winslash = "/")
      addResourcePath("imgplate", p)
      img_parent(p)
      se_name <- input$object
      if (!is.null(se_name) && nzchar(se_name) && !is.null(SEs[[se_name]]) &&
          !is.character(SEs[[se_name]])) {
        S4Vectors::metadata(SEs[[se_name]])$image_path_jpgs <- p
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    output$image_path_jpgs_display <- renderText({
      p <- img_parent()
      if (is.null(p)) "(no folder selected)" else p
    })

    # Match SE Plate_IDs to subdirectories (plate_id substring in dir name)
    img_plate_dirs <- reactive({
      parent <- img_parent()
      se     <- SE()
      req(!is.null(parent), !is.null(se))

      plate_ids  <- unique(colData(se)$Plate_ID)
      subdirs    <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
      subnames   <- basename(subdirs)

      matched <- list()
      for (pid in plate_ids) {
        idx <- which(grepl(pid, subnames, fixed = TRUE))
        if (length(idx) > 0)
          matched[[pid]] <- subdirs[idx[1]]
        else
          # fallback: use parent folder itself if only one plate
          if (length(plate_ids) == 1) matched[[pid]] <- parent
      }
      matched
    })

    observeEvent(img_plate_dirs(), {
      if (img_client_mode()) return()   # client mode drives this instead
      dirs <- img_plate_dirs()
      req(length(dirs) > 0)
      updateSelectInput(session, "plate_img_id",
                        choices  = names(dirs),
                        selected = names(dirs)[1])
    })

    # Client mode: populate plate / channel / class from module metadata
    observeEvent(imgbrowser_data(), {
      d <- imgbrowser_data()
      req(!is.null(d), d$n > 0)

      all_plates <- if (length(d$subdirs) > 0 && nzchar(d$subdirs[1]))
        d$subdirs else "(all)"

      # Filter to only plates that exist in the active SE
      se <- tryCatch(SE(), error = function(e) NULL)
      se_plates <- if (!is.null(se)) unique(colData(se)$Plate_ID) else character(0)

      plates <- if (length(se_plates) > 0 && !identical(all_plates, "(all)")) {
        # Keep subdirs whose name contains a SE Plate_ID as a substring
        matched <- all_plates[vapply(all_plates, function(p)
          any(grepl(p, se_plates, fixed = TRUE) | grepl(se_plates, p, fixed = TRUE)),
          logical(1))]
        if (length(matched) > 0) matched else all_plates   # fallback: show all
      } else {
        all_plates
      }

      updateSelectInput(session, "plate_img_id",
                        choices = plates, selected = plates[1])
      updateSelectInput(session, "img_channel",
                        choices = d$channels,
                        selected = if (length(d$channels)) d$channels[1])
      updateSelectInput(session, "img_class_sel",
                        choices = d$classes,
                        selected = if (length(d$classes)) d$classes[1])
    })

    # JPGs for selected plate folder
    img_files <- reactive({
      dirs <- img_plate_dirs()
      pid  <- input$plate_img_id
      req(length(dirs) > 0, nzchar(pid %||% ""), pid %in% names(dirs))
      list.files(dirs[[pid]],
                 pattern     = "\\.(jpg|jpeg)$",
                 full.names  = TRUE,
                 recursive   = TRUE,
                 ignore.case = TRUE)
    })

    # well/channel/class data.frame for selected plate
    img_well_map <- reactive({
      files <- img_files()
      if (length(files) == 0)
        return(data.frame(well = character(0), channel = character(0),
                          img_class = character(0), file = character(0),
                          stringsAsFactors = FALSE))
      meta <- do.call(rbind, lapply(seq_along(files), function(i) {
        r <- .parse_img_meta(basename(files[i]))
        r$file <- files[i]
        r
      }))
      meta[!is.na(meta$well), , drop = FALSE]
    })

    # Update channel / class dropdowns when server-mode map changes
    observeEvent(img_well_map(), {
      if (img_client_mode()) return()   # handled by imgbrowser_data observer
      df <- img_well_map()
      channels <- if (nrow(df) > 0) sort(unique(df$channel)) else character(0)
      classes  <- if (nrow(df) > 0) sort(unique(df$img_class)) else character(0)
      updateSelectInput(session, "img_channel",   choices = channels,
                        selected = if (length(channels) > 0) channels[1] else NULL)
      updateSelectInput(session, "img_class_sel", choices = classes,
                        selected = if (length(classes) > 0) classes[1] else NULL)
    }, ignoreNULL = TRUE)

    # Filtered map: named vector well -> file for current channel + class
    # (server mode only; client mode uses filterImgMap via img_well_urls)
    img_well_map_filtered <- reactive({
      df      <- img_well_map()
      channel <- input$img_channel
      cls     <- input$img_class_sel
      if (nrow(df) == 0) {
        ids <- as.vector(t(outer(.rows384, sprintf("%02d", 1:24), paste0)))
        return(setNames(rep(NA_character_, length(ids)), ids))
      }
      sub_df <- df
      if (!is.null(channel) && nzchar(channel %||% ""))
        sub_df <- sub_df[sub_df$channel == channel, , drop = FALSE]
      if (!is.null(cls) && nzchar(cls %||% ""))
        sub_df <- sub_df[sub_df$img_class == cls, , drop = FALSE]
      ids <- as.vector(t(outer(.rows384, sprintf("%02d", 1:24), paste0)))
      m   <- setNames(rep(NA_character_, length(ids)), ids)
      # Keep first file per well
      sub_df <- sub_df[!duplicated(sub_df$well), , drop = FALSE]
      m[sub_df$well] <- sub_df$file
      m
    })

    # Unified well -> URL map (blob: URLs in client mode, served paths in server mode)
    img_well_urls <- reactive({
      if (img_client_mode()) {
        d      <- imgbrowser_data()
        pid    <- input$plate_img_id %||% ""
        subdir <- if (nzchar(pid) && pid != "(all)") pid else ""
        filterImgMap(d$map,
                     subdir    = subdir,
                     channel   = input$img_channel   %||% NULL,
                     img_class = input$img_class_sel %||% NULL)
      } else {
        m      <- img_well_map_filtered()
        parent <- img_parent()
        if (is.null(parent)) return(setNames(rep(NA_character_, 384L),
                                             as.vector(t(outer(.rows384, sprintf("%02d", 1:24), paste0)))))
        if (!endsWith(parent, "/")) parent <- paste0(parent, "/")
        vapply(names(m), function(w) {
          fp <- m[[w]]
          if (is.null(fp) || is.na(fp)) return(NA_character_)
          rel <- substring(normalizePath(fp, winslash = "/"), nchar(parent) + 1)
          file.path("imgplate", rel)
        }, character(1))
      }
    })

    # Raw (unfiltered) colData values — used for palette + filter UI
    img_coldata_raw <- reactive({
      se  <- SE()
      var <- input$img_coldata_var
      pid <- input$plate_img_id
      if (is.null(se) || is.null(var) || var == "None" || is.null(pid)) return(NULL)
      cd <- as.data.frame(colData(se))
      cd <- cd[cd$Plate_ID == pid, , drop = FALSE]
      if (nrow(cd) == 0 || !var %in% colnames(cd)) return(NULL)
      setNames(cd[[var]], cd$Well)
    })

    # Filtered values — NAs mark excluded wells
    img_coldata_vals <- reactive({
      vals <- img_coldata_raw()
      if (is.null(vals)) return(NULL)
      if (is.numeric(vals)) {
        rng <- input$img_filter_range
        if (!is.null(rng))
          vals[vals < rng[1] | vals > rng[2]] <- NA_real_
      } else {
        groups <- input$img_filter_groups
        if (!is.null(groups))
          vals[!as.character(vals) %in% groups] <- NA
      }
      vals
    })

    # Filter UI: range slider (numeric) or checkbox group (categorical)
    output$img_filter_ui <- renderUI({
      vals <- img_coldata_raw()
      if (is.null(vals)) return(NULL)
      if (is.numeric(vals)) {
        rng <- range(vals, na.rm = TRUE)
        tagList(
          hr(),
          strong("Filter range"),
          sliderInput("img_filter_range", label = NULL,
                      min = rng[1], max = rng[2],
                      value = c(rng[1], rng[2]),
                      step  = signif(diff(rng) / 100, 2))
        )
      } else {
        lvls <- sort(unique(na.omit(as.character(vals))))
        tagList(
          hr(),
          strong("Filter groups"),
          checkboxGroupInput("img_filter_groups", label = NULL,
                             choices  = lvls,
                             selected = lvls)
        )
      }
    })

    # Legend UI: gradient bar (numeric) or colored squares (categorical)
    output$img_legend_ui <- renderUI({
      vals <- img_coldata_raw()
      if (is.null(vals)) return(NULL)
      if (is.numeric(vals)) {
        rng  <- range(vals, na.rm = TRUE)
        pal  <- viridis::viridis(12)
        grad <- paste0("linear-gradient(to right, ",
                       paste(pal, collapse = ", "), ")")
        tagList(
          hr(),
          strong("Legend"),
          div(style = paste0("background:", grad,
                             "; height:16px; border-radius:4px; margin:6px 0 2px;")),
          div(style = "display:flex; justify-content:space-between; font-size:11px;",
              span(signif(rng[1], 3)), span(signif(rng[2], 3)))
        )
      } else {
        lvls <- sort(unique(na.omit(as.character(vals))))
        pal  <- setNames(rainbow(length(lvls), s = 0.8, v = 0.9), lvls)
        tagList(
          hr(),
          strong("Legend"),
          div(style = "margin-top:6px;",
              lapply(lvls, function(lv) {
                div(style = "display:flex; align-items:center; margin:3px 0;",
                    div(style = paste0("width:14px; height:14px; border-radius:3px;",
                                       " background:", pal[lv],
                                       "; margin-right:7px; flex-shrink:0;")),
                    span(style = "font-size:12px;", lv))
              })
          )
        )
      }
    })

    # Map values → semi-transparent hex colors (names preserved)
    .vals_to_colors <- function(vals, alpha) {
      nms <- names(vals)
      if (is.numeric(vals)) {
        rng <- range(vals, na.rm = TRUE)
        pal <- viridis::viridis(256)
        idx <- if (diff(rng) == 0) rep(128L, length(vals)) else
          round((vals - rng[1]) / diff(rng) * 255) + 1L
        idx[is.na(vals)] <- NA_integer_
        cols <- pal[idx]
        cols[is.na(idx)] <- NA_character_
      } else {
        chars <- as.character(vals)
        lvls  <- unique(na.omit(chars))
        pal   <- setNames(rainbow(length(lvls), s = 0.8, v = 0.9), lvls)
        cols  <- pal[chars]
        cols[is.na(chars) | chars == "NA"] <- NA_character_
      }
      setNames(adjustcolor(cols, alpha.f = alpha), nms)
    }

    output$img_plate_stats <- renderText({
      pid  <- input$plate_img_id
      urls <- img_well_urls()
      n_shown <- sum(!is.na(urls) & nzchar(urls))
      if (img_client_mode()) {
        d     <- imgbrowser_data()
        chans <- if (length(d$channels) > 0) paste(d$channels, collapse = ", ") else "-"
        clses <- if (length(d$classes)  > 0) paste(d$classes,  collapse = ", ") else "-"
        paste0(
          "Plate:       ", pid %||% "-", "\n",
          "JPGs found:  ", d$n,          "\n",
          "Channels:    ", chans,         "\n",
          "Classes:     ", clses,         "\n",
          "Wells shown: ", n_shown, " / 384"
        )
      } else {
        dirs <- img_plate_dirs()
        df   <- img_well_map()
        folder   <- if (!is.null(pid) && pid %in% names(dirs)) basename(dirs[[pid]]) else "-"
        channels <- if (nrow(df) > 0) paste(sort(unique(df$channel)),   collapse = ", ") else "-"
        classes  <- if (nrow(df) > 0) paste(sort(unique(df$img_class)), collapse = ", ") else "-"
        paste0(
          "Folder:      ", folder,              "\n",
          "JPGs found:  ", length(img_files()), "\n",
          "Channels:    ", channels,            "\n",
          "Classes:     ", classes,             "\n",
          "Wells shown: ", n_shown, " / 384"
        )
      }
    })

    output$img_plate_css <- renderUI({
      b <- input$img_brightness %||% 1
      c <- input$img_contrast   %||% 1
      tags$style(HTML(sprintf(
        ".imgplate-well img { filter: brightness(%s) contrast(%s); }", b, c
      )))
    })

    output$img_plate_ui <- renderUI({
      urls   <- img_well_urls()
      vals   <- img_coldata_vals()
      colors <- if (!is.null(vals))
        .vals_to_colors(vals, alpha = input$img_overlay_alpha %||% 0.4) else NULL
      plateGridUI(
        url_map        = urls,
        coldata_colors = colors,
        click_input_id = session$ns("img_well_click"),
        hover_input_id = session$ns("img_well_hover")
      )
    })

    # Hover: show well + color variable + any extra selected colData columns
    output$img_hover_info <- renderText({
      well <- input$img_well_hover
      if (is.null(well) || !nzchar(well)) return("")
      se         <- SE()
      pid        <- input$plate_img_id
      var        <- input$img_coldata_var
      extra_vars <- input$img_hover_vars
      if (!is.null(se) && !is.null(pid)) {
        cd  <- as.data.frame(colData(se))
        row <- cd[cd$Well == well & cd$Plate_ID == pid, , drop = FALSE]
        if (nrow(row) > 0) {
          col_var <- if (!is.null(var) && var != "None" && var %in% colnames(row)) var else "QC"
          txt <- paste0(well, "  ", col_var, ": ", row[[col_var]][1])
          if (length(extra_vars) > 0) {
            extra <- vapply(extra_vars, function(v) {
              if (v %in% colnames(row)) paste0(v, ": ", row[[v]][1]) else ""
            }, character(1))
            extra <- extra[nzchar(extra)]
            if (length(extra) > 0)
              txt <- paste0(txt, "  |  ", paste(extra, collapse = "  "))
          }
          return(txt)
        }
      }
      well
    })

    # Click: enlarge both BF + Fluoro in modal with colData info
    observeEvent(input$img_well_click, {
      well <- input$img_well_click
      if (is.null(well) || !nzchar(well)) return()

      cls <- input$img_class_sel

      if (img_client_mode()) {
        d      <- imgbrowser_data()
        pid    <- input$plate_img_id %||% ""
        subdir <- if (nzchar(pid) && pid != "(all)") pid else ""
        bf_url  <- filterImgMap(d$map, subdir = subdir,
                                channel = "BF", img_class = cls)[well]
        other_ch <- setdiff(d$channels, "BF")
        flu_url  <- if (length(other_ch) > 0)
          filterImgMap(d$map, subdir = subdir,
                       channel = other_ch[1], img_class = cls)[well]
        else NA_character_
        bf_url  <- if (is.null(bf_url)  || is.na(bf_url))  NA_character_ else bf_url
        flu_url <- if (is.null(flu_url) || is.na(flu_url)) NA_character_ else flu_url
      } else {
        parent <- img_parent()
        if (is.null(parent)) return()
        if (!endsWith(parent, "/")) parent <- paste0(parent, "/")
        fp_to_url <- function(fp) {
          if (is.null(fp) || is.na(fp)) return(NA_character_)
          rel <- substring(normalizePath(fp, winslash = "/"), nchar(parent) + 1)
          file.path("imgplate", rel)
        }
        df  <- img_well_map()
        sub <- df[df$well == well &
                    (if (!is.null(cls) && nzchar(cls %||% "")) df$img_class == cls else TRUE),
                  , drop = FALSE]
        bf_row  <- sub[tolower(sub$channel) == "bf", , drop = FALSE]
        flu_row <- sub[tolower(sub$channel) != "bf", , drop = FALSE]
        bf_url  <- if (nrow(bf_row)  > 0) fp_to_url(bf_row$file[1])  else NA_character_
        flu_url <- if (nrow(flu_row) > 0) fp_to_url(flu_row$file[1]) else NA_character_
      }

      if (is.na(bf_url) && is.na(flu_url)) return()

      # Hover info string
      hover_txt  <- well
      se         <- SE()
      pid        <- input$plate_img_id
      var        <- input$img_coldata_var
      extra_vars <- input$img_hover_vars
      if (!is.null(se) && !is.null(pid)) {
        cd  <- as.data.frame(colData(se))
        row <- cd[cd$Well == well & cd$Plate_ID == pid, , drop = FALSE]
        if (nrow(row) > 0) {
          col_var <- if (!is.null(var) && var != "None" && var %in% colnames(row)) var else "QC"
          hover_txt <- paste0(well, "  \u2022  ", col_var, ": ", row[[col_var]][1])
          if (length(extra_vars) > 0) {
            extra <- vapply(extra_vars, function(v)
              if (v %in% colnames(row)) paste0(v, ": ", row[[v]][1]) else "",
              character(1))
            extra <- extra[nzchar(extra)]
            if (length(extra) > 0)
              hover_txt <- paste0(hover_txt, "  |  ", paste(extra, collapse = "  "))
          }
        }
      }

      init_b <- input$img_brightness %||% 1
      init_c <- input$img_contrast   %||% 1
      init_filter <- sprintf("max-width:100%%; max-height:55vh; border-radius:6px; filter:brightness(%s) contrast(%s);",
                             init_b, init_c)

      make_panel <- function(url, label, img_id) {
        if (!is.na(url))
          div(style = "text-align: center;",
              tags$img(src = url, id = img_id, style = init_filter),
              tags$p(style = "font-size: 12px; color: #666; margin-top: 4px;", label))
        else
          div(style = paste0("display:flex; align-items:center; justify-content:center;",
                             " min-height: 180px; background: #2a2a2a;",
                             " border-radius: 6px; color: #888;"),
              tags$p(paste(label, "\u2014 no image")))
      }

      showModal(modalDialog(
        title = paste("Well", well, "\u2014", pid),
        tags$div(style = paste0("font-family: monospace; font-size: 13px;",
                                " margin-bottom: 12px; color: #333;"), hover_txt),
        div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 12px;",
            make_panel(bf_url,  "BF",     "modal_bf_img"),
            make_panel(flu_url, "Fluoro", "modal_flu_img")),
        fluidRow(
          style = "margin-top: 14px; padding: 0 4px;",
          column(6,
            div(style = "display:flex; align-items:center; gap:8px;",
              tags$strong("BF", style = "font-size:12px;"),
              actionButton("modal_bf_reset", icon("rotate-left"), class = "btn-xs btn-default")),
            fluidRow(
              column(6, sliderInput("modal_bf_b", "Brightness", 0.2, 4, init_b, 0.1, width = "100%")),
              column(6, sliderInput("modal_bf_c", "Contrast",   0.2, 4, init_c, 0.1, width = "100%"))
            )
          ),
          column(6,
            div(style = "display:flex; align-items:center; gap:8px;",
              tags$strong("Fluoro", style = "font-size:12px;"),
              actionButton("modal_flu_reset", icon("rotate-left"), class = "btn-xs btn-default")),
            fluidRow(
              column(6, sliderInput("modal_flu_b", "Brightness", 0.2, 4, init_b, 0.1, width = "100%")),
              column(6, sliderInput("modal_flu_c", "Contrast",   0.2, 4, init_c, 0.1, width = "100%"))
            )
          )
        ),
        easyClose = TRUE, size = "l", footer = modalButton("Close")
      ))
    })

    # Modal image brightness / contrast — update CSS filter live via JS
    observe({
      b <- input$modal_bf_b %||% 1
      c <- input$modal_bf_c %||% 1
      shinyjs::runjs(sprintf(
        "var el=document.getElementById('modal_bf_img'); if(el) el.style.filter='brightness(%s) contrast(%s)';",
        b, c))
    })
    observe({
      b <- input$modal_flu_b %||% 1
      c <- input$modal_flu_c %||% 1
      shinyjs::runjs(sprintf(
        "var el=document.getElementById('modal_flu_img'); if(el) el.style.filter='brightness(%s) contrast(%s)';",
        b, c))
    })

    # Reset buttons
    observeEvent(input$img_reset_filter, {
      updateSliderInput(session, "img_brightness", value = 1)
      updateSliderInput(session, "img_contrast",   value = 1)
    })
    observeEvent(input$ann_bf_reset, {
      updateSliderInput(session, "ann_bf_b", value = 1)
      updateSliderInput(session, "ann_bf_c", value = 1)
    })
    observeEvent(input$ann_flu_reset, {
      updateSliderInput(session, "ann_flu_b", value = 1)
      updateSliderInput(session, "ann_flu_c", value = 1)
    })
    observeEvent(input$modal_bf_reset, {
      updateSliderInput(session, "modal_bf_b", value = 1)
      updateSliderInput(session, "modal_bf_c", value = 1)
    })
    observeEvent(input$modal_flu_reset, {
      updateSliderInput(session, "modal_flu_b", value = 1)
      updateSliderInput(session, "modal_flu_c", value = 1)
    })

    # ---------------------------------------------------------------
    # Manual Annotation tab
    # (reuses img_parent() from Image Plate Viewer)
    # ---------------------------------------------------------------

    # All JPGs in the selected folder (across all plates)
    all_img_files <- reactive({
      parent <- img_parent()
      if (is.null(parent)) return(character(0))
      list.files(parent,
                 pattern     = "\\.(jpg|jpeg)$",
                 full.names  = TRUE,
                 recursive   = TRUE,
                 ignore.case = TRUE)
    })

    ann_csv_path <- reactive({
      parent <- img_parent()
      if (is.null(parent)) return(NULL)
      file.path(parent, "manual_classifications.csv")
    })

    ann_rv <- reactiveValues(
      classes  = c("Induced", "Not induced", "Ambiguous"),
      shuffled = character(0),   # wells (not file paths)
      idx      = 0L,
      results  = data.frame(
        timestamp = character(0), well = character(0),
        plate_id  = character(0), img_class = character(0),
        label = character(0), person = character(0),
        stringsAsFactors = FALSE
      )
    )

    # Reset results when client-side folder changes (no CSV to load)
    observeEvent(imgbrowser_data(), {
      if (!img_client_mode()) return()
      ann_rv$results <- ann_rv$results[0, ]
    }, ignoreInit = TRUE)

    # Load existing CSV when server-side folder changes
    observeEvent(img_parent(), {
      cp <- ann_csv_path()
      if (!is.null(cp) && file.exists(cp)) {
        existing <- tryCatch(read.csv(cp, stringsAsFactors = FALSE),
                             error = function(e) NULL)
        if (!is.null(existing) && nrow(existing) > 0) {
          for (col in c("timestamp", "well", "plate_id", "img_class", "label", "person"))
            if (!col %in% colnames(existing)) existing[[col]] <- NA_character_
          ann_rv$results <- existing[, c("timestamp", "well", "plate_id", "img_class", "label", "person")]
        }
      } else {
        ann_rv$results <- ann_rv$results[0, ]
      }
    })

    # Plate dirs available for annotation (mirrors img_plate_dirs but with fallback)
    ann_plate_dirs <- reactive({
      if (img_client_mode()) {
        d <- imgbrowser_data()
        subs <- d$subdirs
        if (length(subs) == 0 || !nzchar(subs[1]))
          return(setNames(list(""), "(all)"))
        setNames(as.list(subs), subs)
      } else {
        tryCatch(img_plate_dirs(), error = function(e) {
          parent <- img_parent()
          if (is.null(parent)) return(list())
          subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
          if (length(subdirs) == 0) return(setNames(list(parent), basename(parent)))
          setNames(as.list(subdirs), basename(subdirs))
        })
      }
    })

    # Helper: which plate_id does a file belong to?
    .ann_file_plate <- function(fp, plate_dirs) {
      fp_norm <- normalizePath(fp, winslash = "/", mustWork = FALSE)
      for (pid in names(plate_dirs)) {
        d <- normalizePath(plate_dirs[[pid]], winslash = "/", mustWork = FALSE)
        if (!endsWith(d, "/")) d <- paste0(d, "/")
        if (startsWith(fp_norm, d)) return(pid)
      }
      NA_character_
    }

    # Populate ann_plate_ids when plate dirs change — all selected by default
    observeEvent(ann_plate_dirs(), {
      pids <- names(ann_plate_dirs())
      updateSelectizeInput(session, "ann_plate_ids",
                           choices = pids, selected = pids)
    }, ignoreNULL = TRUE)

    # Parse all images into a meta data.frame: well / channel / img_class / file / plate_id
    # "file" column holds file paths (server mode) or blob: URLs (client mode)
    ann_img_meta_df <- reactive({
      if (img_client_mode()) {
        d    <- imgbrowser_data()
        keys <- names(d$map)
        if (length(keys) == 0)
          return(data.frame(well = character(0), channel = character(0),
                            img_class = character(0), file = character(0),
                            plate_id  = character(0), stringsAsFactors = FALSE))
        split1  <- strsplit(keys, "|", fixed = TRUE)
        subdirs <- vapply(split1, function(x) if (length(x) >= 1) x[1] else "", character(1))
        wcc     <- vapply(split1, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
        split2  <- strsplit(wcc, ".", fixed = TRUE)
        wells   <- vapply(split2, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
        chans   <- vapply(split2, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
        clses   <- vapply(split2, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))
        df <- data.frame(well = wells, channel = chans, img_class = clses,
                         file = unname(d$map), plate_id = subdirs,
                         stringsAsFactors = FALSE)
        df[!is.na(df$well), , drop = FALSE]
      } else {
        files <- all_img_files()
        if (length(files) == 0)
          return(data.frame(well = character(0), channel = character(0),
                            img_class = character(0), file = character(0),
                            plate_id  = character(0), stringsAsFactors = FALSE))
        meta <- do.call(rbind, lapply(seq_along(files), function(i) {
          r <- .parse_img_meta(basename(files[i]))
          r$file     <- files[i]
          r$plate_id <- NA_character_
          r
        }))
        meta[!is.na(meta$well), , drop = FALSE]
      }
    })

    # Populate ann_img_class dropdown from available img_class values
    observeEvent(ann_img_meta_df(), {
      df  <- ann_img_meta_df()
      cls <- if (nrow(df) > 0) sort(unique(df$img_class)) else character(0)
      updateSelectInput(session, "ann_img_class", choices = cls,
                        selected = if (length(cls) > 0) cls[1] else NULL)
    }, ignoreNULL = TRUE)

    # Helper: composite key "well||plate_id" for unique identification across plates
    .ann_key <- function(well, plate_id) paste0(well, "||", ifelse(is.na(plate_id), "", plate_id))

    # Rebuild queue (stores composite "well||plate_id" keys)
    observeEvent(list(ann_img_meta_df(), input$ann_img_class, input$ann_plate_ids), {
      df   <- ann_img_meta_df()
      cls  <- input$ann_img_class
      pids <- input$ann_plate_ids
      if (nrow(df) == 0 || is.null(cls) || !nzchar(cls %||% "")) {
        ann_rv$shuffled <- character(0); ann_rv$idx <- 0L; return()
      }
      sub <- df[df$img_class == cls, , drop = FALSE]

      # Attach plate_id in server mode (computed from file paths)
      if (!img_client_mode() && nrow(sub) > 0) {
        dirs         <- ann_plate_dirs()
        sub$plate_id <- vapply(sub$file, .ann_file_plate, character(1), plate_dirs = dirs)
      }

      # Filter by selected plates; skip when pids are all empty (flat folder = no structure)
      real_pids <- pids[nzchar(pids %||% "")]
      if (length(real_pids) > 0)
        sub <- sub[!is.na(sub$plate_id) & sub$plate_id %in% real_pids, , drop = FALSE]

      sub       <- sub[!is.na(sub$well), , drop = FALSE]
      available <- unique(.ann_key(sub$well, sub$plate_id))
      available <- available[nzchar(available)]
      if (length(available) == 0) { ann_rv$shuffled <- character(0); ann_rv$idx <- 0L; return() }

      done_keys <- .ann_key(ann_rv$results$well, ann_rv$results$plate_id)[
        !is.na(ann_rv$results$img_class) & ann_rv$results$img_class == cls]
      pool <- setdiff(available, done_keys)
      if (length(pool) == 0) pool <- available
      ann_rv$shuffled <- sample(pool)
      ann_rv$idx      <- 1L
    }, ignoreNULL = FALSE)

    # Current item: composite key "well||plate_id"
    ann_current <- reactive({
      if (length(ann_rv$shuffled) == 0 || ann_rv$idx < 1) return(NA_character_)
      ann_rv$shuffled[[ann_rv$idx]]
    })

    # Split composite key into well + plate_id
    .ann_split_key <- function(key) {
      parts <- strsplit(key, "||", fixed = TRUE)[[1]]
      list(well = parts[1],
           pid  = if (length(parts) >= 2 && nzchar(parts[2])) parts[2] else NA_character_)
    }

    # BF + fluoro file paths (or blob: URLs) for the current well+plate + selected class
    ann_pair <- reactive({
      key <- ann_current()
      cls <- input$ann_img_class
      if (is.na(key) || is.null(cls) || !nzchar(cls %||% ""))
        return(list(bf = NA_character_, fluoro = NA_character_))
      kp      <- .ann_split_key(key)
      df      <- ann_img_meta_df()
      sub     <- df[df$well == kp$well & df$img_class == cls, , drop = FALSE]
      if (!is.na(kp$pid))
        sub <- sub[!is.na(sub$plate_id) & sub$plate_id == kp$pid, , drop = FALSE]
      bf_row  <- sub[tolower(sub$channel) == "bf", , drop = FALSE]
      flu_row <- sub[tolower(sub$channel) != "bf", , drop = FALSE]
      list(
        bf     = if (nrow(bf_row)  > 0) bf_row$file[1]  else NA_character_,
        fluoro = if (nrow(flu_row) > 0) flu_row$file[1] else NA_character_
      )
    })

    ann_advance <- function(remove_current = TRUE) {
      n <- length(ann_rv$shuffled)
      if (n == 0) return()
      if (remove_current) {
        # Labeled: remove current well from queue (it's done)
        idx <- ann_rv$idx
        ann_rv$shuffled <- ann_rv$shuffled[-idx]
        n2 <- length(ann_rv$shuffled)
        ann_rv$idx <- if (n2 == 0L) 0L else if (idx > n2) 1L else idx
      } else {
        # Skipped: keep well in queue, move to next
        ann_rv$idx <- ann_rv$idx %% n + 1L
      }
    }

    # Helper: try to infer plate ID from the file path using SE plate IDs
    .extract_plate_from_path <- function(fp, se) {
      if (is.null(se)) return(NULL)
      pids <- unique(colData(se)$Plate_ID)
      matched <- pids[sapply(pids, function(p) grepl(p, fp, fixed = TRUE))]
      if (length(matched) > 0) matched[1] else NULL
    }

    # Serve image via addResourcePath("imgplate") set in img_parent reactive
    .ann_file_to_url <- function(fp) {
      parent <- img_parent()
      if (is.null(parent)) return(NA_character_)
      if (!endsWith(parent, "/")) parent <- paste0(parent, "/")
      fp <- normalizePath(fp, winslash = "/", mustWork = FALSE)
      if (!startsWith(fp, parent)) return(NA_character_)
      file.path("imgplate", substring(fp, nchar(parent) + 1))
    }

    observe({
      b <- input$ann_bf_b %||% 1
      c <- input$ann_bf_c %||% 1
      shinyjs::runjs(sprintf(
        "var el=document.getElementById('ann_bf_img'); if(el) el.style.filter='brightness(%s) contrast(%s)';",
        b, c))
    })
    observe({
      b <- input$ann_flu_b %||% 1
      c <- input$ann_flu_c %||% 1
      shinyjs::runjs(sprintf(
        "var el=document.getElementById('ann_flu_img'); if(el) el.style.filter='brightness(%s) contrast(%s)';",
        b, c))
    })

    output$ann_image_ui <- renderUI({
      if (length(ann_rv$shuffled) == 0) {
        msg <- if (nrow(ann_img_meta_df()) > 0)
          "All wells in the queue have been annotated."
        else
          "Select an image folder in the Image Plate Viewer tab first."
        return(tags$p(style = "color: #aaa; padding: 20px;", msg))
      }
      key <- ann_current()
      if (is.na(key)) return(tags$p("No well available."))
      kp   <- .ann_split_key(key)
      well <- kp$well

      pair <- ann_pair()

      b_bf  <- input$ann_bf_b  %||% 1
      c_bf  <- input$ann_bf_c  %||% 1
      b_flu <- input$ann_flu_b %||% 1
      c_flu <- input$ann_flu_c %||% 1

      make_panel <- function(fp, label, img_id, b, c) {
        ok  <- !is.null(fp) && !is.na(fp) && nzchar(fp)
        url <- if (!ok) NA_character_
                else if (img_client_mode()) fp
                else .ann_file_to_url(fp)
        if (ok && !is.na(url)) {
          div(style = "text-align: center;",
              tags$img(src = url, id = img_id,
                       style = sprintf("max-width:100%%; max-height:52vh; border-radius:6px; filter:brightness(%s) contrast(%s);", b, c)),
              tags$p(style = "font-size: 12px; color: #666; margin: 4px 0 0;", label))
        } else {
          div(style = paste0("display:flex; align-items:center; justify-content:center;",
                             " min-height:200px; background:#2a2a2a; border-radius:6px; color:#888;"),
              tags$p(paste(label, "\u2014 no image")))
        }
      }

      pid_lbl <- if (!is.na(kp$pid)) paste0(" \u2014 ", kp$pid) else ""
      tagList(
        tags$p(style = "font-size:13px; font-family:monospace; color:#555; margin:4px 8px;",
               paste0(well, pid_lbl)),
        div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 12px; padding: 8px;",
            make_panel(pair$bf,     "BF",     "ann_bf_img", b_bf,  c_bf),
            make_panel(pair$fluoro, "Fluoro", "ann_flu_img", b_flu, c_flu))
      )
    })

    output$ann_class_buttons <- renderUI({
      req(length(ann_rv$classes) > 0)
      div(style = "display: flex; flex-wrap: wrap; gap: 8px;",
          lapply(ann_rv$classes, function(lab) {
            tags$button(lab, class = "btn btn-success",
                        style = "padding:10px 20px; border-radius:8px; font-size:14px;",
                        onclick = sprintf(
                          "Shiny.setInputValue('ann_class_selected', %s, {priority:'event'});",
                          jsonlite::toJSON(lab, auto_unbox = TRUE)))
          })
      )
    })

    observeEvent(input$ann_class_selected, {
      key  <- ann_current(); req(!is.na(key))
      kp   <- .ann_split_key(key)
      well <- kp$well
      pid  <- kp$pid %||% NA_character_
      cls    <- input$ann_img_class %||% NA_character_
      person <- trimws(input$ann_person)

      row <- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        well      = well,
        plate_id  = pid,
        img_class = cls,
        label     = input$ann_class_selected,
        person    = person,
        stringsAsFactors = FALSE
      )
      ann_rv$results <- rbind(ann_rv$results, row)
      if (!img_client_mode()) {
        cp <- ann_csv_path()
        if (!is.null(cp))
          write.table(row, cp, sep = ",", row.names = FALSE,
                      col.names = !file.exists(cp), append = file.exists(cp))
      }
      ann_advance()
    }, ignoreInit = TRUE)

    observeEvent(input$ann_add_class, {
      x <- trimws(input$ann_new_class)
      if (nzchar(x) && !(x %in% ann_rv$classes)) {
        ann_rv$classes <- c(ann_rv$classes, x)
        updateTextInput(session, "ann_new_class", value = "")
      }
    })

    observeEvent(input$ann_clear_classes, { ann_rv$classes <- character(0) })
    observeEvent(input$ann_skip,          { ann_advance(remove_current = FALSE) })

    observeEvent(input$ann_undo, {
      if (nrow(ann_rv$results) == 0) return()
      last       <- ann_rv$results[nrow(ann_rv$results), ]
      last_key   <- .ann_key(last$well, last$plate_id)
      ann_rv$results   <- ann_rv$results[-nrow(ann_rv$results), , drop = FALSE]
      ann_rv$shuffled  <- c(last_key, ann_rv$shuffled)
      ann_rv$idx       <- 1L
      if (!img_client_mode()) {
        cp <- ann_csv_path()
        if (!is.null(cp) && file.exists(cp)) {
          lines <- readLines(cp)
          if (length(lines) > 1) writeLines(lines[-length(lines)], cp)
        }
      }
    })

    observeEvent(input$ann_clear_results, {
      ann_rv$results <- ann_rv$results[0, ]
      if (!img_client_mode()) {
        cp <- ann_csv_path()
        if (!is.null(cp) && file.exists(cp)) file.remove(cp)
      }
    })

    # Helper: read and normalise a CSV into the standard results columns
    .read_ann_csv <- function(path) {
      df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) return(NULL)
      req_cols <- c("timestamp", "well", "plate_id", "img_class", "label", "person")
      for (col in req_cols)
        if (!col %in% colnames(df)) df[[col]] <- NA_character_
      df[, req_cols]
    }

    # Upload CSV → load as current results (resumes annotation)
    observeEvent(input$ann_upload_csv, {
      req(input$ann_upload_csv)
      df <- .read_ann_csv(input$ann_upload_csv$datapath)
      if (!is.null(df)) ann_rv$results <- df
    })

    # Merge & deduplicate: keep the most recent label per (well, plate_id, img_class)
    observeEvent(input$ann_merge_results, {
      df <- ann_rv$results
      if (nrow(df) == 0) return()
      df <- df[order(df$timestamp, decreasing = TRUE), , drop = FALSE]
      key <- paste0(df$well, "||",
                    ifelse(is.na(df$plate_id), "", df$plate_id), "||",
                    ifelse(is.na(df$img_class), "", df$img_class))
      ann_rv$results <- df[!duplicated(key), , drop = FALSE]
    })

    output$ann_progress <- renderText({
      cls       <- input$ann_img_class %||% ""
      df        <- ann_img_meta_df()
      total     <- if (nrow(df) > 0 && nzchar(cls)) {
        sub_cls <- df[!is.na(df$img_class) & df$img_class == cls, , drop = FALSE]
        length(unique(.ann_key(sub_cls$well, sub_cls$plate_id)))
      } else 0L
      labeled   <- sum(!is.na(ann_rv$results$img_class) & ann_rv$results$img_class == cls)
      remaining <- length(ann_rv$shuffled) - ann_rv$idx + 1L
      paste0("Labeled: ", labeled, " / ", total,
             if (remaining > 0) paste0("  |  Queue: ", remaining) else "")
    })

    output$ann_results_preview <- renderTable({
      res  <- ann_rv$results
      cols <- intersect(c("timestamp", "well", "plate_id", "img_class", "label", "person"), colnames(res))
      tail(res[, cols, drop = FALSE], 8)
    })

    output$ann_download_csv <- downloadHandler(
      filename = function() paste0("annotations_", format(Sys.Date(), "%Y%m%d"), ".csv"),
      content  = function(file) write.csv(ann_rv$results, file, row.names = FALSE)
    )

    # ---------------------------------------------------------------
    # Auto Classification pipeline  (tab_img_classify)
    # ---------------------------------------------------------------

    cls_rv <- reactiveValues(
      particles  = NULL,
      filtered   = NULL,
      aggregated = NULL,
      scored     = NULL,
      classified = NULL
    )

    output$cls_load_status  <- renderText({ req(cls_rv$particles);  paste0("Loaded: ", nrow(cls_rv$particles), " rows, ", length(unique(cls_rv$particles$Channel_Name)), " channel(s), ", length(unique(cls_rv$particles$Plate_ID)), " plate(s).") })
    output$cls_filter_status <- renderText({ req(cls_rv$filtered);  n_na <- sum(is.na(cls_rv$filtered$Mean)); paste0("After filter: ", nrow(cls_rv$filtered)," rows, ", n_na, " particles blanked (", round(100*n_na/nrow(cls_rv$filtered)),"%).") })
    output$cls_agg_status    <- renderText({ req(cls_rv$aggregated); paste0("Aggregated: ", nrow(cls_rv$aggregated), " well-channel rows.") })
    output$cls_score_status  <- renderText({ req(cls_rv$scored);     paste0("Scored: ", nrow(cls_rv$scored), " rows. Score range: [", round(min(cls_rv$scored$channel_score, na.rm=TRUE),3), ", ", round(max(cls_rv$scored$channel_score, na.rm=TRUE),3), "].") })
    output$cls_merge_status  <- renderText({ req(input$cls_merge_se); "" })

    # Resolve scale group input → actual column names
    cls_scale_group_vars <- reactive({
      switch(input$cls_scale_groups %||% "channel_plate",
        channel_plate = c("Channel_Name", "Plate_ID"),
        channel       = "Channel_Name",
        channel_well  = c("Channel_Name", "Well")
      )
    })

    # Reactive: scaled preview (no filtering — used for diagnostics only)
    # Always uses scale(center=X, scale=TRUE) to match filterParticles logic.
    cls_scaled_preview <- reactive({
      req(cls_rv$particles)
      df         <- cls_rv$particles
      scale_mode <- input$cls_scale_mode %||% "uncentered"
      gvars      <- cls_scale_group_vars()
      center     <- scale_mode == "centered"
      req("Mean" %in% names(df), "Channel_Name" %in% names(df))
      df %>%
        dplyr::group_by(dplyr::across(dplyr::any_of(gvars))) %>%
        dplyr::mutate(
          Mean_scaled = as.numeric(scale(.data$Mean, center = center, scale = TRUE)),
          Area_scaled = as.numeric(scale(.data$Area, center = center, scale = TRUE))
        ) %>%
        dplyr::ungroup()
    })

    # ── Diagnostic plots ──────────────────────────────────────────────────────

    # Filter method selector — pre-selected based on scale mode
    output$cls_filter_method_ui <- renderUI({
      scale_mode     <- input$cls_scale_mode %||% "uncentered"
      default_method <- if (scale_mode == "centered") "zscore" else "median_ratio"
      cur_method     <- input$cls_filter_method %||% default_method
      mismatch <- (scale_mode == "centered"   && cur_method == "median_ratio") ||
                  (scale_mode == "uncentered" && cur_method == "zscore")
      tagList(
        selectInput("cls_filter_method", "Method:",
          choices  = c("Uncentered SD \u00f7 median/3 (recommended)" = "median_ratio",
                       "Centered z-score (threshold = 0)"             = "zscore"),
          selected = default_method),
        if (mismatch)
          tags$p(style="color:#e67e22; font-size:10px;",
                 "\u26a0 Method and scaling mode mismatch — results may be suboptimal.")
        else NULL
      )
    })

    # Threshold slider — range adapts to the current filter method
    output$cls_filter_threshold_ui <- renderUI({
      method <- input$cls_filter_method %||% "median_ratio"
      if (method == "zscore") {
        tagList(
          sliderInput("cls_filter_threshold", "Threshold (z-score):",
            min=-3, max=3, value=0, step=0.05),
          helpText(style="font-size:10px;",
                   "Particles with z-score below this value are rejected. Default = 0 (below mean)."))
      } else {
        tagList(
          sliderInput("cls_filter_threshold", "Threshold (SD units, NA = auto):",
            min=0, max=3, value=0.33, step=0.01),
          helpText(style="font-size:10px;",
                   "Default auto = median(scaled)/3 per group. Set a fixed value to override."),
          checkboxInput("cls_threshold_auto", "Use auto threshold (median/3 per group)", value=TRUE))
      }
    })

    # Step 1 — ECDF: metric + raw/scaled toggle, with optional faceting
    output$cls_diag_load_ecdf <- renderPlotly({
      metric    <- input$cls_load_metric %||% "Mean"
      view_mode <- input$cls_load_scale  %||% "raw"
      if (view_mode == "scaled") {
        req(cls_scaled_preview())
        df       <- cls_scaled_preview()
        col_name <- paste0(metric, "_scaled")
        req(col_name %in% names(df))
        x_label  <- paste0(metric, " (scaled \u00f7 SD)")
      } else {
        req(cls_rv$particles)
        df       <- cls_rv$particles
        col_name <- metric
        req(col_name %in% names(df))
        x_label  <- paste0(metric, " (raw)")
      }
      df <- df[!is.na(df[[col_name]]), , drop=FALSE]
      df$`.x` <- df[[col_name]]
      meta      <- input$cls_meta_col %||% ""
      use_facet <- nzchar(meta) && meta %in% names(df)
      df$`.facet` <- if (use_facet) df[[meta]] else "All"
      thr <- if (view_mode == "scaled" && metric == "Mean") {
        method <- input$cls_filter_method %||% "median_ratio"
        auto   <- isTRUE(input$cls_threshold_auto)
        thr_v  <- input$cls_filter_threshold %||% NA_real_
        if (!is.na(thr_v) && !auto) thr_v else if (method == "zscore") 0 else NULL
      } else NULL
      p <- ggplot(df, aes(x=.x, color=Channel_Name)) +
        stat_ecdf(geom="step", linewidth=0.8) +
        labs(title=paste0("ECDF — ", metric, if (view_mode=="scaled") " (scaled)" else " (raw)"),
             x=x_label, y="Cumulative fraction", color="Channel") +
        theme_minimal(base_size=11)
      if (!is.null(thr))
        p <- p + geom_vline(xintercept=thr, linetype="dashed", color="red", linewidth=0.8)
      if (use_facet) p <- p + facet_wrap(~.facet, labeller=label_both)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.2))
    })

    # Step 1 — Raw: Histogram of particles per well
    output$cls_diag_load_hist <- renderPlotly({
      req(cls_rv$particles)
      df <- cls_rv$particles
      req("Well" %in% names(df), "Channel_Name" %in% names(df))
      cnt <- as.data.frame(table(Well=df$Well, Channel_Name=df$Channel_Name))
      cnt <- cnt[cnt$Freq > 0, ]
      p <- ggplot(cnt, aes(x=Freq, fill=Channel_Name)) +
        geom_histogram(bins=30, alpha=0.7, position="identity") +
        labs(title="Particles per well", x="# particles", y="# wells", fill="Channel") +
        theme_minimal(base_size=11)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.2))
    })

    # Step 1 — Scaled preview ECDF with threshold line
    output$cls_diag_scaled_preview <- renderPlotly({
      req(cls_scaled_preview())
      df  <- cls_scaled_preview()
      req("Mean_scaled" %in% names(df), "Channel_Name" %in% names(df))
      method   <- input$cls_filter_method %||% "median_ratio"
      auto_thr <- isTRUE(input$cls_threshold_auto)
      thr_v    <- input$cls_filter_threshold %||% NA_real_
      # Show fixed threshold line; skip if auto (per-group) since no single line applies
      thr <- if (!auto_thr && !is.na(thr_v)) thr_v else if (method == "zscore") 0 else NULL
      meta      <- input$cls_meta_col %||% ""
      use_facet <- nzchar(meta) && meta %in% names(df)
      df$`.facet` <- if (use_facet) df[[meta]] else "All"
      df <- df[!is.na(df$Mean_scaled), , drop=FALSE]
      scale_lbl <- if ((input$cls_scale_mode %||% "uncentered") == "centered") "z-score" else "\u00f7 SD"
      p <- ggplot(df, aes(x=Mean_scaled, color=Channel_Name)) +
        stat_ecdf(geom="step", linewidth=0.9) +
        labs(title=paste0("Scaled Mean (", scale_lbl, ") — filter preview"),
             x=paste0("Mean_scaled (", scale_lbl, ")"),
             y="Cumulative fraction", color="Channel") +
        theme_minimal(base_size=11)
      if (!is.null(thr))
        p <- p + geom_vline(xintercept=thr, linetype="dashed", color="red", linewidth=0.9) +
             annotate("text", x=thr, y=0.08, hjust=-0.1, size=3, color="red",
                      label=paste0("threshold = ", round(thr, 3)))
      if (use_facet) p <- p + facet_wrap(~.facet, labeller=label_both)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.2))
    })

    # Step 2 — Threshold preview: scaled ECDF + live threshold line
    output$cls_diag_filter_threshold <- renderPlotly({
      req(cls_scaled_preview())
      df     <- cls_scaled_preview()
      req("Mean_scaled" %in% names(df), "Channel_Name" %in% names(df))
      method   <- input$cls_filter_method %||% "median_ratio"
      auto_thr <- isTRUE(input$cls_threshold_auto)
      thr_v    <- input$cls_filter_threshold %||% NA_real_
      # For auto (median/3 per group): show per-group median/3 lines in annotation
      show_fixed_line <- !auto_thr || method == "zscore"
      thr <- if (!is.na(thr_v) && !auto_thr) thr_v else if (method == "zscore") 0 else NA_real_
      meta      <- input$cls_meta_col %||% ""
      use_facet <- nzchar(meta) && meta %in% names(df)
      df$`.facet` <- if (use_facet) df[[meta]] else "All"
      df <- df[!is.na(df$Mean_scaled), , drop=FALSE]
      p <- ggplot(df, aes(x=Mean_scaled, color=Channel_Name)) +
        stat_ecdf(geom="step", linewidth=0.9) +
        labs(title="Threshold preview — scaled Mean",
             x="Mean_scaled", y="Cumulative fraction", color="Channel") +
        theme_minimal(base_size=11)
      if (show_fixed_line && !is.na(thr)) {
        pct_rej <- round(100 * mean(df$Mean_scaled < thr, na.rm=TRUE), 1)
        p <- p +
          geom_vline(xintercept=thr, linetype="dashed", color="red", linewidth=0.9) +
          annotate("text", x=thr, y=0.06, hjust=-0.1, size=3, color="red",
                   label=paste0("thr=", round(thr, 2), " (", pct_rej, "% rejected)"))
      } else if (auto_thr && method == "median_ratio") {
        # Show per-channel per-group median/3 as vertical segments
        thr_df <- df %>%
          dplyr::group_by(Channel_Name) %>%
          dplyr::summarise(auto_t = median(Mean_scaled, na.rm=TRUE)/3, .groups="drop")
        p <- p + geom_vline(data=thr_df, aes(xintercept=auto_t, color=Channel_Name),
                            linetype="dashed", linewidth=0.8) +
             geom_text(data=thr_df, aes(x=auto_t, label=paste0("median/3=", round(auto_t,2)),
                                         color=Channel_Name), y=0.06, hjust=-0.1, size=2.8)
      }
      if (use_facet) p <- p + facet_wrap(~.facet, labeller=label_both)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.2))
    })

    # Step 2 — Filter: ECDF before vs after (per channel, with optional faceting)
    output$cls_diag_filter_ecdf <- renderPlotly({
      req(cls_rv$particles, cls_rv$filtered)
      meta      <- input$cls_meta_col %||% ""
      keep_cols <- intersect(c("Mean", "Channel_Name", if (nzchar(meta)) meta else NULL),
                             names(cls_rv$particles))
      b <- cls_rv$particles[!is.na(cls_rv$particles$Mean), keep_cols, drop=FALSE]
      a <- cls_rv$filtered[!is.na(cls_rv$filtered$Mean),
                           intersect(keep_cols, names(cls_rv$filtered)), drop=FALSE]
      b$Stage <- "Before"; a$Stage <- "After"
      df <- rbind(b, a)
      use_facet <- nzchar(meta) && meta %in% names(df)
      df$`.facet` <- if (use_facet) df[[meta]] else "All"
      p <- ggplot(df, aes(x=Mean, color=Channel_Name, linetype=Stage)) +
        stat_ecdf(geom="step", linewidth=0.8) +
        labs(title="Raw Mean ECDF: before vs after filter",
             x="Mean intensity", y="Cumulative fraction",
             color="Channel", linetype="Stage") +
        theme_minimal(base_size=11)
      if (use_facet) p <- p + facet_wrap(~.facet, labeller=label_both)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.25))
    })

    # Step 2 — Filter: % particles retained per channel
    output$cls_diag_filter_retain <- renderPlotly({
      req(cls_rv$filtered)
      df <- cls_rv$filtered
      req("Channel_Name" %in% names(df), "Mean" %in% names(df))
      ret <- tapply(!is.na(df$Mean), df$Channel_Name, mean, na.rm=TRUE)
      ret_df <- data.frame(Channel=names(ret), pct=as.numeric(ret)*100, stringsAsFactors=FALSE)
      p <- ggplot(ret_df, aes(x=Channel, y=pct, fill=Channel)) +
        geom_col(alpha=0.85) +
        geom_hline(yintercept=50, linetype="dashed", color="grey50") +
        scale_y_continuous(limits=c(0,100)) +
        labs(title="% particles retained after filter", x="Channel", y="% retained") +
        theme_minimal(base_size=11) + theme(legend.position="none")
      ggplotly(p)
    })

    # Step 3 — Aggregate: boxplot + jitter, metric + raw/scaled toggle, with faceting
    output$cls_diag_agg_box <- renderPlotly({
      req(cls_rv$aggregated)
      df        <- .join_se_meta(cls_rv$aggregated, input$cls_seDataset)
      view_mode <- input$cls_agg_view    %||% "raw"
      metric    <- input$cls_agg_metric  %||% "Mean_agg"
      # If "scaled" selected, map raw metric names to their _z equivalents
      raw_to_z <- c(Mean_agg="Mean_z", Area_agg="Area_z", normArea="normArea_z")
      col_name <- if (view_mode == "scaled" && metric %in% names(raw_to_z)) {
        raw_to_z[[metric]]
      } else {
        metric
      }
      req(col_name %in% names(df), "Channel_Name" %in% names(df))
      df$`.m`     <- df[[col_name]]
      df$`.hover` <- paste0("Well: ",   df$Well     %||% "?",
                            "<br>Plate: ", df$Plate_ID %||% "?",
                            "<br>",  col_name, ": ", round(df[[col_name]], 4))
      meta      <- input$cls_meta_col %||% ""
      use_facet <- nzchar(meta) && meta %in% names(df)
      df$`.facet` <- if (use_facet) df[[meta]] else "All"
      add_hline   <- col_name == "normArea"
      p <- ggplot(df, aes(x=Channel_Name, y=.m, fill=Channel_Name)) +
        geom_boxplot(alpha=0.45, outlier.shape=NA, width=0.5) +
        geom_jitter(aes(text=.hover), width=0.18, size=1.8, alpha=0.65) +
        labs(title=paste(col_name, "per channel"), x="Channel", y=col_name) +
        theme_minimal(base_size=11) + theme(legend.position="none")
      if (add_hline)
        p <- p + geom_hline(yintercept=0.1, linetype="dashed", color="grey40")
      if (use_facet) p <- p + facet_wrap(~.facet, labeller=label_both)
      ggplotly(p, tooltip="text")
    })

    # Step 4 — Score: ECDF of channel_score per channel
    output$cls_diag_score_ecdf <- renderPlotly({
      req(cls_rv$scored)
      df <- cls_rv$scored
      req("channel_score" %in% names(df), "Channel_Name" %in% names(df))
      p <- ggplot(df, aes(x=channel_score, color=Channel_Name)) +
        stat_ecdf(geom="step", linewidth=0.9) +
        labs(title="ECDF — channel_score", x="channel_score", y="Cumulative fraction",
             color="Channel") +
        theme_minimal(base_size=11)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.2))
    })

    # Step 4 — Score: scatter of Mean_z vs normArea_z colored by channel_score
    output$cls_diag_score_scatter <- renderPlotly({
      req(cls_rv$scored)
      df <- cls_rv$scored
      req("Mean_z" %in% names(df), "normArea_z" %in% names(df), "channel_score" %in% names(df))
      p <- ggplot(df, aes(x=Mean_z, y=normArea_z, color=channel_score, shape=Channel_Name,
                          text=paste0("Well: ", Well,
                                      "<br>Channel: ", Channel_Name,
                                      "<br>Score: ",   round(channel_score, 3),
                                      "<br>Mean_z: ",  round(Mean_z, 3),
                                      "<br>normArea_z: ", round(normArea_z, 3)))) +
        geom_point(alpha=0.75, size=2.5) +
        scale_color_viridis_c(option="plasma", name="Score") +
        labs(title="Score components", x="Mean_z", y="normArea_z") +
        theme_minimal(base_size=11)
      ggplotly(p, tooltip="text")
    })

    # Step 4 — Score: ECDF of any selectable metric
    output$cls_diag_metric_ecdf <- renderPlotly({
      req(cls_rv$scored)
      df     <- cls_rv$scored
      metric <- input$cls_score_metric %||% "channel_score"
      req(metric %in% names(df), "Channel_Name" %in% names(df))
      df$`.m` <- df[[metric]]
      df <- df[!is.na(df$.m), , drop=FALSE]
      p <- ggplot(df, aes(x=.m, color=Channel_Name)) +
        stat_ecdf(geom="step", linewidth=0.9) +
        labs(title=paste("ECDF —", metric), x=metric, y="Cumulative fraction", color="Channel") +
        theme_minimal(base_size=11)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.2))
    })

    # ── Helpers for SE metadata joining ──────────────────────────────────────

    # Join SE colData onto any data frame keyed by Well + Plate_ID
    .join_se_meta <- function(df, se_name) {
      if (is.null(se_name) || !nzchar(se_name) || is.null(SEs[[se_name]])) return(df)
      cd <- tryCatch(as.data.frame(SummarizedExperiment::colData(SEs[[se_name]])),
                     error = function(e) NULL)
      if (is.null(cd)) return(df)
      join_by  <- intersect(c("Well", "Plate_ID"), intersect(names(cd), names(df)))
      new_cols <- setdiff(names(cd), names(df))
      if (length(join_by) == 0 || length(new_cols) == 0) return(df)
      dplyr::left_join(df, cd[, c(join_by, new_cols), drop=FALSE], by=join_by)
    }

    # Scored data with SE metadata columns appended
    cls_scored_meta <- reactive({
      .join_se_meta(req(cls_rv$scored), input$cls_seDataset)
    })

    # Classified data pivoted to long format (one row per Well × Channel), + SE metadata
    cls_analysis_long <- reactive({
      req(cls_rv$classified)
      df <- .join_se_meta(cls_rv$classified, input$cls_seDataset)
      score_cols <- setdiff(grep("_score$", names(df), value=TRUE), "max_score")
      if (length(score_cols) == 0) return(df)
      long <- tidyr::pivot_longer(df, cols=dplyr::all_of(score_cols),
                                  names_to="Channel", values_to="Score")
      long$Channel <- sub("_score$", "", long$Channel)
      long
    })

    # Metadata column selector — reads from particles (SE cols joined at load)
    # Plate_ID is always offered; default to it for sensible faceting.
    output$cls_meta_col_ui <- renderUI({
      df <- cls_rv$particles
      if (is.null(df)) {
        return(helpText(style="font-size:10px; color:#888;",
                        "Load particles first. SE metadata is joined at load time."))
      }
      # Per-particle identifiers that should never be used as facet variables
      never_facet <- c("Well", "Channel_Name", "Image_Type", "CorrSel",
                       "Selection", "Image_ID", "Mean", "Area", "IntDen",
                       "Number_of_Particles", "Selection_Area",
                       "Mean_scaled", "Area_scaled", "Row", "Column", "QC")
      char_cols <- names(df)[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
      char_cols <- setdiff(char_cols, never_facet)
      # Plate_ID is per-well and always a useful facet option — ensure it's included
      if ("Plate_ID" %in% names(df)) char_cols <- union("Plate_ID", char_cols)
      if (length(char_cols) == 0) {
        return(helpText(style="font-size:10px; color:#888;",
                        "No grouping columns found. Select an SE before loading to join metadata."))
      }
      default_col <- if ("Plate_ID" %in% char_cols) "Plate_ID" else char_cols[1]
      selectInput("cls_meta_col", "Facet / group by:",
                  choices = c("(none)" = "", char_cols), selected = default_col)
    })

    # Step 5 — ECDF of per-channel scores (from wide classified → long)
    output$cls_diag_score_channel <- renderPlotly({
      req(cls_analysis_long())
      df <- cls_analysis_long()
      req("Score" %in% names(df), "Channel" %in% names(df))
      p <- ggplot(df, aes(x=Score, color=Channel)) +
        stat_ecdf(geom="step", linewidth=0.9) +
        labs(title="Channel score ECDF", x="channel_score",
             y="Cumulative fraction", color="Channel") +
        theme_minimal(base_size=11)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.2))
    })

    # Step 5 — ECDF faceted by selected SE metadata column
    output$cls_diag_score_condition <- renderPlotly({
      req(cls_analysis_long())
      df      <- cls_analysis_long()
      meta    <- input$cls_meta_col
      req(!is.null(meta), nzchar(meta), meta %in% names(df))
      req("Score" %in% names(df), "Channel" %in% names(df))
      df$`.facet` <- df[[meta]]
      p <- ggplot(df, aes(x=Score, color=Channel)) +
        stat_ecdf(geom="step", linewidth=0.9) +
        facet_wrap(~.facet, labeller=label_both) +
        labs(title=paste("Score ECDF — faceted by", meta),
             x="channel_score", y="Cumulative fraction", color="Channel") +
        theme_minimal(base_size=11)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.12))
    })

    # Step 5 — ECDF faceted by Classification, colored by SE metadata (or Channel)
    output$cls_diag_score_cls <- renderPlotly({
      req(cls_analysis_long())
      df      <- cls_analysis_long()
      req("Score" %in% names(df), "Channel" %in% names(df),
          "Classification" %in% names(df))
      meta <- input$cls_meta_col
      use_meta <- !is.null(meta) && nzchar(meta) && meta %in% names(df)
      df$`.color` <- if (use_meta) df[[meta]] else df$Channel
      color_lbl   <- if (use_meta) meta else "Channel"
      p <- ggplot(df, aes(x=Score, color=.color)) +
        stat_ecdf(geom="step", linewidth=0.9) +
        facet_wrap(~Classification, labeller=label_both) +
        labs(title="Score ECDF — faceted by Classification",
             x="channel_score", y="Cumulative fraction", color=color_lbl) +
        theme_minimal(base_size=11)
      ggplotly(p) %>% plotly::layout(legend=list(orientation="h", y=-0.12))
    })

    # Step 5 — Classify: bar chart of classification counts
    output$cls_diag_classify_bar <- renderPlotly({
      req(cls_rv$classified)
      df <- as.data.frame(table(Classification=cls_rv$classified$Classification))
      df <- df[order(df$Freq, decreasing=TRUE), ]
      p <- ggplot(df, aes(x=reorder(Classification, -Freq), y=Freq, fill=Classification,
                          text=paste0(Classification, ": ", Freq, " wells"))) +
        geom_col(alpha=0.85) +
        labs(title="Classification counts", x="Classification", y="# wells") +
        theme_minimal(base_size=11) +
        theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))
      ggplotly(p, tooltip="text")
    })

    # Step 5 — 2D correlation controls (dynamic: populated from classified + SE columns)
    output$cls_corr_controls_2d <- renderUI({
      req(cls_rv$classified)
      df         <- cls_analysis_long()   # has SE meta if available
      score_cols <- grep("_score$|_normArea$", names(cls_rv$classified), value=TRUE)
      if (length(score_cols) < 1) return(helpText("Score columns not available."))
      # Color options: classification + plate + any char SE columns in the joined df
      base_color <- intersect(c("Classification", "Plate_ID"), names(df))
      se_cols    <- if (!is.null(input$cls_meta_col) && nzchar(input$cls_meta_col) &&
                        input$cls_meta_col %in% names(df))
                     input$cls_meta_col else character(0)
      color_opts <- unique(c(base_color, se_cols, score_cols))
      sel_y <- if (length(score_cols) >= 2) score_cols[2] else score_cols[1]
      fluidRow(
        column(4, selectInput("cls_corr_x",     "X axis:", choices=score_cols,
                              selected=score_cols[1])),
        column(4, selectInput("cls_corr_y",     "Y axis:", choices=score_cols,
                              selected=sel_y)),
        column(4, selectInput("cls_corr_color", "Color:",  choices=color_opts,
                              selected=color_opts[1]))
      )
    })

    # Step 5 — 2D score scatter (uses SE-joined classified data)
    output$cls_diag_corr_2d <- renderPlotly({
      req(cls_rv$classified, input$cls_corr_x, input$cls_corr_y)
      df      <- .join_se_meta(cls_rv$classified, input$cls_seDataset)
      x_col   <- input$cls_corr_x
      y_col   <- input$cls_corr_y
      col_col <- input$cls_corr_color %||% "Classification"
      req(x_col %in% names(df), y_col %in% names(df), col_col %in% names(df))
      df$`.x`     <- df[[x_col]]
      df$`.y`     <- df[[y_col]]
      df$`.color` <- df[[col_col]]
      df$`.hover` <- paste0("Well: ", df$Well,
                            "<br>Plate: ", df$Plate_ID,
                            "<br>Class: ", df$Classification)
      p <- ggplot(df, aes(x=.x, y=.y, color=.color, text=.hover)) +
        geom_point(alpha=0.75, size=2.5) +
        labs(title=paste(x_col, "vs", y_col), x=x_col, y=y_col, color=col_col) +
        theme_minimal(base_size=11)
      ggplotly(p, tooltip="text") %>%
        plotly::layout(legend=list(orientation="h", y=-0.25))
    })

    # 3D correlation controls (dynamic, uses SE-joined data for color options)
    output$cls_corr_controls_3d <- renderUI({
      req(cls_rv$classified)
      df       <- .join_se_meta(cls_rv$classified, input$cls_seDataset)
      num_cols <- grep("_score$|_normArea$|^max_score$", names(cls_rv$classified), value=TRUE)
      if (length(num_cols) < 1) return(helpText("Score columns not available."))
      se_cols    <- if (!is.null(input$cls_meta_col) && nzchar(input$cls_meta_col) &&
                        input$cls_meta_col %in% names(df))
                     input$cls_meta_col else character(0)
      color_opts <- unique(c(intersect(c("Classification", "Plate_ID"), names(df)),
                             se_cols, num_cols))
      sel_y <- if (length(num_cols) >= 2) num_cols[2] else num_cols[1]
      sel_z <- if (length(num_cols) >= 3) num_cols[3] else num_cols[1]
      tagList(
        selectInput("cls_3d_x",     "X axis:", choices=num_cols, selected=num_cols[1]),
        selectInput("cls_3d_y",     "Y axis:", choices=num_cols, selected=sel_y),
        selectInput("cls_3d_z",     "Z axis:", choices=num_cols, selected=sel_z),
        selectInput("cls_3d_color", "Color:",  choices=color_opts, selected=color_opts[1]),
        helpText(style="font-size:10px;",
                 "Tip: drag to rotate, scroll to zoom, double-click legend to isolate a group.")
      )
    })

    # 3D score scatter (plotly native, SE-joined)
    output$cls_diag_corr_3d <- renderPlotly({
      req(cls_rv$classified, input$cls_3d_x, input$cls_3d_y, input$cls_3d_z)
      df      <- .join_se_meta(cls_rv$classified, input$cls_seDataset)
      x_col   <- input$cls_3d_x
      y_col   <- input$cls_3d_y
      z_col   <- input$cls_3d_z
      col_col <- input$cls_3d_color %||% "Classification"
      req(all(c(x_col, y_col, z_col) %in% names(df)), col_col %in% names(df))
      color_vec <- df[[col_col]]
      is_numeric_color <- is.numeric(color_vec)
      hover_txt <- paste0("Well: ",   df$Well,
                          "<br>Plate: ", df$Plate_ID,
                          "<br>Class: ", df$Classification,
                          "<br>", x_col, ": ", round(df[[x_col]], 3),
                          "<br>", y_col, ": ", round(df[[y_col]], 3),
                          "<br>", z_col, ": ", round(df[[z_col]], 3))
      if (is_numeric_color) {
        plotly::plot_ly(df,
          x=~df[[x_col]], y=~df[[y_col]], z=~df[[z_col]],
          color=~color_vec,
          colors=viridis::viridis(20),
          type="scatter3d", mode="markers",
          marker=list(size=4, opacity=0.8),
          text=hover_txt, hoverinfo="text") %>%
          plotly::layout(
            scene=list(
              xaxis=list(title=x_col),
              yaxis=list(title=y_col),
              zaxis=list(title=z_col)
            ),
            coloraxis=list(colorbar=list(title=col_col))
          )
      } else {
        plotly::plot_ly(df,
          x=~df[[x_col]], y=~df[[y_col]], z=~df[[z_col]],
          color=~as.factor(color_vec),
          type="scatter3d", mode="markers",
          marker=list(size=4, opacity=0.8),
          text=hover_txt, hoverinfo="text") %>%
          plotly::layout(
            scene=list(
              xaxis=list(title=x_col),
              yaxis=list(title=y_col),
              zaxis=list(title=z_col)
            ),
            legend=list(title=list(text=col_col))
          )
      }
    })

    observeEvent(input$cls_load, {
      tryCatch({
        # Always reset downstream steps
        cls_rv$filtered   <- NULL
        cls_rv$aggregated <- NULL
        cls_rv$scored     <- NULL
        cls_rv$classified <- NULL

        if (is.null(cls_rv$particles)) {
          req(input$img_fileDB)
          withProgress(message = "Loading particles", value = 0.3, {
            l_files <- input$img_fileDB$datapath
            df <- prepareImgDF(l_files, analysis = "pa",
                               aggregate = FALSE, cleanNames = TRUE)
            if ("Image_Type" %in% names(df))
              df <- subset(df, Image_Type == "fluor")
            cls_rv$particles <- df
            incProgress(1)
          })
        }

        # Join SE metadata immediately so all downstream plots can facet by it
        se_name <- input$cls_seDataset
        df      <- cls_rv$particles
        if (!is.null(se_name) && nzchar(se_name) && !is.null(SEs[[se_name]])) {
          df <- .join_se_meta(df, se_name)
        }
        cls_rv$particles <- df

        showNotification(
          paste0("Ready: ", nrow(cls_rv$particles),
                 " particles. Downstream steps reset."),
          type = "message", duration = 4)

      }, error = function(e) showModal(modalDialog(title = "Load error",
          tags$pre(conditionMessage(e)), easyClose = TRUE)))
    })

    observeEvent(input$cls_filter, {
      req(cls_rv$particles)
      tryCatch({
        method   <- input$cls_filter_method %||% "median_ratio"
        auto_thr <- isTRUE(input$cls_threshold_auto)
        thr_v    <- input$cls_filter_threshold %||% NA_real_
        # NULL = use method default (auto median/3 for median_ratio, 0 for zscore)
        thr   <- if (!auto_thr && !is.na(thr_v)) thr_v else NULL
        gvars <- cls_scale_group_vars()
        cls_rv$filtered   <- filterParticles(cls_rv$particles,
                                              method     = method,
                                              threshold  = thr,
                                              group_vars = gvars)
        cls_rv$aggregated <- NULL
        cls_rv$scored     <- NULL
        cls_rv$classified <- NULL
      }, error=function(e) showModal(modalDialog(title="Filter error", tags$pre(conditionMessage(e)), easyClose=TRUE)))
    })

    observeEvent(input$cls_aggregate_btn, {
      src <- cls_rv$filtered %||% cls_rv$particles
      req(src)
      tryCatch({
        agg_fun      <- input$cls_agg_fun %||% "mean"
        scale_choice <- input$cls_agg_scale_within %||% "channel"
        scale_within <- switch(scale_choice,
          channel       = "Channel_Name",
          channel_plate = c("Channel_Name", "Plate_ID"),
          none          = NULL)
        scale_center <- isTRUE(input$cls_agg_scale_center)
        cls_rv$aggregated <- aggregateParticles(src,
                                                agg_fun      = agg_fun,
                                                scale_within = scale_within,
                                                scale_center = scale_center)
        cls_rv$scored     <- NULL
        cls_rv$classified <- NULL
      }, error=function(e) showModal(modalDialog(title="Aggregate error", tags$pre(conditionMessage(e)), easyClose=TRUE)))
    })

    observeEvent(input$cls_score, {
      req(cls_rv$aggregated)
      tryCatch({
        wts    <- c(Mean     = input$cls_w_mean     %||% 1,
                    Area     = input$cls_w_area     %||% 1,
                    normArea = input$cls_w_normarea %||% 1)
        center <- isTRUE(input$cls_score_center)
        cls_rv$scored     <- scoreParticles(cls_rv$aggregated, weights=wts, center=center)
        cls_rv$classified <- NULL
      }, error=function(e) showModal(modalDialog(title="Score error", tags$pre(conditionMessage(e)), easyClose=TRUE)))
    })

    # Dynamic channel labels UI — one text input per Channel_Name found in scored data
    output$cls_channel_labels_ui <- renderUI({
      src <- cls_rv$scored %||% cls_rv$aggregated
      if (is.null(src) || !"Channel_Name" %in% names(src)) return(NULL)
      channels <- sort(unique(src$Channel_Name))
      tagList(
        tags$strong("Channel display labels (optional):"),
        helpText(style="font-size:11px;", "Map raw Channel_Name to display label (e.g. C1 \u2192 GFP). Leave blank to keep raw name."),
        lapply(channels, function(ch) {
          fluidRow(
            column(4, tags$p(style="margin-top:7px; font-family:monospace;", ch)),
            column(8, textInput(paste0("cls_ch_label_", gsub("[^A-Za-z0-9]","_",ch)), NULL,
                                placeholder=paste0("label for ", ch), width="100%"))
          )
        })
      )
    })

    observeEvent(input$cls_classify, {
      req(cls_rv$scored)
      tryCatch({
        src  <- cls_rv$scored
        channels <- sort(unique(src$Channel_Name))
        ch_labels <- NULL
        for (ch in channels) {
          lbl <- input[[paste0("cls_ch_label_", gsub("[^A-Za-z0-9]","_",ch))]]
          if (!is.null(lbl) && nzchar(trimws(lbl)))
            ch_labels <- c(ch_labels, setNames(trimws(lbl), ch))
        }
        cls_rv$classified <- classifyWells(src,
                                            delta          = input$cls_delta    %||% 0.5,
                                            min_area       = input$cls_min_area %||% 0.1,
                                            channel_labels = if (length(ch_labels)>0) ch_labels else NULL)
      }, error=function(e) showModal(modalDialog(title="Classify error", tags$pre(conditionMessage(e)), easyClose=TRUE)))
    })

    output$cls_classify_preview <- renderTable({
      req(cls_rv$classified)
      df <- as.data.frame(table(Classification = cls_rv$classified$Classification))
      df[order(df$Freq, decreasing=TRUE), ]
    }, striped=TRUE, bordered=TRUE, spacing="s")

    observeEvent(input$cls_merge_se, {
      req(cls_rv$classified, input$cls_seDataset)
      tryCatch({
        SEname <- input$cls_seDataset
        req(!is.null(SEs[[SEname]]))
        col_nm <- if (input$cls_merge_mode == "nested") {
          nm <- trimws(input$cls_col_name)
          if (!nzchar(nm)) "img_classification" else nm
        } else NULL
        SEs[[SEname]] <- mergeClassificationToSE(SEs[[SEname]], cls_rv$classified, col_name=col_nm)
        SEinit(SEs[[SEname]])
        updateSelectInput(session, "object", selected=SEname, choices=union(names(objects), names(SEs)))
        showNotification(
          paste0("Classification merged into SE '", SEname, "'",
                 if (!is.null(col_nm)) paste0(" as nested column '", col_nm, "'") else " (flat columns)"),
          type="message", duration=5)
      }, error=function(e) showModal(modalDialog(title="Merge error", tags$pre(conditionMessage(e)), easyClose=TRUE)))
    })

    # ---------------------------------------------------------------
    # Save manual annotations to SE colData  (ann_update_se button)
    # ---------------------------------------------------------------

    observeEvent(input$ann_update_se, {
      se_name <- input$object
      req(!is.null(se_name), nzchar(se_name), !is.null(SEs[[se_name]]))
      se <- SEs[[se_name]]
      df <- ann_rv$results
      if (nrow(df) == 0) {
        showNotification("No annotations to save.", type="warning"); return()
      }
      # Deduplicate: keep latest per (well, plate_id, img_class)
      df <- df[order(df$timestamp, decreasing=TRUE), , drop=FALSE]
      key <- paste0(df$well, "||", ifelse(is.na(df$plate_id),"",df$plate_id), "||",
                    ifelse(is.na(df$img_class),"",df$img_class))
      df <- df[!duplicated(key), , drop=FALSE]

      # Column name: manual_ann_<person> or manual_ann
      person  <- trimws(input$ann_person %||% "")
      col_nm  <- if (nzchar(person)) paste0("manual_ann_", make.names(person)) else "manual_ann"

      cd      <- as.data.frame(SummarizedExperiment::colData(se))
      names(df)[names(df)=="well"]    <- "Well"
      names(df)[names(df)=="plate_id"]<- "Plate_ID"
      joined  <- dplyr::left_join(cd[,c("Well","Plate_ID")], df,
                                  by=c("Well","Plate_ID"))
      new_cols <- setdiff(names(joined), c("Well","Plate_ID"))
      SummarizedExperiment::colData(se)[[col_nm]] <-
        S4Vectors::DataFrame(joined[, new_cols, drop=FALSE])
      SEs[[se_name]] <- se
      SEinit(SEs[[se_name]])
      showNotification(paste0("Saved to colData column '", col_nm, "' (", nrow(df), " annotated wells)."),
                       type="message", duration=6)
    })

    if(is.null(logins)) waiter_hide()
  }
}
