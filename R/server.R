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
#' @import shiny ggplot2 SummarizedExperiment sechm waiter plotly
#' @importFrom shinydashboard updateTabItems
#' @importFrom shinyjs showElement hideElement
#' @importFrom DT datatable renderDT
#' @importFrom ComplexHeatmap draw
#' @importFrom S4Vectors metadata
#' @importFrom shinyauthr loginServer
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
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
      updateSelectInput(session, "assay_id1", choices=assayNames(x),
                        selected=assayNames(x)[1])
      updateSelectInput(session, "color_group1", choices=colnames(as.data.frame(colData(x))),
                        selected=NULL)
      updateSelectizeInput(session, "group_by_meta1",
                           choices=colnames(rowData(x)))



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

          df <- do.call(rbind, lapply(event$key, function(k) {
            parts <- strsplit(k, ",\\s*")[[1]]
            data.frame(well = parts[1], plate_id = parts[2], stringsAsFactors = FALSE)
          }))

          plate_ids <- unique(df$plate_id)
          if (length(plate_ids) != 1) return(NULL)  # Only allow one plate per event

          plate_id <- plate_ids
          current_wells <- get_wells(plate_id)

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
            current_wells <- unique(df$well)
          }

          set_wells(plate_id, current_wells)

          #if (source_id == "cluster_plot") {
          #  updateSelectizeInput(session, "plate_id3", selected = plate_id)
          #}

          update_all_select_inputs()
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

        if(!is.null(SE())){

          se <- SE()[,SE()$Plate_ID == input$plate_id1]

          if(input$assay_id1 %in% assayNames(se)){
            assayNames <- input$assay_id1
          }

          if(is.null(input$color_group1)){
            color_group <- input$group_by_meta1
          }else{
            color_group <- input$color_group1
          }


        if(length(get_wells(input$plate_id1)) > 0){
            se <- se[,se$Well %in% get_wells(input$plate_id1)]

        }

         p <-  plotAssayVSSweeps(se, assayList = assayNames,
                            rowCol = input$group_by_meta1, colorGroup = color_group)

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

    if(is.null(logins)) waiter_hide()
  }
}
