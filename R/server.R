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
#' @importFrom jpeg readJPEG
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

        # Auto-load image folder stored in SE metadata
        folder <- tryCatch(S4Vectors::metadata(x)$image_path_jpgs, error = function(e) NULL)
        if (!is.null(folder) && nzchar(folder) && dir.exists(folder)) {
          addResourcePath("imgplate", folder)
          img_parent(folder)
        }

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

  output$plate_view_col <- renderPlotly({
    # Create grid
    rows <- LETTERS[1:16]
    cols <- sprintf("%02d", 1:24)
    grid <- expand.grid(Row = rows, Col = cols)
    grid$Well <- paste0(grid$Row, grid$Col)
    grid$Selected <-  TRUE#grid$Well %in% selected_wells$data$well

    # Plot it
    # p <- ggplot(grid, aes(x = Col, y = Row)) +
    #   geom_tile(aes(fill = Selected), color = "grey50") +
    #   scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "white")) +
    #   scale_y_discrete(limits = rev) +
    #   theme_void() +
    #   theme(legend.position = "none",
    #         panel.border = element_rect(color = "black", fill = NA))

    # if(length(get_wells(input$plate_id4)) > 0 & !input$all_plates2){
    #   grid$Selected <- ifelse(grid$Well %in% get_wells(input$plate_id4), TRUE, FALSE)
    #   grid$Plate_ID <-input$plate_id4
    # }else{
    #   grid$Selected <- TRUE
    #   grid$Plate_ID <- "All_Plates"
    # }
    #melted.dat$is_selected <- ifelse(melted.dat$Well %in% input$selected_well, TRUE, FALSE)
    grid$Selected <- TRUE
    grid$key_combined <- paste(grid$Well, grid$Plate_ID, sep = ", ")

    coldata <- as.data.frame(colData(SE()))
    coldata <- coldata[coldata$Plate_ID %in% input$plate_id4,]

    if(input$condition == "" | !is.null(input$condition)){
      condition <- "Well"

    }else{
      condition <- input$condition
    }
    p <- ggplot(grid, aes(x = as.numeric(Col), y = Row, key = key_combined, fill=coldata[[condition]])) +
      geom_tile() +
      scale_x_continuous(breaks = 1:24) +
      scale_y_discrete(limits = rev) +
      geom_text(aes(label = paste(Row, Col, sep="")), color = "white") +
      theme_minimal() + labs(fill=input$condition, xlab="Column")

    ggplotly(p, source = "condition_plot")


  })


  observeEvent(input$createCondition, {
    tryCatch({
      se <- SE()
    if(!is.null(input$newCondition) | input$newCondition == ""){

      se[[input$newCondition]] <- "init"
      SEs[[input$object]] <- se
      SEinit(SEs[[input$object]])

      coldat <- colnames(as.data.frame(colData(se))[ , !purrr::map_lgl(as.data.frame(colData(se)), is.numeric)])

      updateSelectInput(session, "condition", choices=coldat, selected = input$newCondition)
    }
      }, error=function(e){
        print(conditionMessage(e))
        print(traceback())
        showModal(modalDialog(easyClose=TRUE, title="Error with upload",
                              "The file was not recognized. Are you sure that it is a R .rds file?",
                              tags$pre(e)))
      })


  })


  # Toggle visibility
  observeEvent(input$toggle_plate, {
    toggle("mini_plate_plot")
  })


  ## Download-Tab

  output$downloadRDS <- downloadHandler(
    filename = function() {
      #paste(input$download_table, ' Edited Table.csv', sep='')
      paste(input$rdsObject, ".rds", sep="")
    },
    content = function(file) {
      saveRDS(SEs[[input$rdsObject]], file)
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

    .coords_to_well <- function(x, y) {
      col_num <- floor(x) + 1L
      row_idx <- 16L - floor(y)
      if (col_num < 1 || col_num > 24 || row_idx < 1 || row_idx > 16)
        return(NA_character_)
      paste0(.rows384[row_idx], sprintf("%02d", col_num))
    }

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
      dirs <- img_plate_dirs()
      req(length(dirs) > 0)
      updateSelectInput(session, "plate_img_id",
                        choices  = names(dirs),
                        selected = names(dirs)[1])
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

    # Update channel / class dropdowns when map changes
    observeEvent(img_well_map(), {
      df <- img_well_map()
      channels <- if (nrow(df) > 0) sort(unique(df$channel)) else character(0)
      classes  <- if (nrow(df) > 0) sort(unique(df$img_class)) else character(0)
      updateSelectInput(session, "img_channel",   choices = channels,
                        selected = if (length(channels) > 0) channels[1] else NULL)
      updateSelectInput(session, "img_class_sel", choices = classes,
                        selected = if (length(classes) > 0) classes[1] else NULL)
    }, ignoreNULL = TRUE)

    # Filtered map: named vector well -> file for current channel + class
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

    # Load jpegs for current plate/channel/class only
    img_loaded <- reactive({
      m    <- img_well_map_filtered()
      fps  <- na.omit(unname(m))
      imgs <- list()
      for (fp in fps)
        imgs[[fp]] <- tryCatch(jpeg::readJPEG(fp), error = function(e) NULL)
      imgs
    })

    output$img_plate_stats <- renderText({
      dirs <- img_plate_dirs()
      pid  <- input$plate_img_id
      df   <- img_well_map()
      m    <- img_well_map_filtered()
      folder <- if (!is.null(pid) && pid %in% names(dirs))
        basename(dirs[[pid]]) else "—"
      channels <- if (nrow(df) > 0) paste(sort(unique(df$channel)), collapse = ", ") else "—"
      classes  <- if (nrow(df) > 0) paste(sort(unique(df$img_class)), collapse = ", ") else "—"
      paste0(
        "Folder:      ", folder,              "\n",
        "JPGs found:  ", length(img_files()), "\n",
        "Channels:    ", channels,            "\n",
        "Classes:     ", classes,             "\n",
        "Wells shown: ", sum(!is.na(m)), " / 384"
      )
    })

    output$img_plate_plot <- renderPlot({
      m    <- img_well_map_filtered()
      imgs <- img_loaded()

      par(mar = c(0, 0, 0, 0), bg = "#111111")
      plot(NULL, xlim = c(0, 24), ylim = c(0, 16),
           xlab = "", ylab = "", axes = FALSE,
           xaxs = "i", yaxs = "i")

      all_wells <- as.vector(t(outer(.rows384, sprintf("%02d", 1:24), paste0)))

      # Pass 1: images / empty cells
      for (well in all_wells) {
        row_idx <- which(.rows384 == substr(well, 1, 1))
        col_num <- as.integer(substr(well, 2, 3))
        x0 <- col_num - 1; x1 <- col_num
        y0 <- 16 - row_idx; y1 <- y0 + 1

        fp  <- m[[well]]
        img <- if (!is.null(fp) && !is.na(fp)) imgs[[fp]] else NULL
        if (!is.null(img))
          rasterImage(img, x0, y0, x1, y1, interpolate = TRUE)
        else
          rect(x0, y0, x1, y1, col = "#2a2a2a", border = NA)
      }

      # Pass 2: colData overlay
      vals <- img_coldata_vals()
      if (!is.null(vals) && length(vals) > 0) {
        cols <- .vals_to_colors(vals, alpha = input$img_overlay_alpha)
        for (well in names(vals)) {
          col <- cols[well]
          if (is.na(col)) next
          row_idx <- which(.rows384 == substr(well, 1, 1))
          col_num <- as.integer(substr(well, 2, 3))
          rect(col_num - 1, 16 - row_idx, col_num, 17 - row_idx,
               col = col, border = NA)
        }
      }
    },
    width  = function() session$clientData$output_img_plate_plot_width,
    height = function() session$clientData$output_img_plate_plot_height,
    bg     = "#111111")

    # Hover: show well + color variable + any extra selected colData columns
    output$img_hover_info <- renderText({
      h <- input$img_plate_hover
      if (is.null(h)) return("")
      well <- .coords_to_well(h$x, h$y)
      if (is.na(well)) return("")
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
    observeEvent(input$img_plate_click, {
      cl   <- input$img_plate_click
      well <- .coords_to_well(cl$x, cl$y)
      if (is.na(well)) return()

      parent <- img_parent()
      if (is.null(parent)) return()
      if (!endsWith(parent, "/")) parent <- paste0(parent, "/")

      fp_to_url <- function(fp) {
        if (is.null(fp) || is.na(fp)) return(NA_character_)
        rel <- substring(normalizePath(fp, winslash = "/"), nchar(parent) + 1)
        file.path("imgplate", rel)
      }

      # Find BF + fluoro for this well + current class selection
      df  <- img_well_map()
      cls <- input$img_class_sel
      sub <- df[df$well == well &
                  (if (!is.null(cls) && nzchar(cls %||% "")) df$img_class == cls else TRUE),
                , drop = FALSE]
      bf_row  <- sub[tolower(sub$channel) == "bf", , drop = FALSE]
      flu_row <- sub[tolower(sub$channel) != "bf", , drop = FALSE]
      bf_url  <- if (nrow(bf_row)  > 0) fp_to_url(bf_row$file[1])  else NA_character_
      flu_url <- if (nrow(flu_row) > 0) fp_to_url(flu_row$file[1]) else NA_character_
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

      make_panel <- function(url, label) {
        if (!is.na(url))
          div(style = "text-align: center;",
              tags$img(src = url,
                       style = "max-width: 100%; max-height: 55vh; border-radius: 6px;"),
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
            make_panel(bf_url,  "BF"),
            make_panel(flu_url, "Fluoro")),
        easyClose = TRUE, size = "l", footer = modalButton("Close")
      ))
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

    # Load existing CSV when folder changes
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
      tryCatch(img_plate_dirs(), error = function(e) {
        parent <- img_parent()
        if (is.null(parent)) return(list())
        subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
        if (length(subdirs) == 0) return(setNames(list(parent), basename(parent)))
        setNames(as.list(subdirs), basename(subdirs))
      })
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

    # Parse all files into a meta data.frame: well / channel / img_class / file
    ann_img_meta_df <- reactive({
      files <- all_img_files()
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

    # Populate ann_img_class dropdown from available img_class values
    observeEvent(ann_img_meta_df(), {
      df  <- ann_img_meta_df()
      cls <- if (nrow(df) > 0) sort(unique(df$img_class)) else character(0)
      updateSelectInput(session, "ann_img_class", choices = cls,
                        selected = if (length(cls) > 0) cls[1] else NULL)
    }, ignoreNULL = TRUE)

    # Rebuild well queue when files, class, or plate selection changes
    observeEvent(list(ann_img_meta_df(), input$ann_img_class, input$ann_plate_ids), {
      df   <- ann_img_meta_df()
      cls  <- input$ann_img_class
      pids <- input$ann_plate_ids
      if (nrow(df) == 0 || is.null(cls) || !nzchar(cls %||% "")) {
        ann_rv$shuffled <- character(0); ann_rv$idx <- 0L; return()
      }
      sub <- df[df$img_class == cls, , drop = FALSE]
      # Filter by selected plates when a selection exists
      if (!is.null(pids) && length(pids) > 0) {
        dirs         <- ann_plate_dirs()
        plate_of_file <- vapply(sub$file, .ann_file_plate, character(1),
                                plate_dirs = dirs)
        sub <- sub[!is.na(plate_of_file) & plate_of_file %in% pids, , drop = FALSE]
      }
      available <- unique(sub$well[!is.na(sub$well)])
      if (length(available) == 0) { ann_rv$shuffled <- character(0); ann_rv$idx <- 0L; return() }
      done_wells <- ann_rv$results$well[
        !is.na(ann_rv$results$img_class) & ann_rv$results$img_class == cls]
      pool <- setdiff(available, done_wells)
      if (length(pool) == 0) pool <- available
      ann_rv$shuffled <- sample(pool)
      ann_rv$idx      <- 1L
    }, ignoreNULL = FALSE)

    # Current item is a well ID (e.g. "A01")
    ann_current <- reactive({
      if (length(ann_rv$shuffled) == 0 || ann_rv$idx < 1) return(NA_character_)
      ann_rv$shuffled[[ann_rv$idx]]
    })

    # BF + fluoro file paths for the current well + selected class
    ann_pair <- reactive({
      well <- ann_current()
      cls  <- input$ann_img_class
      if (is.na(well) || is.null(cls) || !nzchar(cls %||% ""))
        return(list(bf = NA_character_, fluoro = NA_character_))
      df      <- ann_img_meta_df()
      sub     <- df[df$well == well & df$img_class == cls, , drop = FALSE]
      bf_row  <- sub[tolower(sub$channel) == "bf", , drop = FALSE]
      flu_row <- sub[tolower(sub$channel) != "bf", , drop = FALSE]
      list(
        bf     = if (nrow(bf_row)  > 0) bf_row$file[1]  else NA_character_,
        fluoro = if (nrow(flu_row) > 0) flu_row$file[1] else NA_character_
      )
    })

    ann_advance <- function() {
      n <- length(ann_rv$shuffled)
      if (n == 0) return()
      nxt <- ann_rv$idx %% n + 1L
      if (nxt == 1L) ann_rv$shuffled <- sample(ann_rv$shuffled)
      ann_rv$idx <- nxt
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

    output$ann_image_ui <- renderUI({
      if (length(ann_rv$shuffled) == 0)
        return(tags$p(style = "color: #aaa; padding: 20px;",
                      "Select a folder in the Image Plate Viewer tab first."))
      well <- ann_current()
      if (is.na(well)) return(tags$p("No well available."))

      pair <- ann_pair()

      make_panel <- function(fp, label) {
        ok  <- !is.null(fp) && !is.na(fp)
        url <- if (ok) .ann_file_to_url(fp) else NA_character_
        if (ok && !is.na(url)) {
          div(style = "text-align: center;",
              tags$img(src = url,
                       style = "max-width: 100%; max-height: 52vh; border-radius: 6px;"),
              tags$p(style = "font-size: 12px; color: #666; margin: 4px 0 0;", label))
        } else {
          div(style = paste0("display:flex; align-items:center; justify-content:center;",
                             " min-height:200px; background:#2a2a2a; border-radius:6px; color:#888;"),
              tags$p(paste(label, "\u2014 no image")))
        }
      }

      div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 12px; padding: 8px;",
          make_panel(pair$bf,     "BF"),
          make_panel(pair$fluoro, "Fluoro"))
    })

    output$ann_class_buttons <- renderUI({
      req(length(ann_rv$classes) > 0)
      div(style = "display: grid; gap: 8px;",
          lapply(ann_rv$classes, function(lab) {
            tags$button(lab, class = "btn btn-success btn-block",
                        style = "width:100%; padding:10px; border-radius:8px; font-size:14px;",
                        onclick = sprintf(
                          "Shiny.setInputValue('ann_class_selected', %s, {priority:'event'});",
                          jsonlite::toJSON(lab, auto_unbox = TRUE)))
          })
      )
    })

    observeEvent(input$ann_class_selected, {
      well    <- ann_current(); req(!is.na(well))
      cls     <- input$ann_img_class %||% NA_character_
      person  <- trimws(input$ann_person)
      se_name <- input$object
      se      <- if (!is.null(se_name) && nzchar(se_name)) SEs[[se_name]] else NULL

      # Get any file for this well+class to infer plate ID from directory structure
      df     <- ann_img_meta_df()
      sub    <- df[df$well == well & !is.na(df$img_class) & df$img_class == cls, , drop = FALSE]
      fp_pid <- if (nrow(sub) > 0) sub$file[1] else NA_character_
      pid    <- if (!is.na(fp_pid))
        .ann_file_plate(fp_pid, ann_plate_dirs()) else NULL

      row <- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        well      = well,
        plate_id  = pid %||% NA_character_,
        img_class = cls,
        label     = input$ann_class_selected,
        person    = person,
        stringsAsFactors = FALSE
      )
      ann_rv$results <- rbind(ann_rv$results, row)
      cp <- ann_csv_path()
      if (!is.null(cp))
        write.table(row, cp, sep = ",", row.names = FALSE,
                    col.names = !file.exists(cp), append = file.exists(cp))
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
    observeEvent(input$ann_skip,          { ann_advance() })

    observeEvent(input$ann_undo, {
      if (nrow(ann_rv$results) == 0) return()
      last_well        <- ann_rv$results$well[nrow(ann_rv$results)]
      ann_rv$results   <- ann_rv$results[-nrow(ann_rv$results), , drop = FALSE]
      ann_rv$shuffled  <- c(last_well, ann_rv$shuffled)
      ann_rv$idx       <- 1L
      cp <- ann_csv_path()
      if (!is.null(cp) && file.exists(cp)) {
        lines <- readLines(cp)
        if (length(lines) > 1) writeLines(lines[-length(lines)], cp)
      }
    })

    observeEvent(input$ann_clear_results, {
      ann_rv$results <- ann_rv$results[0, ]
      cp <- ann_csv_path()
      if (!is.null(cp) && file.exists(cp)) file.remove(cp)
    })

    output$ann_progress <- renderText({
      cls       <- input$ann_img_class %||% ""
      df        <- ann_img_meta_df()
      total     <- if (nrow(df) > 0 && nzchar(cls))
        length(unique(df$well[!is.na(df$img_class) & df$img_class == cls])) else 0L
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

    if(is.null(logins)) waiter_hide()
  }
}
