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
#' @import shiny ggplot2 SummarizedExperiment sechm waiter
#' @importFrom plotly ggplotly renderPlotly event_data
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

    SEinit <- function(x){
      if(!is.null(logins)) req(credentials()$user_auth)
      if(is.null(assayNames(x)))
        assayNames(x) <- paste0("assay",1:length(assays(x)))
      if(ncol(rowData(x))==0) rowData(x)$name <- row.names(x)
      updateSelectizeInput(session, "gene_input", selected=isolate(selGene()),
                           choices=sort(unique(row.names(x))), server=TRUE)
      colvars <- colnames(colData(x))
      updateSelectInput(session, "input_hm_samples", choices=colnames(x),
                        selected=colnames(x))
      updateSelectInput(session, "assay_input", choices=assayNames(x),
                        selected=getDef(x, "assay", rev(assayNames(x))))
      updateSelectInput(session, "assay_input2", choices=assayNames(x),
                        selected=getDef(x, "assay", rev(assayNames(x))))
      updateSelectizeInput(session, "hm_order", choices=colvars)
      updateSelectizeInput(session, "hm_anno", choices=colvars,
                           selected=getDef(x, c("groupvar","colvar")))
      updateSelectizeInput(session, "hm_gaps", choices=colvars,
                           selected=getDef(x, "gridvar"))
      updateSelectInput(session, "select_groupvar", choices=colvars,
                        selected=getDef(x, "groupvar"))
      updateSelectInput(session, "select_colorvar", choices=colvars,
                        selected=getDef(x, "colvar"))
      updateSelectInput(session, "select_gridvars", choices=colvars,
                        selected=getDef(x, "gridvar"))
      x <- mergeFlists(x)
      if(length(sefl <- names(metadata(x)$feature.lists))==0){
        updateSelectInput(session, "genelist_input", choices="")
      }else{
        updateSelectInput(session, "genelist_input", choices=sefl)
      }
      if(!is.null(getDef(x, "assay")))
        updateCheckboxInput(session, "hm_scale",
                            value=!grepl("FC$",getDef(x, "assay")))
      x
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
          showModal(modalDialog(easyClose=TRUE, title="Error with upload",
                                "The file was not recognized. Are you sure that it is a R .rds file?",
                                tags$pre(e)))
        })
    })

    observeEvent(input$fileEphys, {
      tryCatch({
        if(!is.null(input$fileEphys)){
          #x <- readRDS(input$fileEphys$datapath)

          l_files <- ifelse(length(input$fileEphys$datapath) > 1,
                            as.list(input$fileEphys$datapath),
                            input$fileEphys$datapath)


          withProgress(message = 'Loading Excel-Files into SE', value = 0, {
          incProgress(0.5, detail = "This may take a while..")
          x <- prepareSE(l_files)

          print(x)
          if(is(x,"SingleCellExperiment")){
            SEname <- input$se_id
            SEs[[SEname]] <- x
            SEinit(SEs[[SEname]])
            incProgress(0.75, detail = "Updating UI")
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))
          }
          incProgress(1, detail = "SE loaded")
          })
          }else{
            stop("The object is not a SummarizedExperiment!")
          }
        }, error=function(e){
          showModal(modalDialog(easyClose=TRUE, title="Error with upload",
                                "The file was not recognized. Are you sure that it is an .excel file?",
                                tags$pre(e)))
        })
    })

    observeEvent(input$mergeSE, {
      tryCatch({
        if(!is.null(input$fileEphys) & !is.null(input$fileDB)){
          #x <- readRDS(input$fileEphys$datapath)
          l_files <- ifelse(length(input$fileDB$datapath) > 1,
                            as.list(input$fileDB$datapath),
                            input$fileDB$datapath)

          withProgress(message = 'Loading Imaging Results', value = 0, {
            incProgress(0.5, detail = "This may take a while..")

            req(input$tabletype)
            for(tabletype in input$tabletype){
              df_img <- prepareImgDF(l_files, analysis = tabletype)
              SEname <- input$se_id
              SEs[[SEname]] <-  mergeSEandImg(SEs[[SEname]], df_img,
                                              tableType = tabletype)

            }

            #se <- SEs[[SEname]]
            #se <- mergeSEandImg(se, df_img)
            #SEs[[SEname]] <- se
            SEinit(SEs[[SEname]])
            incProgress(0.75, detail = "Updating UI")
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))

            incProgress(1, detail = "SE updated")
          })
        }else{
          stop("The object is not a SummarizedExperiment!")
        }
      }, error=function(e){
        showModal(modalDialog(easyClose=TRUE, title="Error with upload",
                              "The file was not recognized. Are you sure that it is an .excel file?",
                              tags$pre(e)))
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

          for(i in seq_along(input$datasets)){
            dataset <- input$datasets[i]
            SEname <- selected_label[i]
            if(!(SEname %in% names(SEs))){
            data(list = dataset, package = "ephacRTools", envir = environment())

            SEs[[SEname]] <- get(dataset, envir = environment())
            SEinit(SEs[[SEname]])
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))
            }
          }

        }}, error=function(e){
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
            tags$p(paste("A SummarizedExperiment with ", ncol(SE()), "samples and ",
                         length(unique(SE()$Plate_ID)), " with electrophysiology plates: \n",
                         paste(unique(SE()$Plate_ID), "\n"))),
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
    ### BEGIN DEAs

    DEAs <- reactive({
      if(is.null(SE())) return(list())
      RD <- rowData(SE())
      deas <- grep("^DEA\\.", colnames(RD), value=TRUE)
      deas <- deas[unlist(sapply(deas, FUN=function(x) is.data.frame(RD[[x]]) ||
                                   is(RD[[x]], "DFrame") ))]
      if(length(deas)==0) return(list())
      lapply(setNames(deas, gsub("^DEA\\.","",deas)), FUN=function(x) RD[[x]])
    })

    output$menu_DEA <- renderUI({
      if(is.null(DEAs()) || length(DEAs())==0)
        return(menuItem("DEA results", tabName="tab_dea",
                        badgeLabel="N/A", badgeColor="red"))
      menuItem("DEA results", tabName="tab_dea", badgeLabel=length(DEAs()),
               badgeColor="aqua")
    })

    output$dea_input <- renderUI({
      hideElement("dea_box")
      if(length(DEAs())==0)
        return(tags$h5(icon("exclamation-triangle"),
                       "The object contains no DEA."))
      showElement("dea_box")
      choices <- sapply(names(DEAs()), FUN=function(x){
        if(!is.null(a <- attr(DEAs()[[x]], "description"))) return(a)
        x
      })
      choices <- setNames(names(DEAs()), choices)
      options = list( render = I(paste("{
        item:", DEAselRender, ",
        option:", DEAselRender, "}"))
      )
      selectizeInput("dea", label="Differential expression analysis",
                     choices=choices, options=options)
    })

    DEA <- reactive({
      if(is.null(input$dea) || input$dea=="" || is.null(DEAs()[[input$dea]]))
        return(NULL)
      .homogenizeDEA(DEAs()[[input$dea]])
    })

    observeEvent(input$object, {
      if(!is.null(input$object) && input$object != "" && input$main_tabs=="tab_about"){
        if(!is.null(previous_sel()))
          updateTabItems(session, "main_tabs", selected="tab_object")
        previous_sel(input$object)
      }
    })

    observeEvent(input$dea_geneFilt, {
      if(length(g <- input$dea_table_rows_all)>0){
        g <- row.names(DEA())[head(g,maxPlot)]
        updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
        updateTabItems(session, "main_tabs", selected="tab_heatmap")
      }
    })
    observeEvent(input$dea_geneSel, {
      if(length(g <- input$dea_table_rows_selected)>0){
        g <- row.names(DEA())[head(g,maxPlot)]
        updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
        updateTabItems(session, "main_tabs", selected="tab_heatmap")
      }
    })

    output$dea_overview <- renderText({
      if(is.null(dea <- DEA())) return(NULL)
      paste("A differential expression analysis across", sum(!is.na(dea$FDR)),
            "features, ", sum(dea$FDR<0.05, na.rm=TRUE), "of which are at FDR<0.05.")
    })

    output$dea_pvalues <- renderPlot({
      if(is.null(dea <- DEA()) || is.null(dea$PValue)) return(NULL)
      hist(dea$PValue, xlab="Unadjusted p-values", main="")
    })

    output$dea_table <- renderDT({
      if(is.null(dea <- DEA())) return(NULL)
      datatable( dround(as.data.frame(dea), digits=3, roundGreaterThan1=TRUE),
                 filter="top", class="compact", extensions=c("ColReorder"),
                 options=list( pageLength=30, dom = "fltBip" ) )
    })

    output$dea_volcano <- renderPlotly({
      if(is.null(dea <- DEA())) return(NULL)
      if(is.null(dea$logFC)) return(NULL)
      dea <- head(dea, 2000)
      dea <- dround(as.data.frame(dea), digits=3, roundGreaterThan1=TRUE)
      dea$logFC <- round(dea$logFC,2)
      if(is.null(dea$MeanExpr)){
        if(!is.null(dea$logCPM)) dea$MeanExpr <- dea$logCPM
        if(!is.null(dea$baseMean)) dea$MeanExpr <- dea$baseMean
      }
      dea$feature <- row.names(dea)
      if(is.null(dea$MeanExpr)){
        p <- ggplot(dea, aes(logFC, -log10(FDR), Feature=feature, FDR=FDR))
      }else{
        p <- ggplot(dea, aes(logFC, -log10(FDR), Feature=feature, FDR=FDR,
                             MeanExpr=MeanExpr, colour=MeanExpr))
      }
      p <- p + geom_vline(xintercept=0, linetype="dashed") + geom_point() +
        theme_classic()
      plotlyObs$resume()
      ggplotly(p, tooltip=c("Feature","logFC","PValue","FDR"), source="volcano")
    })

    plotlyObs <- observeEvent(event_data("plotly_click", "volcano",
                                         priority="event"), suspended=TRUE, {
                                           if(is.null(dea <- DEA())) return(NULL)
                                           event <- event_data("plotly_click", "volcano")
                                           if(!is.list(event) || is.null(event$pointNumber) && event$pointNumber>=0)
                                             return(NULL)
                                           g <- row.names(dea)[as.integer(event$pointNumber+1)]
                                           selGene(g)
                                           updateTabItems(session, "main_tabs", "tab_gene")
                                         })

    output$dea_download <- downloadHandler(
      filename = function() {
        if(is.null(DEA())) return(NULL)
        paste0(make.names(input$dea),".csv")
      },
      content = function(con) {
        if(is.null(dea <- DEA())) return(NULL)
        write.csv(dea, con)
      }
    )

    ### END DEAs
    ############
    ### BEGIN EAs

    EAs <- reactive({
      if(is.null(SE()) || is.null(eas <- metadata(SE())$EA) || !is.list(eas) ||
         length(eas)==0) return(list())
      if(is.list(eas[[1]]) && !is.data.frame(eas[[1]]) &&
         !is(eas[[1]], "DFrame")) eas <- unlist(eas, recursive=FALSE)
      eas
    })

    output$menu_Enrichments <- renderUI({
      if(length(EAs())==0) return(menuItem("Enrichments", tabName="tab_ea",
                                           badgeLabel="N/A", badgeColor="red"))
      menuItem("Enrichments", tabName="tab_ea", badgeLabel=length(EAs()),
               badgeColor="aqua")
    })

    output$ea_input <- renderUI({
      hideElement("ea_box")
      if(length(EAs())==0)
        return(tags$h5(icon("exclamation-triangle"),
                       "The object contains no enrichment analysis"))
      showElement("ea_box")
      choices <- sapply(names(EAs()), FUN=function(x){
        if(!is.null(a <- attr(EAs()[[x]], "description"))) return(a)
        x
      })
      choices <- setNames(names(EAs()), choices)
      options = list( render = I(paste("{
        item:", DEAselRender, ",
        option:", DEAselRender, "}"))
      )
      selectizeInput("ea", label="Enrichment analysis", choices=choices,
                     options=options)
    })

    EA <- reactive({
      if(is.null(input$ea) || input$ea=="" || is.null(ea <- EAs()[[input$ea]]))
        return(NULL)
      ea
    })

    output$ea_table <- renderDT({
      if(is.null(ea <- EA())) return(NULL)
      if(!isTRUE(input$ea_showGenes)) ea$genes <- ea$leadingEdge <- NULL
      datatable( dround(as.data.frame(ea), digits=3, roundGreaterThan1=TRUE),
                 filter="top", class="compact", extensions=c("ColReorder"),
                 options=list( pageLength=30, dom = "fltBip", colReorder=TRUE ) )
    })

    observeEvent(input$ea_geneSel, {
      if(is.null(ea <- EA()) || is.null(ea$genes)) return(NULL)
      if(length(RN <- input$ea_table_rows_selected)>0){
        g <- head(unique(unlist(.getWordsFromString(ea[RN,"genes"]))), maxPlot)
        updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
        updateTabItems(session, "main_tabs", selected="tab_heatmap")
      }
    })

    ### END EAs
    ############
    ### BEGIN HEATMAP

    selGenes <- reactive({
      g <- gsub(",|\n|\r|;,"," ",input$input_genes)
      g <- strsplit(g," ",fixed=T)[[1]]
      g <- unique(g[which(g!="")])
      grepGene(row.names(SE()), g)
    })


    output$heatmap <- renderPlot({
      validate( need(!is.null(SE()), "No SummarizedExperiment loaded.") )
      g <- selGenes()
      validate( need(length(g)>0,
                     message=paste("No gene selected. Select one in the ",
                                   "dropdown list above. You can type the first few letters in the box ",
                                   "and select from the matching suggestions.")) )
      if(length(g)>2 && input$hm_clusterRow){
        srow <- 1:ncol(SE())
      }else{
        srow <- NULL
      }
      se <- SE()[g, input$input_hm_samples]
      validate(need(ncol(se)>1,
                    message="At least two samples must be selected."))
      o <- input$hm_order
      if(is.null(o)) o <- c()
      for(f in rev(o)) se <- se[,order(colData(se)[[f]])]
      breaks <- NULL
      if(grepl("logFC|log2FC|scaledLFC|zscore", input$assay_input2,
               ignore.case=TRUE) || input$hm_scale)
        breaks <- 1-(input$hm_breaks/100)
      draw(sechm(se, g, do.scale=input$hm_scale,
                 assayName=input$assay_input2, show_colnames=input$hm_showColnames,
                 sortRowsOn=srow, top_annotation=input$hm_anno, gaps_at=input$hm_gaps,
                 cluster_cols=input$hm_clusterCol, cluster_rows=FALSE,
                 breaks=breaks),
           merge_legends=TRUE, padding=unit(c(0.2,0.2,1.5,0.2), "cm"))

    })

    ### END HEATMAP
    ############




    ############
    ### START GENE TAB

    selGene <- reactiveVal()
    observeEvent(input$gene_input, {
      selGene(input$gene_input)
    })

    output$gene_plot <- renderPlot({
      d <- tryCatch(meltSE(SE(), selGene(), rowDat.columns = NA),
                    error=function(x){ NULL })
      validate( need(!is.null(d) && nrow(d)>0,
                     message=paste("No gene selected. Select one in the ",
                                   "dropdown list above. You can type the first few letters in the box ",
                                   "and select from the matching suggestions.")) )
      d <- d[d$sample %in% input$input_hm_samples,]
      validate( need(nrow(d)>0, message="No sample selected!") )
      gr <- input$select_groupvar
      if(input$asfactor) d[[gr]] <- factor(d[[gr]])
      p <- ggplot(d, aes_string(input$select_groupvar, input$assay_input,
                                colour=input$select_colorvar,
                                fill=input$select_colorvar))
      if(grepl("logFC|log2FC|scaledLFC", input$assay_input)){
        p <- p + geom_hline(yintercept=0, linetype="dashed", colour="grey")
      }
      if(input$plottype_input=="violin plot"){
        p <- p + geom_violin()
      }else{
        p <- p + geom_boxplot(outlier.shape = NA)
      }
      if(input$select_plotpoints){
        p <- p + geom_beeswarm(dodge.width=1, cex=1.4, size=2)
      }else{
        p <- p + stat_summary()
      }
      p <- p + theme_classic() + ggtitle(input$gene_input) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      if(!is.null(input$select_gridvars)){
        form <- as.formula(paste0("~",paste(input$select_gridvars, collapse="+")))
        if(input$select_freeaxis){
          p <- p + facet_wrap(form, scales="free_y")
        }else{
          p <- p + facet_wrap(form)
        }
      }

      if(input$select_colorvar %in% names(metadata(SE())$anno_colors)){
        cols <- metadata(SE())$anno_colors[[input$select_colorvar]]
        cols2 <- sapply(cols, maketrans)
        p <- p + scale_color_manual(values=cols) +
          scale_fill_manual(values=cols2)
      }

      p + ggtitle(selGene())
    })

    output$gene_dea_table <- renderTable({
      if(is.null(selGene())) return(NULL)
      if(is.null(deas <- DEAs()) || length(deas)==0) return(NULL)
      d <- dplyr::bind_rows(lapply(deas, FUN=function(x){
        .homogenizeDEA(x[selGene(),])
      }), .id="Comparison")
      d <- d[!is.na(d$FDR),]
      d$LR <- NULL
      if(!is.null(d$logCPM)) d$meanExpr <- NULL
      if(!is.null(d$logFC)) d <- d[!is.na(d$logFC),]
      if(nrow(d)==0) return(NULL)
      dea_desc <- sapply(deas, FUN=function(x){
        a <- attr(x, which="description")
        if(is.null(a)) return("")
        as.character(a[1])
      })
      dea_desc[sapply(dea_desc,is.null)] <- ""
      if(any(dea_desc != "")){
        d$test_description <- dea_desc[d$Comparison]
      }
      row.names(d) <- NULL
      d
    })

    output$gene_inList <- renderText({
      if(is.null(g <- selGene()) || g=="" || length(flists())==0)
        return("")
      if(grepl(".+\\.[a-zA-Z].+", g)) g <- gsub("^[^.]+\\.","",g)
      g <- strsplit(g,"/")[[1]]
      x <- which(sapply(flists(), FUN=function(x){ any(g %in% x) }))
      if(length(x)==0)
        return("This feature is included in none of the registered lists.")
      out <- "This feature is included in the following list(s):"
      paste(c(out, names(flists())[x]), collapse="\n")
    })


    ### END GENE TAB
    ############
    ### START feature.lists TAB

    output$genelist_size <- renderText({
      if(is.null(input$genelist_input) ||
         is.null(gl <- flists()[[input$genelist_input]])) return(NULL)
      length(gl)
    })
    output$genelist_out <- renderText({
      if(is.null(input$genelist_input) ||
         is.null(gl <- flists()[[input$genelist_input]])) return(NULL)
      paste(gl, collapse=", ")
    })
    observeEvent(input$btn_importGenelist, {
      if(is.null(input$genelist_input) ||
         is.null(g <- flists()[[input$genelist_input]])) return(NULL)
      g <- head(g,maxPlot)
      updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
      updateTabItems(session, "main_tabs", selected="tab_heatmap")
    })


    ### END feature.lists TAB
    ############


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
