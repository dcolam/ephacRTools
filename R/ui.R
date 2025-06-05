#' tiny Summarized Experiment Viewer
#'
#' @param objects A named list of (paths to)
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} objects
#' @param title The title of the app (displayed in the header)
#' @param waiterContent Optional content of the loading mask; should be a
#' `tagList`, NULL to use default, or FALSE to disable the waiter
#' @param about Optional content of the introduction page (NULL to disable
#' intro page)
#' @param skin The dashboard skin color, passed to
#'   \code{\link[shinydashboard]{dashboardPage}}.
#' @param uploadMaxSize The maximum upload size. Set to zero to disable upload.
#' @param logins An optional dataframe containing possible logins. Must contain
#'   the columns "user" and "password_hash" (sodium-encoded). Not providing the
#'   argument disables login.
#' @param ... Passed to \code{\link{tinySEV.server}}
#'
#' @return Launches a shiny app
#' @import shiny
#' @export
tinySEV <- function(objects=NULL, title="tinySEV", waiterContent=NULL,
                    about=NULL, skin="blue", uploadMaxSize=1000*1024^2,
                    logins=NULL, ...){
  shinyApp(tinySEV.ui(title, waiterContent, about, skin=skin,
                      hasLogin=!is.null(logins)),
           tinySEV.server(objects, uploadMaxSize, logins=logins, ...))
}


#' tinySEV.ui
#'
#' @param title The title of the app (displayed in the header)
#' @param waiterContent Optional content of the loading mask; should be a
#' `tagList`, NULL to use default, or FALSE to disable the waiter
#' @param about Optional content of the introduction page (NULL to disable
#' intro page)
#' @param skin The dashboard skin color, passed to
#'   \code{\link[shinydashboard]{dashboardPage}}.
#' @param hasLogin Logical; whether login is required (credentials must also be
#'   provided to the server function). Default FALSE.
#'
#' @return a shiny UI
#' @export
#' @import shiny shinydashboard shinyjqui waiter
#' @importFrom shinycssloaders withSpinner
#' @importFrom plotly plotlyOutput
#' @importFrom shinyjs useShinyjs
#' @importFrom DT DTOutput
#' @importFrom shinyauthr loginUI
tinySEV.ui <- function(title="tinySEV", waiterContent=NULL, about=NULL,
                       skin="blue", hasLogin=FALSE){
  if(is.null(waiterContent) || isTRUE(waiterContent)){
    waiterContent <- tagList(
      tags$h3("Please wait while the application is initialized..."), spin_1())
  }
  if(isFALSE(waiterContent)){
    waiterContent <- NULL
  }else{
    waiterContent <- waiter_show_on_load(html=waiterContent)
  }
  if(hasLogin) waiterContent <- tagList(shinyauthr::loginUI("login"))
  aboutMenu <- NULL
  if(!is.null(about)) aboutMenu <- menuItem("About", tabName="tab_about")

  shinyUI( dashboardPage(skin=skin,
                         dashboardHeader(title=title,
                                         tags$li(class="dropdown",
                                                 actionLink("quickStart", label="Quick start", icon=icon("question")),
                                                 style="float: right; width: 112px;"),
                                         tags$li(class = "dropdown",
                                                 tags$div("Dataset: ", style="margin: 10px; font-weight: bold; font-size: 16px; color: #fff; float: left;"),
                                                 div(
                                                   selectInput("object", label=NULL, choices=c(), selectize=FALSE, width = "100%"),
                                                   style= "width: 80%; margin: 5px auto; display: block; float: left;"),
                                                 div(" "),
                                                 style="width: 80%; display: block;")
                         ),
                         dashboardSidebar(collapsed=hasLogin, disable=hasLogin,
                                          sidebarMenu(id="main_tabs", aboutMenu,
                                                      .modify_stop_propagation(
                                                        menuItem("Prepare Object", startExpanded=TRUE,
                                                                 menuSubItem("Overview", tabName="tab_object"),
                                                                 menuItemOutput("uploadMenu"),
                                                                 menuSubItem("Samples", tabName="tab_samples"),
                                                                 menuSubItem("Features", tabName="tab_features")
                                                                 )),
                                                      menuItemOutput("menu_DEA"),
                                                      menuItemOutput("menu_Enrichments"),
                                                      hr(),
                                                      menuItem("Subset samples", tabName="tab_hm_samples"),
                                                      menuItem("Plot gene", tabName="tab_gene"),
                                                      .modify_stop_propagation(menuItem("Heatmap", startExpanded=TRUE,
                                                                                        menuSubItem("Selection of genes", tabName="tab_hm_genes"),
                                                                                        menuSubItem("Plot heatmap", tabName="tab_heatmap")
                                                      )),
                                                      menuItemOutput("menu_genelist"),
                                                      tags$li(class="shinydashboard-menu-output pkgversion",
                                                              tags$span(paste0("tinySEV v",
                                                                               as.character(packageVersion("tinySEV")))))
                                          )
                         ),
                         dashboardBody(
                           tags$head(tags$style(HTML("
        .sidebar-menu li.treeview, .sidebar-menu li.treeview:hover a{
        	background-color: #2c3b41;
        }
        .navbar-custom-menu {
          width: 90%;
          display: block;
        }
        .navbar-custom-menu .navbar-nav {
          width: 100%;
          display: block;
        }
        li.pkgversion {
          margin-top: 50px;
          margin-left: 15px;
          color: #b8c7ce;
        }
      "))),
                           use_waiter(), useShinyjs(), waiterContent,
                           tabItems(
                             tabItem("tab_object", withSpinner(uiOutput("objOverview"))),
                             tabItem("tab_fileinput",
                                     box(width=6,
                                         tags$p("You may upload your own SummarizedExperiment (SE) object
                     saved as a R .rds file. Once uploaded, it will be added to
                     the list of available objects (in the dropdown list on the
                     top left). For instructions on how to optimally prepare
                     the object, ", actionLink("help_SE", "click here"), "."),
                                         fileInput("file", "Choose SE .rds file", multiple=FALSE,
                                                   accept=c(".rds",".RDS", ".rda"))
                                         ),
                                     box(width=6,
                                         tags$p("Load one or more of the sample Datasets bundled with the ephacRTools package. For a description of the single datasets ", actionLink("help_SE", "click here"), "."),
                                         selectInput("datasets", label= "Choose a pre-bundled Dataset:",
                                                      choices = list("Human Adrenal Glands" = "se_hAG",
                                                                     "Primary Neurons"  = "se_pn",
                                                                     "iPSC-Tricultures" = "se_iN",
                                                                     "ROMK" = "se_romk"), multiple = T),
                                         actionButton("dataset_button", label = "Load pre-bundled datasets")
                                         ),
                                     box(width = 12,
                                         tags$p("You may upload one or more Excel-files generated directly by DataControl and one or more imaging results (.db) from Cluster_Analysis. For guidelines ",
                                                actionLink("help_SE", "click here"), "."),

                                         textInput("se_id", "Name your Dataset:", value = "Custom Dataset"),

                                         fluidRow(
                                           column(
                                             width = 6,
                                             fileInput("fileEphys", "Ephys Excel File (.xlsx)", multiple = TRUE,
                                                       accept = ".xlsx")
                                           ),
                                           column(
                                             width = 6,
                                             fileInput("fileDB", "Imaging Results (.db)", multiple = TRUE,
                                                       accept = ".db")
                                           )
                                         ),

                                         box(
                                           title = "Extra Options",
                                           collapsible = TRUE, collapsed = TRUE, width = 12,
                                           selectInput("tabletype", label= "Choose Imaging Table type:",
                                                       choices = list("Particle Analysis Table" = "pa",
                                                                      "Colocalization Table" = "coloc"),
                                                       multiple = T,selected = "pa"),
                                           uiOutput("optionalControls")
                                         ),

                                         actionButton("mergeSE", label = "Connect Ephys and Imaging Results")
                                     )
                             ),
                             tabItem("tab_samples",
                                     box(width=12, tags$div(style="width: 100%; overflow-x: scroll;",
                                                            withSpinner(DTOutput("samples"))))),
                             tabItem("tab_features",
                                     box(width=12, tags$div(style="width: 100%; overflow-x: scroll;",
                                                            withSpinner(DTOutput("features"))))),
                             tabItem("tab_dea", uiOutput("dea_input"),
                                     tabBox(id="dea_box", width=12,
                                            tabPanel("Overview", textOutput("dea_overview"),
                                                     withSpinner(plotOutput("dea_pvalues", width="400px", height="300px"))),
                                            tabPanel("Table", withSpinner(DTOutput("dea_table")),
                                                     actionButton("dea_geneFilt", "Transfer filtered genes to the heatmap"),
                                                     actionButton("dea_geneSel", "Transfer selected rows to the heatmap"),
                                                     downloadLink('dea_download', label="Download whole table")),
                                            tabPanel("Volcano", withSpinner(plotlyOutput("dea_volcano")),
                                                     tags$p("You may click on a gene to view it."))
                                     )
                             ),
                             tabItem("tab_ea",
                                     box(width=12, id="ea_box", fluidRow(
                                       column(8, uiOutput("ea_input")),
                                       column(4, checkboxInput("ea_showGenes", "Show genes", value=FALSE))
                                     ),
                                     tags$div(style="width: 100%; overflow-x: scroll; font-size: 80%",
                                              withSpinner(DTOutput("ea_table"))),
                                     actionButton("ea_geneSel", "Transfer genes from selected rows to the heatmap"))),
                             tabItem("tab_gene",
                                     column(4, selectizeInput("gene_input", "Select Gene", choices=c(), multiple=FALSE)),
                                     box(width=8, title="Options", collapsible=TRUE, collapsed=TRUE,
                                         selectInput("assay_input", choices=c(), multiple=FALSE,
                                                     tags$span("Assay ", actionLink("help_gassay", "[?]"))),
                                         column(6,
                                                selectInput("plottype_input", "Type of Plot",
                                                            choices=c("violin plot","box plot"), multiple=FALSE),
                                                checkboxInput('select_plotpoints','Plot Points', value=TRUE),
                                                checkboxInput('select_logaxis','Logarithmic Axis', value=FALSE)),
                                         column(6,
                                                selectInput("select_groupvar", choices=c(), multiple=FALSE,
                                                            tags$span("Group by ", actionLink("help_ggroup", "[?]"))),
                                                checkboxInput('asfactor','As factor', value=TRUE),
                                                selectInput("select_colorvar", "Color by", choices=c(), multiple=FALSE),
                                                selectizeInput("select_gridvars", choices=c(),
                                                               options = list(maxItems = 2), multiple=TRUE,
                                                               tags$span("Grid by ", actionLink("help_ggrid", "[?]"))),
                                                checkboxInput('select_freeaxis',value=TRUE,
                                                              tags$span("Free axes ", actionLink("help_gfreeaxes", "[?]"))))
                                     ),
                                     box(width=12, withSpinner(shinyjqui::jqui_resizable(plotOutput("gene_plot")))),
                                     box(width=8, withSpinner(tableOutput("gene_dea_table"))),
                                     column(width=4, shiny::verbatimTextOutput("gene_inList"))
                             ),
                             tabItem("tab_hm_genes", box(width=12, title="Select genes to plot",
                                                         textAreaInput('input_genes','Genes to plot', width="90%", rows=10,
                                                                       placeholder="Enter genes symbols separated by commas, spaces, or line breaks..."),
                                                         tags$p("If your the row names of the object are dot-separated IDs,
                 such as 'ensemblID.symbol' (you can view this in the 'Features' tab),
                 you may also enter just the genes symbols and the corresponding
                 rows will be fetched."),
                                                         tags$p(tags$strong("Important:"), "Note that the number of input genes
                 is capped to ", textOutput("maxGenes", inline=TRUE)))
                             ),
                             tabItem("tab_hm_samples", box(width=12, title="Select samples",
                                                           selectInput("input_hm_samples", multiple=TRUE, selectize=FALSE,
                                                                       label="Select samples to include", choices=c(),
                                                                       size=15)
                             )
                             ),
                             tabItem("tab_heatmap",
                                     box( width=12, title="Heatmap parameters", collapsible=TRUE,
                                          column(4, selectInput("assay_input2", choices=c(), multiple=FALSE,
                                                                tags$span("Assay ", actionLink("help_hmassay", "[?]"))),
                                                 checkboxInput('hm_scale',
                                                               tags$span("Scale rows ", actionLink("help_hmscale", "[?]"))),
                                                 sliderInput("hm_breaks", min=0, max=2, value=0.5, step=0.25,
                                                             tags$span("Colorscale trim ", actionLink("help_hmtrim", "[?]")))
                                          ),
                                          column(4, selectizeInput('hm_anno', "Column annotation", choices=c(), multiple=TRUE),
                                                 selectizeInput('hm_gaps', "Gaps at", choices=c(), multiple=TRUE)
                                          ),
                                          column(4, selectizeInput('hm_order', "Column ordering", choices=c(), multiple=TRUE),
                                                 checkboxInput('hm_clusterCol','Cluster columns', value=FALSE),
                                                 checkboxInput('hm_showColnames','Column names', value=FALSE),
                                                 checkboxInput('hm_clusterRow','Sort rows', value=TRUE)
                                          )
                                     ),
                                     box(width=12, title="Heatmap",
                                         withSpinner(shinyjqui::jqui_resizable(plotOutput("heatmap", height=600))))
                             ),
                             tabItem("tab_genelists",
                                     box(width=12, title="Gene lists",
                                         fluidRow(
                                           column(width=6,
                                                  selectInput("genelist_input", choices=c(), multiple=FALSE,
                                                              tags$span("Select gene list ",
                                                                        actionLink("help_genelists", "[?]")))),
                                           column(width=2,tags$strong("Size: "),textOutput("genelist_size")),
                                           column(width=4, actionButton("btn_importGenelist",
                                                                        label="Use as heatmap input"))
                                         ),
                                         verbatimTextOutput("genelist_out")
                                     )
                             ),
                             tabItem("tab_about", about)
                           ), tags$div(style="clear: both;")
                         )))
}
