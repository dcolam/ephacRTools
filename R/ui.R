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
#' @import shinyjqui
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
#' @import shiny shinydashboard shinyjqui waiter shinyjs
#' @importFrom shinycssloaders withSpinner
#' @importFrom plotly plotlyOutput
#' @importFrom shinyjs useShinyjs
#' @importFrom DT DTOutput
#' @importFrom shinyauthr loginUI
tinySEV.ui <- function(title="tinySEV", waiterContent=NULL, about=NULL,
                       skin="blue", hasLogin=FALSE){
  if(is.null(waiterContent) || isTRUE(waiterContent)){
    waiterContent <- tagList(
      tags$h3("Please wait while the application is initialized..."), waiter::spin_1())
  }
  if(isFALSE(waiterContent)){
    waiterContent <- NULL
  }else{
    waiterContent <- waiter::waiter_show_on_load(html=waiterContent)
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
                                                                 menuSubItem("Column Data", tabName="tab_samples"),
                                                                 menuSubItem("Sweeps", tabName="tab_features")
                                                                 )),
                                                      .modify_stop_propagation(
                                                        menuItem("Customize Object Groups", startExpanded=TRUE,
                                                                 menuSubItem("Define Conditions", tabName="tab_coldata"),
                                                                 menuSubItem("Define Sweeps", tabName="tab_rowdata"),
                                                                 menuSubItem("Change Assays", tabName="tab_assays"),
                                                                 menuSubItem("Filtering", tabName="tab_assays")
                                                        )),
                                                      .modify_stop_propagation(
                                                        menuItem("Plotting", startExpanded=TRUE,
                                                                 menuSubItem("Plate Overview", tabName="tab_plate"),
                                                                 menuSubItem("Plot Sweeps", tabName="tab_sweeps"),
                                                                 menuSubItem("Image Plate Viewer", tabName="tab_imgplate")
                                                        )),


                                                      hr(),
                                                      menuItem("Manual Annotation", tabName = "tab_annotation",
                                                               icon = icon("tags")),
                                                      menuItem("Clustering", tabName = "tab_cluster"),
                                                      menuItem("Export", tabName="tab_export"),
                                                      tags$li(class="shinydashboard-menu-output pkgversion",
                                                              tags$span(paste0("ephacRTools v",
                                                                               as.character(packageVersion("ephacRTools"))))),
                                                      tags$li(
                                                        class = "shinydashboard-menu-output memory-monitor",
                                                        style = "padding-left: 15px; padding-top: 5px; font-size: 13px; color: #888;",
                                                        fluidRow(
                                                          column(width = 6, style = "padding-left: 15px; padding-top: 5px; font-size: 13px; color: #888;",
                                                                 icon("microchip"),
                                                                 textOutput("memoryUsage", inline = TRUE)),
                                                          column(width = 3, actionButton("refreshMem", label = NULL, icon = icon("refresh"),
                                                                                         style = "padding: 2px 6px; font-size: 6px;"))
                                                        )
                                                      )
                                          )
                         ),
                         dashboardBody(
                           useShinyjs(),

                           # absolutePanel(id = "plate_viewer", class = "panel panel-default",
                           #               fixed = TRUE,
                           #               draggable = TRUE,
                           #               top = "auto", left = "auto", right = 20, bottom = 20,
                           #               width = 150, height = 100,
                           #               style = "z-index: 1000; overflow: hidden;",
                           #
                           #               div(style = "display: flex; justify-content: space-between; align-items: center;",
                           #
                           #                   actionButton("toggle_plate", label = NULL, icon = icon("window-minimize"), style = "padding: 2px 6px; font-size: 12px;")
                           #               ),
                           #               plotOutput("mini_plate_plot", height = "250px")
                           # ),
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
                                     fluidRow(
                                     box(width=6,title = "Load RDS",
                                         tags$p("You may upload your own SummarizedExperiment (SE) object
                     saved as a R .rds file. Once uploaded, it will be added to
                     the list of available objects (in the dropdown list on the
                     top left). For instructions on how to optimally prepare
                     the object, ", actionLink("help_SE", "click here"), "."),
                                         fileInput("file", "Choose SE .rds file", multiple=FALSE,
                                                   accept=c(".rds",".RDS", ".rda"))
                                         ),
                                     box(width=6,title = "Load Prebundlet Dataset",
                                         tags$p("Load one or more of the sample Datasets bundled with the ephacRTools package. For a description of the single datasets ", actionLink("help_SE", "click here"), "."),
                                         selectInput("datasets", label= "Choose a pre-bundled Dataset:",
                                                      choices = list("Human Adrenal Glands" = "se_hAG",
                                                                     "Primary Neurons"  = "se_pn",
                                                                     "iPSC-Tricultures" = "se_iN",
                                                                     "ROMK" = "se_romk"), multiple = T),
                                         actionButton("dataset_button", label = "Load pre-bundled datasets")
                                         )
                                     ),


                                      fluidRow(
                                     box(
                                       width = 6,title = "Load Excel-File",
                                       tags$p("You may upload one or more Excel files generated directly by DataControl. Make sure to follow the guidelines."),


                                       fluidRow(
                                         column(
                                           width = 6,
                                           textInput("se_id", "Name your Dataset:", value = "Custom Dataset")
                                         ),
                                         column(
                                           width = 6,
                                           fileInput("fileEphys", "Ephys Excel File (.xlsx)", multiple = TRUE, accept = ".xlsx")
                                         )
                                       ),


                                       tableOutput("importXl")
                                       # actionButton("loadEphys", label = "Load Excel into SE")
                                     ),
                                     box(
                                       width = 6,title = "Load Imaging Results",
                                       tags$p("You may upload imaging result `.db` files and link them to a selected SummarizedExperiment (SE) dataset. Make sure the data types match."),

                                       # Row with selectizeInput and fileInput side-by-side
                                       fluidRow(
                                         column(
                                           width = 6,
                                           selectizeInput(
                                             "seDataset",
                                             label = "Select SE",
                                             choices = c(),
                                             multiple = FALSE
                                           )
                                         ),
                                         column(
                                           width = 6,
                                           fileInput(
                                             "fileDB",
                                             label = "Imaging Results (.db)",
                                             multiple = TRUE,
                                             accept = ".db"
                                           )
                                         )
                                       ),

                                       # Table type selector
                                       fluidRow(
                                         column(
                                           width = 12,
                                           selectInput(
                                             "tabletype",
                                             label = "Choose Imaging Table type:",
                                             choices = list(
                                               "Particle Analysis Table" = "pa",
                                               "Colocalization Table" = "coloc"
                                             ),
                                             multiple = TRUE,
                                             selected = "pa"
                                           ),
                                           tableOutput("importDB")

                                         )
                                       ),

                                       # Optional UI controls
                                       fluidRow(
                                         column(
                                           width = 12,
                                           uiOutput("optionalControls")
                                         )
                                       ),

                                       # Action button
                                       fluidRow(
                                         column(
                                           width = 12,
                                           actionButton("mergeSE", label = "Connect Ephys and Imaging Results")
                                         )
                                       )
                                     )
                                      )


                             ),
                             tabItem("tab_samples",
                                     box(width=12, tags$div(style="width: 100%; overflow-x: scroll;",
                                                            withSpinner(DTOutput("samples"))))),
                             tabItem("tab_features",
                                     box(width=12, tags$div(style="width: 100%; overflow-x: scroll;",
                                                            withSpinner(DTOutput("features"))))),
                             tabItem("tab_plate",
                                     box(width = 12,

                                         # Row for plot
                                         fluidRow(
                                           column(width = 12,
                                                  withSpinner(plotlyOutput("plate_view"))
                                           )
                                         ),

                                         # Row for main input controls
                                         fluidRow(
                                           column(width = 4,
                                                  selectInput("plate_id", label = "Plate ID", choices = c()),

                                                  # NEW: Well selector
                                                  selectizeInput("selected_well", "Selected Well", choices = c(),
                                                              multiple = TRUE),

                                                  # Add spacing between selectInput and button
                                                  div(style = "margin-top: 10px;"),
                                                  actionButton("reset_well", "Reset Well Selection"),
                                                  div(style = "margin-top: 50px;")
                                           ),
                                           column(width = 4,
                                                  selectInput("assay_id", label = "Assay", choices = c()),
                                                  radioButtons("assay_option", "Assay Display Mode:",
                                                               choices = c("Raw" = "raw", "Log10" = "log10", "Z-score" = "scale"),
                                                               inline = TRUE)
                                           ),
                                           column(width = 4,
                                                  selectInput("sweep_id", label = "Sweep", choices = c()),

                                                  # Group & Aggregate Box
                                                  box(
                                                    title = "Group & Aggregate Options",
                                                    width = 12,
                                                    collapsible = TRUE,
                                                    collapsed = TRUE,
                                                    selectInput("sweep_group", "Group Sweeps (aggregate):",
                                                                choices = c(), multiple = TRUE),
                                                    radioButtons("agg_method", "Aggregation method:",
                                                                 choices = c("Mean" = "mean", "Median" = "median", "Sum" = "sum"),
                                                                 inline = TRUE),
                                                    div(style = "text-align: center; font-weight: bold; margin: 10px 0;", "— or —"),
                                                    selectInput("group_by_meta", "Group by metadata column:", choices = c())
                                                  )
                                           )
                                         ),

                                         # Row for custom grouping box — separated from main inputs
                                         fluidRow(
                                           box(
                                             width = 12,
                                             title = "Create Custom Groupings:",
                                             collapsible = TRUE,
                                             collapsed = TRUE,
                                             column(width = 4,
                                                    sliderInput("selected_slider", "Select a Range:", min = 0, max = 1, value = c(0, 1))
                                             )
                                           )
                                         )

                                     )
                             ),
                             tabItem("tab_imgplate",
                                     fluidRow(
                                       column(width = 3,
                                         box(width = 12, title = "Image Folder",
                                           localImgBrowserUI("imgbrowser"),
                                           helpText(style = "font-size:11px;",
                                             "Works locally and when deployed. ",
                                             "Filenames must follow: ",
                                             tags$code("..._Well-site_Channel_Class_crop.jpg")),
                                           hr(),
                                           tags$small(style = "color:#aaa;",
                                             "Local / server only:"),
                                           shinyFiles::shinyDirButton(
                                             "img_dir_btn",
                                             label       = "Browse server folder…",
                                             title       = "Select the folder containing plate images",
                                             class       = "btn-default btn-block btn-sm"
                                           ),
                                           div(style = "font-size: 11px; font-family: monospace; color: #666; margin-top: 4px; word-break: break-all;",
                                               textOutput("img_folder_display")),
                                           hr(),
                                           selectInput("plate_img_id", "Matched Plate ID", choices = c()),
                                           selectInput("img_channel",   "Channel:",         choices = c()),
                                           selectInput("img_class_sel", "Class:",            choices = c()),
                                           hr(),
                                           selectInput("img_coldata_var", "Color overlay column:",
                                                       choices = "None", selected = "None"),
                                           selectizeInput("img_hover_vars", "Extra hover info:",
                                                          choices  = c(),
                                                          multiple = TRUE,
                                                          options  = list(placeholder = "None (add columns…)")),
                                           sliderInput("img_overlay_alpha", "Overlay opacity:",
                                                       min = 0, max = 1, value = 0.4, step = 0.05),
                                           uiOutput("img_filter_ui"),
                                           uiOutput("img_legend_ui"),
                                           hr(),
                                           verbatimTextOutput("img_plate_stats")
                                         )
                                       ),
                                       column(width = 9,
                                         box(width = 12, style = "padding: 0; overflow: hidden;",
                                           div(style = paste0(
                                             "padding: 5px 12px; background: #f8f8f8;",
                                             " border-bottom: 1px solid #e0e0e0;",
                                             " font-family: monospace; font-size: 13px;",
                                             " min-height: 28px; color: #333;"),
                                             textOutput("img_hover_info", inline = TRUE)
                                           ),
                                           div(style = "overflow:auto;",
                                             uiOutput("img_plate_ui"))
                                         )
                                       )
                                     )
                             ),
                             tabItem("tab_sweeps",
                                     box(width = 12,

                                         # Row for plot
                                         fluidRow(
                                           column(width = 12,
                                                  withSpinner(plotlyOutput("sweep_view"))
                                           )
                                         ),
                                         # Row for input controls
                                         fluidRow(
                                           column(width = 4,
                                                  selectInput("plate_id1", label = "Plate ID", choices = c()),

                                                  checkboxInput("all_plates", label = "Show All Plates"),

                                                  # NEW: Well selector that updates on click
                                                  selectizeInput("selected_well1", "Selected Well",
                                                              choices = c(),
                                                              multiple = T),

                                                  # NEW: Reset button
                                                  actionButton("reset_well", "Reset Well Selection")
                                           ),
                                           column(width = 4,
                                                  selectInput("assay_id1", label = "Assay", choices = c(), multiple=T),
                                                  radioButtons("assay_option1", "Assay Display Mode:",
                                                               choices = c("Raw" = "raw", "Log10" = "log10", "Z-score" = "scale"),
                                                               inline = TRUE),
                                                  textInput("x.labels", label="X-Axis Labeling", value = "Sweep"),
                                                  textInput("y.labels", label="Y-Axis Labeling", value = "Assay Value")
                                           ),
                                           column(width = 4,
                                                  # Primary sweep selection
                                                  selectInput("group_by_meta1", "Choose X-Axis Variable",
                                                              choices = c()),
                                                  selectInput("color_group1", "Color by:",
                                                              choices = c(), multiple = F),
                                                  selectInput("facet_group1", "Facet Grouping:",
                                                              choices = c(), multiple = F),

                                                  radioButtons("facets", label="Facetting",
                                                               choices=c("Grid" = "grid", "Wrap" = "wrap", "No Facet" = "nothing"),
                                                               inline=T, selected = "grid"),
                                                  checkboxInput("invertFacet", label = "Invert Facetting:", value = F)


                                           )
                                         )

                                     )
                             ),



                             tabItem("tab_coldata",


                                              box(width = 6,
                                                  fluidRow(
                                                  column(width=4,
                                                  selectInput("plate_id4", label = "Plate ID", choices = c())
                                                  ),
                                                  column(width=4,
                                                  selectInput("condition", label = "Highlight Condition:", choices = c(),
                                                              selected=FALSE)
                                                  ),
                                                  column(width=4,
                                                  checkboxInput("all_plates2", label = "Same for All Plates", value = T))
                                                  ),
                                                  withSpinner(plotlyOutput("plate_view_col", height = "400px")),
                                                  fluidRow(
                                                    column(width=6,
                                                           textInput("newCondition", label = "Name the New Condition Group:", value= "Ex. Genotype"),
                                                           actionButton("createCondition", label="Submit")
                                                    ),
                                                    column(width=6,
                                                           selectInput("subGroup", label = "Choose the Wells:", choices = c()),
                                                           textInput("newCondition", label = "Name the selected Group:", value = "Ex. Control"),
                                                           actionButton("createGroup", label="Submit")

                                                    )
                                                  )




                                              ),


                                     box(width = 6,
                                         withSpinner(DTOutput("features_col"))

                                         )


                             ),


                             tabItem("tab_annotation",
                                     fluidRow(
                                       column(width = 9,
                                         box(width = 12, title = "Current Image (BF  |  Fluoro)",
                                           uiOutput("ann_image_ui"),
                                           div(style = "margin-top: 12px;",
                                               actionButton("ann_skip", "Skip / Next", class = "btn-default"),
                                               actionButton("ann_undo", "Undo last",   class = "btn-warning"),
                                               span(style = "margin-left: 12px; color: #888; font-size: 12px;",
                                                    textOutput("ann_progress", inline = TRUE)))
                                         )
                                       ),
                                       column(width = 3,
                                         box(width = 12, title = "Annotation",
                                           helpText("Select a folder in the Image Plate Viewer tab. Note: auto-save CSV requires the server folder picker; client-side folders use Download CSV."),
                                           selectInput("ann_img_class", "Class to annotate:",
                                                       choices = c(), width = "100%"),
                                           selectizeInput("ann_plate_ids", "Plates to include:",
                                                          choices  = c(),
                                                          multiple = TRUE,
                                                          options  = list(placeholder = "All plates")),
                                           textInput("ann_person", "Your name:",
                                                     placeholder = "Annotator name"),
                                           hr(),
                                           textInput("ann_new_class", NULL,
                                                     placeholder = "New label…"),
                                           actionButton("ann_add_class",     "Add label",    class = "btn-primary btn-sm"),
                                           actionButton("ann_clear_classes", "Clear labels", class = "btn-danger btn-sm"),
                                           hr(),
                                           uiOutput("ann_class_buttons"),
                                           hr(),
                                           tags$h5("Results"),
                                           downloadButton("ann_download_csv", "Download CSV"),
                                           actionButton("ann_clear_results", "Clear results", class = "btn-danger btn-sm"),
                                           div(style = "margin-top: 10px;",
                                               tableOutput("ann_results_preview"))
                                         )
                                       )
                                     )
                             ),
                             tabItem("tab_cluster",

                                     box(width = 12, title = "Set up Clustering:", collapsible = T, collapsed = T,

                                        fluidRow(

                                          column( width=4,

                                                 selectInput("clusterAssay", label = "Use Assay:", choices = c(), multiple=T),
                                                 textInput("suffix_cluster", label = "Add Suffix:", value=".1"),
                                                 actionButton("clustering", label="Perform Clustering")
                                          ),
                                          column( width=4,
                                                 selectInput("clusterColData", label = "Use Column Data:", choices = c(), multiple=T)
                                                 ),
                                          column( width=4,
                                                 sliderInput("k_cluster", label= "Number of Clusters", min=1, max=20, value=3, step = 1),

                                                 )
                                        )
                                         ),

                                     tabBox(
                                       width = 12,
                                       title = "Cluster Plots",
                                       tabPanel("t-SNE",
                                                withSpinner(uiOutput("cluster_tsne_ui")),
                                                selectInput("clustercolor1", "Color your points:", choices = c(), multiple = TRUE)
                                       ),
                                       tabPanel("UMAP",
                                                withSpinner(uiOutput("cluster_umap_ui")),
                                                selectInput("clustercolor2", "Color your points:", choices = c(),  multiple=T)),
                                       tabPanel("PCA",
                                                withSpinner(uiOutput("cluster_pca_ui")),
                                                selectInput("clustercolor3", "Color your points:", choices = c(),  multiple=T))
                                     ),

                                      box(width=12,

                                          fluidRow(
                                            column(width=4,
                                                   selectInput("plate_id3", "Select Plate:",
                                                               choices = c()),
                                                selectizeInput("clusteredwells", "Selected Wells:",
                                                               choices = c(),
                                                               multiple =T),
                                                # Add spacing between selectInput and button
                                                div(style = "margin-top: 10px;"),
                                                actionButton("reset_well", "Reset Well Selection"),
                                                div(style = "margin-top: 50px;")

                                                   )
                                          )

                                          )



                             ),

                          tabItem("tab_export",


                                  box("Export as RDS-Object", width=6,
                                      column(width=12,
                                      fluidRow(



                                      selectInput("rdsObject", choices = c(), label = "Choose the Dataset"),
                                      downloadButton("downloadRDS", "Download")
                                        )
                                      )
                                      )

                                  ),




                             tabItem("tab_about", about)
                           ), tags$div(style="clear: both;")
                         )))
}
