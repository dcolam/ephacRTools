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
                                                        menuItem("Plotting", startExpanded=TRUE,
                                                                 menuSubItem("Plate Overview", tabName="tab_plate"),
                                                                 menuSubItem("Plot Sweeps", tabName="tab_sweeps"),
                                                                 menuSubItem("Show Images", tabName="tab_images")
                                                        )),

                                                      hr(),
                                                      menuItem("Clustering", tabName = "tab_cluster"),
                                                      menuItem("Define Groupings", tabName="tab_groupings"),
                                                      menuItem("Export", tabName="tab_export"),
                                                      tags$li(class="shinydashboard-menu-output pkgversion",
                                                              tags$span(paste0("ephacRTools v",
                                                                               as.character(packageVersion("ephacRTools")))))
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
                                                       accept = ".xlsx"),
                                             actionButton("loadEphys", label = "Load Excel into SE")
                                           ),
                                           column(
                                             width = 6,
                                             fileInput("fileDB", "Imaging Results (.db)", multiple = TRUE,
                                                       accept = ".db"),
                                             actionButton("loadDB", label = "Connect Ephys and Imaging Results")
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
                             tabItem("tab_images",
                                     box(width = 12,

                                         # Row for plot
                                         fluidRow(
                                           column(width = 12,
                                                  withSpinner(plotlyOutput("sweep_view1"))
                                           )
                                         ),
                                         fluidRow(
                                           column(width = 4,
                                                  selectInput("plate_id2", label = "Plate ID", choices = c()),

                                                  # NEW: Well selector that updates on click
                                                  selectizeInput("selected_well3", "Selected Well",
                                                                 choices = c()),

                                                  # NEW: Reset button
                                                  actionButton("reset_well", "Reset Well Selection")
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

                                                  # NEW: Well selector that updates on click
                                                  selectizeInput("selected_well1", "Selected Well",
                                                              choices = c(),
                                                              multiple = T),

                                                  # NEW: Reset button
                                                  actionButton("reset_well", "Reset Well Selection")
                                           ),
                                           column(width = 4,
                                                  selectInput("assay_id1", label = "Assay", choices = c()),
                                                  radioButtons("assay_option1", "Assay Display Mode:",
                                                               choices = c("Raw" = "raw", "Log10" = "log10", "Z-score" = "scale"),
                                                               inline = TRUE)
                                           ),
                                           column(width = 4,
                                                  # Primary sweep selection
                                                  selectInput("group_by_meta1", "Group by metadata column:",
                                                              choices = c()),
                                                  selectInput("color_group1", "Color by:",
                                                              choices = c(), multiple = TRUE)

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



                             tabItem("tab_about", about)
                           ), tags$div(style="clear: both;")
                         )))
}
