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
                                                        menuItem("Import", startExpanded=TRUE,
                                                                 menuSubItem("Overview",           tabName="tab_object"),
                                                                 menuSubItem("Load SE",            tabName="tab_load_se"),
                                                                 menuSubItem("Voltage Clamp → SE", tabName="tab_vc_import"),
                                                                 menuSubItem("Combine → MAE", tabName="tab_mae")
                                                                 )),
                                                      .modify_stop_propagation(
                                                        menuItem("Current Clamp → SE",
                                                                 icon = icon("bolt"), startExpanded=TRUE,
                                                                 menuSubItem("1 — CSV → Parquet",    tabName="tab_cc_csv"),
                                                                 menuSubItem("2 — Detection & Preview",   tabName="tab_cc_detect"),
                                                                 menuSubItem("3 — Build SE",             tabName="tab_cc_build")
                                                        )),
                                                      .modify_stop_propagation(
                                                        menuItem("Inspect & Edit", startExpanded=TRUE,
                                                                 menuSubItem("Column Data",       tabName="tab_samples"),
                                                                 menuSubItem("Sweep Data",        tabName="tab_features"),
                                                                 menuSubItem("Define Conditions", tabName="tab_coldata"),
                                                                 menuSubItem("Define Sweeps",     tabName="tab_rowdata"),
                                                                 menuSubItem("Change Assays",     tabName="tab_assays"),
                                                                 menuSubItem("Filter Wells",      tabName="tab_filter_wells")
                                                        )),
                                                      .modify_stop_propagation(
                                                        menuItem("Visualize", startExpanded=TRUE,
                                                                 menuSubItem("Plate Overview", tabName="tab_plate"),
                                                                 menuSubItem("Plot Sweeps",    tabName="tab_sweeps")
                                                        )),
                                                      .modify_stop_propagation(
                                                        menuItem("Image Analysis", startExpanded=TRUE,
                                                                 menuSubItem("Import Data",         tabName="tab_img_import"),
                                                                 menuSubItem("Auto Classification", tabName="tab_img_classify"),
                                                                 menuSubItem("Manual Scoring",      tabName="tab_annotation"),
                                                                 menuSubItem("Plate Viewer",        tabName="tab_imgplate")
                                                        )),

                                                      hr(),
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
                           tags$head(tags$script(HTML("
                             var _tinyWarnClose = false;

                             // RStudio's viewer pane runs inside Qt WebEngine.
                             // It fires beforeunload but cannot show the native
                             // confirmation dialog, so e.preventDefault() blocks
                             // closing permanently. Skip the warning there.
                             var _inRStudioViewer = /QtWebEngine/i.test(navigator.userAgent);

                             window.addEventListener('beforeunload', function(e) {
                               if (_tinyWarnClose && !_inRStudioViewer) {
                                 e.preventDefault();
                                 e.returnValue = '';
                               }
                             });

                             Shiny.addCustomMessageHandler('setWarnBeforeClose', function(val) {
                               _tinyWarnClose = val;
                             });
                           "))),

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
                             tabItem("tab_load_se",
                                     fluidRow(
                                       box(width=6, title = "Load SE from File",
                                         tags$p("Upload a SummarizedExperiment (SE) object saved as an .rds file.
                                                 Once uploaded it will appear in the dataset selector at the top.
                                                 For preparation instructions, ", actionLink("help_SE", "click here"), "."),
                                         fileInput("file", "Choose SE .rds file", multiple=FALSE,
                                                   accept=c(".rds",".RDS",".rda"))
                                       ),
                                       box(width=6, title = "Load Pre-bundled Dataset",
                                         tags$p("Load one of the sample datasets bundled with ephacRTools.
                                                 For dataset descriptions, ", actionLink("help_SE2", "click here"), "."),
                                         selectInput("datasets", label = "Choose a pre-bundled dataset:",
                                                     choices = list("Human Adrenal Glands" = "se_hAG",
                                                                    "Primary Neurons"       = "se_pn",
                                                                    "iPSC-Tricultures"      = "se_iN",
                                                                    "ROMK"                  = "se_romk"),
                                                     multiple = TRUE),
                                         actionButton("dataset_button", label = "Load pre-bundled datasets")
                                       )
                                     )
                             ),
                             tabItem("tab_vc_import",
                                     fluidRow(
                                       box(width=8, title = "Voltage Clamp → SE (DataControl Excel)",
                                         tags$p("Upload one or more DataControl Excel files to build a
                                                 SingleCellExperiment. Make sure to follow the formatting guidelines."),
                                         fluidRow(
                                           column(width=6,
                                                  textInput("se_id", "Name your dataset:", value = "Custom Dataset")),
                                           column(width=6,
                                                  fileInput("fileEphys", "DataControl Excel (.xlsx)",
                                                            multiple = TRUE, accept = ".xlsx"))
                                         ),
                                         tableOutput("importXl")
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
                                         box(width = 12, title = "View Controls",
                                           helpText(style="font-size:11px; color:#888;",
                                             "Select image folder in ", tags$b("Image Analysis \u2192 Import Data"), "."),
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
                                           uiOutput("img_plate_css"),
                                           div(style = "overflow:auto;",
                                             uiOutput("img_plate_ui")),
                                           div(style = "padding: 8px 16px 4px;",
                                             fluidRow(
                                               column(6, sliderInput("img_brightness", "Brightness:",
                                                         min = 0.2, max = 4, value = 1, step = 0.1,
                                                         width = "100%")),
                                               column(6, sliderInput("img_contrast", "Contrast:",
                                                         min = 0.2, max = 4, value = 1, step = 0.1,
                                                         width = "100%"))
                                             ),
                                             actionButton("img_reset_filter", icon("rotate-left"),
                                                          class = "btn-xs btn-default",
                                                          style = "margin-top:2px;")
                                           )
                                         )
                                       )
                                     )
                             ),
                             tabItem("tab_img_import",
                               fluidRow(
                                 box(width = 7, title = "Load Imaging Database (.db)",
                                   tags$p(style="color:#666;",
                                     "Upload Cluster Analysis SQLite databases (.db), then click ",
                                     tags$strong("Scan"), " to inspect the database structure ",
                                     "and configure the import below."),
                                   fluidRow(
                                     column(7,
                                       fileInput("img_fileDB", "Imaging Results (.db)",
                                                 multiple=TRUE, accept=".db")
                                     ),
                                     column(5,
                                       selectInput("img_tabletype", "Table type:",
                                         choices  = list("Particle Analysis"="pa",
                                                         "Colocalization"="coloc"),
                                         multiple = FALSE, selected="pa"),
                                       br(),
                                       actionButton("img_scan_btn", "Scan Databases",
                                                    class="btn-info btn-block",
                                                    icon=icon("magnifying-glass"))
                                     )
                                   ),
                                   tableOutput("img_importDB_preview"),
                                   hr(),
                                   uiOutput("img_config_ui")
                                 ),
                                 box(width = 5, title = "Image Folder (JPG Thumbnails)",
                                   tags$p("Select the folder containing cropped JPG thumbnails exported by Cluster Analysis.",
                                          "Filenames must follow: ", tags$code("..._Well-site_Channel_Class_crop.jpg")),
                                   localImgBrowserUI("imgbrowser"),
                                   helpText(style="font-size:11px; color:#aaa;",
                                     "Works locally and when deployed (client-side). ",
                                     "For server/local use, also browse by folder:"),
                                   shinyFiles::shinyDirButton(
                                     "img_dir_btn", label="Browse server folder\u2026",
                                     title="Select the folder containing plate images",
                                     class="btn-default btn-block btn-sm"),
                                   div(style="font-size:11px;font-family:monospace;color:#666;margin-top:4px;word-break:break-all;",
                                       textOutput("img_folder_display"))
                                 )
                               ),
                               fluidRow(
                                 box(width = 5, title = "Load Pre-prepared Particle Table (.rds)",
                                   status = "warning", solidHeader = FALSE,
                                   tags$p(style="color:#666;",
                                     "Load a previously exported particle data frame (.rds) directly ",
                                     "into the classification pipeline, skipping the database import."),
                                   fileInput("img_rds_file", "Particle table (.rds)", accept=".rds"),
                                   actionButton("img_load_rds", "Load", class="btn-warning btn-block"),
                                   verbatimTextOutput("img_rds_status")
                                 )
                               )
                             ),
                             tabItem("tab_img_classify",
                               fluidRow(
                                 box(width=12, title="Automatic Cell Classification from Particle Data",
                                     tags$p(style="color:#666;",
                                            "Run the pipeline step by step. Each step uses the output of the previous. ",
                                            "Diagnostic plots update after each step."))
                               ),

                               # ── Step 1: Load ──────────────────────────────────────────────────────────
                               fluidRow(
                                 box(width=4, title="1 \u2014 Load", collapsible=TRUE,
                                   tags$p(style="color:#666; font-size:11px;",
                                          "Load particles from the Import Data tab. SE metadata will be joined ",
                                          "immediately and is available for faceting in all plots. ",
                                          tags$em("Reload always resets all downstream steps.")),
                                   selectizeInput("cls_seDataset", "SE (for step 6 & faceting):", choices=c(), multiple=FALSE),
                                   hr(),
                                   selectInput("cls_scale_mode", "Scaling mode (for filter step):",
                                     choices=c(
                                       "Uncentered \u00f7 median (recommended)"  = "uncentered",
                                       "Centered z-score"                         = "centered"),
                                     selected="uncentered"),
                                   helpText(style="font-size:10px;",
                                            "Uncentered: Mean \u00f7 median(Mean) per group. ",
                                            "Centered: (Mean \u2212 \u03bc) \u00f7 \u03c3 per group. ",
                                            "Affects the preview and links to the recommended filter method."),
                                   hr(),
                                   selectInput("cls_scale_groups", "Scale within:",
                                     choices=c(
                                       "Channel \u00d7 Plate"  = "channel_plate",
                                       "Channel only"          = "channel",
                                       "Channel \u00d7 Well"   = "channel_well"),
                                     selected="channel_plate"),
                                   hr(),
                                   uiOutput("cls_meta_col_ui"),
                                   hr(),
                                   actionButton("cls_load", "Check / Reload", class="btn-primary btn-block"),
                                   hr(),
                                   verbatimTextOutput("cls_load_status")
                                 ),
                                 box(width=8, title="Raw Particles: Distribution", collapsible=TRUE,
                                   fluidRow(
                                     column(4,
                                       selectInput("cls_load_metric", "Metric:",
                                         choices=c("Mean", "Area"), selected="Mean")),
                                     column(4,
                                       radioButtons("cls_load_scale", "View:",
                                         choices=c("Raw"="raw", "Scaled (÷SD)"="scaled"),
                                         selected="raw", inline=TRUE))
                                   ),
                                   tabsetPanel(
                                     tabPanel("ECDF",
                                       br(),
                                       helpText(style="font-size:10px; color:#888;",
                                                "Scaled view uses scale(center=F,scale=T) — same scaling as the filter step. ",
                                                "Dashed line = current Step 2 threshold (left = would be rejected)."),
                                       withSpinner(plotlyOutput("cls_diag_load_ecdf",  height="280px"))),
                                     tabPanel("Particles / well",
                                       br(),
                                       withSpinner(plotlyOutput("cls_diag_load_hist",  height="280px"))),
                                     tabPanel("Scaled preview (Mean only)",
                                       br(),
                                       helpText(style="font-size:10px; color:#888;",
                                                "Scaled Mean with dashed line at the current Step 2 threshold. ",
                                                "Particles LEFT of the line would be rejected."),
                                       withSpinner(plotlyOutput("cls_diag_scaled_preview", height="260px")))
                                   )
                                 )
                               ),

                               # ── Step 2: Filter ────────────────────────────────────────────────────────
                               fluidRow(
                                 box(width=4, title="2 \u2014 Filter", collapsible=TRUE,
                                   tags$p(style="color:#666; font-size:11px;",
                                          "Optional — aggregation can skip directly to Step 3 using raw particles. ",
                                          "The method is pre-selected based on your scaling mode."),
                                   uiOutput("cls_filter_method_ui"),
                                   uiOutput("cls_filter_threshold_ui"),
                                   helpText(style="font-size:10px;",
                                            "Leave at default (NA) to use the method default. ",
                                            "z-score default = 0, uncentered SD default = median(scaled)/3 per group. ",
                                            "Check the Threshold preview tab for a live view of the cut."),
                                   actionButton("cls_filter", "Filter", class="btn-primary btn-block"),
                                   hr(),
                                   verbatimTextOutput("cls_filter_status")
                                 ),
                                 box(width=8, title="Filter Diagnostics", collapsible=TRUE,
                                   tabsetPanel(
                                     tabPanel("Threshold preview",
                                       br(),
                                       helpText(style="font-size:10px; color:#888;",
                                                "Scaled distribution with dashed line at current threshold. ",
                                                "Change the threshold input and this updates live."),
                                       withSpinner(plotlyOutput("cls_diag_filter_threshold", height="270px"))),
                                     tabPanel("Before vs After (raw)",
                                       br(),
                                       withSpinner(plotlyOutput("cls_diag_filter_ecdf",      height="280px"))),
                                     tabPanel("Retention %",
                                       br(),
                                       withSpinner(plotlyOutput("cls_diag_filter_retain",    height="280px")))
                                   )
                                 )
                               ),

                               # ── Step 3: Aggregate ─────────────────────────────────────────────────────
                               fluidRow(
                                 box(width=4, title="3 \u2014 Aggregate + Scale", collapsible=TRUE,
                                   tags$p(style="color:#666; font-size:11px;",
                                          "Step 3a: summarise to one row per Channel \u00d7 Plate \u00d7 Well. ",
                                          "Step 3b: scale the aggregated metrics within groups (for scoring). ",
                                          "Uses filtered particles if available, otherwise raw."),
                                   selectInput("cls_agg_fun", "Aggregation function:",
                                     choices=c(
                                       "Mean (recommended)" = "mean",
                                       "Median"             = "median",
                                       "Sum"                = "sum"),
                                     selected="mean"),
                                   helpText(style="font-size:10px;",
                                            "Applied to Mean and Area. normArea always = sum(Area)/mean(Selection_Area)."),
                                   hr(),
                                   selectInput("cls_agg_scale_within", "Scale within (Step 3b):",
                                     choices=c(
                                       "Channel only (recommended)" = "channel",
                                       "Channel \u00d7 Plate"       = "channel_plate",
                                       "None (skip scaling)"        = "none"),
                                     selected="channel"),
                                   checkboxInput("cls_agg_scale_center",
                                     "Center when scaling (z-score)?", value=FALSE),
                                   helpText(style="font-size:10px;",
                                            "Uncentered (\u00f7SD) is recommended: keeps negative wells near 0."),
                                   actionButton("cls_aggregate_btn", "Aggregate", class="btn-primary btn-block"),
                                   hr(),
                                   verbatimTextOutput("cls_agg_status")
                                 ),
                                 box(width=8, title="Aggregated: Per-channel Distribution", collapsible=TRUE,
                                   fluidRow(
                                     column(5,
                                       selectInput("cls_agg_metric", "Metric:",
                                         choices=c("Mean_agg", "Area_agg", "normArea", "n_particles",
                                                   "Mean_z",   "Area_z",   "normArea_z"),
                                         selected="Mean_agg")),
                                     column(4,
                                       radioButtons("cls_agg_view", "Show:",
                                         choices=c("Raw"="raw", "Scaled"="scaled"),
                                         selected="raw", inline=TRUE))
                                   ),
                                   withSpinner(plotlyOutput("cls_diag_agg_box", height="290px"))
                                 )
                               ),

                               # ── Step 4: Score ─────────────────────────────────────────────────────────
                               fluidRow(
                                 box(width=4, title="4 \u2014 Score", collapsible=TRUE,
                                   tags$p(style="color:#666; font-size:11px;",
                                          "Scales each metric within Channel \u00d7 Plate, then combines with weights."),
                                   numericInput("cls_w_mean",    "Weight: Mean",     value=1, min=0, step=0.5),
                                   numericInput("cls_w_area",    "Weight: Area",     value=1, min=0, step=0.5),
                                   numericInput("cls_w_normarea","Weight: normArea", value=1, min=0, step=0.5),
                                   hr(),
                                   checkboxInput("cls_score_center", "Center scores (z-score)?", value=FALSE),
                                   conditionalPanel("input.cls_score_center == true",
                                     tags$div(
                                       style="background:#fff3cd; border:1px solid #ffc107; padding:7px; border-radius:4px; font-size:10px; margin-bottom:6px;",
                                       tags$strong("\u26a0 Warning:"),
                                       " Centered z-score scoring can assign near-average scores to empty/negative ",
                                       "wells (those zeroed out by aggregation). Only use if you know what you are doing. ",
                                       tags$strong("Uncentered (default) is recommended.")
                                     )
                                   ),
                                   actionButton("cls_score", "Score", class="btn-primary btn-block"),
                                   hr(),
                                   verbatimTextOutput("cls_score_status")
                                 ),
                                 box(width=8, title="Score Diagnostics", collapsible=TRUE,
                                   tabsetPanel(
                                     tabPanel("Score ECDF",
                                       br(),
                                       withSpinner(plotlyOutput("cls_diag_score_ecdf", height="300px"))),
                                     tabPanel("Metric ECDFs",
                                       br(),
                                       fluidRow(column(5,
                                         selectInput("cls_score_metric", "Metric:",
                                           choices=c("channel_score",
                                                     "Mean_agg", "Area_agg", "normArea",
                                                     "Mean_z",   "Area_z",   "normArea_z"),
                                           selected="channel_score"))),
                                       withSpinner(plotlyOutput("cls_diag_metric_ecdf", height="270px"))),
                                     tabPanel("Scatter",
                                       br(),
                                       withSpinner(plotlyOutput("cls_diag_score_scatter", height="300px")))
                                   )
                                 )
                               ),

                               # ── Step 5: Classify ──────────────────────────────────────────────────────
                               fluidRow(
                                 box(width=4, title="5 \u2014 Classify", collapsible=TRUE,
                                   sliderInput("cls_delta",    "Delta:", min=0.1, max=2,   value=0.5, step=0.1),
                                   sliderInput("cls_min_area", "Min area:", min=0, max=1, value=0.1, step=0.05),
                                   uiOutput("cls_channel_labels_ui"),
                                   actionButton("cls_classify", "Classify", class="btn-primary btn-block"),
                                   hr(),
                                   tableOutput("cls_classify_preview")
                                 ),
                                 box(width=8, title="Classification Diagnostics", collapsible=TRUE,
                                   tabsetPanel(
                                     tabPanel("Counts",
                                       br(),
                                       withSpinner(plotlyOutput("cls_diag_classify_bar", height="290px"))),
                                     tabPanel("Score | Channel",
                                       br(),
                                       withSpinner(plotlyOutput("cls_diag_score_channel", height="300px"))),
                                     tabPanel("Score | Condition",
                                       br(),
                                       helpText(style="font-size:10px; color:#888;",
                                                "ECDF of channel scores faceted by the SE column selected above."),
                                       withSpinner(plotlyOutput("cls_diag_score_condition", height="290px"))),
                                     tabPanel("Score | Classification",
                                       br(),
                                       helpText(style="font-size:10px; color:#888;",
                                                "ECDF per channel, faceted by Classification. Color = SE column above (or Channel if not set)."),
                                       withSpinner(plotlyOutput("cls_diag_score_cls", height="290px"))),
                                     tabPanel("2D Correlation",
                                       br(),
                                       uiOutput("cls_corr_controls_2d"),
                                       withSpinner(plotlyOutput("cls_diag_corr_2d", height="250px")))
                                   )
                                 )
                               ),

                               # ── 3D Score Correlation ──────────────────────────────────────────────────
                               fluidRow(
                                 box(width=12, title="3D Score Correlation Explorer", collapsible=TRUE,
                                   fluidRow(
                                     column(3,
                                       tags$p(style="color:#666; font-size:11px;",
                                              "Select axes and color variable from the scored / classified columns. ",
                                              "Rotate the plot by dragging."),
                                       uiOutput("cls_corr_controls_3d")
                                     ),
                                     column(9, withSpinner(plotlyOutput("cls_diag_corr_3d", height="480px")))
                                   )
                                 )
                               ),

                               # ── Step 6: Add to SE ─────────────────────────────────────────────────────
                               fluidRow(
                                 box(width=12, title="6 \u2014 Add to SE", collapsible=TRUE,
                                   tags$p(style="color:#666; font-size:11px;",
                                          "Joins result to colData by Well + Plate_ID."),
                                   fluidRow(
                                     column(4,
                                       radioButtons("cls_merge_mode", "Storage:",
                                         choices=c("Flat"="flat", "Nested DataFrame"="nested"),
                                         selected="flat"),
                                       conditionalPanel("input.cls_merge_mode == 'nested'",
                                         textInput("cls_col_name", "Column name:", value="img_classification")
                                       ),
                                       actionButton("cls_merge_se", "Add to SE", class="btn-success btn-block")
                                     ),
                                     column(8, verbatimTextOutput("cls_merge_status"))
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



                             # ── Define Conditions ─────────────────────────────────────────────────
                             tabItem("tab_coldata",
                               fluidRow(
                                 box(width = 5, title = "Rule Builder",
                                   helpText("Define a new colData column by mapping Plate_ID + Column ranges to label values."),
                                   fluidRow(
                                     column(8, textInput("cond_new_col", "New column name:", placeholder = "e.g. Condition")),
                                     column(4, br(), actionButton("cond_add_rule", "Add Rule", icon = icon("plus"), class = "btn-info btn-sm"))
                                   ),
                                   uiOutput("cond_rules_ui"),
                                   hr(),
                                   fluidRow(
                                     column(6, actionButton("cond_apply", "Apply to SE", icon = icon("check"), class = "btn-success")),
                                     column(6, actionButton("cond_clear_rules", "Clear Rules", icon = icon("trash"), class = "btn-danger btn-sm"))
                                   ),
                                   br(),
                                   helpText(style = "font-size:11px; color:#888;",
                                     "Each rule assigns a label to wells matching the selected Plate_ID(s) and Column range. Rules are applied top-to-bottom; later rules override earlier ones.")
                                 ),
                                 box(width = 7, title = "Plate Grid Preview",
                                   fluidRow(
                                     column(6, selectInput("plate_id4", "Plate to preview:", choices = c())),
                                     column(6, selectInput("cond_preview_col", "Color by:", choices = c()))
                                   ),
                                   withSpinner(plotlyOutput("plate_view_col", height = "320px")),
                                   hr(),
                                   withSpinner(DTOutput("cond_coldata_preview"))
                                 )
                               )
                             ),

                             # ── Define Sweeps ──────────────────────────────────────────────────────
                             tabItem("tab_rowdata",
                               fluidRow(
                                 box(width = 5, title = "Recode Column Values",
                                   selectInput("row_col_select", "Column to recode:", choices = c()),
                                   uiOutput("row_recode_ui"),
                                   fluidRow(
                                     column(6, actionButton("row_add_recode", "Add mapping", icon = icon("plus"), class = "btn-info btn-sm")),
                                     column(6, actionButton("row_apply_recode", "Apply recoding", icon = icon("check"), class = "btn-success btn-sm"))
                                   ),
                                   hr(),
                                   tags$strong("Create logical column (LP step)"),
                                   helpText(style = "font-size:11px;", "Creates a TRUE/FALSE column where the selected column equals the chosen value."),
                                   fluidRow(
                                     column(4, selectInput("lp_col", "Column:", choices = c())),
                                     column(4, selectInput("lp_val", "Value:", choices = c())),
                                     column(4, textInput("lp_new_col", "New col name:", value = "LP"))
                                   ),
                                   actionButton("lp_apply", "Create logical column", class = "btn-primary btn-sm")
                                 ),
                                 box(width = 7, title = "Current rowData",
                                   withSpinner(DTOutput("row_data_table"))
                                 )
                               )
                             ),

                             # ── Change Assays ──────────────────────────────────────────────────────
                             tabItem("tab_assays",
                               tabsetPanel(
                                 tabPanel("Aggregate (colAG)",
                                   br(),
                                   fluidRow(
                                     column(5,
                                       selectInput("ag_assays", "Assays to aggregate:", choices = c(), multiple = TRUE),
                                       tags$strong("Sweep selection:"),
                                       radioButtons("ag_sweep_mode", NULL,
                                         choices = c("All sweeps" = "all", "Row index range" = "range", "Use logical column" = "logical"),
                                         selected = "all", inline = TRUE),
                                       uiOutput("ag_sweep_ui"),
                                       actionButton("ag_run", "Run colAG", icon = icon("play"), class = "btn-primary"),
                                       br(), br(),
                                       verbatimTextOutput("ag_status")
                                     ),
                                     column(7,
                                       tags$strong("Result columns added to colData:"),
                                       br(), br(),
                                       withSpinner(DTOutput("ag_result_preview"))
                                     )
                                   )
                                 ),
                                 tabPanel("Transform Columns",
                                   br(),
                                   fluidRow(
                                     column(5,
                                       selectInput("transform_col", "Column to transform:", choices = c()),
                                       selectInput("transform_fn", "Transformation:",
                                         choices = c("× 1000" = "*1000", "÷ 1000" = "/1000",
                                                     "× 1e9" = "*1e9", "÷ 1e9" = "/1e9",
                                                     "log1p" = "log1p", "exp" = "exp",
                                                     "abs" = "abs", "negate" = "negate",
                                                     "Custom expression" = "custom")),
                                       conditionalPanel("input.transform_fn == 'custom'",
                                         textInput("transform_custom_expr", "R expression (use .x for column):",
                                                   placeholder = "e.g. log1p(.x * 1000)")
                                       ),
                                       textInput("transform_new_col", "Save as (blank = overwrite):",
                                                 placeholder = "e.g. Minima_mean_nA"),
                                       actionButton("transform_apply", "Apply transform", icon = icon("check"), class = "btn-success"),
                                       br(), br(),
                                       verbatimTextOutput("transform_status")
                                     ),
                                     column(7,
                                       tags$strong("Preview (before / after):"),
                                       br(), br(),
                                       withSpinner(plotlyOutput("transform_preview", height = "350px"))
                                     )
                                   )
                                 ),
                                 tabPanel("Dimensionality Reduction",
                                   br(),
                                   fluidRow(
                                     column(5,
                                       selectInput("dr_assays", "Assays:", choices = c(), multiple = TRUE),
                                       selectInput("dr_colnames", "Numeric colData columns:", choices = c(), multiple = TRUE),
                                       numericInput("dr_k", "k clusters:", value = 3, min = 1, max = 30),
                                       selectInput("dr_scaling", "Scaling:",
                                         choices = c("Within assay" = "within", "Global" = "global", "None" = "none"),
                                         selected = "within"),
                                       checkboxInput("dr_center", "Center before scaling", value = FALSE),
                                       actionButton("dr_run", "Run reducedDim.Cellwise", icon = icon("play"), class = "btn-primary"),
                                       br(), br(),
                                       verbatimTextOutput("dr_status")
                                     ),
                                     column(7,
                                       tags$strong("PCA / Cluster preview:"),
                                       br(), br(),
                                       withSpinner(plotlyOutput("dr_preview", height = "380px")),
                                       br(),
                                       selectInput("dr_color_col", "Color by:", choices = c())
                                     )
                                   )
                                 )
                               )
                             ),

                             # ── Filter Wells ───────────────────────────────────────────────────────
                             tabItem("tab_filter_wells",
                               fluidRow(
                                 box(width = 5, title = "Filter Wells",
                                   helpText("Subset the SE to only keep wells matching all active filters."),
                                   selectInput("fw_col", "Filter by column:", choices = c()),
                                   uiOutput("fw_val_ui"),
                                   radioButtons("fw_mode", "Keep wells where value is:",
                                     choices = c("IN selection" = "keep", "NOT IN selection" = "exclude"),
                                     selected = "keep", inline = TRUE),
                                   actionButton("fw_add", "Add filter", icon = icon("plus"), class = "btn-info btn-sm"),
                                   br(), br(),
                                   uiOutput("fw_active_filters_ui"),
                                   hr(),
                                   textOutput("fw_preview_n"),
                                   br(),
                                   actionButton("fw_apply", "Apply filter to SE", icon = icon("scissors"), class = "btn-danger"),
                                   br(), br(),
                                   verbatimTextOutput("fw_status")
                                 ),
                                 box(width = 7, title = "Preview",
                                   withSpinner(plotlyOutput("fw_plate_preview", height = "340px")),
                                   br(),
                                   withSpinner(DTOutput("fw_coldata_preview"))
                                 )
                               )
                             ),


                             tabItem("tab_annotation",
                                     fluidRow(
                                       column(width = 9,
                                         box(width = 12, title = "Current Image (BF  |  Fluoro)",
                                           uiOutput("ann_image_ui"),
                                           div(style = "margin-top: 16px;",
                                             uiOutput("ann_class_buttons")
                                           ),
                                           div(style = "margin-top: 12px; display: flex; align-items: center; gap: 8px;",
                                               actionButton("ann_skip", "Skip / Next", class = "btn-default"),
                                               actionButton("ann_undo", "Undo last",   class = "btn-warning"),
                                               span(style = "margin-left: 4px; color: #888; font-size: 12px;",
                                                    textOutput("ann_progress", inline = TRUE))),
                                           hr(),
                                           fluidRow(
                                             column(6,
                                               div(style = "display:flex; align-items:center; gap:8px;",
                                                 tags$strong("BF", style = "font-size:12px;"),
                                                 actionButton("ann_bf_reset", icon("rotate-left"),
                                                              class = "btn-xs btn-default")),
                                               fluidRow(
                                                 column(6, sliderInput("ann_bf_b", "Brightness",
                                                           0.2, 4, 1, 0.1, width = "100%")),
                                                 column(6, sliderInput("ann_bf_c", "Contrast",
                                                           0.2, 4, 1, 0.1, width = "100%"))
                                               )
                                             ),
                                             column(6,
                                               div(style = "display:flex; align-items:center; gap:8px;",
                                                 tags$strong("Fluoro", style = "font-size:12px;"),
                                                 actionButton("ann_flu_reset", icon("rotate-left"),
                                                              class = "btn-xs btn-default")),
                                               fluidRow(
                                                 column(6, sliderInput("ann_flu_b", "Brightness",
                                                           0.2, 4, 1, 0.1, width = "100%")),
                                                 column(6, sliderInput("ann_flu_c", "Contrast",
                                                           0.2, 4, 1, 0.1, width = "100%"))
                                               )
                                             )
                                           ),
                                           hr(),
                                           fluidRow(
                                             column(12,
                                               tags$h5("Resume from CSV"),
                                               fileInput("ann_upload_csv", NULL,
                                                         accept = ".csv",
                                                         placeholder = "Upload previous CSV…",
                                                         buttonLabel = "Browse…",
                                                         width = "100%"),
                                               helpText(style = "font-size:11px; margin-top:-10px;",
                                                 "Loads results and skips already-labeled wells."),
                                               hr(),
                                               tags$h5("Results"),
                                               div(style = "display:flex; gap:8px; flex-wrap:wrap; margin-bottom:10px;",
                                                 downloadButton("ann_download_csv", "Download CSV"),
                                                 actionButton("ann_merge_results", "Merge & deduplicate",
                                                              class = "btn-info btn-sm"),
                                                 actionButton("ann_clear_results", "Clear results",
                                                              class = "btn-danger btn-sm")
                                               ),
                                               div(style="margin-bottom:10px;",
                                                 actionButton("ann_update_se", icon("save"), "Save to SE colData",
                                                              class="btn-success btn-sm"),
                                                 helpText(style="font-size:11px; margin-top:4px;",
                                                   "Writes merged annotations to colData(SE)[[\"manual_ann_<YourName>\"]]. ",
                                                   "Only triggered on click.")
                                               ),
                                               tableOutput("ann_results_preview")
                                             )
                                           )
                                         )
                                       ),
                                       column(width = 3,
                                         box(width = 12, title = "Setup",
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
                                           actionButton("ann_clear_classes", "Clear labels", class = "btn-danger btn-sm")
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
                                                   selectInput("plate_id3", "Select Plate:", choices = c()),
                                                   selectizeInput("clusteredwells", "Selected Wells:",
                                                                  choices = c(), multiple = TRUE),
                                                   div(style = "margin-top: 10px;"),
                                                   actionButton("reset_well", "Reset Well Selection"),
                                                   div(style = "margin-top: 50px;")
                                                   )
                                          )
                                      ),

                                     # ── Supervised LDA ──────────────────────────────────────────────────
                                     box(width = 12,
                                       title = "Supervised Classification (LDA)",
                                       collapsible = TRUE, collapsed = TRUE,
                                       tabsetPanel(type = "tabs",

                                         # ── Setup tab ──────────────────────────────────────────────────
                                         tabPanel("Setup",
                                           br(),
                                           fluidRow(
                                             column(4,
                                               tags$strong("Features"),
                                               selectizeInput("lda_assays", "Assays:",
                                                 choices = c(), multiple = TRUE,
                                                 options = list(placeholder = "Select assays…")),
                                               radioButtons("lda_assay_agg", "Assay aggregation:",
                                                 choices = c("All sweeps (IV curve)" = "all",
                                                             "Mean across sweeps"    = "mean",
                                                             "SD across sweeps"      = "sd"),
                                                 selected = "all", inline = FALSE),
                                               selectizeInput("lda_coldata_feats", "colData columns:",
                                                 choices = c(), multiple = TRUE,
                                                 options = list(placeholder = "Optional numeric columns…")),
                                               tags$strong("Scaling pipeline"),
                                               checkboxGroupInput("lda_scale_steps",
                                                 label = NULL,
                                                 choices = c(
                                                   "Assay-type normalisation (preserves IV shape)" = "assay",
                                                   "Within-plate centering (removes batch drift)"  = "plate",
                                                   "Per-feature z-score (LDA stability)"           = "feature"),
                                                 selected = c("assay", "plate", "feature"))
                                             ),
                                             column(4,
                                               tags$strong("Training"),
                                               selectInput("lda_label_col", "Class label column:",
                                                 choices = c()),
                                               selectInput("lda_method", "LDA method:",
                                                 choices = c(
                                                   "Global (one model)" = "global",
                                                   "Per-plate"          = "per_plate",
                                                   "PCA → LDA"          = "pca"),
                                                 selected = "global"),
                                               conditionalPanel("input.lda_method == 'pca'",
                                                 sliderInput("lda_pca_var", "PCA variance to retain:",
                                                   min = 0.5, max = 0.999, value = 0.95, step = 0.01)
                                               ),
                                               br(),
                                               actionButton("lda_fit_btn", "Fit LDA",
                                                 icon = icon("play"), class = "btn-primary btn-block"),
                                               br(),
                                               verbatimTextOutput("lda_fit_status")
                                             ),
                                             column(4,
                                               tags$strong("Predict & write to SE"),
                                               textInput("lda_out_prefix", "Output column prefix:",
                                                 value = "lda"),
                                               actionButton("lda_predict_btn", "Predict all wells",
                                                 icon = icon("wand-magic-sparkles"),
                                                 class = "btn-success btn-block"),
                                               br(),
                                               actionButton("lda_cv_btn", "Run LOPO cross-validation",
                                                 icon = icon("rotate"), class = "btn-default btn-block"),
                                               br(),
                                               uiOutput("lda_cv_status_ui")
                                             )
                                           )
                                         ),

                                         # ── Results tab ────────────────────────────────────────────────
                                         tabPanel("Results",
                                           br(),
                                           fluidRow(
                                             column(7,
                                               tags$strong("LD Projection"),
                                               fluidRow(
                                                 column(5, selectInput("lda_ld_x", "X axis:", choices = c())),
                                                 column(5, selectInput("lda_ld_y", "Y axis:", choices = c())),
                                                 column(2, br(),
                                                   actionButton("lda_scatter_refresh", icon("refresh"),
                                                     class = "btn-xs btn-default"))
                                               ),
                                               withSpinner(plotlyOutput("lda_scatter_plot", height = "360px"))
                                             ),
                                             column(5,
                                               tags$strong("Feature importance"),
                                               fluidRow(
                                                 column(7, selectInput("lda_imp_ld", "LD axis:",
                                                   choices = c(), multiple = FALSE)),
                                                 column(5, numericInput("lda_imp_top", "Top N:",
                                                   value = 20, min = 5, max = 100, step = 5))
                                               ),
                                               withSpinner(plotlyOutput("lda_importance_plot", height = "360px"))
                                             )
                                           ),
                                           fluidRow(
                                             column(12,
                                               tags$strong("Confusion matrix"),
                                               fluidRow(
                                                 column(4, selectInput("lda_conf_label", "True label column:",
                                                   choices = c())),
                                                 column(8, verbatimTextOutput("lda_confusion_out"))
                                               )
                                             )
                                           )
                                         ),

                                         # ── Validation tab ──────────────────────────────────────────────
                                         tabPanel("Validation",
                                           br(),
                                           helpText(style = "color:#666;",
                                             "Leave-one-plate-out cross-validation. ",
                                             "Click 'Run LOPO cross-validation' in the Setup tab first."),
                                           fluidRow(
                                             column(6, withSpinner(plotlyOutput("lda_cv_plot", height = "300px"))),
                                             column(6, withSpinner(tableOutput("lda_cv_table")))
                                           ),
                                           conditionalPanel("input.lda_method == 'pca'",
                                             hr(),
                                             tags$strong("PCA scree (retained PCs highlighted)"),
                                             withSpinner(plotlyOutput("lda_pca_scree", height = "240px"))
                                           )
                                         )
                                       )
                                     )


                             ),

                          tabItem("tab_cc_csv",
                                  ccPreviewUI_csv("cc_preview")
                          ),
                          tabItem("tab_cc_detect",
                                  ccPreviewUI_detect("cc_preview")
                          ),
                          tabItem("tab_cc_build",
                                  ccPreviewUI_build("cc_preview")
                          ),
                          tabItem("tab_mae",
                            fluidRow(
                              box(width=8, title = "Combine SEs into MultiAssayExperiment",
                                tags$p(style="color:#666;",
                                  "Select SingleCellExperiment objects to combine into a MultiAssayExperiment. ",
                                  "Each must have ", tags$code("well_id"), " and ", tags$code("plate_id"),
                                  " in its colData. Current Clamp SEs built in this app already have these. ",
                                  "For Voltage Clamp SEs, ", tags$code("Well"), " / ", tags$code("Plate_ID"),
                                  " are used automatically if the lowercase versions are absent."),
                                uiOutput("mae_se_list_ui"),
                                hr(),
                                textInput("mae_name", "MAE object name:",
                                          placeholder = "e.g. CC_VC_combined"),
                                actionButton("build_mae_btn", "Build MAE",
                                             icon  = icon("layer-group"),
                                             class = "btn-primary btn-block"),
                                uiOutput("mae_status_ui")
                              ),
                              box(width=4, title = "About MultiAssayExperiment",
                                helpText(style="font-size:11px; color:#666;",
                                  "A MultiAssayExperiment (MAE) holds multiple assay types ",
                                  "(e.g., voltage clamp + current clamp) linked by well identity. ",
                                  "Once built, the MAE appears in the dataset selector at the top. ",
                                  "Use downstream packages (MOFA2, etc.) for multi-modal analysis.")
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
