.getHelp <- function(topic){
  switch(topic,
         general=modalDialog(title="Quick start", easyClose=TRUE, tags$p(
           "The ", tags$em("ephacRTools")," is an R-package that facilitates the analysis of
           High-throughput Automated patch-clamp experiments (HT-APC) into standardized SummarizedExperiment objects.
           Furthermore, we included the integration for imaging results, coined as correlative imaging. This allows
           us to determine the identity of the patched cell."),
           footer=paste("ephacRTools version", as.character(packageVersion("ephacRTools")))
         ),
         genelists=modalDialog(title="Gene lists", easyClose=TRUE,
                               "Gene lists are simply pre-loaded lists of genes/features, which can be
      used as selection for plotting. List membership will additionally be
      flagged in the 'Plot gene' tab."),
         assay=modalDialog(title="Assays", easyClose=TRUE,
                           "The 'assay' represents the type of values to plot. These could for
      instance be raw un-normalized data (e.g. counts, intensity), normalized
      signals (e.g. tpm, logcpm, log-normalized intensity), variance-stabilized
      or corrected data, or relative signals like log-foldchanges relative to
      a reference condition (pre-defined in the object)."),
         group=modalDialog(title="Grouping", easyClose=TRUE,
                           "'Group by' determines the varialbe on the basis of which the points
      are grouped together. Typically, 'group_by' will be the same as 'Color by'."),
         grid=modalDialog(title="Grid/faceting", easyClose=TRUE,
                          "'Grid by' can be used to split a plot into subplots showing different
      according to this variable. When doing so, the 'Free axes' input determines
      whether each subplot is allowed to have its own axis and limits, or whether
      they should have the same."),
         scale=modalDialog(title="Scaling", easyClose=TRUE,
                           "'Scale rows' determines whether to scale the data by rows before plotting.
      If enabled, each row is centered around its mean and scaled by unit
      variance (z-scores)."),
         scaletrim=modalDialog(title="Colorscale trimming", easyClose=TRUE,
                               "In heatmaps of highly heteroscedastic data such as foldchanges, isolated
      outlier values can cause most of the colorscale to span a range with very
      few data points (e.g. a very dark heatmap with a single very bright data
      point). To circumvent this, it is common to trim out of the extremes
      values before mapping to the colorscale. For instance, choosing a trimming
      of 1% means that the colorscale will be based on the remaining 99% of the
      data, and values outside this range will be displayed as if it was at the
      extreme of the range.", tags$br(),
                               "Note that this parameter is only used for symmetrical scales, such as
      log-foldchanges or scaled data."),
         SE=modalDialog(title="Preparing a SummarizedExperiment object", easyClose=TRUE,
                        tags$ul("See the ", tags$a(
                            href=paste0("https://bioconductor.org/packages/release/bioc/vignettes/",
                                        "SummarizedExperiment/inst/doc/SummarizedExperiment.html"),
                            "SummarizedExperiment documentation", target="_blank"),
                            " for a general introduction to SummarizedExperiment (SE) objects. Here,
                            you can find a brief introduction on how to initialize an SE object using the ephacRTools package:",
                          tags$li(tags$b("Prepare SE object"),
                                  "In R, you can use the", tags$code("prepareSE"), " function to load an Excel file
                                  created by DataControl. The Excel-file should contain Online Analyses outpu values per sweep.
                                  We additionally recommend to add the Nanion Barcode ID (which becomes Plate_ID), the QC as well as
                                  QC metrices such as Series, Seal and Capacitance.",
                                  tags$code("prepareSE"), " will detect the single sweeps and all numeric columns which are then
                                  automatically transformed into", tags$code("rowData(se)"), " and ", tags$code("assays(se)"), ". On
                                  the other hand, Plate_ID, QC and Well information will be loaded into ", tags$code("colData(se)"),".
                                  You can add and manipulate whichever table asociated with the se",
                                  tags$pre("se <- prepare('path/to/excel')\nrowData(se)$LiquidPeriod <- paste('LP',rep(1:nrow(se)), sep='')\ncolData(se)$Condition <- ifelse(as.numeric(se$Column) %% 2 == 0,
                                           'WT',
                                           'KO')\nassays(se)$CurrenDensity <- assays(se)$Current / assays(se)$Capacitance")
                                  ),
                          tags$li(tags$b("Adding Image Data: "),
                                  "We integrated the ability to import SQLite databases outputed by the Fiji ", tags$a(
                                    href="https://github.com/dcolam/Cluster-Analysis-Plugin", "Cluster Analysis Plugin.", target="_blank"),
                          "Refer to our original publication for optimal setup of the imaging pipeline,
                          as well as the segmentation of the patching pore. To connect the two datasets you can follow these steps: ",
                          tags$pre("df_img <- prepareImgDF('path/to/database.db', scale_num = TRUE)\nse <- mergeSEandImg(se, df_img)"),
                          "You can import both single channels and colocalized channels.
                          In order for the function to work, you will need a column with Well information of the image and ideally a Plate_ID.
                          If your database does not have that, you can also add it after importing the .db into a dataframe with ", tags$code("prepareImgDF"),
                          ". Your se object is now ready for saving and importing into this ShinyApp.",
                          tags$pre("saveRDS(se, 'path/to/save/se.rds)'\nephacRTools::tinySEV(list('My SE' = se)")
                                  )
                        )
         ),
         modalDialog(title="Unknown topic", easyClose=TRUE,
                     "No help is currently available on this topic.")
  )
}
