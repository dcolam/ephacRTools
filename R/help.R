.getHelp <- function(topic){
  switch(topic,
         general=modalDialog(title="Quick start", easyClose=TRUE, tags$p(
           "The ", tags$em("tiny SummarizedExperiment Viewer")," (tinySEV) is
      organized around SummarizedExperiment (SE) objects, which contain all data
      relative to an experiment. You can select among the different SE objects
      available using the drowndown list at the top-left of the dashboard, just
      below the application's title. If the list is empty, it means that the app
      was not preloaded with objects, but you can upload your own object in the
      'Upload object' tab on the left (assuming the feature hasn't been
      disabled)."),
           tags$p("When selecting an object, the 'DEA' and 'Enrichments' tabs on the
        left will be updated with the number of differential expression analyses
        (DEAs) and enrichment analyses included in the object (if none are
        included, 'N/A' will be shown. These can be browsed as tables in the
        respective tabs; in addition, DEA results can also be visualized as
        volcano plots."),
           tags$p("Independently of whether DEA or enrichment results are available,
        the expression of single genes or features can be plotted in the
        'Plot gene' tab."),
           tags$p("The heatmap functionalities instead allow you to visualize sets
        of genes. To this end, first enter your genes of interest in the 'Genes'
        tab, and then go to the 'Heatmap' tab to visualize them."),
           tags$p("Finally, note that gene selections can be transferred form the
        DEA/Enrichments tabs into the gene/heatmap tabs!"),
           footer=paste("tinySEV version", as.character(packageVersion("tinySEV")))
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
                        tags$ul(
                          tags$li("See the ", tags$a(
                            href=paste0("https://bioconductor.org/packages/release/bioc/vignettes/",
                                        "SummarizedExperiment/inst/doc/SummarizedExperiment.html"),
                            "SummarizedExperiment documentation", target="_blank"),
                            " for a general introduction to SummarizedExperiment objects."),
                          tags$li(tags$b("Differential expression analyses (DEAs): "),
                                  "To be available to tinySEV, differential expression analyses
        should be placed in the ", tags$code("rowData"), " of the object, in a
        column prefixed with 'DEA.'. For example, assuming that ",
                                  tags$code("dea"), " is your dataframe of differential expression results
        (with feature names as row names), you can simply do:",
                                  tags$pre("rowData(se)$DEA.myAnalysis <- dea[row.names(se),]"),
                                  "DEAs are accepted in formats of common differential expression packages
        (e.g. edgeR, limma, DESeq2)."),
                          tags$li(tags$b("Enrichment analysis: "),
                                  "Enrichment analysis results should be stored as a named list of
        data.frames in the metadata as follows:",
                                  tags$pre("metadata(se)$EA <- list( analysis1=df1, analysis2=df2 )"),
                                  "The dataframes will be displayed as such and require no special
        formatting. However, if a column named 'genes' is present (containing
        comma-separated lists of genes), it will be used to support extra
        functionalities."),
                          tags$li(tags$b("Metadata: "), "Additional meta-data information can be
        saved as text in the 'title', 'name', 'source', or 'description' slots,
        e.g. ",
                                  tags$pre(
                                    'metadata(se)$description <- "Whatever description you wish displayed..."')),
                          tags$li(tags$b("Default plotting parameters: "),
                                  "Some default plotting parameters can be stored in the ",
                                  tags$code('default_view'), "metadata slot. For example, if the object
        contains an assay 'logcpm' which we would like to be selected by
        default, we can save this information in the object as follows:",
                                  tags$pre("metadata(se)$default_view <- list(assay='logcpm')"),
                                  "The following 'default_view' elements are recognized: assay, groupvar,
        colvar, and gridvar."),
                          tags$li(tags$b("Annotation colors: "), "Instead of the randomly-generated
        annotation colors, colors for specific annotation variables can be
        specified via the", tags$code("anno_colors"), " metadata slot, for
        example:", tags$pre(
          'metadata(se)anno_colors <- list(genotype=c(WT="grey", mutant="red"))'),
          "See the ", tags$a("sechm documentation", target="_blank",
                             href=paste0("https://bioconductor.org/packages/release/bioc/",
                                         "vignettes/sechm/inst/doc/sechm")),
          " for more information.")
                        )
         ),
         modalDialog(title="Unknown topic", easyClose=TRUE,
                     "No help is currently available on this topic.")
  )
}
