install.packages("RInno")
library(RInno)

# One-time: download Inno Setup
install_inno()

# Create the installer config
create_app(
  app_name    = "ephacRTools",
  app_dir     = "app",           # folder with app.R
  pkgs        = c("shiny", "shinydashboard", "plotly", "DT",
                  "SingleCellExperiment", "SummarizedExperiment",
                  "data.table", "dplyr", "ggplot2", "sechm",
                  "waiter", "shinyjs", "shinycssloaders"),
  include_R   = TRUE,            # bundle portable R
  R_version   = "4.5.1"
)
