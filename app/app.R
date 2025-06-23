cat("Launching app\n")
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ephacRTools)
cat("Package loaded\n")
#l_files <- list.files(path = "./" ,pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)
#se_hAG <- ephacRTools::prepareSE(l_files)
#print(se_hAG)
waiterContent <- tagList(
  tags$div(
    style = "text-align:center;",
    tags$img(src = "https://raw.githubusercontent.com/dcolam/ephacRTools/master/inst/figures/ggplot2.svg", height = "100px", style = "margin-bottom: 20px;"),
    waiter::spin_1(),
    tags$h3("Please wait while the application is initialized...")
  )
)

login <- data.frame(user="guest", password_hash="$7$C6..../....KBtKRNulXGoaOUSNltgKqkOukRedn06odEeU1CkM3X/$BrdNlK7/U1rVOnJKy2qRPqtKrWd3zyK6c009oZ68VVB")
data("se_romk")
ephacRTools::tinySEV(title = "ephacRTools",waiterContent =waiterContent, logins = login, objects = list("ROMK"= se_romk))
