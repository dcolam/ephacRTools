## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

devtools::load_all()

l_files <- list.files(path = "data-raw/hAG/Ephys/" ,pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)

prepareMultipleDFs(l_files)
se <- prepareSE("data-raw/hAG/Ephys/NCI_ramp_ATP1A1_080325_18T39421.xlsx")

se <- prepareSE(l_files)

class(l_files)

colData(se)

library(SummarizedExperiment)




