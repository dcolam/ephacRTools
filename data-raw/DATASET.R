## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

devtools::load_all()

l_files <- list.files(path = "data-raw/hAG/Ephys/" ,pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)

prepareMultipleDFs(l_files)


se <- prepareSE("data-raw/iNeurons/IV neurons_14.28.48_18T39265_LC_new_LC.xlsx")

library(SummarizedExperiment)

rowData(se)$Description



