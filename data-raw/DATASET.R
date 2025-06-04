## code to prepare `DATASET` dataset goes here

devtools::load_all()
## Human Adrenal Gland Dataset SE with imaging data
l_files <- list.files(path = "data-raw/hAG/Ephys/" ,pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)
se_hAG <- prepareSE(l_files)

l_filel_filel_files <- list.files(path = "data-raw/hAG/Imaging/" ,pattern = "*.db$", recursive = TRUE, full.names = TRUE)
se_hAG <- mergeSEandImg(se_hAG, df_img)
usethis::use_data(se_hAG, overwrite = TRUE)


## iNeurons dataset with imaging data
l_files <- list.files(path = "data-raw/iNeurons/" ,pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)
se_iN <- prepareSE(l_files)


l_files <- list.files(path = "data-raw/iNeurons/" ,pattern = "*.db$", recursive = TRUE, full.names = TRUE)
df_img <- prepareImgDF(l_files, scale_num = TRUE, analysis = "coloc")

se_iN <- mergeSEandImg(se_iN, df_img, tableType = "coloc")
usethis::use_data(se_iN, overwrite = TRUE)

# ROMK dataset

l_files <- list.files(path = "data-raw/ROMK/" ,pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)
se_romk <- prepareSE(l_files)
usethis::use_data(se_romk, overwrite = TRUE)
