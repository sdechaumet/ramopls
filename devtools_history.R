
require(magrittr)

c(
  "MetStaT",
  "kopls",
  "magrittr",
  "data.table",
  "ggplot2",
  "plyr",
  "pracma",
  "ggrepel"
) %>%
  lapply(., usethis::use_package)

usethis::use_build_ignore("devtools_history.R")

# DATASET 1 -----------------------------------------------------------
require(data.table)
temp <- readRDS("../../BAPTISTE/CODES/_data/681024_AMOPLS.rds")

temp$variablemetadata[, varid := paste0("Var", formatC(1:.N, width = 5, format = "d", flag = "0"))]
temp$variablemetadata[, Group := ifelse(is.na(Compound), F, T)]
setnames(temp$datamatrix, temp$variablemetadata$variableid, temp$variablemetadata$varid)
temp$variablemetadata <- temp$variablemetadata[, .(varid, Group)]
col_id <- c("class", "stage", "trt")
setnames(temp$samplemetadata, c("class", "Stage", "Mouse_treatment"), col_id)

temp$samplemetadata[, (col_id) := lapply(.SD, function(x) {factor(x) %>% as.numeric()}), .SDcols = col_id]
data_sample <- temp
usethis::use_data(data_sample, overwrite = T)

# DATASET 2 -----------------------------------------------------------
library(mixOmics)
data(liver.toxicity)
usethis::use_data(liver.toxicity)
rm(temp, data_sample, col_id, liver.toxicity)
