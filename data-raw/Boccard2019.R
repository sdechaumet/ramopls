## Create Boccard2019 dataset

## 129 compounds are available in the dataset from Ruiz whereas 126 are described in the article.
## There are also 2 duplicated compound :
temp <- data.table::fread("./data-raw/data_Boccard2019_2.txt", header = T, check.names = T)
temp[, sampleid := paste0("ID", formatC(1:.N, width = 2, flag = "0"))]
datamatrix <- temp[, -c("Chemical", "Batch"), keyby = "sampleid", with = F]
setkey(datamatrix, sampleid)
setcolorder(datamatrix)

data_Boccard2019 <- list("datamatrix" = datamatrix,
                         "samplemetadata" = temp[, .(Chemical, Batch), keyby = "sampleid"])

## Export to formatted data.table
usethis::use_data(data_Boccard2019, overwrite = T)

