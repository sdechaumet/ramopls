## Create Ruiz2017 dataset

## 129 compounds are available in the dataset from Ruiz whereas 126 are described in the article.
## There are also 2 duplicated compound :
temp <- data.table::fread("./data-raw/ARO36TMTHILIC+DATA_Jolivet.csv")
dup_cp <- names(temp) %>% {.[which(duplicated(.))]}
dup_cp %>% paste(collapse = "\n") %>% message()

## Since there is no available evidence to choose, the duplicated compound are renamed as _1 and _2
lapply(dup_cp, function(x) {
  # x <- dup_cp[[1]]
  col_names <- grep(x, names(temp), value = T)
  col_names_i <- grep(x, names(temp))
  setnames(temp, col_names_i, paste0(col_names, "_", 1:length(col_names)))
}) %>% invisible()

## Export to formatted data.table
data_Ruiz2017 <- list("datamatrix" = data.table(temp[, -c(2:3)], key = "Sample code"),
                      "samplemetadata" = data.table(temp[, 1:3], key = "Sample code"))
usethis::use_data(data_Ruiz2017, overwrite = T)

## Minimal dataset for exemples
data_exemple <- list("datamatrix" = NULL, "samplemetadata" = NULL)
data_exemple$samplemetadata <- data_Ruiz2017$samplemetadata[Dose != 0.5]
data_exemple$datamatrix <- data_Ruiz2017$datamatrix[`Sample code` %in% data_exemple$samplemetadata$`Sample code`, 1:6]

usethis::use_data(data_exemple, overwrite = T)
