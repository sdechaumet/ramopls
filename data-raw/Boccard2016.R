## code to prepare `Boccard2016` dataset goes here
if (!require("R.matlab")) {
  install.packages("R.matlab")
  require("R.matlab")
}

temp <- readMat("./data-raw/AMOPLS_demo_data.mat")
dt_samples <- data.table("sampleid" = unlist(temp$ObsNames))
data_Boccard2016 <- list("datamatrix" = data.table(dt_samples, temp$AllData),
                         "samplemetadata" = data.table(dt_samples, temp$F),
                         "variablemetadata" = data.table("variableid" = unlist(temp$VarNames))
)

setnames(data_Boccard2016$datamatrix, c("sampleid", data_Boccard2016$variablemetadata[, variableid]))
setnames(data_Boccard2016$samplemetadata, c("sampleid", "Time", "Distance"))

usethis::use_data(data_Boccard2016, overwrite = T)

