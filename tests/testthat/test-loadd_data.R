## Load data (from run_AMOPLS)
test_that("Check load data", {
  factor_names <- c("Exposure time", "Dose")
  t <- lapply(
    list("datatable" = data_exemple$datamatrix,
         "matrix" = as.matrix(data_exemple$datamatrix, rownames = 1)
    ), function(x) {
      test_data <- fun_load_data(x, data_exemple$samplemetadata, factor_names)
      expect_true(length(test_data) == 2)
      expect_true(all(names(test_data) %in% c("dataset", "factors")))
      expect_true(nrow(test_data$dataset) ==  nrow(data_exemple$datamatrix))
      expect_true(ncol(test_data$dataset) ==  ncol(data_exemple$datamatrix) - 1)
      expect_true(nrow(test_data$factors) == nrow(data_exemple$samplemetadata))
      expect_true(ncol(test_data$factors) == length(factor_names))
      expect_true(all(rownames(test_data$dataset) %in% data_exemple$datamatrix[,1][[1]]))
      expect_true(all(rownames(test_data$factors) %in% data_exemple$datamatrix[,1][[1]]))
    })
})

