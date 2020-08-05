test_that("Test balance data (balanced dataset)", {
  factor_names <- c("Exposure time", "Dose")
  test_data <- fun_load_data(data_exemple$datamatrix, data_exemple$samplemetadata, factor_names)
  expect_true(fun_is_balanced(test_data, "Dose", 0))
  test_result <- rAMOPLS:::fun_balance_data(test_data)
  expect_identical(test_data, test_result)
})

test_that("Test balance data (unbalanced dataset)", {
  factor_names <- c("Exposure time", "Dose")
  test_data <- fun_load_data(data_exemple$datamatrix[-10], data_exemple$samplemetadata[-10], factor_names)
  expect_false(fun_is_balanced(test_data, "Dose", 0))
  test_result <- rAMOPLS:::fun_balance_data(test_data)
  expect_true(data.table(test_result$factors)[, .N, by = .(`Exposure time`, Dose)][, unique(N)] %>% length() == 1)
})


