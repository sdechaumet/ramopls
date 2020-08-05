
M <- data.frame(matrix(nrow = 10, ncol = 3, 1:(3*10)))
colnames(M) <- c('Dose','Time','Age')

test_that("Check fun_factor_names", {
  sapply(c('1', '1,2,3', '2,3'), function(x) {expect_vector(fun_factor_names(M,x))})
  expect_equal(fun_factor_names(M,'1'), "Dose")
  expect_equal(fun_factor_names(M,'2'), "Time")
  expect_equal(fun_factor_names(M,'3'), "Age")
  expect_equal(fun_factor_names(M,'23'), "Time x Age")
  expect_equal(fun_factor_names(M,'13'), "Dose x Age")
  expect_equal(fun_factor_names(M,'1,2,23'), c("Dose", "Time", "Time x Age"))
  expect_error(fun_factor_names(M,'123'))
})
