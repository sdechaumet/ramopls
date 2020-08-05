context("Pure Effect factor permutation")

data_test <- s <- list("general" = list("equation.elements" = list(1,2,c(1,2)),
                                        "factors" = data.table("id" = paste0("ID", 1:100),
                                                               "class" = rep(c("classA","classB"), each = 50),
                                                               "grp" = rep(c("grp01","grp02","grp03"), length.out = 100),
                                                               "trt" = sample(c("trtA","trtB", "trtC", "trtD", "trtE"), size = 100, replace = T)) %>% as.matrix(rownames = "id")))
test_that("Permutation on 2 classes and interaction", {
  temp_data <- data_test
  temp_data$general$factors <- temp_data$general$factors[, 1:2]
  temp <- fun_perm_settings(s = temp_data, ee = 1, perm_t = "1")$general$factor
  expect_identical(dim(temp_data$general$factors), dim(temp))
  expect_identical(temp_data$general$factors[, 2], temp[, 2])
  expect_false(identical(temp_data$general$factors[, 1], temp[, 1]))
  expect_false(identical(temp_data$general$factors[, 1], temp[, 2]))

  temp <- fun_perm_settings(s = temp_data, ee = 2, perm_t = "2")$general$factor
  expect_identical(dim(temp_data$general$factors), dim(temp))
  expect_identical(temp_data$general$factors[, 1], temp[, 1])
  expect_false(identical(temp_data$general$factors[, 2], temp[, 2]))
  expect_false(identical(temp_data$general$factors[, 2], temp[, 1]))

  temp <- fun_perm_settings(s = temp_data, ee = 3, perm_t = "12")$general$factor
  expect_identical(dim(temp_data$general$factors), dim(temp))
  expect_identical(temp_data$general$factors[, 1], temp[, 1])
  expect_identical(temp_data$general$factors[, 2], temp[, 2])
})

test_that("Permutation on 3 classes and interaction", {
  temp_data <- data_test
  temp_data$general$factors <- temp_data$general$factors[, 1:3]
  temp <- fun_perm_settings(s = temp_data, ee = 1, perm_t = "1")$general$factor
  expect_identical(dim(temp_data$general$factors), dim(temp))
  expect_identical(temp_data$general$factors[, 2], temp[, 2])
  expect_identical(temp_data$general$factors[, 3], temp[, 3])
  expect_false(identical(temp_data$general$factors[, 1], temp[, 1]))
  expect_false(identical(temp_data$general$factors[, 1], temp[, 2]))
  expect_false(identical(temp_data$general$factors[, 1], temp[, 3]))

  temp <- fun_perm_settings(s = temp_data, ee = 2, perm_t = "2")$general$factor
  expect_identical(dim(temp_data$general$factors), dim(temp))
  expect_identical(temp_data$general$factors[, 1], temp[, 1])
  expect_identical(temp_data$general$factors[, 3], temp[, 3])
  expect_false(identical(temp_data$general$factors[, 2], temp[, 2]))
  expect_false(identical(temp_data$general$factors[, 2], temp[, 1]))
  expect_false(identical(temp_data$general$factors[, 2], temp[, 3]))

  temp <- fun_perm_settings(s = temp_data, ee = 3, perm_t = "12")$general$factor
  expect_identical(dim(temp_data$general$factors), dim(temp))
  expect_identical(temp_data$general$factors[, 1], temp[, 1])
  expect_identical(temp_data$general$factors[, 2], temp[, 2])
  expect_identical(temp_data$general$factors[, 3], temp[, 3])
  expect_false(identical(temp_data$general$factors[, 3], temp[, 1]))
  expect_false(identical(temp_data$general$factors[, 3], temp[, 2]))
})









