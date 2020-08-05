context("Testing run_AMOPLS function")

result <- run_AMOPLS(datamatrix = data_Ruiz2017$datamatrix,
                     samplemetadata = data_Ruiz2017$samplemetadata,
                     factor_names = c("Exposure time", "Dose"),
                     nb_compo_orthos = 1:2,
                     nb_perm = 10)


test_that("run_AMOPLS output str check", {
  expect_equal(length(result), 2)
  expect_equal(names(result[[1]]), c("general", "decompo_ANOVA", "kOPLS", "output"))
})


