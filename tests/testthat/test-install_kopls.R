skip_on_cran()

test_that("kopls intern installation", {
  ## Uninstall kopls package if available
  res <- rAMOPLS::install_kopls()
  expect_true(require(kopls))
})
