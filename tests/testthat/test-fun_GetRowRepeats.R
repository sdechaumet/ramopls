M <- matrix(nrow = 10, ncol = 3, 1:(3*10))
M <- rbind(M, M)
colnames(M) <- c('Dose','Time','Age')

test_that("fun_GetRowRepeats", {
  expect_true(is.list(rAMOPLS:::fun_GetRowRepeats(M)))
  expect_identical(rAMOPLS:::fun_GetRowRepeats(M)[[1]], unique(M))
  expect_equal(length(rAMOPLS:::fun_GetRowRepeats(M)[[2]]), nrow(unique(M)))
})




