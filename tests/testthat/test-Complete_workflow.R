
## Run AMOPLS
# require(rAMOPLS) ; require(testthat)

PermNb <- 10
nb_compo_orthos <- 1:2

AMOPLS_result <- run_AMOPLS(
  data_exemple$datamatrix,
  data_exemple$samplemetadata,
  factor_names = c("Exposure time", "Dose"),
  interaction_level = 1,
  scaling = T,
  nb_perm = PermNb,
  subsampling = NULL,
  nb_compo_orthos = nb_compo_orthos,
  parallel = F,
  debug = F)

AMOPLS_result_parallel <- run_AMOPLS(
  data_exemple$datamatrix,
  data_exemple$samplemetadata,
  factor_names = c("Exposure time", "Dose"),
  interaction_level = 1,
  scaling = T,
  nb_perm = PermNb,
  subsampling = NULL,
  nb_compo_orthos = nb_compo_orthos,
  parallel = T,
  debug = F)

## Unbalanced groups
AMOPLS_result_unbalanced <- run_AMOPLS(
  data_Boccard2016$datamatrix[-5, 1:10],
  data_Boccard2016$samplemetadata[-5, ],
  factor_names = c("Time", "Distance"),
  interaction_level = 1,
  scaling = T,
  nb_perm = PermNb,
  subsampling = 3,
  nb_compo_orthos = nb_compo_orthos,
  parallel = F,
  debug = F)

data_list <- list(
  AMOPLS_result,
  AMOPLS_result_parallel,
  AMOPLS_result_unbalanced
)

test_that("Testing run_AMOPLS structure", {
  t <- lapply(
    data_list, function(x) {
      expect_equal(class(x), "list")
      expect_equal(length(x), max(nb_compo_orthos))
      expect_equal(unique(x[[1]]$output$Permutation_result$PermNb), PermNb)
    }
  )
})

test_that("Testing orthogonal plot", {
  t <- lapply(
    data_list, function(x) {
      temp <- fun_plot_ortho(x)
      expect_true("ggplot" %in% class(temp))
      expect_true(nrow(temp$data) == max(nb_compo_orthos))
      expect_true(all(temp$data$R2Y %in% sapply(x, function(y) {y$output$R2Y})))
      rm(temp)
    }
  )
})

test_that("Check summary", {
  t <- lapply(
    data_list, function(x) {
      temp <- fun_get_summary(x[[1]])
      expect_true("data.table" %in% class(temp))
      expect_true(all(names(x$general$factors_data)[-1] %in% temp$`Effect Name`))
      expect_equal(dim(temp), dim(x[[1]]$output$Summary))
    }
  )
})

test_that("Check optimal score plot", {
  t <- lapply(
    data_list, function(x) {
      temp <- fun_plot_optimal_scores(x[[1]])
      expect_true(all(c("ggplot", "ggarrange") %in% class(temp)))
    }
  )
})

test_that("Check optimal loadings plot", {
  t <- lapply(
    data_list, function(x) {
      temp <- fun_plot_optimal_loadings(x[[1]])
      expect_true(all(c("ggplot", "ggarrange") %in% class(temp)))
    }
  )
})

test_that("Check VIPs", {
  t <- lapply(
    data_list, function(x) {
      sapply(1:5, function(y) {
        temp <- fun_plot_VIPs(x[[1]], VIP_nb = y)
        expect_true(all(c("ggplot") %in% class(temp)))
        expect_true(length(temp$data[, unique(id)]) == y)
      })
    }
  )
})


#### STEP by STEP


test_that("Check pre-processing levels as letter", {
  factor_names <- c("Exposure time", "Dose")
  test_samplemetadata <- data_exemple$samplemetadata[, `Exposure time` := rep(c("A", "B"), each = 6)]
  test_data <- fun_load_data(data_exemple$datamatrix, test_samplemetadata, factor_names)
  test_result <- rAMOPLS:::fun_pre_processing(
    test_data,
    Nb_compo_ortho = 2,
    equation.elements = c("1", "2"),
    scaling = T,
    only.means.matrix = F
  )
  expect_true(class(test_result) == "list")
  expect_true(is.numeric(test_result$general$factors))

})

test_that("Check pre-processing different row numbers", {
  factor_names <- c("Exposure time", "Dose")
  test_data <- fun_load_data(data_exemple$datamatrix, data_exemple$samplemetadata, factor_names)
  test_data$dataset <- test_data$dataset[-10,]
  expect_error(
    rAMOPLS:::fun_pre_processing(
      test_data,
      Nb_compo_ortho = 2,
      equation.elements = c("1", "2"),
      scaling = T,
      only.means.matrix = F
    )
  )
})

test_that("Check pre-processing scaling", {
  factor_names <- c("Exposure time", "Dose")
  test_samplemetadata <- data_exemple$samplemetadata[, `Exposure time` := rep(c("A", "B"), each = 6)]
  test_data <- fun_load_data(data_exemple$datamatrix, test_samplemetadata, factor_names)
  lapply(c(NULL, FALSE, TRUE), function (x) {
    # x <- TRUE
    test_result <- rAMOPLS:::fun_pre_processing(
      test_data,
      Nb_compo_ortho = 2,
      equation.elements = c("1", "2"),
      scaling = x,
      only.means.matrix = F
    )
    expect_true(class(test_result) == "list")
    expect_true(is.numeric(test_result$general$factors))
  })
})





