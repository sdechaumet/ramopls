require("rAMOPLS")
source("./R/functions.R")

result <- run_AMOPLS(
  datamatrix = data_sample$datamatrix[, which(lapply(.SD, function(x) {any(is.na(x))}) == T)] %>% {data_sample$datamatrix[, -., with = F]},
  samplemetadata = data_sample$samplemetadata,
  factor_names = c("class", "stage", "trt"),
  interaction_level = 1,
  scaling = T,
  nb_perm = 100,
  nb_compo_orthos = 1,
  parallel = 3,
  debug = F)



#### TEST 2 EFFECTS
datamatrix <- data_Ruiz2017$datamatrix
samplemetadata <- data_Ruiz2017$samplemetadata
factor_names <- c("Exposure time", "Dose")
interaction_level = 1
scaling = T
nb_perm = 100
subsampling = NULL
nb_compo_orthos = 1
parallel = F
debug = F


result <- run_AMOPLS(datamatrix,
                     samplemetadata,
                     factor_names,
                     interaction_level,
                     scaling,
                     nb_perm,
                     subsampling,
                     nb_compo_orthos,
                     parallel,
                     debug)

### DEV
temp <- liver.toxicity$treatment %>% data.table()
temp[, Animal.Number := paste0("ID", Animal.Number)]
temp_data <- liver.toxicity$gene %>% data.table(keep.rownames = "Animal.number")

result <- run_AMOPLS(
  datamatrix = temp_data,
  samplemetadata = temp,
  factor_names = c("Dose.Group", "Time.Group"),
  interaction_level = 1,
  scaling = T,
  nb_perm = 100,
  nb_compo_orthos = 1,
  parallel = F,
  debug = F)

fun_plot_ortho(result)
fun_get_summary(result[[1]]) %>% knitr::kable(digits = 3, align = 'c')
fun_plot_optimal_scores(result[[1]])
fun_plot_optimal_loadings(result[[1]])
fun_plot_VIPs(result[[1]], "Dose.Group")
result[[1]]$output$VIP
result[[1]]$general$data[1:5,1:5]


## DEV
result <- run_AMOPLS(
  datamatrix = data_Boccard2016$datamatrix,
  samplemetadata = data_Boccard2016$samplemetadata,
  factor_names = c("Time", "Distance"),
  interaction_level = 1,
  scaling = T,
  nb_perm = 100,
  nb_compo_orthos = 1,
  parallel = F,
  debug = F)

result$orthoNb_1$general$Nb_compo_pred


## Test ensembl of anova decompo
ggpubr::ggarrange(
  ggpubr::ggarrange(
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`1`[, -1], graph = F), "ind", habillage = result_optimal$general$factors_data$Chemical),
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`1`[, -1], graph = F), "var"),
    ncol = 2
  ),
  ggpubr::ggarrange(
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`2`[, -1], graph = F), "ind", habillage = result_optimal$general$factors_data$Batch),
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`2`[, -1], graph = F), "var"),
    ncol = 2
  ),
  ggpubr::ggarrange(
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`12`[, -1], graph = F), "ind"),
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`12`[, -1], graph = F), "var"),
    ncol = 2
  ),
  ggpubr::ggarrange(
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`residuals`[, -1], graph = F), "ind"),
    factoextra::fviz(FactoMineR::PCA(result_optimal$decompo_ANOVA$`residuals`[, -1], graph = F), "var"),
    ncol = 2
  ),
  nrow = 4
)

