require(magrittr) ; require(rAMOPLS) ; require(data.table) ; require(ggplot2)


result <- run_AMOPLS(datamatrix = data_Boccard2016$datamatrix,
                     samplemetadata = data_Boccard2016$samplemetadata,
                     factor_names = c("Time", "Distance"),
                     interaction_level = 1,
                     nb_compo_orthos = c(1,2,3),
                     scaling = T,
                     nb_perm = 1000,
                     parallel = 3)

lapply(result, function(x) {
  # x <- result[[1]]
  data.table("R2Y" = x$kOPLS$R2Yhat %>% {.[length(.)]},
             "R2Y_pval" = fun_get_perm(x) %>% {.[, sum(`R2Y p-value`) / .N]})
}) %>%
  rbindlist(use.names = T, idcol = "rn") %>%
  ggplot(., aes(rn, R2Y)) +
  geom_bar(stat = "identity", color = "black", fill = "grey") +
  theme_bw() +
  ylim(0,1) +
  geom_text(aes(label = R2Y_pval), vjust = -0.5) +
  labs(title = "Optimal number of orthogonal component",
       x = "Number of orthogonal component",
       y = "R²Y")

lapply(result, function(x) {fun_AMOPLS_summary(x, 'All') %>% {.[x$general$ee.names]} %>% knitr::kable(digits = 3)})
