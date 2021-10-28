rm(list = ls())
library(grf)
library(randomForestSRC)
library(survival)
library(gbm)
library(glmnet)
#library(bcf)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team")
source("./grf_CSF_simulations/r-package/grf/R/dgps.R")
set.seed(123)

# *** Comparison methods ***
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/grf_CSF_simulations/experiments/csf/comparison_estimators.R")
# estimators = list(estimate_rfsrc_X_W = estimate_rfsrc_X_W,
#                   estimate_rfsrc_XW_W = estimate_rfsrc_XW_W,
#                   estimate_rfsrc_twin = estimate_rfsrc_twin,
#                   estimate_IPCW_grf = estimate_IPCW_grf,
#                   estimate_grf = estimate_grf)

estimators = list(estimate_lasso_slearner = estimate_lasso_slearner,
                  estimate_lasso_tlearner = estimate_lasso_tlearner,
                  estimate_ipcw_wocf_lasso_xlearner = estimate_ipcw_wocf_lasso_xlearner,
                  estimate_ipcw_wcf_lasso_xlearner = estimate_ipcw_wcf_lasso_xlearner,
                  estimate_ipcw_wocf_lasso_flearner = estimate_ipcw_wocf_lasso_flearner,
                  estimate_ipcw_wcf_lasso_flearner = estimate_ipcw_wcf_lasso_flearner,
                  estimate_ipcw_wocf_lasso_rlearner = estimate_ipcw_wocf_lasso_rlearner,
                  estimate_ipcw_wcf_lasso_rlearner = estimate_ipcw_wcf_lasso_rlearner,

                  estimate_gbm_slearner = estimate_gbm_slearner,
                  estimate_gbm_tlearner = estimate_gbm_tlearner,
                  estimate_ipcw_wocf_gbm_xlearner = estimate_ipcw_wocf_gbm_xlearner,
                  estimate_ipcw_wcf_gbm_xlearner = estimate_ipcw_wcf_gbm_xlearner,
                  estimate_ipcw_wocf_gbm_flearner = estimate_ipcw_wocf_gbm_flearner,
                  estimate_ipcw_wcf_gbm_flearner = estimate_ipcw_wcf_gbm_flearner,
                  estimate_ipcw_wocf_gbm_rlearner = estimate_ipcw_wocf_gbm_rlearner,
                  estimate_ipcw_wcf_gbm_rlearner = estimate_ipcw_wcf_gbm_rlearner,
                  
                  estimate_grf_slearner = estimate_grf_slearner,
                  estimate_grf_tlearner = estimate_grf_tlearner,
                  estimate_ipcw_grf_xlearner = estimate_ipcw_grf_xlearner,
                  estimate_ipcw_grf_flearner = estimate_ipcw_grf_flearner,
                  estimate_ipcw_grf_rlearner = estimate_ipcw_grf_rlearner)



# *** Setup ***
out = list()
n.sim = 1 # 200
n.mc = 100000
grid = expand.grid(n = 500, # c(500, 1000, 2000, 5000)
                   p = 5,
                   n.test = 500, # 2000
                   dgp = "type1", # c("type1", "type2", "type3", "type4"), 
                   stringsAsFactors = FALSE)

for (i in 1:nrow(grid)) {
  print(paste("grid", i))
  n = grid$n[i]
  p = grid$p[i]
  n.test = grid$n.test[i]
  dgp = grid$dgp[i]

  for (sim in 1:n.sim) {
    print(paste("sim", sim))
    data = generate_causal_survival_data(n = n, p = p, dgp = dgp, n.mc = 10)
    data.test = generate_causal_survival_data(n = n.test, p = p, dgp = dgp, n.mc = n.mc)
    true.cate = data.test$cate
    true.cate.sign = data.test$cate.sign
    true.catesp = data.test$catesp
    true.catesp.sign = data.test$catesp.sign
    estimator.output = list()
    for (j in 1:length(estimators)) {
      estimator.name = names(estimators)[j]
      print(estimator.name)
      predictions = as.numeric(unlist(estimators[[estimator.name]](data, data.test, times = median(data$Y))))
      correct.classification = sign(predictions) == true.cate.sign
      dfj = data.frame(
        estimator.name = estimator.name,
        mse = mean((predictions - true.cate)^2),
        classif.rate = mean(correct.classification, na.rm = TRUE) # NA: to ignore X1 < 0.3 in DGP 4.
        )
      estimator.output[[j]] = dfj
    }
    df = do.call(rbind, estimator.output)
    df$n = n
    df$p = p
    df$n.test = n.test
    df$dgp = dgp
    df$sim = sim

    out = c(out, list(df))
  }
}
out.df = do.call(rbind, out)
write.csv(out.df, gzfile("./simulation.csv.gz"), row.names = FALSE)
