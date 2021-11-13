rm(list = ls())
library(survival)
library(grf)
library(randomForestSRC)
library(survival)
library(gbm)
library(glmnet)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/risk-vs-hte-R")
source("sprint_parametric_simulation.R")

setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/HTE_experiments")
source("./dgps.R")
source("./comparison_estimators.R")
source("./flearner_surv/Flasso.R")
source("./flearner_surv/Fgbm.R")
source("./flearner_surv/Fgrf.R")
source("./rlearner_surv/utils.R")
source("./rlearner_surv/scoxph.R")
source("./rlearner_surv/slasso_surv.R")
source("./rlearner_surv/rlasso.R")
source("./rlearner_surv/rgbm.R")
source("./rlearner_surv/rgrf.R")

source("./grf-csf-prob/r-package/grf/R/causal_survival_forest.R")
source("./grf-csf-prob/r-package/grf/R/input_utilities.R")
source("./grf-csf-prob/r-package/grf/R/RcppExports.R")

set.seed(123)


# *** Comparison methods ***
# estimators = list(estimate_rfsrc_X_W = estimate_rfsrc_X_W,
#                   estimate_rfsrc_XW_W = estimate_rfsrc_XW_W,
#                   estimate_rfsrc_twin = estimate_rfsrc_twin,
#                   estimate_IPCW_grf = estimate_IPCW_grf,
#                   estimate_grf = estimate_grf)

estimators = list(estimate_coxph_slearner = estimate_coxph_slearner,
                  estimate_coxph_tlearner = estimate_coxph_tlearner,
                  
                  estimate_lasso_slearner = estimate_lasso_slearner,
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
                  estimate_ipcw_grf_rlearner = estimate_ipcw_grf_rlearner,
                  
                  estimate_csf_probs = estimate_csf_probs)

estimators = list(estimate_csf_probs = estimate_csf_probs)
# *** Setup ***
out = list()
n.sim = 200 # 200
n.mc = 100000
grid = expand.grid(n = c(2000, 5000),
                   p = 5,
                   n.test = 2000,
                   dgp = c("type1", "type2", "type3", "type4", "type5", "type6", "type7"), 
                   stringsAsFactors = FALSE)
grid <- grid[c(1:11,13),]
grid[11:12,1] <- rep(9000, 2)  # use N=9000 for SPRINT dgps
i=12

for (i in 1:nrow(grid)) {
  print(paste("grid", i))
  n = grid$n[i]
  p = grid$p[i]
  n.test = grid$n.test[i]
  dgp = grid$dgp[i]

  for (sim in 1:n.sim) {
    print(paste("sim", sim))
    data = generate_causal_survival_data(n = n, p = p, dgp = dgp, n.mc = 10)
    if (dgp == "type6" | dgp == "type7"){
      data.test = generate_causal_survival_data(n = n.test, p = p, dgp = dgp, n.mc = 10)
    }else{
      data.test = generate_causal_survival_data(n = n.test, p = p, dgp = dgp, n.mc = n.mc)
    }
    data$Y = pmax(rep(0.001, length(data$Y)), data$Y)
    true.cate = data.test$cate
    true.cate.sign = data.test$cate.sign
    true.catesp = data.test$catesp
    true.catesp.sign = data.test$catesp.sign
    estimator.output = list()
    for (j in 1:length(estimators)) {
      estimator.name = names(estimators)[j]
      print(estimator.name)
      if (dgp == "type5"){
        predictions = as.numeric(unlist(estimators[[estimator.name]](data, data.test, ps = mean(as.numeric(data$W)), times = median(data$Y))))
      }else if (dgp == "type6"|dgp == "type7"){
        predictions = as.numeric(unlist(estimators[[estimator.name]](data, data.test, times = 365.25*3)))
      }else{
        predictions = as.numeric(unlist(estimators[[estimator.name]](data, data.test, times = median(data$Y))))
      }
      correct.classification = sign(predictions) == true.catesp.sign
      dfj = data.frame(
        estimator.name = estimator.name,
        mse = mean((predictions - true.catesp)^2),
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
