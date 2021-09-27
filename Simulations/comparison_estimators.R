# *** Comparison methods ***

#' Compute E[T | X]
#'
#' @param S.hat The estimated survival curve.
#' @param Y.grid The time values corresponding to S.hat.
#' @return A vector of expected values.
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))

  c(cbind(1, S.hat) %*% grid.diff)
}

# "SRC1"
estimate_rfsrc_X_W = function(data, data.test) {
  df = data.frame(Y=data$Y, D=data$D, cbind(x=data$X, w=data$W))
  ntree = 500
  fit = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df)
  n.test = nrow(data.test$X)
  W1 = rep(1, n.test)
  W0 = rep(0, n.test)
  p1 = predict(fit, data.frame(cbind(data.test$X, w=W1)))
  p0 = predict(fit, data.frame(cbind(data.test$X, w=W0)))

  expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
}

# "SRC2"
estimate_rfsrc_XW_W = function(data, data.test) {
  df = data.frame(Y=data$Y, D=data$D, cbind(x=data$X, xw=data$X*data$W, w=data$W))
  ntree = 500
  fit = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df)
  n.test = nrow(data.test$X)
  W1 = rep(1, n.test)
  W0 = rep(0, n.test)
  p1 = predict(fit, data.frame(cbind(data.test$X, data.test$X*W1, w=W1)))
  p0 = predict(fit, data.frame(cbind(data.test$X, data.test$X*W0, w=W0)))

  expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
}

# "VT"
estimate_rfsrc_twin = function(data, data.test) {
  ntree = 500
  W = data$W
  df1 = data.frame(Y=data$Y[W==1], D=data$D[W==1], data$X[W==1, ])
  df0 = data.frame(Y=data$Y[W==0], D=data$D[W==0], data$X[W==0, ])
  fit1 = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df1)
  fit0 = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df0)
  p1 = predict(fit1, data.frame(data.test$X))
  p0 = predict(fit0, data.frame(data.test$X))

  expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
}

# "CSF"
estimate_grf = function(data, data.test) {
  fit = causal_survival_forest(data$X, data$Y, data$W, data$D)
  p = predict(fit, data.test$X)

  p$predictions
}

# "IPCW" - causal forest
estimate_IPCW_grf = function(data, data.test) {
  sf.censor = survival_forest(data$X, data$Y, 1 - data$D, prediction.type = "Nelson-Aalen",
                              num.trees = 500)
  C.hat = predict(sf.censor)$predictions
  Y.relabeled = sf.censor$Y.relabeled # The event time values relabeled to consecutive integers 0/1 to length(Y.grid)
  C.Y.hat = rep(1, nrow(data$X)) # (for events before the first failure, C.Y.hat is one)
  # Pick out S_C(Yi, X) from the estimated survival curve
  C.Y.hat[Y.relabeled != 0] = C.hat[cbind(1:nrow(data$X), Y.relabeled)]
  sample.weights = 1 / C.Y.hat
  subset = data$D == 1
  cf = causal_forest(data$X[subset, ], data$Y[subset], data$W[subset], sample.weights = sample.weights[subset])
  p = predict(cf, data.test$X)

  p$predictions
}


#---------------------------------------------------------------------------------------------------------------------------#
# Using modified causal_survival_forest simulation code in "grf" package with the following modifications: 
# 1. dgps.R   -- "generate_causal_survival_data" function: the ground truths "cate" and "cate.sign" 
#                 are defined as difference in survival probabilities for all five settings 
# 2. comparison_estimators.R -- added more estimators, e.g., R-, X-, S-, T-learner, PTO, causal BART
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/grf_CSF_simulations/experiments/csf")
files.sources = list.files()
sapply(files.sources, source)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/grf_CSF_simulations/r-package/grf/R/dgps.R")


# "IPCW" - R-learner boost 
# library(devtools); install_github("xnie/rlearner"); library(rlearner)
# Using modified r-learner code with the following modifications: 
# 1. rboost.R -- "y_fit"   (line #57) changed to a binary outcome
#             -- "weights" (line #95) added the ipcw weight 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/rlearner-master/R")
files.sources = list.files()
sapply(files.sources, source)

estimate_ipcw_rlearner <- function(data, data.test){
  # IPCW weights
  sf.censor <- survival_forest(data$X, data$Y, 1 - data$D, num.trees = 500)
  C.hat <- predict(sf.censor)$predictions
  cent <- data$Y
  cent[data$D==0] <- median(data$Y)
  C.index <- rep(NA, length(cent))
  for (h in 1:length(cent)){
    C.index[h] <- which.min(abs(sort(unique(data$Y[data$D==0]))-cent[h]))
  }
  C.Y.hat <- C.hat[cbind(1:length(data$Y), C.index)]
  data$sample.weights <- 1 / C.Y.hat
  
  # Subset of uncensored subjects
  t50 <- median(data$Y)
  tempdat <- data.frame(data$Y, data$D, data$W, data$sample.weights, data$X)
  colnames(tempdat)[1:4] <- c("Y", "D", "W", "ipcw")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y>t50,]          # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$Y > t50] <- 0   # recode the event status for subjects who had events after t50

  # R-learner boosting
  rboost_fit <- rboost(x = as.matrix(binary_data[,5:9]), 
                       w = binary_data$W, 
                       y = binary_data$D,
                      ipcw = binary_data$ipcw)
  rboost_est <- -predict(rboost_fit, as.matrix(data.test$X)) # negative as predicted values are difference in risks
  rboost_est
}
mean((rboost_est-true.cate)^2)
mean((sign(rboost_est)==true.cate.sign))

# PTO 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/causalLearning-master/R")
files.sources = list.files()
sapply(files.sources, source)

estimate_ipcw_pto <- function(data, data.test){
  # IPCW weights
  sf.censor <- survival_forest(data$X, data$Y, 1 - data$D, num.trees = 500)
  C.hat <- predict(sf.censor)$predictions
  cent <- data$Y
  cent[data$D==0] <- median(data$Y)
  C.index <- rep(NA, length(cent))
  for (h in 1:length(cent)){
    C.index[h] <- which.min(abs(sort(unique(data$Y[data$D==0]))-cent[h]))
  }
  C.Y.hat <- C.hat[cbind(1:length(data$Y), C.index)]
  data$sample.weights <- 1 / C.Y.hat
  
  # Propensity score 
  ps_fit <- regression_forest(data$X, data$W, num.trees = 500)
  data$ps_score <- predict(ps_fit)$predictions
  
  # Subset of uncensored subjects
  t50 <- median(data$Y)
  tempdat <- data.frame(data$Y, data$D, data$W, data$sample.weights, data$ps_score, data$X)
  colnames(tempdat)[1:5] <- c("Y", "D", "W", "ipcw", "ps_score")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > t50,]          # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$Y > t50] <- 0     # recode the event status for subjects who had events after t50
  
  pto_fit <- PTOforest(x = as.matrix(binary_data[,6:10]), tx = binary_data$W, y = binary_data$D, 
                       pscore = binary_data$ps_score, weight = binary_data$sample.weights, 
                       min.node.size=100, num.trees=1000)
  pred_pto <- -predict(pto_fit, data.test$X)
  pred_pto
}

# Bayesian causal forest
library(bcf); library(grf)
estimate_ipcw_causalbart <- function(data, data.test){
  # IPCW weights
  sf.censor <- survival_forest(data$X, data$Y, 1 - data$D, num.trees = 500)
  C.hat <- predict(sf.censor)$predictions
  cent <- data$Y
  cent[data$D==0] <- median(data$Y)
  C.index <- rep(NA, length(cent))
  for (h in 1:length(cent)){
    C.index[h] <- which.min(abs(sort(unique(data$Y[data$D==0]))-cent[h]))
  }
  C.Y.hat <- C.hat[cbind(1:length(data$Y), C.index)]
  data$sample.weights <- 1 / C.Y.hat
  
  # Propensity score 
  ps_fit <- regression_forest(data$X, data$W, num.trees = 500)
  data$ps_score <- predict(ps_fit)$predictions
  
  # Subset of uncensored subjects
  t50 <- median(data$Y)
  tempdat <- data.frame(data$Y, data$D, data$W, data$sample.weights, data$ps_score, data$X)
  colnames(tempdat)[1:5] <- c("Y", "D", "W", "ipcw", "ps_score")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > t50,]          # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$Y > t50] <- 0     # recode the event status for subjects who had events after t50
  
  bcf_fit <- bcf(y = binary_data$D, 
                 z = binary_data$W,
                 x_control = as.matrix(binary_data[, 6:10]),
                 x_moderate = as.matrix(binary_data[, 6:10]),
                 pihat = binary_data$ps_score, 
                 w = binary_data$ipcw, 
                 nburn = 1000, 
                 nsim = 1000,
                 n_cores = 1)
  
  summary(bcf_fit) # assess convergence of our run 
  
  ps_pred_fit <- regression_forest(data.test$X, data.test$W, num.trees = 500)
  data.test$ps_pred <- predict(ps_pred_fit)$predictions
  pred_out <- predict(object = bcf_fit,
                     x_predict_control = data.test$X,
                     x_predict_moderate = data.test$X,
                     pi_pred = data.test$ps_pred,
                     z_pred = data.test$W,
                     save_tree_directory = '..')
  pred_bcf <- data.frame(Mean  = -colMeans(bcf_fit$tau),
                         Low95 = -apply(bcf_fit$tau, 2, function(x) quantile(x, 0.025)),
                         Up95  = -apply(bcf_fit$tau, 2, function(x) quantile(x, 0.975)))
  
}

# S-learner 
estimate_slearner <- function(data, data.test){
  # Q: do we want to tune hyperparameters?
  t50 <- median(data$Y)
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- data.frame(data.test$Y, data.test$D, data.test$W, data.test$X)
  colnames(testdat)[1:3] <- c("Y", "D", "W")
  testdat1 <- testdat; testdat1$W <- 1
  testdat0 <- testdat; testdat0$W <- 0
  
  gbmfit <- gbm(Surv(Y, D)~., 
                data = traindat, 
                distribution = "coxph",
                shrinkage = 0.01, 
                n.trees = 2000,   
                interaction.depth = 2, 
                n.minobsinnode = 15, 
                bag.fraction = 0.8,
                cv.folds = 10,
                n.cores = 10) 
  best_iter <- gbm.perf(gbmfit, plot.it = F, method = 'cv')
  
  # return a vector of prediction on n.trees indicating log hazard scale.f(x)
  pred.train <- predict(gbmfit_final, traindat, n.trees = best_iter)
  pred.test1 <- predict(gbmfit_final, testdat1, n.trees = best_iter)
  pred.test0 <- predict(gbmfit_final, testdat0, n.trees = best_iter)
  
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum <- basehaz.gbm(traindat$t_cvds, traindat$cvd, pred.train, t.eval = t50, cumulative = TRUE)
  
  # S(X) at median survival time 
  surf1 <- surf0 <- rep(NA,dim(testdat1)[1])
  for (j in 1:dim(testdat1)[1]){
    surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum)
    surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum)
  }
  
  pred_S_gbm <- surf1 - surf0
  pred_S_gbm
}

# T-learner 
estimate_tlearner <- function(data, data.test){
  # Q: do we want to tune hyperparameters?
  t50 <- median(data$Y)
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  traindat1 <- traindat[traindat$W==1, !colnames(traindat1) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat0) %in% c("W")]
  testdat <- data.frame(data.test$Y, data.test$D, data.test$W, data.test$X)
  colnames(testdat)[1:3] <- c("Y", "D", "W")
  testdat <- testdat[, !colnames(testdat) %in% c("W")]
  
  # Model for W = 1
  gbm_fit1 <- gbm(Surv(Y, D)~.,
                  data = traindat1,
                  distribution = "coxph",
                  shrinkage = 0.01,
                  n.trees = 2000,
                  interaction.depth = 2,
                  n.minobsinnode = 15,
                  bag.fraction = 0.8,
                  cv.folds = 10,
                  n.cores = 10)
  best_iter <- gbm.perf(gbm_fit1, plot.it = F, method = 'cv')
  
  # return a vector of prediction on n.trees indicating log hazard scale.f(x)
  pred.train1 <- predict(gbm_fit1, traindat1, n.trees = best_iter)
  pred.test1  <- predict(gbm_fit1, testdat, n.trees = best_iter)
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum1 <- basehaz.gbm(traindat1$t_cvds, traindat1$cvd, pred.train1, t.eval = t50, cumulative = TRUE)
  
  # Model for W = 0
  gbm_fit0 <- gbm(Surv(Y, D)~.,
                  data = traindat0,
                  distribution = "coxph",
                  shrinkage = 0.001,
                  n.trees = 10000,
                  interaction.depth = 2,
                  n.minobsinnode = 15,
                  bag.fraction = 0.8,
                  cv.folds = 10,
                  n.cores = 10)
  best_iter <- gbm.perf(gbm_fit0, plot.it = F, method = 'cv')
  
  # return a vector of prediction on n.trees indicating log hazard scale.f(x)
  pred.train0 <- predict(gbm_fit0, traindat0, n.trees = best_iter)
  pred.test0  <- predict(gbm_fit0, testdat, n.trees = best_iter)
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum0 <- basehaz.gbm(traindat0$t_cvds, traindat0$cvd, pred.train0, t.eval = t50, cumulative = TRUE)
  
  surf1 <- surf0 <- rep(NA,dim(testdat)[1])
  for (j in 1:dim(testdat)[1]){
    surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
    surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
  }
  
  pred_T_gbm <- surf1 - surf0
  pred_T_gbm
}

# "IPCW" - X-learner gradient boosting
estimate_ipcw_xlearner <- function(data, data.test){
  
  # fit model on W==1 
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  t50 <- median(data$Y[data$D==1])
  gbm_fit1 <- gbm(Surv(Y, D) ~ ., 
                  data = train1,
                  distribution = "coxph",
                  cv.folds = 10)
  best.iter1 <- gbm.perf(gbm_fit1, plot.it = F, method = 'cv')
  pred.train1 <- predict(gbm_fit1, train1, n.trees = best.iter1)
  pred.test1  <- predict(gbm_fit1, testdat, n.trees = best.iter1)
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum1 <- basehaz.gbm(train1$Y, train1$D, pred.train1, t.eval = t50, cumulative = TRUE)
  
  # fit model on W==0
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  gbm_fit0 <- gbm(Surv(Y, D) ~ ., 
                  data = train0,
                  distribution = "coxph",
                  cv.folds = 10)
  best.iter0 <- gbm.perf(gbm_fit0, plot.it = F, method = 'cv')
  pred.train0 <- predict(gbm_fit0, train0, n.trees = best.iter0)
  pred.test0  <- predict(gbm_fit0, testdat, n.trees = best.iter0)
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum0 <- basehaz.gbm(train0$Y, train0$D, pred.train0, t.eval = t50, cumulative = TRUE)
  
  # T-learner estimates of 1-S(t) at the median survival time
  surf1 <- surf0 <- rep(NA,dim(testdat)[1])
  for (j in 1:dim(testdat)[1]){
    surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
    surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
  }
  traindat$Tgbm1 <- 1-surf1
  traindat$Tgbm0 <- 1-surf0
  
  # X-learner
  delta <- ifelse(data$D==1|data$Y>t50,1,0)
  binary_data <- traindat[delta==1,]                           # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$Y > t50] <- 0   # recode the event status for subjects who had events after t50
  binary_data$Y[binary_data$D==1 & binary_data$Y > t50] <- t50 # recode the survival times for subjects who had events after t50
  
  # IPCW weights
  sf.censor <- survival_forest(binary_data[,4:8], binary_data$Y, 1 - binary_data$D, num.trees = 500)
  C.hat <- predict(sf.censor)$predictions
  Y.relabeled <- sf.censor$Y.relabeled 
  C.Y.hat <- rep(1, nrow(binary_data)) 
  C.Y.hat[Y.relabeled != 0] <- C.hat[cbind(1:nrow(binary_data), Y.relabeled)]
  sample.weights <- 1 / C.Y.hat
  
  XLfit1 <- gbm(Y ~ ., 
                distribution = "gaussian",
                data = data.frame(Y = binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                                  binary_data[binary_data$W==1,4:8]),
                weights = sample.weights[binary_data$W==1])
  XLtau1 <- predict(XLfit1, data.frame(data.test$X))
  
  XLfit0 <- gbm(Y ~ ., 
                distribution = "gaussian",
                data = data.frame(Y = binary_data$Tgbm1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                                  binary_data[binary_data$W==0,4:8]),
                weights = sample.weights[binary_data$W==0])
  XLtau0 <- predict(XLfit0, data.frame(data.test$X))
  
  # propensity score 
  ps <- regression_forest(data.test$X, data.test$W, num.trees = 500)
  cate <- XLtau1 * ps$predictions + XLtau0 * (1 - ps$predictions)
  cate
}











