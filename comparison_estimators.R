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
estimate_rfsrc_X_W = function(data, data.test, times = NULL) {
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
estimate_rfsrc_XW_W = function(data, data.test, times = NULL) {
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
estimate_rfsrc_twin = function(data, data.test, times = NULL) {
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
estimate_grf = function(data, data.test, times = NULL) {
  fit = causal_survival_forest(data$X, data$Y, data$W, data$D)
  p = predict(fit, data.test$X)

  p$predictions
}

# "IPCW" - causal forest
estimate_IPCW_grf = function(data, data.test, times = NULL) {
  sf.censor = survival_forest(cbind(data$W, data$X),
                              data$Y,
                              1 - data$D,
                              prediction.type = "Nelson-Aalen")
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

setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/rlearner")
files.sources = list.files()
sapply(files.sources, source)

setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/flearner")
files.sources = list.files()
sapply(files.sources, source)

base_surv <- function(fit, Y, D, x, lambda){
  data <- data.frame(t_event=Y, event=D, x)
  tab <- data.frame(table(data[data$event == 1, "t_event"])) 
  y <- as.numeric(as.character(sort(unique(tab[,1]))))
  d <- tab[,2]  # number of events at each unique time                               
  
  betaHat <- as.vector((fit$glmnet.fit$beta)[,fit$lambda==lambda])
  h0 <- rep(NA, length(y))
  for(l in 1:length(y)){
    h0[l] <- d[l] / sum(exp(x[data$t_event >= y[l], rownames(fit$glmnet.fit$beta)] %*% betaHat))
  }    
  
  S0 <- exp(-cumsum(h0))
  outcome <- data.frame(time=y,survival=S0)
  outcome
}
pred_surv <- function(fit, S0, x, times, lambda){
  link <- predict(fit$glmnet.fit,x,type = "link")[,fit$lambda==lambda]
  colnames(link) <- NULL
  
  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }
  
  surv <- S0_t^exp(link)
  surv
}
pred_surv_preval <- function(fit, S0, times, lambda){
  link <- fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == lambda] 
  colnames(link) <- NULL
  
  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }
  
  surv <- S0_t^exp(link)
  surv
}

# S-learner 
estimate_coxph_slearner <- function(data, data.test, times){
  
  scoxph_fit <- scoxph(x = data$X,
                       w = data$W,
                       y = data$Y, 
                       D = data$D, 
                       times = times)
  
  pred_S_coxph <- predict(scoxph_fit, newx = data.test$X, times = times)
}

estimate_lasso_slearner <- function(data, data.test, nfolds = 10, alpha = 1, times = NULL){
  
  slasso_fit <- slasso_surv(x = data$X,
                            w = data$W,
                            y = data$Y, 
                            D = data$D, 
                            times = times)
  
  pred_S_lasso <- predict(slasso_fit, newx = data.test$X, times = times)
  pred_S_lasso
}

estimate_gbm_slearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- data.frame(data.test$Y, data.test$D, data.test$W, data.test$X)
  colnames(testdat)[1:3] <- c("Y", "D", "W")
  testdat1 <- testdat; testdat1$W <- 1
  testdat0 <- testdat; testdat0$W <- 0
  
  gbmfit <- gbm(Surv(Y, D)~., 
                data = traindat, 
                distribution = "coxph",
                shrinkage = shrinkage, 
                n.trees = n.trees,   
                cv.folds = cv.folds, 
                n.cores = detectCores()) 
  best_iter <- gbm.perf(gbmfit, plot.it = F, method = 'cv')
  
  # return a vector of prediction on n.trees indicating log hazard scale.f(x)
  pred.train <- predict(gbmfit, traindat, n.trees = best_iter)
  pred.test1 <- predict(gbmfit, testdat1, n.trees = best_iter)
  pred.test0 <- predict(gbmfit, testdat0, n.trees = best_iter)
  
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum <- basehaz.gbm(traindat$Y, traindat$D, pred.train, t.eval = times, cumulative = TRUE)
  
  # S(X) at median survival time 
  surf1 <- surf0 <- rep(NA,dim(testdat1)[1])
  for (j in 1:dim(testdat1)[1]){
    surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum)
    surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum)
  }
  
  pred_S_gbm <- surf1 - surf0
  pred_S_gbm
}

estimate_grf_slearner <- function(data, data.test, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- data.frame(data.test$Y, data.test$D, data.test$W, data.test$X)
  colnames(testdat)[1:3] <- c("Y", "D", "W")
  testdat1 <- testdat; testdat1$W <- 1
  testdat0 <- testdat; testdat0$W <- 0
  
  grffit <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                            testdat$Y,
                            testdat$D)
  surf1 <- predict(grffit, as.matrix(testdat1[,3:dim(testdat1)[2]]))$predictions[, which.min(abs(grffit$failure.times-times))]
  surf0 <- predict(grffit, as.matrix(testdat0[,3:dim(testdat0)[2]]))$predictions[, which.min(abs(grffit$failure.times-times))]
  pred_S_grf <- surf1 - surf0
  pred_S_grf
}


# T-learner 
estimate_coxph_tlearner <- function(data, data.test, times){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  
  # Model for W = 1
  coxph_fit1 <- coxph(Surv(Y, D) ~., data = traindat1)
  bh_dat <- basehaz(coxph_fit1)
  bh <- bh_dat[which.min(abs(bh_dat$time - times)),]$hazard
  est_r1 <- predict(coxph_fit1, newdata = data.frame(data.test$X), type="risk")
  surf1 <- exp(-bh)^est_r1
  
  # Model for W = 0
  coxph_fit0 <- coxph(Surv(Y, D) ~., data = traindat0)
  bh_dat <- basehaz(coxph_fit0)
  bh <- bh_dat[which.min(abs(bh_dat$time - times)),]$hazard
  est_r0 <- predict(coxph_fit0, newdata = data.frame(data.test$X), type="risk")
  surf0 <- exp(-bh)^est_r0
  
  pred_T_coxph <- surf1 - surf0
  pred_T_coxph
}

estimate_lasso_tlearner <- function(data, data.test, nfolds = 10, alpha = 1, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  
  # Model for W = 1
  foldid <- sample(rep(seq(nfolds), length = length(traindat1$Y)))
  lasso_fit1 <- cv.glmnet(as.matrix(traindat1[,3:dim(traindat1)[2]]),
                          Surv(traindat1$Y, traindat1$D),
                          family = "cox",
                          alpha = alpha,
                          foldid = foldid)
  
  S0 <- base_surv(fit = lasso_fit1, 
                  Y = traindat1$Y,
                  D = traindat1$D,
                  x = as.matrix(traindat1[,3:dim(traindat1)[2]]),
                  lambda = lasso_fit1$lambda.min)
  
  surf1 <- pred_surv(fit = lasso_fit1,
                     S0 = S0,
                     x = data.test$X,
                     times = times,
                     lambda = lasso_fit1$lambda.min)
  
  # Model for W = 0
  foldid <- sample(rep(seq(nfolds), length = length(traindat0$Y)))
  lasso_fit0 <- cv.glmnet(as.matrix(traindat0[,3:dim(traindat0)[2]]),
                          Surv(traindat0$Y, traindat0$D),
                          family = "cox",
                          alpha = alpha,
                          foldid = foldid)
  
  S0 <- base_surv(fit = lasso_fit0, 
                  Y = traindat0$Y,
                  D = traindat0$D,
                  x = as.matrix(traindat0[,3:dim(traindat0)[2]]),
                  lambda = lasso_fit0$lambda.min)
  
  surf0 <- pred_surv(fit = lasso_fit0,
                     S0 = S0,
                     x = data.test$X,
                     times = times,
                     lambda = lasso_fit0$lambda.min)
  
  pred_T_lasso <- surf1 - surf0
  pred_T_lasso
}

estimate_gbm_tlearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  
  # Model for W = 1
  gbm_fit1 <- gbm(Surv(Y, D)~.,
                  data = traindat1,
                  distribution = "coxph",
                  shrinkage = shrinkage, 
                  n.trees = n.trees,   
                  cv.folds = cv.folds, 
                  n.cores = detectCores())
  best_iter <- gbm.perf(gbm_fit1, plot.it = F, method = 'cv')
  
  # return a vector of prediction on n.trees indicating log hazard scale.f(x)
  pred.train1 <- predict(gbm_fit1, traindat1, n.trees = best_iter)
  pred.test1  <- predict(gbm_fit1, data.frame(data.test$X), n.trees = best_iter)
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum1 <- basehaz.gbm(traindat1$Y, traindat1$D, pred.train1, t.eval = times, cumulative = TRUE)
  
  # Model for W = 0
  gbm_fit0 <- gbm(Surv(Y, D)~.,
                  data = traindat0,
                  distribution = "coxph",
                  shrinkage = shrinkage, 
                  n.trees = n.trees,   
                  cv.folds = cv.folds)
  best_iter <- gbm.perf(gbm_fit0, plot.it = F, method = 'cv')
  
  # return a vector of prediction on n.trees indicating log hazard scale.f(x)
  pred.train0 <- predict(gbm_fit0, traindat0, n.trees = best_iter)
  pred.test0  <- predict(gbm_fit0, data.frame(data.test$X), n.trees = best_iter)
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum0 <- basehaz.gbm(traindat0$Y, traindat0$D, pred.train0, t.eval = times, cumulative = TRUE)
  
  surf1 <- surf0 <- rep(NA,dim(data.frame(data.test$X))[1])
  for (j in 1:dim(data.frame(data.test$X))[1]){
    surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
    surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
  }
  
  pred_T_gbm <- surf1 - surf0
  pred_T_gbm
}

estimate_grf_tlearner <- function(data, data.test, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  traindat1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  traindat0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  
  # Model for W = 1
  grffit1 <- survival_forest(as.matrix(traindat1[,3:dim(traindat1)[2]]),
                            traindat1$Y,
                            traindat1$D)
  surf1 <- predict(grffit1, data.test$X)$predictions[, which.min(abs(grffit1$failure.times-times))]

  # Model for W = 0
  grffit0 <- survival_forest(as.matrix(traindat0[,3:dim(traindat0)[2]]),
                            traindat0$Y,
                            traindat0$D)
  surf0 <- predict(grffit0, data.test$X)$predictions[, which.min(abs(grffit0$failure.times-times))]
  
  pred_T_gbm <- surf1 - surf0
  pred_T_gbm
}


# "IPCW" - X-learner gradient boosting
estimate_ipcw_wocf_lasso_xlearner <- function(data, data.test, nfolds = 10, alpha = 1, times = NULL){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]   
  
  # fit model on W==1
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  foldid <- sample(rep(seq(nfolds), length = length(train1$Y)))
  lasso_fit1 <- cv.glmnet(as.matrix(train1[,3:dim(train1)[2]]),
                          Surv(train1$Y, train1$D),
                          family = "cox",
                          alpha = alpha,
                          foldid = foldid)
  
  S0 <- base_surv(fit = lasso_fit1, 
                  Y = train1$Y,
                  D = train1$D,
                  x = as.matrix(train1[,3:dim(train1)[2]]),
                  lambda = lasso_fit1$lambda.min)
  
  surf1 <- pred_surv(fit = lasso_fit1,
                     S0 = S0,
                     x = as.matrix(testdat[,3:dim(testdat)[2]]),
                     times = times,
                     lambda = lasso_fit1$lambda.min)
  
  # fit model on W==0
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  foldid <- sample(rep(seq(nfolds), length = length(train0$Y)))
  lasso_fit0 <- cv.glmnet(as.matrix(train0[,3:dim(train0)[2]]),
                          Surv(train0$Y, train0$D),
                          family = "cox",
                          alpha = alpha,
                          foldid = foldid)
  
  S0 <- base_surv(fit = lasso_fit0, 
                  Y = train0$Y,
                  D = train0$D,
                  x = as.matrix(train0[,3:dim(train0)[2]]),
                  lambda = lasso_fit0$lambda.min)
  
  surf0 <- pred_surv(fit = lasso_fit0,
                     S0 = S0,
                     x = as.matrix(testdat[,3:dim(testdat)[2]]),
                     times = times,
                     lambda = lasso_fit0$lambda.min)
  
  Tgbm1 <- 1-surf1
  Tgbm0 <- 1-surf0
  
  # IPCW weights
  foldid <- sample(rep(seq(nfolds), length = length(traindat$Y)))
  c_fit <- cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]), 
                     Surv(traindat$Y, 1-traindat$D),
                     family = "cox",
                     foldid = foldid,
                     alpha = alpha)
  S0 <- base_surv(c_fit, traindat$Y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_fit$lambda.min)
  cent <- rep(times, length(traindat$Y))
  cent[traindat$D==1] <- traindat$Y[traindat$D==1]
  c_hat <- pred_surv(c_fit, S0, as.matrix(traindat[,3:dim(traindat)[2]]), times = cent, lambda = c_fit$lambda.min)
  
  # Propensity score
  w_fit <- cv.glmnet(data$X, 
                     data$W,
                     foldid = foldid,
                     alpha = 1)
  ps <- as.vector(predict(w_fit, newx = data$X, s = w_fit$lambda.min))
  
  weight <- (1 / c_hat) * (1 / ps)
  
  # X-learner
  tempdat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tgbm0, Tgbm1)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          # remove subjects who got censored before the time of interest times
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     # recode the event status for subjects who had events after times
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  foldid <- sample(rep(seq(nfolds), length = length(binary_data[binary_data$W==1,]$Y)))
  XLfit1 <- cv.glmnet(as.matrix(binary_data[binary_data$W==1,5:9]), 
                      binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                      weights = binary_data$weight[binary_data$W==1],
                      foldid = foldid,
                      alpha = alpha)
  XLtau1 <- as.vector(-predict(XLfit1, data.test$X, s = XLfit1$lambda.min))
  
  foldid <- sample(rep(seq(nfolds), length = length(binary_data[binary_data$W==0,]$Y)))
  XLfit0 <- cv.glmnet(as.matrix(binary_data[binary_data$W==0,5:9]), 
                      binary_data$Tgbm1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                      weights = binary_data$weight[binary_data$W==0],
                      foldid = foldid,
                      alpha = alpha)
  XLtau0 <- as.vector(-predict(XLfit0, data.test$X, s = XLfit0$lambda.min))
  
  # propensity score
  ps.test <- as.vector(predict(w_fit, newx = data.test$X, s = w_fit$lambda.min))
  pred_X_lasso <- XLtau1 * ps.test + XLtau0 * (1 - ps.test)
  as.vector(pred_X_lasso)
}

estimate_ipcw_wcf_lasso_xlearner <- function(data, data.test, nfolds = 10, alpha = 1, times = NULL){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]   
  
  # fit model on W==1 (cross-fitted using 'preval' in glmnet)
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  foldid <- sample(rep(seq(nfolds), length = length(train1$Y)))
  lasso_fit1 <- cv.glmnet(as.matrix(train1[,3:dim(train1)[2]]),
                          Surv(train1$Y, train1$D),
                          family = "cox",
                          alpha = alpha,
                          keep = TRUE,
                          foldid = foldid)

  lambda_1_min <- lasso_fit1$lambda[which.min(lasso_fit1$cvm[!is.na(colSums(lasso_fit1$fit.preval))])]
  S0 <- base_surv(lasso_fit1, train1$Y, train1$D, as.matrix(train1[, 3:dim(train1)[2]]), lambda = lambda_1_min)
  surf1 <- rep(NA, length(traindat$Y))
  surf1[traindat$W==1] <- pred_surv_preval(lasso_fit1, S0, times = times, lambda = lambda_1_min)
  surf1[traindat$W==0] <- pred_surv(fit = lasso_fit1,
                                    S0 = S0,
                                    x = as.matrix(traindat[traindat$W==0, 4:dim(traindat)[2]]),
                                    times = times,
                                    lambda = lasso_fit1$lambda.min)

  # fit model on W==0 (cross-fitted using 'preval' in glmnet)
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  foldid <- sample(rep(seq(nfolds), length = length(train0$Y)))
  lasso_fit0 <- cv.glmnet(as.matrix(train0[,3:dim(train0)[2]]),
                          Surv(train0$Y, train0$D),
                          family = "cox",
                          alpha = alpha,
                          keep = TRUE,
                          foldid = foldid)
  
  lambda_0_min <- lasso_fit0$lambda[which.min(lasso_fit0$cvm[!is.na(colSums(lasso_fit0$fit.preval))])]
  S0 <- base_surv(lasso_fit0, train0$Y, train0$D, as.matrix(train0[, 3:dim(train0)[2]]), lambda = lambda_0_min)
  surf0 <- rep(NA, length(traindat$Y))
  surf0[traindat$W==0] <- pred_surv_preval(lasso_fit0, S0, times = times, lambda = lambda_0_min)
  surf0[traindat$W==1] <- pred_surv(fit = lasso_fit0,
                                    S0 = S0,
                                    x = as.matrix(traindat[traindat$W==1, 4:dim(traindat)[2]]),
                                    times = times,
                                    lambda = lasso_fit0$lambda.min)
  
  Tgbm1 <- 1-surf1
  Tgbm0 <- 1-surf0
  
  # IPCW weights (cross-fitted using 'preval' in glmnet)
  foldid <- sample(rep(seq(nfolds), length = length(traindat$Y)))
  c_fit <- cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]), 
                     Surv(traindat$Y, 1-traindat$D),
                     family = "cox",
                     foldid = foldid,
                     keep = TRUE,
                     alpha = alpha)
  
  c_lambda_min <- c_fit$lambda[which.min(c_fit$cvm[!is.na(colSums(c_fit$fit.preval))])]
  S0 <- base_surv(c_fit, traindat$Y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_lambda_min)
  cent <- rep(times, length(traindat$Y))
  cent[traindat$D==1] <- traindat$Y[traindat$D==1]
  c_hat <- pred_surv_preval(c_fit, S0, times = cent, lambda = c_lambda_min)
  
  # Propensity score (cross-fitted using 'preval' in glmnet)
  w_fit <- cv.glmnet(as.matrix(traindat[, 4:dim(traindat)[2]]),
                     traindat$W,
                     foldid = foldid,
                     keep = TRUE,
                     alpha = alpha)
  
  w_lambda_min <- w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
  ps <- w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
  
  weight <- (1 / c_hat)*(1 / ps)
  
  # X-learner
  tempdat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tgbm0, Tgbm1)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]         
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  foldid <- sample(rep(seq(nfolds), length = length(binary_data[binary_data$W==1,]$Y)))
  XLfit1 <- cv.glmnet(as.matrix(binary_data[binary_data$W==1, 5:9]), 
                      binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                      weights = binary_data$weight[binary_data$W==1],
                      foldid = foldid,
                      alpha = alpha)
  XLtau1 <- as.vector(-predict(XLfit1, data.test$X, s = XLfit1$lambda.min))
  
  foldid <- sample(rep(seq(nfolds), length = length(binary_data[binary_data$W==0,]$Y)))
  XLfit0 <- cv.glmnet(as.matrix(binary_data[binary_data$W==0, 5:9]), 
                      binary_data$Tgbm1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                      weights = binary_data$weight[binary_data$W==0],
                      foldid = foldid,
                      alpha = alpha)
  XLtau0 <- as.vector(-predict(XLfit0, data.test$X, s = XLfit0$lambda.min))
  
  # propensity score
  ps.test <- as.vector(predict(w_fit, newx = data.test$X, s = w_fit$lambda.min))
  pred_X_lasso <- XLtau1 * ps.test + XLtau0 * (1 - ps.test)
  as.vector(pred_X_lasso)
}

estimate_ipcw_wocf_gbm_xlearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times = NULL){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]   
  
  # fit model on W==1
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  gbm_fit1 <- gbm(Surv(Y, D) ~ .,
                  data = train1,
                  distribution = "coxph",
                  shrinkage = shrinkage, 
                  n.trees = n.trees,   
                  cv.folds = cv.folds, 
                  n.cores = detectCores())
  best.iter1 <- gbm.perf(gbm_fit1, plot.it = F, method = 'cv')
  pred.train1 <- predict(gbm_fit1, train1, n.trees = best.iter1)
  pred.test1  <- predict(gbm_fit1, testdat, n.trees = best.iter1)
  basehaz.cum1 <- basehaz.gbm(train1$Y, train1$D, pred.train1, t.eval = times, cumulative = TRUE)
  
  # fit model on W==0
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  gbm_fit0 <- gbm(Surv(Y, D) ~ .,
                  data = train0,
                  distribution = "coxph",
                  shrinkage = shrinkage, 
                  n.trees = n.trees,   
                  cv.folds = cv.folds, 
                  n.cores = detectCores())
  best.iter0 <- gbm.perf(gbm_fit0, plot.it = F, method = 'cv')
  pred.train0 <- predict(gbm_fit0, train0, n.trees = best.iter0)
  pred.test0  <- predict(gbm_fit0, testdat, n.trees = best.iter0)
  basehaz.cum0 <- basehaz.gbm(train0$Y, train0$D, pred.train0, t.eval = times, cumulative = TRUE)
  
  # T-learner estimates of 1-S(t) at the median survival time
  surf1 <- surf0 <- rep(NA,dim(testdat)[1])
  for (j in 1:dim(testdat)[1]){
    surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
    surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
  }
  Tgbm1 <- 1-surf1
  Tgbm0 <- 1-surf0
  
  # IPCW weights
  c_fit <- gbm(Surv(Y, 1-D)~.,
               data = traindat, 
               distribution = "coxph",
               shrinkage = shrinkage, 
               n.trees = n.trees,   
               cv.folds = cv.folds, 
               n.cores = detectCores()) 
  best_iter <- gbm.perf(c_fit, plot.it = F, method = 'cv')
  
  # baseline hazard
  pred.train <- predict(c_fit, traindat, n.trees = best_iter)
  cent <- traindat$Y
  cent[traindat$D==0] <- times
  basehaz.cum <- rep(NA, length(cent))
  for (z in 1:length(cent)){
    basehaz.cum[z] <- basehaz.gbm(traindat$Y, 1-traindat$D, pred.train, t.eval = cent[z], cumulative = TRUE)
  }
  pred.test <- predict(c_fit, traindat[, 3:dim(traindat)[2]], n.trees = best_iter)
  c_hat <- exp(-exp(pred.test)*basehaz.cum)
  
  # Propensity score 
  tmpdat <- traindat[, 3:dim(traindat)[2]]
  w_fit <- gbm(W ~ ., 
               data = tmpdat, 
               shrinkage = shrinkage, 
               n.trees = n.trees,   
               cv.folds = cv.folds, 
               n.cores = detectCores()) 
  best_iter_ps <- gbm.perf(w_fit, plot.it = F, method = 'cv')
  log_odds <- predict(w_fit, tmpdat[, 2:dim(tmpdat)[2]], n.trees = best_iter_ps)
  ps <- 1/(1 + exp(-log_odds))
  
  weight <- (1 / c_hat) * (1 / ps)
  
  # X-learner
  tempdat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tgbm0, Tgbm1)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  XLfit1 <- gbm(Y ~ .,
                distribution = "gaussian",
                data = data.frame(Y = binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                                  binary_data[binary_data$W==1, 5:9]),
                weights = binary_data$weight[binary_data$W==1],
                shrinkage = shrinkage, 
                n.trees = n.trees,   
                cv.folds = cv.folds, 
                n.cores = detectCores())
  best_iter_tau1 <- gbm.perf(XLfit1, plot.it = F, method = 'cv')
  XLtau1 <- -predict(XLfit1, data.frame(data.test$X), n.trees = best_iter_tau1)
  
  XLfit0 <- gbm(Y ~ .,
                distribution = "gaussian",
                data = data.frame(Y = binary_data$Tgbm1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                                  binary_data[binary_data$W==0, 5:9]),
                weights = binary_data$weight[binary_data$W==0],
                shrinkage = shrinkage, 
                n.trees = n.trees,   
                cv.folds = cv.folds, 
                n.cores = detectCores())
  best_iter_tau0 <- gbm.perf(XLfit0, plot.it = F, method = 'cv')
  XLtau0 <- -predict(XLfit0, data.frame(data.test$X), n.trees=best_iter_tau0)
  
  # propensity score
  log_odds <- predict(w_fit, data.frame(data.test$X), n.trees = best_iter_ps)
  ps_score <- 1/(1 + exp(-log_odds))
  pred_X_gbm <- XLtau1 * ps_score + XLtau0 * (1 - ps_score)
  pred_X_gbm
}

estimate_ipcw_wcf_gbm_xlearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times = NULL){
  
  tempdat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  foldid <- sample(rep(seq(nfolds), length = length(tempdat$W)))
  Tgbm1 <- Tgbm0 <- rep(NA, length(tempdat$W))
  
  # fit model on W==1
  for (k in 1:nfolds){
    testdat  <- tempdat[foldid==k, !colnames(tempdat) %in% c("W")]
    traindat <- tempdat[!foldid==k, ]
    train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
    
    gbm_fit1 <- gbm(Surv(Y, D) ~ .,
                    data = train1,
                    distribution = "coxph",
                    shrinkage = shrinkage, 
                    n.trees = n.trees,   
                    cv.folds = cv.folds, 
                    n.cores = detectCores())
    best.iter1 <- gbm.perf(gbm_fit1, plot.it = F, method = 'cv')
    pred.train1 <- predict(gbm_fit1, train1, n.trees = best.iter1)
    pred.test1  <- predict(gbm_fit1, testdat, n.trees = best.iter1)
    basehaz.cum1 <- basehaz.gbm(train1$Y, train1$D, pred.train1, t.eval = times, cumulative = TRUE)
    surf1 <- rep(NA,dim(testdat)[1])
    for (j in 1:dim(testdat)[1]){
      surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
    }
    
    Tgbm1[foldid==k] <- 1 - surf1
  }
  
  # fit model on W==0
  for (k in 1:nfolds){
    testdat  <- tempdat[foldid==k, !colnames(tempdat) %in% c("W")]
    traindat <- tempdat[!foldid==k, ]
    train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
    
    gbm_fit0 <- gbm(Surv(Y, D) ~ .,
                    data = train0,
                    distribution = "coxph",
                    shrinkage = shrinkage, 
                    n.trees = n.trees,   
                    cv.folds = cv.folds, 
                    n.cores = detectCores())
    best.iter0 <- gbm.perf(gbm_fit0, plot.it = F, method = 'cv')
    pred.train0 <- predict(gbm_fit0, train0, n.trees = best.iter0)
    pred.test0  <- predict(gbm_fit0, testdat, n.trees = best.iter0)
    basehaz.cum0 <- basehaz.gbm(train0$Y, train0$D, pred.train0, t.eval = times, cumulative = TRUE)
    surf0 <- rep(NA,dim(testdat)[1])
    for (j in 1:dim(testdat)[1]){
      surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
    }
    
    Tgbm0[foldid==k] <- 1 - surf0
  }  
  
  # IPCW weights
  c_hat <- rep(NA, length(tempdat$W))
  for (k in 1:nfolds){
    testdat <- tempdat[foldid==k, ]
    traindat <- tempdat[!foldid==k, ]
    c_fit <- gbm(Surv(Y, 1-D)~.,
                 data = traindat, 
                 distribution = "coxph",
                 shrinkage = shrinkage, 
                 n.trees = n.trees,   
                 cv.folds = cv.folds, 
                 n.cores = detectCores()) 
    best_iter <- gbm.perf(c_fit, plot.it = F, method = 'cv')
    
    # baseline hazard
    pred.train <- predict(c_fit, traindat, n.trees = best_iter)
    cent <- testdat$Y
    cent[testdat$D==0] <- times
    basehaz.cum <- rep(NA, length(cent))
    for (z in 1:length(cent)){
      basehaz.cum[z] <- basehaz.gbm(traindat$Y, 1-traindat$D, pred.train, t.eval = cent[z], cumulative = TRUE)
    }
    pred.test <- predict(c_fit, testdat[,3:dim(testdat)[2]], n.trees = best_iter)
    c_hat[foldid==k] <- exp(-exp(pred.test)*basehaz.cum)
  }
  
  # Propensity score
  ps <- rep(NA, length(tempdat$W))
  for (k in 1:nfolds){
    testdat <- tempdat[foldid==k, ]
    traindat <- tempdat[!foldid==k, ]
    w_dat <- traindat[, 3:dim(traindat)[2]]
    w_fit <- gbm(W ~ ., 
                 data = w_dat, 
                 distribution = "bernoulli",
                 shrinkage = shrinkage, 
                 n.trees = n.trees,   
                 cv.folds = cv.folds, 
                 n.cores = detectCores()) 
    best_iter <- gbm.perf(w_fit, plot.it = F, method = 'cv')
    log_odds <- predict(w_fit, testdat[, 4:dim(testdat)[2]], n.trees = best_iter)
    ps[foldid==k] <- 1/(1 + exp(-log_odds))
  }
  
  weight <- (1 / c_hat) * (1 / ps)
  
  # X-learner
  xldat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tgbm0, Tgbm1)
  colnames(xldat)[1:3] <- c("Y", "D", "W")
  binary_data <- xldat[xldat$D==1|xldat$Y > times,]         
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  XLfit1 <- gbm(Y ~ .,
                distribution = "gaussian",
                data = data.frame(Y = binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                                  binary_data[binary_data$W==1,5:9]),
                weights = binary_data$weight[binary_data$W==1],
                shrinkage = shrinkage, 
                n.trees = n.trees,   
                cv.folds = cv.folds, 
                n.cores = detectCores())
  best_iter_tau1 <- gbm.perf(XLfit1, plot.it = F, method = 'cv')
  XLtau1 <- -predict(XLfit1, data.frame(data.test$X), n.trees = best_iter_tau1)
  
  XLfit0 <- gbm(Y ~ .,
                distribution = "gaussian",
                data = data.frame(Y = binary_data$Tgbm1[binary_data$W==0] - binary_data$D[binary_data$W==0],
                                  binary_data[binary_data$W==0,5:9]),
                weights = binary_data$weight[binary_data$W==0],
                shrinkage = shrinkage, 
                n.trees = n.trees,   
                cv.folds = cv.folds, 
                n.cores = detectCores())
  best_iter_tau0 <- gbm.perf(XLfit0, plot.it = F, method = 'cv')
  XLtau0 <- -predict(XLfit0, data.frame(data.test$X), n.trees = best_iter_tau0)
  
  # propensity score
  tmpdat <- tempdat[, 3:dim(tempdat)[2]]
  w_fit_test <- gbm(W ~ ., 
                    distribution = "bernoulli",
                    data = tmpdat, 
                    shrinkage = shrinkage, 
                    n.trees = n.trees,   
                    cv.folds = cv.folds, 
                    n.cores = detectCores()) 
  best_iter_ps <- gbm.perf(w_fit_test, plot.it = F, method = 'cv')
  log_odds <- predict(w_fit_test, data.frame(data.test$X), n.trees = best_iter_ps)
  ps_score <- 1/(1 + exp(-log_odds))
  pred_X_gbm <- XLtau1 * ps_score + XLtau0 * (1 - ps_score)
  pred_X_gbm
}

estimate_ipcw_grf_xlearner <- function(data, data.test, times = NULL){
  
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]   
  
  # fit model on W==1
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
  grffit1 <- survival_forest(as.matrix(train1[,3:dim(train1)[2]]),
                            train1$Y,
                            train1$D)
  surf1 <- rep(NA, length(traindat$Y))
  surf1[traindat$W==1] <- grffit1$predictions[, which.min(abs(grffit1$failure.times-times))]
  surf1[traindat$W==0]<- predict(grffit1, as.matrix(traindat[traindat$W==0, 4:dim(traindat)[2]]))$predictions[, which.min(abs(grffit1$failure.times-times))]
  
  # fit model on W==0
  train0 <- traindat[traindat$W==0, !colnames(traindat) %in% c("W")]
  grffit0 <- survival_forest(as.matrix(train0[,3:dim(train0)[2]]),
                             train0$Y,
                             train0$D)
  surf0 <- rep(NA, length(traindat$Y))
  surf0[traindat$W==0] <- grffit0$predictions[, which.min(abs(grffit0$failure.times-times))]
  surf0[traindat$W==1] <- predict(grffit0, as.matrix(traindat[traindat$W==1, 4:dim(traindat)[2]]))$predictions[, which.min(abs(grffit0$failure.times-times))]
  
  Tgbm1 <- 1-surf1
  Tgbm0 <- 1-surf0
  
  # IPCW weights
  sf.censor <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                               traindat$Y,
                               1 - traindat$D,
                               prediction.type = "Nelson-Aalen")
  C.hat <- sf.censor$predictions
  cent <- traindat$Y
  cent[traindat$D==0] <- times
  C.index <- rep(NA, length(cent))
  for (h in 1:length(cent)){
    C.index[h] <- which.min(abs(sort(unique(traindat$Y[traindat$D==0]))-cent[h]))
  }
  C.Y.hat <- C.hat[cbind(1:length(traindat$Y), C.index)]
  ipcw <- 1 / C.Y.hat
  
  # Propensity score 
  tmpdat <- traindat[, 3:dim(traindat)[2]]
  psfit <- regression_forest(as.matrix(tmpdat[,2:dim(tmpdat)[2]]), 
                             tmpdat$W, 
                             tune.parameters = "all") 
  ps <- psfit$predictions
  
  weight <- ipcw/ps  # censoring weight * treatment weight
    
  # X-learner
  tempdat <- data.frame(data$Y, data$D, data$W, weight, data$X, Tgbm0, Tgbm1)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  XLfit1 <- regression_forest(as.matrix(binary_data[binary_data$W==1, 5:9]),
                              binary_data$D[binary_data$W==1] - binary_data$Tgbm0[binary_data$W==1],
                              sample.weights = binary_data$weight[binary_data$W==1], 
                              tune.parameters = "all")
  XLtau1 <- -predict(XLfit1, data.frame(data.test$X))
  
  XLfit0 <- regression_forest(as.matrix(binary_data[binary_data$W==0, 5:9]),
                              binary_data$D[binary_data$W==0] - binary_data$Tgbm0[binary_data$W==0],
                              sample.weights = binary_data$weight[binary_data$W==0], 
                              tune.parameters = "all")
  XLtau0 <- -predict(XLfit0, data.frame(data.test$X))

  ps_test <- predict(psfit, data.test$X)$predictions
  pred_X_grf <- XLtau1 * ps_test + XLtau0 * (1 - ps_test)
  pred_X_grf
}


# "IPCW" R-learner 
estimate_ipcw_wocf_lasso_rlearner <- function(data, data.test, nfolds = 10, alpha = 1, times=NULL){
  rlasso_fit <- rlasso(x = data$X, 
                       w = data$W, 
                       y = data$Y, 
                       D = data$D,
                       alpha = alpha,
                       k_folds = nfolds,
                       cf = FALSE,
                       times = times)
  rlasso_est <- predict(object = rlasso_fit, newx = as.matrix(data.test$X)) 
  as.vector(rlasso_est)
}

estimate_ipcw_wcf_lasso_rlearner <- function(data, data.test, nfolds = 10, alpha = 1, times=NULL){
  rlasso_fit <- rlasso(x = data$X, 
                       w = data$W, 
                       y = data$Y, 
                       D = data$D, 
                       alpha = alpha,
                       k_folds = nfolds,
                       times = times)
  rlasso_est <- predict(object = rlasso_fit, newx = as.matrix(data.test$X)) 
  as.vector(rlasso_est)
}

estimate_ipcw_wocf_gbm_rlearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times=NULL){
  rgbm_fit <- rgbm(x = data$X, 
                   w = data$W, 
                   y = data$Y, 
                   D = data$D,
                   shrinkage = shrinkage, 
                   n.trees = n.trees, 
                   cv.folds = cv.folds,
                   cf = FALSE,
                   times = times)
  rgbm_est <- predict(object = rgbm_fit, newx = data.test$X) 
  rgbm_est
}

estimate_ipcw_wcf_gbm_rlearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times=NULL){
  rgbm_fit <- rgbm(x = data$X, 
                   w = data$W, 
                   y = data$Y, 
                   D = data$D, 
                   shrinkage = shrinkage, 
                   n.trees = n.trees, 
                   cv.folds = cv.folds,
                   times = times)
  rgbm_est <- predict(object = rgbm_fit, newx = data.test$X) 
  rgbm_est
}

estimate_ipcw_grf_rlearner <- function(data, data.test, times=NULL){
  rgrf_fit <- rgrf(x = data$X, 
                   w = data$W, 
                   y = data$Y, 
                   D = data$D, 
                   times = times)
  rgrf_est <- predict(object = rgrf_fit, data.test$X) 
  rgrf_est
}


# "IPCW" F-learner
estimate_ipcw_wocf_lasso_flearner <- function(data, data.test, nfolds = 10, alpha = 1, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  foldid <- sample(rep(seq(nfolds), length = length(data$Y)))
  
  # IPCW weights
  c_fit <- cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]), 
                     Surv(traindat$Y, 1 - traindat$D),
                     family = "cox",
                     foldid = foldid,
                     alpha = alpha)
  
  S0 <- base_surv(c_fit, traindat$Y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_fit$lambda.min)
  cent <- rep(times, length(traindat$Y))
  cent[traindat$D==1] <- traindat$Y[traindat$D==1]
  c_hat <- pred_surv(c_fit, S0, as.matrix(traindat[,3:dim(traindat)[2]]), times = cent, lambda = c_fit$lambda.min)
  ipcw <- 1 / c_hat
  
  # Propensity score
  w_fit <- cv.glmnet(as.matrix(traindat[,4:dim(traindat)[2]]), 
                     traindat$W,
                     foldid = foldid,
                     alpha = alpha)
  
  ps_score <- as.vector(predict(w_fit, newx = as.matrix(traindat[,4:dim(traindat)[2]]), s = w_fit$lambda.min))
  
  # Subset of uncensored subjects
  tempdat <- data.frame(data$Y, data$D, data$W, ps_score, ipcw, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]         
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0    
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  flasso_fit <- Flasso(x = as.matrix(binary_data[,6:10]), 
                       tx = binary_data$W, 
                       y = binary_data$D, 
                       pscore = binary_data$ps_score, 
                       weight = binary_data$ipcw,
                       alpha = alpha,
                       nfolds = nfolds)
  pred_flasso <- as.vector(-predict(flasso_fit, data.test$X))
  pred_flasso
}

estimate_ipcw_wcf_lasso_flearner <- function(data, data.test, nfolds = 10, alpha = 1, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  foldid <- sample(rep(seq(nfolds), length = length(data$Y)))
  
  # IPCW weights
  c_fit <- cv.glmnet(as.matrix(traindat[,3:dim(traindat)[2]]), 
                     Surv(traindat$Y, 1-traindat$D),
                     family = "cox",
                     foldid = foldid,
                     keep = TRUE,
                     alpha = alpha)
  
  c_lambda_min <- c_fit$lambda[which.min(c_fit$cvm[!is.na(colSums(c_fit$fit.preval))])]
  S0 <- base_surv(c_fit, traindat$Y, 1-traindat$D, as.matrix(traindat[,3:dim(traindat)[2]]), lambda = c_lambda_min)
  cent <- rep(times, length(traindat$Y))
  cent[traindat$D==1] <- traindat$Y[traindat$D==1]
  c_hat <- pred_surv_preval(c_fit, S0, times = cent, lambda = c_lambda_min)
  ipcw <- 1 / c_hat
  
  # Propensity score
  w_fit <- cv.glmnet(as.matrix(traindat[,4:dim(traindat)[2]]), 
                     traindat$W,
                     foldid = foldid,
                     keep = TRUE,
                     alpha = alpha)
  
  w_lambda_min <- w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
  ps_score <- w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
  
  # Subset of uncensored subjects
  tempdat <- data.frame(data$Y, data$D, data$W, ps_score, ipcw, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0    
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  flasso_fit <- Flasso(x = as.matrix(binary_data[,6:10]), 
                       tx = binary_data$W, 
                       y = binary_data$D, 
                       pscore = binary_data$ps_score, 
                       weight = binary_data$ipcw, 
                       alpha = alpha,
                       nfolds = nfolds)
  pred_flasso <- as.vector(-predict(flasso_fit, data.test$X))
  pred_flasso
}

estimate_ipcw_wocf_gbm_flearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  
  # IPCW weights
  gbmfit <- gbm(Surv(Y, 1-D)~., 
                data = traindat, 
                distribution = "coxph",
                shrinkage = shrinkage, 
                n.trees = n.trees, 
                cv.folds = cv.folds, 
                n.cores = detectCores()) 
  best_iter <- gbm.perf(gbmfit, plot.it = F, method = 'cv')
  pred.train <- predict(gbmfit, traindat, n.trees = best_iter)
  cent <- data$Y
  cent[data$D==0] <- times
  basehaz.cum <- rep(NA, length(cent))
  for (k in 1:length(cent)){
    basehaz.cum[k] <- basehaz.gbm(traindat$Y, 1-traindat$D, pred.train, t.eval = cent[k], cumulative = TRUE)
  }
  c_hat <- exp(-exp(pred.train)*basehaz.cum)
  ipcw <- 1 / c_hat
  
  # Propensity score 
  tmpdat <- traindat[, 3:dim(traindat)[2]]
  w_fit <- gbm(W ~ ., 
               data = tmpdat, 
               distribution = "bernoulli",
               shrinkage = shrinkage, 
               n.trees = n.trees, 
               cv.folds = cv.folds, 
               n.cores = detectCores()) 
  best_iter <- gbm.perf(w_fit, plot.it = F, method = 'cv')
  log_odds <- predict(w_fit, tmpdat, n.trees = best_iter)
  ps_score <- 1/(1 + exp(-log_odds))
  
  # Subset of uncensored subjects
  tempdat <- data.frame(data$Y, data$D, data$W, ipcw, ps_score, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
                  
  fgbm_fit <- Fgbm(x = as.matrix(binary_data[,6:10]), 
                   tx = binary_data$W, 
                   y = binary_data$D, 
                   pscore = binary_data$ps_score, 
                   weight = binary_data$ipcw,
                   shrinkage = shrinkage, 
                   n.trees = n.trees, 
                   cv.folds = cv.folds)
  pred_fgbm <- -predict(fgbm_fit, newx = data.frame(data.test$X))
  pred_fgbm
}

estimate_ipcw_wcf_gbm_flearner <- function(data, data.test, nfolds = 10, cv.folds = 5, shrinkage = 0.005, n.trees = 2000, times = NULL){
  tempdat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  foldid <- sample(rep(seq(nfolds), length = length(tempdat$Y)))
  
  # IPCW weights
  c_hat <- rep(NA, length(tempdat$W))
  for (k in 1:nfolds){
    testdat <- tempdat[foldid==k, ]
    traindat <- tempdat[!foldid==k, ]
    c_fit <- gbm(Surv(Y, 1-D)~.,
                 data = traindat, 
                 distribution = "coxph",
                 shrinkage = shrinkage, 
                 n.trees = n.trees, 
                 cv.folds = cv.folds) 
    best_iter <- gbm.perf(c_fit, plot.it = F, method = 'cv')
    pred.train <- predict(c_fit, traindat, n.trees = best_iter)
    cent <- testdat$Y
    cent[testdat$D==0] <- times
    basehaz.cum <- rep(NA, length(cent))
    for (z in 1:length(cent)){
      basehaz.cum[z] <- basehaz.gbm(traindat$Y, 1-traindat$D, pred.train, t.eval = cent[z], cumulative = TRUE)
    }
    pred.test <- predict(c_fit, testdat, n.trees = best_iter)
    c_hat[foldid==k] <- exp(-exp(pred.test)*basehaz.cum)
  }
  ipcw <- 1 / c_hat

  # Propensity score 
  ps_score <- rep(NA, length(tempdat$W))
  for (k in 1:nfolds){
    testdat <- tempdat[foldid==k, ]
    traindat <- tempdat[!foldid==k, ]
    w_dat <- traindat[, 3:dim(traindat)[2]]
    w_fit <- gbm(W ~ ., 
                 data = w_dat, 
                 distribution = "bernoulli",
                 shrinkage = shrinkage, 
                 n.trees = n.trees, 
                 cv.folds = cv.folds, 
                 n.cores = detectCores()) 
    best_iter <- gbm.perf(w_fit, plot.it = F, method = 'cv')
    log_odds <- predict(w_fit, testdat[, 4:dim(testdat)[2]], n.trees = best_iter)
    ps_score[foldid==k] <- 1/(1 + exp(-log_odds))
  }

  # Subset of uncensored subjects
  tmpdat <- data.frame(data$Y, data$D, data$W, ipcw, ps_score, data$X)
  colnames(tmpdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tmpdat[tmpdat$D==1|tmpdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  fgbm_fit <- Fgbm(x = as.matrix(binary_data[,6:10]), 
                   tx = binary_data$W, 
                   y = binary_data$D, 
                   pscore = binary_data$ps_score, 
                   weight = binary_data$ipcw,
                   shrinkage = shrinkage, 
                   n.trees = n.trees, 
                   cv.folds = cv.folds)
  pred_fgbm <- -predict(fgbm_fit, data.frame(data.test$X))
  pred_fgbm
}

estimate_ipcw_grf_flearner <- function(data, data.test, times = NULL){
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  
  # IPCW weights
  sf.censor <- survival_forest(as.matrix(traindat[,3:dim(traindat)[2]]),
                               data$Y,
                               1 - data$D,
                               prediction.type = "Nelson-Aalen")
  C.hat <- sf.censor$predictions
  cent <- data$Y
  cent[data$D==0] <- times
  C.index <- rep(NA, length(cent))
  for (h in 1:length(cent)){
    C.index[h] <- which.min(abs(sort(unique(data$Y[data$D==0]))-cent[h]))
  }
  C.Y.hat <- C.hat[cbind(1:length(data$Y), C.index)]
  ipcw <- 1 / C.Y.hat
  
  # Propensity score 
  psfit <- regression_forest(as.matrix(traindat[,4:dim(traindat)[2]]), 
                             traindat$W, 
                             tune.parameters = "all") 
  ps <- psfit$predictions
  
  # Subset of uncensored subjects
  tempdat <- data.frame(data$Y, data$D, data$W, ps, ipcw, data$X)
  colnames(tempdat)[1:3] <- c("Y", "D", "W")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     
  binary_data <- binary_data[complete.cases(binary_data), ]
  
  fgrf_fit <- Fgrf(x = as.matrix(binary_data[,6:10]), 
                   tx = binary_data$W, 
                   y = binary_data$D, 
                   pscore = binary_data$ps, 
                   weight = binary_data$ipcw)
  pred_fgrf <- -predict(fgrf_fit, data.test$X)
  pred_fgrf
}



# Bayesian causal forest
estimate_ipcw_causalbart <- function(data, data.test, times = NULL){
  # IPCW weights
  sf.censor <- survival_forest(cbind(data$W, data$X),
                               data$Y,
                               1 - data$D,
                               prediction.type = "Nelson-Aalen")
  C.hat <- predict(sf.censor)$predictions
  cent <- data$Y
  cent[data$D==0] <- times
  C.index <- rep(NA, length(cent))
  for (h in 1:length(cent)){
    C.index[h] <- which.min(abs(sort(unique(data$Y[data$D==0]))-cent[h]))
  }
  C.Y.hat <- C.hat[cbind(1:length(data$Y), C.index)]
  data$sample.weights <- 1 / C.Y.hat
  
  # Propensity score 
  ps_fit <- regression_forest(data$X, data$W, tune.parameters = "all")
  data$ps_score <- predict(ps_fit)$predictions
  
  # Subset of uncensored subjects
  tempdat <- data.frame(data$Y, data$D, data$W, data$sample.weights, data$ps_score, data$X)
  colnames(tempdat)[1:5] <- c("Y", "D", "W", "ipcw", "ps_score")
  binary_data <- tempdat[tempdat$D==1|tempdat$Y > times,]          # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$Y > times] <- 0     # recode the event status for subjects who had events after t50
  
  bcf_fit <- bcf(y = binary_data$D, 
                 z = binary_data$W,
                 x_control = as.matrix(binary_data[, 6:10]),
                 x_moderate = as.matrix(binary_data[, 6:10]),
                 pihat = binary_data$ps_score, 
                 w = binary_data$ipcw, 
                 nburn = 10, 
                 nsim = 10,
                 n_cores = 2,
                 update_interval = 1)
  
  data.test$ps_pred <- predict(ps_fit, newdata = data.test$X)$predictions
  pred_out <- predict(object = bcf_fit,
                      x_predict_control = data.test$X,
                      x_predict_moderate = data.test$X,
                      pi_pred = data.test$ps_pred,
                      z_pred = data.test$W,
                      n_cores = 2,
                      save_tree_directory = '..')
  pred_bcf <- data.frame(Mean  = -colMeans(bcf_fit$tau),
                         Low95 = -apply(bcf_fit$tau, 2, function(x) quantile(x, 0.025)),
                         Up95  = -apply(bcf_fit$tau, 2, function(x) quantile(x, 0.975)))
  pred_bcf
}

