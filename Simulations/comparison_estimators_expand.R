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

# "IPCW" - X-learner gradient boosting
t50 <- median(data$Y)
gbm_fit <- function(data, data.test){
  
  # fit model on W==1 
  traindat <- data.frame(data$Y, data$D, data$W, data$X)
  colnames(traindat)[1:3] <- c("Y", "D", "W")
  testdat <- traindat[, !colnames(traindat) %in% c("W")]
  train1 <- traindat[traindat$W==1, !colnames(traindat) %in% c("W")]
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
  sf.censor <- survival_forest(binary_data[,4:8], binary_data$Y, 1 - binary_data$D, prediction.type = "Nelson-Aalen", num.trees = 500)
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
  XLtau1 * ps$predictions + XLtau0 * (1 - ps$predictions)
}












