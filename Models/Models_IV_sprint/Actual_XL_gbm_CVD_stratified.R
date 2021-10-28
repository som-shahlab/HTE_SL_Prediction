##---------------------- Apply GBM method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021  -------------------##
library(survival);library(gbm) 
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
source("eval_funs.R")
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")

# Combine baseline variables and outcomes into one data frame 
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8580 57

# set counterfacutal treatment levels 
Xvar <- tmpdata[,2:dim(blvars)[2]]

Xvar1 <- Xvar
Xvar1$INTENSIVE <- 1

Xvar0 <- Xvar
Xvar0$INTENSIVE <- 0

#ith fold
if(length(args <- commandArgs(T))>0){
  caseid <- as.integer(args[[1]])
  message("running for parameter set ", caseid)
}
i <- caseid+1

#-------------------- Benefit model ----------------------------#
tmpdata_cvd <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar)
names(tmpdata_cvd)[1:2] <- c("t_cvds","cvd")
tmpdata_cvd1 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar1)
names(tmpdata_cvd1)[1:2] <- c("t_cvds","cvd")
tmpdata_cvd0 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar0)
names(tmpdata_cvd0)[1:2] <- c("t_cvds","cvd")

t50 <- floor(3.26*365.25)

# # 10-fold CV to get the out-of-bag estimates
# # Randomly shuffle the data
# set.seed(seed_num)
# index <- sample(nrow(tmpdata_cvd))
# tmpdata_cvd <- tmpdata_cvd[index,]
# tmpdata_cvd1 <- tmpdata_cvd1[index,]
# tmpdata_cvd0 <- tmpdata_cvd0[index,]
# folds <- cut(seq(1,nrow(tmpdata_cvd)),breaks=10,labels=FALSE)

# 10-fold CV to get the out-of-bag estimates (stratified to ensure equal number of events in each fold)
nfold <- 10
rownames(tmpdata_cvd) <- as.character(1:dim(tmpdata_cvd)[1])
rownames(tmpdata_cvd1) <- as.character(1:dim(tmpdata_cvd1)[1])
rownames(tmpdata_cvd0) <- as.character(1:dim(tmpdata_cvd0)[1])
tmpdata_cvd_event <- tmpdata_cvd[tmpdata_cvd$cvd==1,]
tmpdata_cvd_noevent <- tmpdata_cvd[tmpdata_cvd$cvd==0,]
folds_event <- cut(seq(1,nrow(tmpdata_cvd_event)),breaks=nfold,labels=FALSE)
folds_noevent <- cut(seq(1,nrow(tmpdata_cvd_noevent)),breaks=nfold,labels=FALSE)

# Sgbm_risk_treated_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# Sgbm_risk_control_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# Tgbm_risk_treated_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# Tgbm_risk_control_cvd <- rep(NA,dim(tmpdata_cvd)[1])


### T-Learner
# testIndexes <- which(folds==i,arr.ind=TRUE)
# traindat <- tmpdata_cvd[-testIndexes, ]
# traindat1 <- traindat[traindat$INTENSIVE==1,c(1:2,4:dim(traindat)[2])]
# traindat0 <- traindat[traindat$INTENSIVE==0,c(1:2,4:dim(traindat)[2])]
# testdat  <- tmpdata_cvd[testIndexes, c(1:2,4:dim(tmpdata_cvd)[2])]

testIndexes_event <- which(folds_event==i,arr.ind=TRUE)
testIndexes_noevent <- which(folds_noevent==i,arr.ind=TRUE)
testIndexes <- as.numeric(c(rownames(tmpdata_cvd_event[testIndexes_event,]), rownames(tmpdata_cvd_noevent[testIndexes_noevent,])))
traindat <- tmpdata_cvd[-testIndexes, ]
traindat1 <- traindat[traindat$INTENSIVE==1,c(1:2,4:dim(traindat)[2])]
traindat0 <- traindat[traindat$INTENSIVE==0,c(1:2,4:dim(traindat)[2])]
testdat  <- tmpdata_cvd[testIndexes,c(1:2,4:dim(tmpdata_cvd)[2])]

## Train using treated obs only 
# Tuning hyperparameters
maxdepths <- c(1,2,3)
minobs <- c(5,15,25)
subsamples <- c(0.7,0.8,0.9)
CV_errors_treated <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      #set.seed(1)
      gbm_tune_treated <- gbm(Surv(t_cvds,cvd)~., data=traindat1,
                              distribution = "coxph",
                              shrinkage = 0.005,
                              n.trees = 3000,
                              interaction.depth = interaction.depth,
                              n.minobsinnode = n.minobsinnode,
                              bag.fraction = bag.fraction,
                              cv.folds = 5,
                              n.cores = 5)
      best.iter <- gbm.perf(gbm_tune_treated, plot.it = F, method = 'cv')
      CV_errors_treated <- rbind(CV_errors_treated, c(interaction.depth, n.minobsinnode, bag.fraction, best.iter, gbm_tune_treated$cv.error[best.iter]))
      print(paste0("T1 model - interaction.depth=",interaction.depth," n.minobsinnode=",n.minobsinnode," bag.fraction=",bag.fraction))
    }
  }
}
CV_errors_treated <- data.frame(CV_errors_treated)
names(CV_errors_treated) <- c("interaction.depth", "n.minobsinnode", "bag.fraction", "best_iter", "deviance_err")
best_params_treated <- CV_errors_treated[which.min(CV_errors_treated$deviance_err),]

# Final tuned model for treated arm
set.seed(1)
gbm_final_treated <- gbm(Surv(t_cvds,cvd)~., data=traindat1,
                         distribution = "coxph",
                         shrinkage = 0.001,
                         n.trees = 10000,
                         interaction.depth = best_params_treated$interaction.depth,
                         n.minobsinnode = best_params_treated$n.minobsinnode,
                         bag.fraction = best_params_treated$bag.fraction,
                         cv.folds = 10,
                         n.cores = 10)
best_params_treated$best_iter <- gbm.perf(gbm_final_treated, plot.it = F, method = 'cv')

# return a vector of prediction on n.trees indicating log hazard scale.f(x)
pred.train1 <- predict(gbm_final_treated, traindat1, n.trees = best_params_treated$best_iter)
pred.test1  <- predict(gbm_final_treated, testdat, n.trees = best_params_treated$best_iter)
# Estimate the cumulative baseline hazard function using training data
basehaz.cum1 <- basehaz.gbm(traindat1$t_cvds, traindat1$cvd, pred.train1, t.eval = t50, cumulative = TRUE)


## Train using control obs only 
# tune hyperparameter
CV_errors_control <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      # set.seed(1)
      gbm_tune_control <- gbm(Surv(t_cvds,cvd)~., data=traindat0,
                              distribution = "coxph",
                              shrinkage = 0.005,
                              n.trees = 3000,
                              interaction.depth = interaction.depth,
                              n.minobsinnode = n.minobsinnode,
                              bag.fraction = bag.fraction,
                              cv.folds = 5,
                              n.cores = 5)
      best.iter <- gbm.perf(gbm_tune_control, plot.it = F, method = 'cv')
      CV_errors_control <- rbind(CV_errors_control, c(interaction.depth, n.minobsinnode, bag.fraction, best.iter, gbm_tune_control$cv.error[best.iter]))
      print(paste0("T0 model - interaction.depth=",interaction.depth," n.minobsinnode=",n.minobsinnode," bag.fraction=",bag.fraction))
    }
  }
}
CV_errors_control <- data.frame(CV_errors_control)
names(CV_errors_control) <- c("interaction.depth", "n.minobsinnode", "bag.fraction", "best_iter", "deviance_err")
best_params_control <- CV_errors_control[which.min(CV_errors_control$deviance_err),]

# Final tuned model for control arm
set.seed(1)
gbm_final_control <- gbm(Surv(t_cvds,cvd)~., data=traindat0,
                         distribution = "coxph",
                         shrinkage = 0.001,
                         n.trees = 10000,
                         interaction.depth = best_params_control$interaction.depth,
                         n.minobsinnode = best_params_control$n.minobsinnode,
                         bag.fraction = best_params_control$bag.fraction,
                         cv.folds = 10,
                         n.cores = 10)
best_params_control$best_iter <- gbm.perf(gbm_final_control, plot.it = F, method = 'cv')

# return a vector of prediction on n.trees indicating log hazard scale.f(x)
pred.train0 <- predict(gbm_final_control, traindat0, n.trees = best_params_control$best_iter)
pred.test0  <- predict(gbm_final_control, testdat, n.trees = best_params_control$best_iter)
# Estimate the cumulative baseline hazard function using training data
basehaz.cum0 <- basehaz.gbm(traindat0$t_cvds, traindat0$cvd, pred.train0, t.eval = t50, cumulative = TRUE)

# S(X) at median survival time
surf1 <- surf0 <- rep(NA,dim(testdat)[1])
for (j in 1:dim(testdat)[1]){
  surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
  surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
}
Tgbm_risk_treated_cvd <- 1-surf1
Tgbm_risk_control_cvd <- 1-surf0


### S-Learner
traindat <- tmpdata_cvd[-testIndexes,]
testdat1 <- tmpdata_cvd1[testIndexes,]
testdat0 <- tmpdata_cvd0[testIndexes,]

# For a given shrinkage and n.trees, tune the other 3 params
validerrors <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      # set.seed(1)
      gbmfit2 <- gbm(Surv(t_cvds,cvd)~., data=traindat,
                     distribution = "coxph",
                     shrinkage = 0.005,
                     n.trees = 10000,
                     interaction.depth = interaction.depth,
                     n.minobsinnode = n.minobsinnode,
                     bag.fraction = bag.fraction,
                     cv.folds = 5,
                     n.cores = 5)
      best.iter <- gbm.perf(gbmfit2, plot.it = F, method = 'cv')
      validerrors <- rbind(validerrors, c(interaction.depth, n.minobsinnode, bag.fraction, best.iter, gbmfit2$cv.error[best.iter]))
      print(paste0("S model-interaction.depth=",interaction.depth," n.minobsinnode=",n.minobsinnode," bag.fraction=",bag.fraction))
    }
  }
}
validerrors <- data.frame(validerrors)
names(validerrors) <- c("interaction.depth", "n.minobsinnode", "bag.fraction", "best_iter", "deviance_err")
best_params_S <- validerrors[which.min(validerrors$deviance_err),]

# Final tuned S-learner model 
set.seed(1)
gbmfit_final <- gbm(Surv(t_cvds,cvd)~., data=traindat, 
                    distribution = "coxph",
                    shrinkage = 0.001, 
                    n.trees = 20000,   
                    interaction.depth = best_params_S$interaction.depth, 
                    n.minobsinnode = best_params_S$n.minobsinnode, 
                    bag.fraction = best_params_S$bag.fraction,
                    cv.folds = 10,
                    n.cores = 10) 
best_params_S$best_iter <- gbm.perf(gbmfit_final, plot.it = F, method = 'cv')

# return a vector of prediction on n.trees indicating log hazard scale.f(x)
pred.train <- predict(gbmfit_final, traindat, n.trees = best_params_S$best_iter)
pred.test1 <- predict(gbmfit_final, testdat1, n.trees = best_params_S$best_iter)
pred.test0 <- predict(gbmfit_final, testdat0, n.trees = best_params_S$best_iter)

# Estimate the cumulative baseline hazard function using training data
basehaz.cum <- basehaz.gbm(traindat$t_cvds, traindat$cvd, pred.train, t.eval = t50, cumulative = TRUE)

# S(X) at median survival time 
surf1 <- surf0 <- rep(NA,dim(testdat1)[1])
for (j in 1:dim(testdat1)[1]){
  surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum)
  surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum)
}

# Sgbm
Sgbm_risk_treated_cvd <- 1-surf1
Sgbm_risk_control_cvd <- 1-surf0

# save results 
risk_gbm_cvd <- data.frame(testIndexes, Sgbm_risk_control_cvd, Sgbm_risk_treated_cvd, Tgbm_risk_control_cvd, Tgbm_risk_treated_cvd)
write.csv(risk_gbm_cvd, paste0("./cvd_tuned_results/stratified_split_gbm_cvd",i,".csv"), row.names = F)

# best_params <- data.frame(t(best_params_control), t(best_params_treated), t(best_params_S))
# write.csv(best_params, paste0("./cvd_tuned_results/gbm_best_params_cvd",i,".csv"), row.names = F)
