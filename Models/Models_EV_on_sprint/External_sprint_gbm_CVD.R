##---------------------- Apply GBM method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021  -------------------##
library(survival);library(gbm) 
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
source("eval_funs.R")

# ACCORD as training data 
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_train <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_train <- accord_train[complete.cases(accord_train),]; dim(accord_train) # 4448 43
accord_train$t_cvds <- ceiling(accord_train$t_cvds)
accord_train$t_saes <- ceiling(accord_train$t_saes)
Xvar <- accord_train[,2:dim(accord_blvars)[2]]

# SPRINT as testing data 
sprint_blvars <- read.csv("blvars.csv")
sprint_outcome <- read.csv("outcome.csv")
sprint_test <- merge(sprint_blvars, sprint_outcome, by="MASKID",all=T)
sprint_test <- sprint_test[complete.cases(sprint_test),]; dim(sprint_test) # 8969 43

# Set counterfactual treatment levels 
testXvar <- sprint_test[,2:dim(sprint_blvars)[2]]
testXvar1 <- testXvar; testXvar1$INTENSIVE <- 1
testXvar0 <- testXvar; testXvar0$INTENSIVE <- 0

#-------------------- Benefit model ----------------------------#
accord_cvd <- data.frame(accord_train$t_cvds, accord_train$cvd, Xvar)
names(accord_cvd)[1:2] <- c("t_cvds","cvd")
sprint_cvd <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar)
names(sprint_cvd)[1:2] <- c("t_cvds","cvd")
sprint_cvd1 <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar1)
names(sprint_cvd1)[1:2] <- c("t_cvds","cvd")
sprint_cvd0 <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar0)
names(sprint_cvd0)[1:2] <- c("t_cvds","cvd")

t50 <- floor(3.26*365.25)

# T-Learner
accord_cvd1 <- accord_cvd[accord_cvd$INTENSIVE==1,c(1:2,4:dim(accord_cvd)[2])]
accord_cvd0 <- accord_cvd[accord_cvd$INTENSIVE==0,c(1:2,4:dim(accord_cvd)[2])]

sprint_cvdtr1 <- sprint_cvd[sprint_cvd$INTENSIVE==1,c(1:2,4:dim(sprint_cvd)[2])]
sprint_cvdtr0 <- sprint_cvd[sprint_cvd$INTENSIVE==0,c(1:2,4:dim(sprint_cvd)[2])]

# Tuning hyperparameters
maxdepths <- c(1,2,3)
minobs <- c(5,15,25)
subsamples <- c(0.7,0.8,0.9)

#--------------------------- train using treated obs only ------------------------------#
CV_errors_treated <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      gbm_tune_treated <- gbm(Surv(t_cvds,cvd)~., data=accord_cvd1,
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
gbm_final_treated <- gbm(Surv(t_cvds,cvd)~., data=accord_cvd1,
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
pred.test1  <- predict(gbm_final_treated, sprint_cvd[,4:dim(sprint_cvd)[2]], n.trees = best_params_treated$best_iter)
# Estimate the cumulative baseline hazard function using treated patients in ACCORD as the training data
sprint.train1 <- predict(gbm_final_treated, sprint_cvdtr1[,3:dim(sprint_cvdtr1)[2]], n.trees = best_params_treated$best_iter)
basehaz.cum1 <- basehaz.gbm(sprint_cvdtr1$t_cvds, sprint_cvdtr1$cvd, sprint.train1, t.eval = t50, cumulative = T)


#--------------------------- train using control obs only ------------------------------##
# tune hyperparameter
CV_errors_control <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      gbm_tune_control <- gbm(Surv(t_cvds,cvd)~., data=accord_cvd0,
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
gbm_final_control <- gbm(Surv(t_cvds,cvd)~., data=accord_cvd0,
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
pred.test0  <- predict(gbm_final_control, sprint_cvd[,4:dim(sprint_cvd)[2]], n.trees = best_params_control$best_iter)
# Estimate the cumulative baseline hazard function using control patients in ACCORD as the training data
sprint.train0 <- predict(gbm_final_control, sprint_cvdtr0[,3:dim(sprint_cvdtr0)[2]], n.trees = best_params_control$best_iter)
basehaz.cum0 <- basehaz.gbm(sprint_cvdtr0$t_cvds, sprint_cvdtr0$cvd, sprint.train0, t.eval = t50, cumulative = T)


# S(X) at median survival time
surf1 <- surf0 <- rep(NA,dim(sprint_cvd)[1])
for (j in 1:dim(sprint_cvd)[1]){
  surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
  surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
}
Tgbm_risk_treated_cvd <- 1-surf1
Tgbm_risk_control_cvd <- 1-surf0


### S-Learner
# For a given shrinkage and n.trees, tune the other 3 params
validerrors <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      gbmfit2 <- gbm(Surv(t_cvds,cvd)~., data=accord_cvd,
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
gbmfit_final <- gbm(Surv(t_cvds,cvd)~., data=accord_cvd, 
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
pred.test1 <- predict(gbmfit_final, sprint_cvd1[,3:dim(sprint_cvd1)[2]], n.trees = best_params_S$best_iter)
pred.test0 <- predict(gbmfit_final, sprint_cvd0[,3:dim(sprint_cvd0)[2]], n.trees = best_params_S$best_iter)

# Estimate the cumulative baseline hazard function using ACCORD as the training data
sprint.train <- predict(gbmfit_final, sprint_cvd[,3:dim(sprint_cvd)[2]], n.trees = best_params_S$best_iter)
basehaz.cum <- basehaz.gbm(sprint_cvd$t_cvds, sprint_cvd$cvd, sprint.train, t.eval = t50, cumulative = T)

# S(X) at median survival time 
surf1 <- surf0 <- rep(NA,dim(sprint_cvd1)[1])
for (j in 1:dim(sprint_cvd1)[1]){
  surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum)
  surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum)
}

# Sgbm
Sgbm_risk_treated_cvd <- 1-surf1
Sgbm_risk_control_cvd <- 1-surf0

# save results 
risk_gbm_cvd <- data.frame(Sgbm_risk_control_cvd, Sgbm_risk_treated_cvd, Tgbm_risk_control_cvd, Tgbm_risk_treated_cvd)
write.csv(risk_gbm_cvd, paste0("./cvd_tuned_results/external_sprint_gbm_cvd.csv"), row.names = F)

