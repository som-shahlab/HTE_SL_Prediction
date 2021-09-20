##---------------------- Apply GBM method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021  -------------------##
library(survival);library(gbm) 
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
source("eval_funs.R")

# SPRINT as training data 
sprint_blvars <- read.csv("blvars.csv")
sprint_outcome <- read.csv("outcome.csv")
sprint_train <- merge(sprint_blvars, sprint_outcome, by="MASKID",all=T)
sprint_train <- sprint_train[complete.cases(sprint_train),]; dim(sprint_train) # 8969 43
Xvar <- sprint_train[,2:dim(sprint_blvars)[2]]

# ACCORD as testing data 
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_test <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_test <- accord_test[complete.cases(accord_test),]; dim(accord_test) # 4448 43
accord_test$t_cvds <- ceiling(accord_test$t_cvds)
accord_test$t_saes <- ceiling(accord_test$t_saes)

# Set counterfacutal treatment levels 
testXvar <- accord_test[,2:dim(accord_blvars)[2]]
testXvar1 <- testXvar; testXvar1$INTENSIVE <- 1
testXvar0 <- testXvar; testXvar0$INTENSIVE <- 0

#-------------------- Risk model ----------------------------#
sprint_sae <- data.frame(sprint_train$t_saes, sprint_train$sae, Xvar)
names(sprint_sae)[1:2] <- c("t_saes","sae")
accord_sae <- data.frame(accord_test$t_saes, accord_test$sae, testXvar)
names(accord_sae)[1:2] <- c("t_saes","sae")
accord_sae1 <- data.frame(accord_test$t_saes, accord_test$sae, testXvar1)
names(accord_sae1)[1:2] <- c("t_saes","sae")
accord_sae0 <- data.frame(accord_test$t_saes, accord_test$sae, testXvar0)
names(accord_sae0)[1:2] <- c("t_saes","sae")

t50 <- floor(3.26*365.25)

# T-Learner
sprint_sae1 <- sprint_sae[sprint_sae$INTENSIVE==1,c(1:2,4:dim(sprint_sae)[2])]
sprint_sae0 <- sprint_sae[sprint_sae$INTENSIVE==0,c(1:2,4:dim(sprint_sae)[2])]

accord_saetr1 <- accord_sae[accord_sae$INTENSIVE==1,c(1:2,4:dim(accord_sae)[2])]
accord_saetr0 <- accord_sae[accord_sae$INTENSIVE==0,c(1:2,4:dim(accord_sae)[2])]

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
      gbm_tune_treated <- gbm(Surv(t_saes,sae)~., data=sprint_sae1,
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
gbm_final_treated <- gbm(Surv(t_saes,sae)~., data=sprint_sae1,
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
#pred.train1 <- predict(gbm_final_treated, sprint_sae1, n.trees = best_params_treated$best_iter)
pred.test1  <- predict(gbm_final_treated, accord_sae[,4:dim(accord_sae)[2]], n.trees = best_params_treated$best_iter)
# Estimate the cumulative baseline hazard function using treated patients in ACCORD as the training data
accord.train1 <- predict(gbm_final_treated, accord_saetr1[,3:dim(accord_saetr1)[2]], n.trees = best_params_treated$best_iter)
basehaz.cum1 <- basehaz.gbm(accord_saetr1$t_saes, accord_saetr1$sae, accord.train1, t.eval = t50, cumulative = T)


#--------------------------- train using control obs only ------------------------------##
# tune hyperparameter
CV_errors_control <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      gbm_tune_control <- gbm(Surv(t_saes,sae)~., data=sprint_sae0,
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
gbm_final_control <- gbm(Surv(t_saes,sae)~., data=sprint_sae0,
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
#pred.train0 <- predict(gbm_final_control, sprint_sae0, n.trees = best_params_control$best_iter)
pred.test0  <- predict(gbm_final_control, accord_sae[,4:dim(accord_sae)[2]], n.trees = best_params_control$best_iter)
# Estimate the cumulative baseline hazard function using control patients in ACCORD as the training data
accord.train0 <- predict(gbm_final_control, accord_saetr0[,3:dim(accord_saetr0)[2]], n.trees = best_params_control$best_iter)
basehaz.cum0 <- basehaz.gbm(accord_saetr0$t_saes, accord_saetr0$sae, accord.train0, t.eval = t50, cumulative = T)


# S(X) at median survival time
surf1 <- surf0 <- rep(NA,dim(accord_sae)[1])
for (j in 1:dim(accord_sae)[1]){
  surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum1)
  surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum0)
}
Tgbm_risk_treated_sae <- 1-surf1
Tgbm_risk_control_sae <- 1-surf0


### S-Learner
# For a given shrinkage and n.trees, tune the other 3 params
validerrors <- NULL
for (l in 1:length(maxdepths)){
  interaction.depth <- maxdepths[l]
  for (j in 1:length(minobs)){
    n.minobsinnode <- minobs[j]
    for (k in 1:length(subsamples)){
      bag.fraction <- subsamples[k]
      gbmfit2 <- gbm(Surv(t_saes,sae)~., data=sprint_sae,
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
gbmfit_final <- gbm(Surv(t_saes,sae)~., data=sprint_sae, 
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
pred.test1 <- predict(gbmfit_final, accord_sae1[,3:dim(accord_sae1)[2]], n.trees = best_params_S$best_iter)
pred.test0 <- predict(gbmfit_final, accord_sae0[,3:dim(accord_sae0)[2]], n.trees = best_params_S$best_iter)

# Estimate the cumulative baseline hazard function using ACCORD as the training data
accord.train <- predict(gbmfit_final, accord_sae[,3:dim(accord_sae)[2]], n.trees = best_params_S$best_iter)
basehaz.cum <- basehaz.gbm(accord_sae$t_saes, accord_sae$sae, accord.train, t.eval = t50, cumulative = T)

# S(X) at median survival time 
surf1 <- surf0 <- rep(NA,dim(accord_sae1)[1])
for (j in 1:dim(accord_sae1)[1]){
  surf1[j] <- exp(-exp(pred.test1[j])*basehaz.cum)
  surf0[j] <- exp(-exp(pred.test0[j])*basehaz.cum)
}

# Sgbm
Sgbm_risk_treated_sae <- 1-surf1
Sgbm_risk_control_sae <- 1-surf0

# save results 
risk_gbm_sae <- data.frame(Sgbm_risk_control_sae, Sgbm_risk_treated_sae, Tgbm_risk_control_sae, Tgbm_risk_treated_sae)
write.csv(risk_gbm_sae, paste0("./sae_tuned_results/external_gbm_sae.csv"), row.names = F)

best_params <- data.frame(t(best_params_control), t(best_params_treated), t(best_params_S))
write.csv(best_params, paste0("./sae_tuned_results/external_gbm_best_params_sae.csv"), row.names = F)

