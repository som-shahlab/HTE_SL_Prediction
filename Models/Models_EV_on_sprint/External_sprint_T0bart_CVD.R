##---------------------- Apply BART method to the actual SPRINT data -------------------##
##----------------------     Crystal Xu               02/23/2021     -------------------##
library(survival);library(BART)
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
#source("eval_funs.R")

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

t50 <- floor(3.26*12)

# T-Learner
accord_cvd0 <- accord_cvd[accord_cvd$INTENSIVE==0, c(1:2,4:dim(accord_cvd)[2])]
bartSurv0 <- mc.surv.bart(x.train=accord_cvd0[,3:dim(accord_cvd0)[2]],
                          times=ceiling(accord_cvd0$t_cvds/30),
                          delta=accord_cvd0$cvd,
                          ntree=50,
                          x.test=sprint_cvd[,4:dim(sprint_cvd)[2]],
                          seed=99, 
                          mc.cores=4)
surv_ests <- t(matrix(bartSurv0$surv.test.mean, nrow = length(bartSurv0$times)))
surv_est0 <- surv_ests[,which.min(abs(bartSurv0$times-t50))]
Tbart_risk_control_cvd <- 1-surv_est0

# save results 
risk_bart_cvd <- data.frame(Tbart_risk_control_cvd)
write.csv(risk_bart_cvd, paste0("./cvd_tuned_results/external_sprint_T0bart_cvd.csv"), row.names = F)

