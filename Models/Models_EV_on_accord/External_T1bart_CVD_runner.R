##---------------------- Apply BART method to the actual SPRINT data -------------------##
##----------------------     Crystal Xu               02/23/2021     -------------------##
library(survival);library(BART)
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
#source("eval_funs.R")

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

#-------------------- Benefit model ----------------------------#
sprint_cvd <- data.frame(sprint_train$t_cvds, sprint_train$cvd, Xvar)
names(sprint_cvd)[1:2] <- c("t_cvds","cvd")
accord_cvd <- data.frame(accord_test$t_cvds, accord_test$cvd, testXvar)
names(accord_cvd)[1:2] <- c("t_cvds","cvd")
accord_cvd1 <- data.frame(accord_test$t_cvds, accord_test$cvd, testXvar1)
names(accord_cvd1)[1:2] <- c("t_cvds","cvd")
accord_cvd0 <- data.frame(accord_test$t_cvds, accord_test$cvd, testXvar0)
names(accord_cvd0)[1:2] <- c("t_cvds","cvd")

t50 <- floor(3.26*12)

# T-Learner
sprint_cvd1 <- sprint_cvd[sprint_cvd$INTENSIVE==1, c(1:2,4:dim(sprint_cvd)[2])]
bartSurv1 <- mc.surv.bart(x.train=sprint_cvd1[,3:dim(sprint_cvd1)[2]],
                          times=ceiling(sprint_cvd1$t_cvds/30),
                          delta=sprint_cvd1$cvd,
                          ntree=50,
                          x.test=accord_cvd[,4:dim(accord_cvd)[2]],
                          seed=99, 
                          mc.cores=4)
surv_ests <- t(matrix(bartSurv1$surv.test.mean, nrow = length(bartSurv1$times)))
surv_est1 <- surv_ests[,which.min(abs(bartSurv1$times-t50))]
Tbart_risk_treated_cvd <- 1-surv_est1

# save results 
risk_bart_cvd <- data.frame(Tbart_risk_treated_cvd)
write.csv(risk_bart_cvd, paste0("./cvd_tuned_results/external_T1bart_cvd.csv"), row.names = F)

# Sgbm_predinc_cvd <- Sbart_risk_treated_cvd*accord_cvd$INTENSIVE + Sbart_risk_control_cvd*(1-accord_cvd$INTENSIVE)
# Spredgrp <- as.numeric(cut2(Sgbm_predinc_cvd,g=10))
# GND_S <- GND.calib(pred=Sgbm_predinc_cvd, tvar=ceiling(accord_cvd$t_cvds/30), out=accord_cvd$cvd,
#                    cens.t=t50, groups=Spredgrp, adm.cens=t50,name="")

