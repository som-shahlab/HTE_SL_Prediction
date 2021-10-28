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

#-------------------- Risk model ----------------------------#
sprint_sae <- data.frame(sprint_train$t_saes, sprint_train$sae, Xvar)
names(sprint_sae)[1:2] <- c("t_saes","sae")
accord_sae <- data.frame(accord_test$t_saes, accord_test$sae, testXvar)
names(accord_sae)[1:2] <- c("t_saes","sae")
accord_sae1 <- data.frame(accord_test$t_saes, accord_test$sae, testXvar1)
names(accord_sae1)[1:2] <- c("t_saes","sae")
accord_sae0 <- data.frame(accord_test$t_saes, accord_test$sae, testXvar0)
names(accord_sae0)[1:2] <- c("t_saes","sae")

t50 <- floor(3.26*12)

# T-Learner
sprint_sae1 <- sprint_sae[sprint_sae$INTENSIVE==1, c(1:2,4:dim(sprint_sae)[2])]
bartSurv1 <- mc.surv.bart(x.train=sprint_sae1[,3:dim(sprint_sae1)[2]],
                          times=ceiling(sprint_sae1$t_saes/30),
                          delta=sprint_sae1$sae,
                          ntree=50,
                          x.test=accord_sae[,4:dim(accord_sae)[2]],
                          seed=99, 
                          mc.cores=4)
surv_ests <- t(matrix(bartSurv1$surv.test.mean, nrow = length(bartSurv1$times)))
surv_est1 <- surv_ests[,which.min(abs(bartSurv1$times-t50))]
Tbart_risk_treated_sae <- 1-surv_est1

# save results 
risk_bart_sae <- data.frame(Tbart_risk_treated_sae)
write.csv(risk_bart_sae, paste0("./sae_tuned_results/external_T1bart_sae.csv"), row.names = F)

# Sgbm_predinc_sae <- Sbart_risk_treated_sae*accord_sae$INTENSIVE + Sbart_risk_control_sae*(1-accord_sae$INTENSIVE)
# Spredgrp <- as.numeric(cut2(Sgbm_predinc_sae,g=10))
# GND_S <- GND.calib(pred=Sgbm_predinc_sae, tvar=ceiling(accord_sae$t_saes/30), out=accord_sae$sae,
#                    cens.t=t50, groups=Spredgrp, adm.cens=t50,name="")

