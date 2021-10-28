##---------------------- Apply Survival Forest method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021              -------------------##
library(survival);library(grf);library(pec)
source("eval_funs.R")
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 

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

# Tuning hyperparameters
alphas <- c(0.001,0.005,0.01)
nodesizes <- c(15,20,25)
ntrees <- c(200,500,1000)
mtrys <- c(21,24,27)
pet1 <- pet0 <- NULL
for (j in 1:length(ntrees)){
  ntree <- ntrees[j]
  for (k in 1:length(mtrys)){
    mtry <- mtrys[k]
    for (l in 1:length(nodesizes)){
      nodesize <- nodesizes[l]
      for (h in 1: length(alphas)){
        alpha <- alphas[h]
        
        # train using treated obs only
        eventtimes <- sort(unique(accord_cvd1$t_cvds[accord_cvd1$cvd==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit1 <- survival_forest(X = accord_cvd1[,3:dim(accord_cvd1)[2]],
                                   Y = accord_cvd1$t_cvds,
                                   D = accord_cvd1$cvd,
                                   failure.times = failure.times,
                                   num.trees = ntree,
                                   mtry = mtry,
                                   sample.fraction = 0.8,
                                   min.node.size = nodesize,
                                   alpha = alpha,
                                   honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                                   seed = 99)
        surf1 <- 1-sf_fit1$predictions[, which.min(abs(sf_fit1$failure.times-t50))]
        
        # train using control obs only
        eventtimes <- sort(unique(accord_cvd0$t_cvds[accord_cvd0$cvd==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit0 <- survival_forest(X = accord_cvd0[,3:dim(accord_cvd0)[2]],
                                   Y = accord_cvd0$t_cvds,
                                   D = accord_cvd0$cvd,
                                   failure.times = failure.times,
                                   num.trees = ntree,
                                   mtry = mtry,
                                   sample.fraction = 0.8,
                                   min.node.size = nodesize,
                                   alpha = alpha,
                                   honesty.prune.leaves = F,
                                   seed = 99)
        surf0 <- 1-sf_fit0$predictions[, which.min(abs(sf_fit0$failure.times-t50))]

        # Calibration
        Spredgrp <- as.numeric(cut2(surf1,g=5))
        GND_T1 <- GND.calib.simp(pred=surf1, tvar=accord_cvd1$t_cvds, out=accord_cvd1$cvd,
                                 cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("T1-learner: ",c(ntree, mtry, nodesize, alpha, GND_T1[3])))
        pet1 <- rbind(pet1, c(ntree, mtry, nodesize, alpha, GND_T1[3]))

        # Calibration
        Spredgrp <- as.numeric(cut2(surf0,g=5))
        GND_T0 <- GND.calib.simp(pred=surf0, tvar=accord_cvd0$t_cvds, out=accord_cvd0$cvd,
                                 cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("T0-learner: ",c(ntree, mtry, nodesize, alpha, GND_T0[3])))
        pet0 <- rbind(pet0, c(ntree, mtry, nodesize, alpha, GND_T0[3]))
      }
    }
  }
}
pet1 <- data.frame(pet1); colnames(pet1)[dim(pet1)[2]] <- "error"
opt_params1 <- as.list(pet1[which.max(pet1$error),])
eventtimes <- sort(unique(accord_cvd1$t_cvds[accord_cvd1$cvd==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit1 <- survival_forest(X = accord_cvd1[,3:dim(accord_cvd1)[2]],
                           Y = accord_cvd1$t_cvds,
                           D = accord_cvd1$cvd,
                           failure.times = failure.times,
                           num.trees = opt_params1[[1]],
                           mtry = opt_params1[[2]],
                           sample.fraction = 1,
                           min.node.size = opt_params1[[3]],
                           alpha = opt_params1[[4]],
                           honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                           seed = 99)
surf1 <- predict(sf_fit1, newdata = sprint_cvd[,4:dim(sprint_cvd)[2]])
Tsf_risk_treated_cvd <- 1-surf1$predictions[, which.min(abs(surf1$failure.times-t50))]


pet0 <- data.frame(pet0); colnames(pet0)[dim(pet0)[2]] <- "error"
opt_params0 <- as.list(pet0[which.max(pet0$error),])
eventtimes <- sort(unique(accord_cvd0$t_cvds[accord_cvd0$cvd==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit0 <- survival_forest(X = accord_cvd0[,3:dim(accord_cvd0)[2]],
                           Y = accord_cvd0$t_cvds,
                           D = accord_cvd0$cvd,
                           failure.times = failure.times,
                           num.trees = opt_params0[[1]],
                           mtry = opt_params0[[2]],
                           sample.fraction = 1,
                           min.node.size = opt_params0[[3]],
                           alpha = opt_params0[[4]],
                           honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                           seed = 99)
surf0 <- predict(sf_fit0, newdata = sprint_cvd[,4:dim(sprint_cvd)[2]])
Tsf_risk_control_cvd <- 1-surf0$predictions[, which.min(abs(surf0$failure.times-t50))]


## S-Learner
# Tuning hyperparameters
pet <- NULL
for (j in 1:length(ntrees)){
  ntree <- ntrees[j]
  for (k in 1:length(mtrys)){
    mtry <- mtrys[k]
    for (l in 1:length(nodesizes)){
      nodesize <- nodesizes[l]
      for (h in 1: length(alphas)){
        alpha <- alphas[h]
        
        eventtimes <- sort(unique(accord_cvd$t_cvds[accord_cvd$cvd==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit <- survival_forest(X = accord_cvd[,3:dim(accord_cvd)[2]],
                                  Y = accord_cvd$t_cvds,
                                  D = accord_cvd$cvd,
                                  failure.times = failure.times,
                                  sample.fraction = 0.8,
                                  num.trees = ntree,
                                  mtry = mtry,
                                  min.node.size = nodesize,
                                  alpha = alpha,
                                  honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                                  seed = 99)
        surf <- 1-sf_fit$predictions[, which.min(abs(sf_fit$failure.times-t50))]

        # Calibration
        Spredgrp <- as.numeric(cut2(surf,g=5))
        GND_S <- GND.calib.simp(pred=surf, tvar=accord_cvd$t_cvds, out=accord_cvd$cvd,
                           cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("S-learner: ", c(ntree, mtry, nodesize, alpha, GND_S[3])))
        pet <- rbind(pet, c(ntree, mtry, nodesize, alpha, GND_S[3]))
      }
    }
  }
}
pet <- data.frame(pet); colnames(pet)[dim(pet)[2]] <- "error"
opt_params <- as.list(pet[which.max(pet$error),])
eventtimes <- sort(unique(accord_cvd$t_cvds[accord_cvd$cvd==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit <- survival_forest(X = accord_cvd[,3:dim(accord_cvd)[2]],
                          Y = accord_cvd$t_cvds,
                          D = accord_cvd$cvd,
                          failure.times = failure.times,
                          sample.fraction = 1,
                          num.trees = opt_params[[1]],
                          mtry = opt_params[[2]],
                          min.node.size = opt_params[[3]],
                          alpha = opt_params[[4]],
                          honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                          seed = 99)
surf1 <- predict(sf_fit, newdata = sprint_cvd1[,3:dim(sprint_cvd1)[2]])
surf0 <- predict(sf_fit, newdata = sprint_cvd0[,3:dim(sprint_cvd0)[2]])
Ssf_risk_treated_cvd <- 1-surf1$predictions[, which.min(abs(surf1$failure.times-t50))]
Ssf_risk_control_cvd <- 1-surf0$predictions[, which.min(abs(surf0$failure.times-t50))]

risk_sf_cvd <- data.frame(Ssf_risk_control_cvd,Ssf_risk_treated_cvd,Tsf_risk_control_cvd,Tsf_risk_treated_cvd)
write.csv(risk_sf_cvd, paste0("./cvd_tuned_results/external_sprint_sf_cvd.csv"), row.names = F)

