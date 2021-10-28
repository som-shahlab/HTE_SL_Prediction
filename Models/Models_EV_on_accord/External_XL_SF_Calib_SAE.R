##---------------------- Apply Survival Forest method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021              -------------------##
library(survival);library(grf);library(pec)
source("eval_funs.R")
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 

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
        eventtimes <- sort(unique(sprint_sae1$t_saes[sprint_sae1$sae==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit1 <- survival_forest(X = sprint_sae1[,3:dim(sprint_sae1)[2]],
                                   Y = sprint_sae1$t_saes,
                                   D = sprint_sae1$sae,
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
        eventtimes <- sort(unique(sprint_sae0$t_saes[sprint_sae0$sae==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit0 <- survival_forest(X = sprint_sae0[,3:dim(sprint_sae0)[2]],
                                   Y = sprint_sae0$t_saes,
                                   D = sprint_sae0$sae,
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
        GND_T1 <- GND.calib.simp(pred=surf1, tvar=sprint_sae1$t_saes, out=sprint_sae1$sae,
                           cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("T1-learner: ",c(ntree, mtry, nodesize, alpha, GND_T1[3])))
        pet1 <- rbind(pet1, c(ntree, mtry, nodesize, alpha, GND_T1[3]))

        # Calibration
        Spredgrp <- as.numeric(cut2(surf0,g=5))
        GND_T0 <- GND.calib.simp(pred=surf0, tvar=sprint_sae0$t_saes, out=sprint_sae0$sae,
                           cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("T0-learner: ",c(ntree, mtry, nodesize, alpha, GND_T0[3])))
        pet0 <- rbind(pet0, c(ntree, mtry, nodesize, alpha, GND_T0[3]))
      }
    }
  }
}
pet1 <- data.frame(pet1); colnames(pet1)[dim(pet1)[2]] <- "error"
opt_params1 <- as.list(pet1[which.max(pet1$error),])
eventtimes <- sort(unique(sprint_sae1$t_saes[sprint_sae1$sae==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit1 <- survival_forest(X = sprint_sae1[,3:dim(sprint_sae1)[2]],
                           Y = sprint_sae1$t_saes,
                           D = sprint_sae1$sae,
                           failure.times = failure.times,
                           num.trees = opt_params1[[1]],
                           mtry = opt_params1[[2]],
                           sample.fraction = 1,
                           min.node.size = opt_params1[[3]],
                           alpha = opt_params1[[4]],
                           honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                           seed = 99)
surf1 <- predict(sf_fit1, newdata = accord_sae[,4:dim(accord_sae)[2]])
Tsf_risk_treated_sae <- 1-surf1$predictions[, which.min(abs(surf1$failure.times-t50))]


pet0 <- data.frame(pet0); colnames(pet0)[dim(pet0)[2]] <- "error"
opt_params0 <- as.list(pet0[which.max(pet0$error),])
eventtimes <- sort(unique(sprint_sae0$t_saes[sprint_sae0$sae==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit0 <- survival_forest(X = sprint_sae0[,3:dim(sprint_sae0)[2]],
                           Y = sprint_sae0$t_saes,
                           D = sprint_sae0$sae,
                           failure.times = failure.times,
                           num.trees = opt_params0[[1]],
                           mtry = opt_params0[[2]],
                           sample.fraction = 1,
                           min.node.size = opt_params0[[3]],
                           alpha = opt_params0[[4]],
                           honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                           seed = 99)
surf0 <- predict(sf_fit0, newdata = accord_sae[,4:dim(accord_sae)[2]])
Tsf_risk_control_sae <- 1-surf0$predictions[, which.min(abs(surf0$failure.times-t50))]


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
        
        eventtimes <- sort(unique(sprint_sae$t_saes[sprint_sae$sae==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit <- survival_forest(X = sprint_sae[,3:dim(sprint_sae)[2]],
                                  Y = sprint_sae$t_saes,
                                  D = sprint_sae$sae,
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
        GND_S <- GND.calib.simp(pred=surf, tvar=sprint_sae$t_saes, out=sprint_sae$sae,
                           cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("S-learner: ", c(ntree, mtry, nodesize, alpha, GND_S[3])))
        pet <- rbind(pet, c(ntree, mtry, nodesize, alpha, GND_S[3]))
      }
    }
  }
}
pet <- data.frame(pet); colnames(pet)[dim(pet)[2]] <- "error"
opt_params <- as.list(pet[which.max(pet$error),])
eventtimes <- sort(unique(sprint_sae$t_saes[sprint_sae$sae==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit <- survival_forest(X = sprint_sae[,3:dim(sprint_sae)[2]],
                          Y = sprint_sae$t_saes,
                          D = sprint_sae$sae,
                          failure.times = failure.times,
                          sample.fraction = 1,
                          num.trees = opt_params[[1]],
                          mtry = opt_params[[2]],
                          min.node.size = opt_params[[3]],
                          alpha = opt_params[[4]],
                          honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                          seed = 99)
surf1 <- predict(sf_fit, newdata = accord_sae1[,3:dim(accord_sae1)[2]])
surf0 <- predict(sf_fit, newdata = accord_sae0[,3:dim(accord_sae0)[2]])
Ssf_risk_treated_sae <- 1-surf1$predictions[, which.min(abs(surf1$failure.times-t50))]
Ssf_risk_control_sae <- 1-surf0$predictions[, which.min(abs(surf0$failure.times-t50))]

risk_sf_sae <- data.frame(Ssf_risk_control_sae,Ssf_risk_treated_sae,Tsf_risk_control_sae,Tsf_risk_treated_sae)
write.csv(risk_sf_sae, paste0("./sae_tuned_results/external_sf_sae.csv"), row.names = F)
best_params <- data.frame(unlist(opt_params1), unlist(opt_params0), unlist(opt_params))
write.csv(best_params, paste0("./sae_tuned_results/external_sf_best_params_sae.csv"), row.names = F)

# sum_sae <- data.frame(Ssf_risk_control_sae,Ssf_risk_treated_sae)
# S_predinc_sae <- sum_sae[,2]*accord_sae$INTENSIVE + sum_sae[,1]*(1-accord_sae$INTENSIVE)
# Spredgrp <- as.numeric(cut2(S_predinc_sae,g=10))
# GND_Ssae <- GND.calib(pred=S_predinc_sae, tvar=accord_sae$t_saes, out=accord_sae$sae,
#                       cens.t=t50, groups=Spredgrp, adm.cens=t50, name="")

