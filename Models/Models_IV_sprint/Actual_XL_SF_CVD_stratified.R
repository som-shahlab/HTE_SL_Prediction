##---------------------- Apply Survival Forest method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021              -------------------##
library(survival);library(grf);library(pec)
source("eval_funs.R")
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
# Combine baseline variables and outcomes into one data frame 
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969 43

# set counterfacutal treatment levels 
Xvar <- tmpdata[,2:dim(blvars)[2]]
Xvar1 <- Xvar
Xvar1$INTENSIVE <- 1
Xvar0 <- Xvar
Xvar0$INTENSIVE <- 0

# ith fold 
if(length(args <- commandArgs(T))>0){
  caseid <- as.integer(args[[1]])
  message("running for parameter set ", caseid)
}
i <- caseid+1
t50 <- floor(365.25*3.26)

#-------------------- Benefit model ----------------------------#
tmpdata_cvd <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar)
names(tmpdata_cvd)[1:2] <- c("t_cvds","cvd")
tmpdata_cvd1 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar1)
names(tmpdata_cvd1)[1:2] <- c("t_cvds","cvd")
tmpdata_cvd0 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar0)
names(tmpdata_cvd0)[1:2] <- c("t_cvds","cvd")

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

# Ssf_risk_treated_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# Ssf_risk_control_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# Tsf_risk_treated_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# Tsf_risk_control_cvd <- rep(NA,dim(tmpdata_cvd)[1])

# testIndexes <- which(folds==i,arr.ind=TRUE)
# traindat <- tmpdata_cvd[-testIndexes,]
# traindat1 <- traindat[traindat$INTENSIVE==1,c(1:2,4:dim(traindat)[2])]
# traindat0 <- traindat[traindat$INTENSIVE==0,c(1:2,4:dim(traindat)[2])]
# testdat  <- tmpdata_cvd[testIndexes,c(1:2,4:dim(tmpdata_cvd)[2])]

testIndexes_event <- which(folds_event==i,arr.ind=TRUE)
testIndexes_noevent <- which(folds_noevent==i,arr.ind=TRUE)
testIndexes <- as.numeric(c(rownames(tmpdata_cvd_event[testIndexes_event,]), rownames(tmpdata_cvd_noevent[testIndexes_noevent,])))
traindat <- tmpdata_cvd[-testIndexes, ]
traindat1 <- traindat[traindat$INTENSIVE==1,c(1:2,4:dim(traindat)[2])]
traindat0 <- traindat[traindat$INTENSIVE==0,c(1:2,4:dim(traindat)[2])]
testdat  <- tmpdata_cvd[testIndexes,c(1:2,4:dim(tmpdata_cvd)[2])]

# T-Learner
# Tuning hyperparameters
alphas <- c(0.001, 0.005, 0.01)
nodesizes <- c(5,15,25)
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
        eventtimes <- sort(unique(traindat1$t_cvds[traindat1$cvd==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit1 <- survival_forest(X = traindat1[,3:dim(traindat1)[2]],
                                   Y = traindat1$t_cvds,
                                   D = traindat1$cvd,
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
        eventtimes <- sort(unique(traindat0$t_cvds[traindat0$cvd==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit0 <- survival_forest(X = traindat0[,3:dim(traindat0)[2]],
                                   Y = traindat0$t_cvds,
                                   D = traindat0$cvd,
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
        GND_T1 <- GND.calib.simp(pred=surf1, tvar=traindat1$t_cvds, out=traindat1$cvd,
                                 cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("T1-learner: ",c(ntree, mtry, nodesize, alpha, GND_T1[3])))
        pet1 <- rbind(pet1, c(ntree, mtry, nodesize, alpha, GND_T1[3]))
        
        # Calibration
        Spredgrp <- as.numeric(cut2(surf0,g=5))
        GND_T0 <- GND.calib.simp(pred=surf0, tvar=traindat0$t_cvds, out=traindat0$cvd,
                                 cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("T0-learner: ",c(ntree, mtry, nodesize, alpha, GND_T0[3])))
        pet0 <- rbind(pet0, c(ntree, mtry, nodesize, alpha, GND_T0[3]))
      }
    }
  }
}
# Final treated model 
pet1 <- data.frame(pet1); colnames(pet1)[dim(pet1)[2]] <- "error"
opt_params1 <- as.list(pet1[which.max(pet1$error),])
eventtimes <- sort(unique(traindat1$t_cvds[traindat1$cvd==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit1 <- survival_forest(X = traindat1[,3:dim(traindat1)[2]],
                           Y = traindat1$t_cvds,
                           D = traindat1$cvd,
                           failure.times = failure.times,
                           num.trees = opt_params1[[1]],
                           mtry = opt_params1[[2]],
                           sample.fraction = 1,
                           min.node.size = opt_params1[[3]],
                           alpha = opt_params1[[4]],
                           honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                           seed = 99)
surf1 <- predict(sf_fit1, newdata = testdat[,3:dim(testdat)[2]])
Tsf_risk_treated_cvd <- 1-surf1$predictions[, which.min(abs(surf1$failure.times-t50))]

# Final control model 
pet0 <- data.frame(pet0); colnames(pet0)[dim(pet0)[2]] <- "error"
opt_params0 <- as.list(pet0[which.max(pet0$error),])
eventtimes <- sort(unique(traindat0$t_cvds[traindat0$cvd==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit0 <- survival_forest(X = traindat0[,3:dim(traindat0)[2]],
                           Y = traindat0$t_cvds,
                           D = traindat0$cvd,
                           failure.times = failure.times,
                           num.trees = opt_params0[[1]],
                           mtry = opt_params0[[2]],
                           sample.fraction = 1,
                           min.node.size = opt_params0[[3]],
                           alpha = opt_params0[[4]],
                           honesty.prune.leaves = F, # improve performance on small/marginally powered data, but requires more trees
                           seed = 99)
surf0 <- predict(sf_fit0, newdata = testdat[,3:dim(testdat)[2]])
Tsf_risk_control_cvd <- 1-surf0$predictions[, which.min(abs(surf0$failure.times-t50))]


# S-Learner
traindat <- tmpdata_cvd[-testIndexes,]
testdat1 <- tmpdata_cvd1[testIndexes, ]
testdat0 <- tmpdata_cvd0[testIndexes, ]

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
        eventtimes <- sort(unique(traindat$t_cvds[traindat$cvd==1]))
        nn <- round(length(eventtimes)/100,0)
        failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
        sf_fit <- survival_forest(X = traindat[,3:dim(traindat)[2]],
                                  Y = traindat$t_cvds,
                                  D = traindat$cvd,
                                  failure.times = failure.times,
                                  sample.fraction = 0.8,
                                  num.trees = ntree,
                                  mtry = mtry,
                                  min.node.size = nodesize,
                                  alpha = alpha,
                                  honesty.prune.leaves = F, 
                                  seed = 99)
        surf <- 1-sf_fit$predictions[, which.min(abs(sf_fit$failure.times-t50))]
        
        # Calibration
        Spredgrp <- as.numeric(cut2(surf,g=5))
        GND_S <- GND.calib.simp(pred=surf, tvar=traindat$t_cvds, out=traindat$cvd,
                                cens.t=t50, groups=Spredgrp, adm.cens=t50)
        print(paste0("S-learner: ", c(ntree, mtry, nodesize, alpha, GND_S[3])))
        pet <- rbind(pet, c(ntree, mtry, nodesize, alpha, GND_S[3]))
      }
    }
  }
}

# Final S-learner model 
pet <- data.frame(pet); colnames(pet)[dim(pet)[2]] <- "error"
opt_params <- as.list(pet[which.max(pet$error),])
eventtimes <- sort(unique(traindat$t_cvds[traindat$cvd==1]))
nn <- round(length(eventtimes)/100,0)
failure.times <- eventtimes[seq(1,length(eventtimes),nn)]
sf_fit <- survival_forest(X = traindat[,3:dim(traindat)[2]],
                          Y = traindat$t_cvds,
                          D = traindat$cvd,
                          failure.times = failure.times,
                          sample.fraction = 1,
                          num.trees = opt_params[[1]],
                          mtry = opt_params[[2]],
                          min.node.size = opt_params[[3]],
                          alpha = opt_params[[4]],
                          honesty.prune.leaves = F, 
                          seed = 99)
surf1 <- predict(sf_fit, newdata = testdat1[,3:dim(testdat1)[2]])
surf0 <- predict(sf_fit, newdata = testdat0[,3:dim(testdat0)[2]])
Ssf_risk_treated_cvd <- 1-surf1$predictions[, which.min(abs(surf1$failure.times-t50))]
Ssf_risk_control_cvd <- 1-surf0$predictions[, which.min(abs(surf0$failure.times-t50))]

risk_sf_cvd <- data.frame(testIndexes, Ssf_risk_control_cvd,Ssf_risk_treated_cvd,Tsf_risk_control_cvd,Tsf_risk_treated_cvd)
write.csv(risk_sf_cvd, paste0("./cvd_tuned_results/stratified_new_split_sf_cvd",i,".csv"), row.names = F)
# best_params <- data.frame(unlist(opt_params1), unlist(opt_params0), unlist(opt_params))
# write.csv(best_params, paste0("./cvd_tuned_results/sf_best_params_cvd",seed_num,".csv"), row.names = F)


