##---------------------- Apply BART method to the actual SPRINT data -------------------##
##----------------------     Crystal Xu               02/23/2021     -------------------##
library(survival);library(BART)
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
source("eval_funs.R")
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")

# Combine baseline variables and outcomes into one data frame 
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969

# Set counterfacutal treatment levels
Xvar <- tmpdata[, 2:dim(blvars)[2]]
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

#-------------------- Benefit model ----------------------------#
tmpdata_cvd <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar)
names(tmpdata_cvd)[1:2] <- c("t_cvds","cvd")
tmpdata_cvd1 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar1)
names(tmpdata_cvd1)[1:2] <- c("t_cvds","cvd")
tmpdata_cvd0 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar0)
names(tmpdata_cvd0)[1:2] <- c("t_cvds","cvd")

# BART Method 
# nfold <- 10 
# n <- round(dim(tmpdata_cvd)[1]/nfold,0)
t50 <- 3.26*12

# 10-fold CV to get the out-of-bag estimates (stratified to ensure equal number of events in each fold)
nfold <- 10
rownames(tmpdata_cvd) <- as.character(1:dim(tmpdata_cvd)[1])
rownames(tmpdata_cvd1) <- as.character(1:dim(tmpdata_cvd1)[1])
rownames(tmpdata_cvd0) <- as.character(1:dim(tmpdata_cvd0)[1])
tmpdata_cvd_event <- tmpdata_cvd[tmpdata_cvd$cvd==1,]
tmpdata_cvd_noevent <- tmpdata_cvd[tmpdata_cvd$cvd==0,]
folds_event <- cut(seq(1,nrow(tmpdata_cvd_event)),breaks=nfold,labels=FALSE)
folds_noevent <- cut(seq(1,nrow(tmpdata_cvd_noevent)),breaks=nfold,labels=FALSE)

testIndexes_event <- which(folds_event==i,arr.ind=TRUE)
testIndexes_noevent <- which(folds_noevent==i,arr.ind=TRUE)
testIndexes <- as.numeric(c(rownames(tmpdata_cvd_event[testIndexes_event,]), rownames(tmpdata_cvd_noevent[testIndexes_noevent,])))
traindat <- tmpdata_cvd[-testIndexes, ]
testdat1 <- tmpdata_cvd1[testIndexes,]
testdat0 <- tmpdata_cvd0[testIndexes,]


# S-Learner
# if(i==nfold){
#   traindat <- tmpdata_cvd[-((1+(i-1)*n):(dim(tmpdata_cvd)[1])),]
#   testdat <- tmpdata_cvd[(1+(i-1)*n):(dim(tmpdata_cvd)[1]),]
#   testdat1 <- tmpdata_cvd1[(1+(i-1)*n):(dim(tmpdata_cvd)[1]),]
#   testdat0 <- tmpdata_cvd0[(1+(i-1)*n):(dim(tmpdata_cvd)[1]),]
# }else{
#   traindat <- tmpdata_cvd[-((1+(i-1)*n):(i*n)),]
#   testdat <- tmpdata_cvd[(1+(i-1)*n):(i*n),]
#   testdat1 <- tmpdata_cvd1[(1+(i-1)*n):(i*n),]
#   testdat0 <- tmpdata_cvd0[(1+(i-1)*n):(i*n),]
# }
x.test.both <- rbind(testdat1,testdat0)

# # Tuning hyperparameters
# ntrees <- c(20, 50, 100)
# ntune <- 5
# n2 <- floor(dim(traindat)[1]/ntune)
# pet <- NULL
# for (j in 1:length(ntrees)){
#   ntree <- ntrees[j]
#   surv_est <- NULL
#   for (m in 1:ntune){
#     if(m==ntune){
#       train_dat <- traindat[-((1+(m-1)*n2):dim(traindat)[1]),]
#       valid_dat <- traindat[((1+(m-1)*n2):dim(traindat)[1]),]
#     }else{
#       train_dat <- traindat[-((1+(m-1)*n2):(m*n2)),]
#       valid_dat <- traindat[((1+(m-1)*n2):(m*n2)),]
#     }
#     bartSurv <- mc.surv.bart(x.train=train_dat[,3:dim(train_dat)[2]],
#                              times=ceiling(train_dat$t_cvds/30),
#                              delta=train_dat$cvd,
#                              ntree=ntree,
#                              x.test = valid_dat[,3:dim(valid_dat)[2]],
#                              seed = 99, 
#                              mc.cores = 4)
#     surv_ests <- t(matrix(bartSurv$surv.test.mean, nrow = length(bartSurv$times)))
#     surv_est <- 1-c(surv_est, surv_ests[,which.min(abs(bartSurv$times-t50))])
#   }
#   # Calibration
#   Spredgrp <- as.numeric(cut2(surv_est,g=10))
#   GND_S <- GND.calib(pred=surv_est, tvar=ceiling(traindat$t_cvds/30), out=traindat$cvd,
#                      cens.t=t50, groups=Spredgrp, adm.cens=t50)
#   print(c(ntree, ndpost, GND_S[3]))
#   pet <- rbind(pet, c(ntree, ndpost, GND_S[3]))
# }
# sum_pe <- data.frame(pet)
# colnames(sum_pe)[3] <- "ave_pe"
# opt_params <- as.list(sum_pe[which.max(sum_pe$ave_pe),1:2])

bartSurv <- mc.surv.bart(x.train=traindat[,3:dim(traindat)[2]],
                         times=ceiling(traindat$t_cvds/30),
                         delta=traindat$cvd,
                         ntree=50,
                         x.test = x.test.both[,3:dim(testdat1)[2]],
                         seed=99, 
                         mc.cores=4)
surv_ests <- t(matrix(bartSurv$surv.test.mean, nrow = length(bartSurv$times)))
surv_est1 <- surv_ests[1:dim(testdat1)[1],which.min(abs(sort(bartSurv$times)-t50))]
surv_est0 <- surv_ests[-(1:dim(testdat1)[1]),which.min(abs(sort(bartSurv$times)-t50))]
Sbart_risk_treated_cvd <- 1-surv_est1
Sbart_risk_control_cvd <- 1-surv_est0

# variable importance
vimp <- bartSurv$varcount.mean[-1]
varimp <- data.frame(names(vimp),vimp)
colnames(varimp) <- c("varnames","imp")

# save results 
risk_bart_cvd <- data.frame(testIndexes, Sbart_risk_control_cvd,Sbart_risk_treated_cvd)
write.csv(risk_bart_cvd, paste0("./cvd_tuned_results/stratified_Sbart_cvd_replica",i,".csv"), row.names = F)
write.csv(varimp, paste0("./cvd_tuned_results/Sbart_varimp",i,".csv"), row.names = F)
# best_params <- data.frame(unlist(opt_params))
# write.csv(best_params, paste0("./cvd_tuned_results/bart_best_params_cvd",i,".csv"), row.names = F)

# summary(Sbart_risk_treated_cvd-Sbart_risk_control_cvd)
# Sgbm_predinc_cvd <- Sbart_risk_treated_cvd*testdat$INTENSIVE + Sbart_risk_control_cvd*(1-testdat$INTENSIVE)
# Spredgrp <- as.numeric(cut2(Sgbm_predinc_cvd,g=3))
# GND_S <- GND.calib(pred=Sgbm_predinc_cvd, tvar=ceiling(testdat$t_cvds/30), out=testdat$cvd,
#                    cens.t=t50, groups=Spredgrp, adm.cens=t50)
# GND_S_pvalue <- round(GND_S[3:5],4);GND_S_pvalue

