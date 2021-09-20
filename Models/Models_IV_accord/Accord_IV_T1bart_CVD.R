##---------------------- Apply BART method to the actual SPRINT data -------------------##
##----------------------     Crystal Xu               02/23/2021     -------------------##
library(survival);library(BART)
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
source("eval_funs.R")
blvars <- read.csv("ACCORD_blvars.csv")
outcome <- read.csv("ACCORD_outcomes.csv")
tmpdat <- merge(blvars, outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 4448 43
tmpdata$t_cvds <- ceiling(tmpdata$t_cvds)
tmpdata$t_saes <- ceiling(tmpdata$t_saes)

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
traindat1 <- traindat[traindat$INTENSIVE==1, c(1:2,4:dim(traindat)[2])]
testdat <- tmpdata_cvd[testIndexes, c(1:2,4:dim(tmpdata_cvd)[2])]

# T-Learner
# if(i==nfold){
#   traindat <- tmpdata_cvd[-((1+(i-1)*n):dim(tmpdata_cvd)[1]),]
#   traindat1 <- traindat[traindat$INTENSIVE==1, c(1:2,4:dim(traindat)[2])]
#   testdat  <- tmpdata_cvd[((1+(i-1)*n):dim(tmpdata_cvd)[1]), c(1:2,4:dim(tmpdata_cvd)[2])]
# }else{
#   traindat <- tmpdata_cvd[-((1+(i-1)*n):(i*n)),]
#   traindat1 <- traindat[traindat$INTENSIVE==1, c(1:2,4:dim(traindat)[2])]
#   testdat  <- tmpdata_cvd[((1+(i-1)*n):(i*n)), c(1:2,4:dim(tmpdata_cvd)[2])]
# }

# # 5-fold CV for tuning hyperparameters
# ntrees <- c(50,100,200)
# pr_sds <- c(1,2,3)
# ntune <- 5
# n1 <- floor(dim(traindat1)[1]/ntune)
# n0 <- floor(dim(traindat0)[1]/ntune)

bartSurv1 <- mc.surv.bart(x.train=traindat1[,3:dim(traindat1)[2]],
                          times=ceiling(traindat1$t_cvds/30),
                          delta=traindat1$cvd,
                          ntree=50,
                          x.test=testdat[,3:dim(testdat)[2]],
                          seed=99, 
                          mc.cores=4)
surv_ests <- t(matrix(bartSurv1$surv.test.mean, nrow = length(bartSurv1$times)))
surv_est1 <- surv_ests[,which.min(abs(bartSurv1$times-t50))]
Tbart_risk_treated_cvd <- 1-surv_est1


# save results 
risk_bart_cvd <- data.frame(testIndexes, Tbart_risk_treated_cvd)
write.csv(risk_bart_cvd, paste0("./cvd_tuned_results/accord_IV_T1bart_cvd",i,".csv"), row.names = F)
