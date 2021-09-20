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
# 
#----------------------------------- train using treated obs only --------------------------------#
# petable1 <- list(NA,ntune)
# for (m in 1:ntune){
#   print(paste0("OOB estimates fold ",i," hyperparameter tuning fold ",m, " -- T-learner"))
#   if(m==ntune){
#     train_dat1 <- traindat1[-((1+(m-1)*n1):dim(traindat1)[1]),]
#     valid_dat1 <- traindat1[((1+(m-1)*n1):dim(traindat1)[1]),]
#     train_dat0 <- traindat0[-((1+(m-1)*n0):dim(traindat0)[1]),]
#     valid_dat0 <- traindat0[((1+(m-1)*n0):dim(traindat0)[1]),]
#   }else{
#     train_dat1 <- traindat1[-((1+(m-1)*n1):(m*n1)),]
#     valid_dat1 <- traindat1[((1+(m-1)*n1):(m*n1)),]
#     train_dat0 <- traindat0[-((1+(m-1)*n0):(m*n0)),]
#     valid_dat0 <- traindat0[((1+(m-1)*n0):(m*n0)),]
#   }
#   
#   pet1 <- NULL
#   for (j in 1:length(ntrees)){
#     ntree <- ntrees[j]
#     for (k in 1:length(pr_sds)){
#       sd <- pr_sds[k]
# 
#       # train using treated obs only
#       bartSurv1 <- surv.bart(x.train=train_dat1[,3:dim(train_dat1)[2]],
#                              times=train_dat1$t_cvds,
#                              delta=train_dat1$cvd,
#                              ntree=ntree,
#                              k=sd,
#                              power=2,
#                              base=0.95,
#                              x.test=valid_dat1[,3:dim(valid_dat1)[2]],
#                              ndpost=100, nskip=50, seed=99)
#       surv_ests <- matrix(bartSurv1$surv.test.mean, ncol = length(bartSurv1$times))
#       surv_est1 <- surv_ests[,which.min(abs(bartSurv1$times-t50))]
# 
#       # negative log-likelihood
#       D1 <- ifelse(valid_dat1$t_cvds>t50,0,valid_dat1$cvd)
#       NLL1 <- -1*mean((1-D1)*log(surv_est1)+D1*log(1-surv_est1))
#       pet1 <- rbind(pet1, c(ntree, sd, NLL1))
# 
#       # # IBS
#       # PredError1 <- pec(object=bartSurv1,
#       #                   formula = as.formula("Surv(t_cvds, cvd)~1"),
#       #                   cens.model="marginal",
#       #                   start=min(valid_dat1$t_cvds),
#       #                   maxtime = t50,
#       #                   exact = F,
#       #                   exactness = 100,
#       #                   reference = F,
#       #                   data=valid_dat1)
#       # PredError1$AppErr[is.na(PredError1$AppErr)]<-0
#       # pe1 <- crps(PredError1, times=t50, start=min(valid_dat1$t_cvds))  # OOB prediction error based on Graf's integrated BS (1999)
#     }
#   }
#   petable1[[m]] <- pet1
# }

# # Final model - treated obs only
# sum_pe1 <- data.frame(petable1[[1]],petable1[[2]][,3],petable1[[3]][,3],petable1[[4]][,3],petable1[[5]][,3])
# sum_pe1$ave_pe <- rowMeans(sum_pe1[,3:7])
# opt_params1 <- as.list(sum_pe1[which.min(sum_pe1$ave_pe),1:2])
# bartSurv1 <- surv.bart(x.train=traindat1[,3:dim(traindat1)[2]],
#                        times=traindat1$t_cvds,
#                        delta=traindat1$cvd,
#                        events = traindat1$t_cvds[traindat1$cvd==1],
#                        ntree=opt_params1[[1]],
#                        k=opt_params1[[2]],
#                        power=2,
#                        base=0.95,
#                        x.test=testdat[,3:dim(testdat)[2]],
#                        ndpost=100, nskip=100, seed=99)

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
write.csv(risk_bart_cvd, paste0("./cvd_tuned_results/stratified__T1bart_cvd",i,".csv"), row.names = F)
