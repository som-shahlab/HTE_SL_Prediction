
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

#-------------------- Risk model ----------------------------#
tmpdata_sae <- data.frame(tmpdata$t_saes, tmpdata$sae, Xvar)
names(tmpdata_sae)[1:2] <- c("t_saes","sae")
tmpdata_sae1 <- data.frame(tmpdata$t_saes, tmpdata$sae, Xvar1)
names(tmpdata_sae1)[1:2] <- c("t_saes","sae")
tmpdata_sae0 <- data.frame(tmpdata$t_saes, tmpdata$sae, Xvar0)
names(tmpdata_sae0)[1:2] <- c("t_saes","sae")

# BART Method 
# nfold <- 10 
# n <- round(dim(tmpdata_sae)[1]/nfold,0)
t50 <- 3.26*12

# 10-fold CV to get the out-of-bag estimates (stratified to ensure equal number of events in each fold)
nfold <- 10
rownames(tmpdata_sae) <- as.character(1:dim(tmpdata_sae)[1])
rownames(tmpdata_sae1) <- as.character(1:dim(tmpdata_sae1)[1])
rownames(tmpdata_sae0) <- as.character(1:dim(tmpdata_sae0)[1])
tmpdata_sae_event <- tmpdata_sae[tmpdata_sae$sae==1,]
tmpdata_sae_noevent <- tmpdata_sae[tmpdata_sae$sae==0,]
folds_event <- cut(seq(1,nrow(tmpdata_sae_event)),breaks=nfold,labels=FALSE)
folds_noevent <- cut(seq(1,nrow(tmpdata_sae_noevent)),breaks=nfold,labels=FALSE)

testIndexes_event <- which(folds_event==i,arr.ind=TRUE)
testIndexes_noevent <- which(folds_noevent==i,arr.ind=TRUE)
testIndexes <- as.numeric(c(rownames(tmpdata_sae_event[testIndexes_event,]), rownames(tmpdata_sae_noevent[testIndexes_noevent,])))
traindat <- tmpdata_sae[-testIndexes, ]
traindat0 <- traindat[traindat$INTENSIVE==0, c(1:2,4:dim(traindat)[2])]
testdat <- tmpdata_sae[testIndexes, c(1:2,4:dim(tmpdata_sae)[2])]



# T-Learner
# if(i==nfold){
#   traindat <- tmpdata_sae[-((1+(i-1)*n):dim(tmpdata_sae)[1]),]
#   traindat0 <- traindat[traindat$INTENSIVE==0, c(1:2,4:dim(traindat)[2])]
#   testdat  <- tmpdata_sae[((1+(i-1)*n):dim(tmpdata_sae)[1]), c(1:2,4:dim(tmpdata_sae)[2])]
# }else{
#   traindat <- tmpdata_sae[-((1+(i-1)*n):(i*n)),]
#   traindat0 <- traindat[traindat$INTENSIVE==0, c(1:2,4:dim(traindat)[2])]
#   testdat  <- tmpdata_sae[((1+(i-1)*n):(i*n)), c(1:2,4:dim(tmpdata_sae)[2])]
# }

# # 5-fold CV for tuning hyperparameters
# ntrees <- c(50,100,200)
# pr_sds <- c(2)
# ntune <- 5
# n1 <- floor(dim(traindat1)[1]/ntune)
# n0 <- floor(dim(traindat0)[1]/ntune)
# 
# #----------------------------------- train using treated obs only --------------------------------#
# petable0 <- list(NA,ntune)
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
#   pet0 <- NULL
#   for (j in 1:length(ntrees)){
#     ntree <- ntrees[j]
#     for (k in 1:length(pr_sds)){
#       sd <- pr_sds[k]
# 
#       # train using control obs only
#       bartSurv0 <- surv.bart(x.train=train_dat0[,3:dim(train_dat0)[2]],
#                              times=train_dat0$t_saes,
#                              delta=train_dat0$sae,
#                              ntree=ntree,
#                              k=sd,
#                              power=2,
#                              base=0.95,
#                              x.test=valid_dat0[,3:dim(valid_dat0)[2]],
#                              ndpost=100, nskip=50, seed=99)
#       surv_ests <- matrix(bartSurv0$surv.test.mean, ncol = length(bartSurv0$times))
#       surv_est0 <- surv_ests[,which.min(abs(bartSurv0$times-t50))]
# 
#       # # negative log-likelihood
#       # D0 <- ifelse(valid_dat0$t_saes>t50,0,valid_dat0$sae)
#       # NLL0 <- -1*mean((1-D0)*log(surv_est0)+D0*log(1-surv_est0))
#       # pet0 <- rbind(pet0, c(ntree, sd, NLL0))
# 
#       # Calibration
#               Spredgrp <- as.numeric(cut2(surv_est,g=10))
#               GND_S <- GND.calib(pred=surv_est, tvar=traindat$t_saes, out=traindat$sae,
#                                  cens.t=t50, groups=Spredgrp, adm.cens=t50)
#               print(c(ntree, sd, base, ndpost, GND_S[3]))
#               pet <- rbind(pet, c(ntree, sd, base, ndpost, GND_S[3]))
# 
#     }
#   }
#   petable0[[m]] <- pet0
# }

# # Final model - control obs only
# sum_pe0 <- data.frame(petable0[[1]],petable0[[2]][,3],petable0[[3]][,3],petable0[[4]][,3],petable0[[5]][,3])
# sum_pe0$ave_pe <- rowMeans(sum_pe0[,3:7])
# opt_params0 <- as.list(sum_pe0[which.min(sum_pe0$ave_pe),1:2])
# bartSurv0 <- surv.bart(x.train=traindat0[,3:dim(traindat0)[2]],
#                        times=traindat0$t_saes,
#                        delta=traindat0$sae,
#                        events = traindat0$t_saes[traindat0$sae==1],
#                        ntree=opt_params0[[1]],
#                        k=opt_params0[[2]],
#                        power=2,
#                        base=0.95,
#                        x.test=testdat[,3:dim(testdat)[2]],
#                        ndpost=100, nskip=100, seed=99)

bartSurv0 <- mc.surv.bart(x.train=traindat0[,3:dim(traindat0)[2]],
                       times=ceiling(traindat0$t_saes/30),
                       delta=traindat0$sae,
                       ntree=50,
                       x.test=testdat[,3:dim(testdat)[2]],
                       seed=99, 
                       mc.cores=4)
surv_ests <- t(matrix(bartSurv0$surv.test.mean, nrow = length(bartSurv0$times)))
surv_est0 <- surv_ests[,which.min(abs(bartSurv0$times-t50))]
Tbart_risk_control_sae <- 1-surv_est0

# save results 
risk_bart_sae <- data.frame(testIndexes, Tbart_risk_control_sae)
write.csv(risk_bart_sae, paste0("./sae_tuned_results/stratified_T0bart_sae",i,".csv"), row.names = F)



# # SBART CATE 
# Sbart_RD_sae <- Sbart_risk_treated_sae - Sbart_risk_control_sae
# 
# # TBART CATE 
# Tbart_RD_sae <- Tbart_risk_treated_sae - Tbart_risk_control_sae
# 
# # X-learner + BART CATE
# tmpdata_sae$Tbart_risk_control_sae <- Tbart_risk_control_sae
# tmpdata_sae$Tbart_risk_treated_sae <- Tbart_risk_treated_sae
# tmpdata_sae_treated <- tmpdata_sae[tmpdata_sae$INTENSIVE==1, ] 
# tmpdata_sae_control <- tmpdata_sae[tmpdata_sae$INTENSIVE==0, ] 
# tmpdata_sae_treated$ITE_treated <- tmpdata_sae_treated$sae - tmpdata_sae_treated$Tbart_risk_control_sae
# tmpdata_sae_control$ITE_control <- tmpdata_sae_control$Tbart_risk_treated_sae - tmpdata_sae_control$sae
# 
# # using a linear regression as the base learner
# newdat <- tmpdata_sae[4:56]
# tmpdata_sae_treated <- tmpdata_sae_treated[,c(4:56,59)]
# fit1 <- lm(ITE_treated ~ ., data = tmpdata_sae_treated)
# tau1 <- predict(fit1, newdata=newdat, type='response')
# tmpdata_sae_control <- tmpdata_sae_control[,c(4:56,59)]
# fit0 <- lm(ITE_control ~ ., data = tmpdata_sae_control)
# tau0 <- predict(fit0, newdata=newdat, type='response')
# 
# # XL BART CATE
# XLbart_RD_sae <- (tau1+tau0)/2
# sum_bart_sae <- data.frame(Sbart_risk_treated_sae, Sbart_risk_control_sae, Sbart_RD_sae,
#                            Tbart_risk_treated_sae, Tbart_risk_control_sae, Tbart_RD_sae,
#                            tau1, tau0, XLbart_RD_sae)
# write.csv(sum_bart_sae, "sum_bart_sae.csv", row.names = F)




