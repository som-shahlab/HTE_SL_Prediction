
##----------------------        X-learner for HTE Estimates      -------------------##
##----------------------       Crystal Xu         03/15/2021     -------------------##
library(survival);library(ggplot2);library(gtsummary);library(boot);library(survcomp);library(nricens); library(png)
library(gridExtra);library(grid);library(lattice); library(riskRegression); library(ranger)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969
t50 <- floor(365.25*3.26)

##------------------------------------- Internal Validation ----------------------------------------##
setwd("../Analyses Stanford Team/Analysis Results/Tuned Results") 
cvdpath <- "bart_cvd_tuned_results/stratified_sum_cvd_bart"
saepath <- "bart_sae_tuned_results/stratified_sum_sae_bart"
cvdpath <- "deepsurv_cvd_tuned_results/deepsurv_tuned_cvd_sum"
saepath <- "deepsurv_sae_tuned_results/deepsurv_tuned_sae_sum"
cvdpath <- "gbm_cvd_tuned_results/stratified_split_gbm_cvd_sum"
saepath <- "gbm_sae_tuned_results/stratified_split_gbm_sae_sum"
cvdpath <- "sf_cvd_tuned_results/stratified_split_sf_cvd_sum"
saepath <- "sf_sae_tuned_results/stratified_split_sf_sae_sum"

# load estimates
sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)
sum_sae <- read.csv(paste0(saepath,".csv"));head(sum_sae)

# X-learner CATE
tmpdata$T_risk_control_cvd <- sum_cvd[,3]
tmpdata$T_risk_treated_cvd <- sum_cvd[,4]
tmpdata$T_risk_control_sae <- sum_sae[,3]
tmpdata$T_risk_treated_sae <- sum_sae[,4]

# Indicator for not being censored/had event by time t50
tmpdata$uncensored_cvd <- ifelse(tmpdata$cvd==1|tmpdata$t_cvds>t50,1,0)
tmpdata$uncensored_sae <- ifelse(tmpdata$sae==1|tmpdata$t_saes>t50,1,0)

# Compute IPCW (out-of-sample)
vars_all <- colnames(tmpdata)[2:39]
nfold <- 10
n <- round(dim(tmpdata)[1]/nfold,0)
cvd_ipcw <- sae_ipcw <- NULL
for (i in 1:nfold){
  if(i==nfold){
    traindat <- tmpdata[-((1+(i-1)*n):(dim(tmpdata)[1])),]
    testdat  <- tmpdata[((1+(i-1)*n):(dim(tmpdata)[1])),]
  }else{
    traindat <- tmpdata[-((1+(i-1)*n):(i*n)),]
    testdat  <- tmpdata[((1+(i-1)*n):(i*n)),]
  }

  Wdata <- traindat[,colnames(traindat) %in% c(vars_all, "t_cvds","cvd")]
  cvdwfit <- coxph(Surv(t_cvds, 1-cvd) ~., data = Wdata)
  bh_dat <- basehaz(cvdwfit)
  bh <- bh_dat[bh_dat$time==t50,]$hazard
  est_r <- predict(cvdwfit, newdata = testdat[,vars_all], type="risk")
  cvd_ipcw <- c(cvd_ipcw, 1/exp(-bh)^est_r)

  Wdata <- traindat[,colnames(traindat) %in% c(vars_all, "t_saes","sae")]
  saewfit <- coxph(Surv(t_saes, 1-sae) ~., data = Wdata)
  bh_dat <- basehaz(saewfit)
  bh <- bh_dat[bh_dat$time==t50,]$hazard
  est_r <- predict(saewfit, newdata = testdat[,vars_all], type="risk")
  sae_ipcw <- c(sae_ipcw, 1/exp(-bh)^est_r)
}

# Trim extreme weights
cvd_ipcw[cvd_ipcw > 10] <- 10
sae_ipcw[sae_ipcw > 10] <- 10
tmpdata$cvd_ipcw <- cvd_ipcw
tmpdata$sae_ipcw <- sae_ipcw

### Using an RF + IPCW to estimate the CATE function
vars <- vars_all[-which(vars_all=="INTENSIVE")]

# X-learner 10-fold CV
nfold <- 10
n <- round(dim(tmpdata)[1]/nfold,0)
XL_RD_sae <- XL_RD_cvd <- NULL
cvd1_imp <- cvd0_imp <- sae1_imp <- sae0_imp <- matrix(NA, length(vars), nfold)
for (i in 1:nfold){
  
  print(i)
  
  if(i==nfold){
    traindat <- tmpdata[-((1+(i-1)*n):(dim(tmpdata)[1])),]
    testdat  <- tmpdata[((1+(i-1)*n):(dim(tmpdata)[1])),]
  }else{
    traindat <- tmpdata[-((1+(i-1)*n):(i*n)),]
    testdat  <- tmpdata[((1+(i-1)*n):(i*n)),]
  }
  
  # CVD
  traindatUC <- traindat[traindat$uncensored_cvd==1, ]
  traindatUC$cvd <- ifelse(traindatUC$t_cvds>t50, 0, traindatUC$cvd)
  
  fit1 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$cvd[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_cvd[traindatUC$INTENSIVE==1],
                                          traindatUC[traindatUC$INTENSIVE==1,vars]),
                 num.trees = 1000, mtry = 5, min.node.size = 1, regularization.factor = 0.5,
                 case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
  tau1 <- predict(fit1, testdat[,vars])$predictions
  cvd1_imp[,i] <- fit1$variable.importance
  
  fit0 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_cvd[traindatUC$INTENSIVE==0] - traindatUC$cvd[traindatUC$INTENSIVE==0],
                                          traindatUC[traindatUC$INTENSIVE==0,vars]),
                 num.trees = 1000, mtry = 5, min.node.size = 1, regularization.factor = 0.5,
                 case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
  tau0 <- predict(fit0, testdat[,vars])$predictions
  cvd0_imp[,i] <- fit0$variable.importance
  
  XL_RD_cvd <- c(XL_RD_cvd, (tau1+tau0)/2)
  
  # SAE
  traindatUC <- traindat[traindat$uncensored_sae==1, ]
  traindatUC$sae <-  ifelse(traindatUC$t_saes>t50, 0, traindatUC$sae)
  
  fit1 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$sae[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_sae[traindatUC$INTENSIVE==1],
                                          traindatUC[traindatUC$INTENSIVE==1,vars]),
                 num.trees = 1000, mtry = 5, min.node.size = 1, regularization.factor = 0.5,
                 case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
  tau1 <- predict(fit1, testdat[,vars])$predictions
  sae1_imp[,i] <- fit1$variable.importance
  
  fit0 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_sae[traindatUC$INTENSIVE==0] - traindatUC$sae[traindatUC$INTENSIVE==0],
                                          traindatUC[traindatUC$INTENSIVE==0,vars]),
                 num.trees = 1000, mtry = 5, min.node.size = 1, regularization.factor = 0.5,
                 case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
  tau0 <- predict(fit0, testdat[,vars])$predictions
  sae0_imp[,i] <- fit0$variable.importance
  
  XL_RD_sae <- c(XL_RD_sae, (tau1+tau0)/2)
  index <- c(index,testIndexes)
}

# Evaluate 
cvd_S_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = XL_RD_cvd_ordered)
CB_Scvd <- cforbenefit(cvd_S_RD);CB_Scvd
calib_Scvd <- HTEcalib(data = cvd_S_RD, times = t50, name = "CVD Risk Reduction, X-learner BART");calib_Scvd

sae_S_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_saes, D = tmpdata$sae, pred_benefit = XL_RD_sae_ordered)
CB_Ssae <- cforbenefit(sae_S_RD);CB_Ssae
calib_Ssae <- HTEcalib(data = sae_S_RD, times = t50, name = "SAE Risk Increase, X-learner BART");calib_Ssae

# Save results
sum_cvd$XL_RD_cvd <- XL_RD_cvd
sum_sae$XL_RD_sae <- XL_RD_sae
write.csv(sum_cvd,paste0(cvdpath,".csv"), row.names = F)
write.csv(sum_sae,paste0(saepath,".csv"), row.names = F)

cvd_imp <- (rowMeans(cvd1_imp) + rowMeans(cvd0_imp))/2
cvd_imp_dat <- data.frame(vars, cvd_imp)
cvd_imp_dat <- cvd_imp_dat[order(-cvd_imp_dat$cvd_imp),]

sae_imp <- (rowMeans(sae1_imp) + rowMeans(sae0_imp))/2
sae_imp_dat <- data.frame(vars, sae_imp)
sae_imp_dat <- sae_imp_dat[order(-sae_imp_dat$sae_imp),]

sum_imp <- data.frame(cvd_imp_dat, sae_imp_dat)
colnames(sum_imp) <- c("cvd_vars", "cvd_imp","sae_vars", "sae_imp")
write.csv(sum_imp,"IV_imp_sf.csv", row.names = F)


