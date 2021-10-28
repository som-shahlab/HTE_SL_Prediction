
##----------------------        X-learner for HTE Estimates      -------------------##
##----------------------       Crystal Xu         03/15/2021     -------------------##
library(survival);library(ggplot2);library(gtsummary);library(boot);library(survcomp);library(nricens); library(png)
library(gridExtra);library(grid);library(lattice); library(riskRegression); library(ranger)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_test <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
tmpdata <- accord_test[complete.cases(accord_test),]; dim(tmpdata) # 4448 43
tmpdata$t_cvds <- ceiling(tmpdata$t_cvds)
tmpdata$t_saes <- ceiling(tmpdata$t_saes)
t50 <- floor(365.25*3.26)

# SPRINT as testing data
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
sprint_test <- merge(blvars,outcome, by="MASKID",all=T)
sprint_test <- sprint_test[complete.cases(sprint_test),]; dim(sprint_test) # 8969
testXvar <- sprint_test[,2:dim(blvars)[2]]

##------------------------------------- Internal Validation ----------------------------------------##
setwd("../Analyses Stanford Team/Analysis Results/ACCORD_internal_results") 
cvdpath <- "accord_IV_TL_deepsurv_cvd"
cvdpath <- "accord_IV_TL_gbm_cvd"
cvdpath <- "accord_IV_TL_sf_cvd"
cvdpath <- "accord_IV_T1bart_cvd"

# Load estimates
#sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)
sum_cvd <- NULL
for (i in 1:10){
  dat <- read.csv(paste0(cvdpath,i,".csv"))
  sum_cvd <- rbind(sum_cvd,dat)
}
sum_cvd_T0 <- sum_cvd[order(sum_cvd$testIndexes),2]; head(sum_cvd_T0)
sum_cvd_T1 <- sum_cvd[order(sum_cvd$testIndexes),2]; head(sum_cvd_T1)

# X-learner CATE
tmpdata$T_risk_control_cvd <- sum_cvd_T0
tmpdata$T_risk_treated_cvd <- sum_cvd_T1

# Indicator for not being censored/had event by time t50
tmpdata$uncensored_cvd <- ifelse(tmpdata$cvd==1|tmpdata$t_cvds>t50,1,0)

# Compute IPCW 
vars_all <- colnames(tmpdata)[2:39]
nfold <- 10
n <- round(dim(tmpdata)[1]/nfold,0)
cvd_ipcw <- NULL
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
  bh <- bh_dat[bh_dat$time==bh_dat$time[which.min(abs(bh_dat$time-t50))],]$hazard
  est_r <- predict(cvdwfit, newdata = testdat[,vars_all], type="risk")
  cvd_ipcw <- c(cvd_ipcw, 1/exp(-bh)^est_r)
}

# Trim extreme weights
cvd_ipcw[cvd_ipcw > 10] <- 10
tmpdata$cvd_ipcw <- cvd_ipcw

## Using an RF + IPCW to estimate the CATE function
# X-learner
vars <- vars_all[-which(vars_all=="INTENSIVE")]
traindatUC <- tmpdata[tmpdata$uncensored_cvd==1, ]
traindatUC$cvd <- ifelse(traindatUC$t_cvds>t50, 0, traindatUC$cvd)

fit1_cvd <- ranger(Y ~ ., data = data.frame(Y = traindatUC$cvd[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_cvd[traindatUC$INTENSIVE==1],
                                            traindatUC[traindatUC$INTENSIVE==1,vars]),
                   num.trees = 1000, mtry = 5, min.node.size = 1, regularization.factor = 0.5,
                   case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
tau1 <- predict(fit1_cvd, sprint_test[,vars])$predictions
cvd1_imp <- fit1_cvd$variable.importance

fit0_cvd <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_cvd[traindatUC$INTENSIVE==0] - traindatUC$cvd[traindatUC$INTENSIVE==0],
                                            traindatUC[traindatUC$INTENSIVE==0,vars]),
                   num.trees = 1000, mtry = 5, min.node.size = 1, regularization.factor = 0.5,
                   case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
tau0 <- predict(fit0_cvd, sprint_test[,vars])$predictions
cvd0_imp <- fit0_cvd$variable.importance
XL_RD_cvd <- (tau1+tau0)/2

# Evaluate 
cvd_S_RD <- data.frame(W = sprint_test$INTENSIVE, Y = sprint_test$t_cvds, D = sprint_test$cvd, pred_benefit = XL_RD_cvd)
CB_Scvd <- cforbenefit(cvd_S_RD);CB_Scvd
calib_Scvd <- HTEcalib(data = cvd_S_RD, times = t50, name = "");calib_Scvd

# Save results
respath <- "C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/SPRINT_external_results"
write.csv(XL_RD_cvd, paste0(respath,"/XL_bart_CVD.csv"), row.names = F)

# Summarize variable importance
cvd_imp <- (cvd1_imp + cvd0_imp)/2
cvd_imp_dat <- data.frame(vars, cvd_imp)
cvd_imp_dat <- cvd_imp_dat[order(-cvd_imp_dat$cvd_imp),]
sum_imp <- data.frame(cvd_imp_dat)
colnames(sum_imp) <- c("cvd_vars", "cvd_imp")
write.csv(sum_imp, paste0(respath,"/XL_imp_bart.csv"), row.names = F)




