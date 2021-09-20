
##----------------------        X-learner for HTE Estimates      -------------------##
##----------------------       Crystal Xu         03/15/2021     -------------------##
library(survival);library(ggplot2);library(gtsummary);library(boot);library(survcomp);library(nricens); library(png)
library(gridExtra);library(grid);library(lattice); library(riskRegression); library(ranger)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
# blvars <- read.csv("blvars.csv")
# outcome <- read.csv("outcome.csv")
# tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
# tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969
t50 <- floor(365.25*3.26)

# ACCORD as testing data
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_test <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_test <- accord_test[complete.cases(accord_test),]; dim(accord_test) # 4448 43
accord_test$t_cvds <- ceiling(accord_test$t_cvds)
accord_test$t_saes <- ceiling(accord_test$t_saes)
testXvar <- accord_test[,2:dim(accord_blvars)[2]]
tmpdata <- accord_test

# For combined results 
# tmpdata <- read.csv("sprint_accord.csv")

##------------------------------------- Internal Validation ----------------------------------------##
setwd("../Analyses Stanford Team/Analysis Results/External Results") 
cvdpath <- "bart_cvd_tuned_results/stratified_sum_cvd_bart"
saepath <- "bart_sae_tuned_results/stratified_sum_sae_bart"

cvdpath <- "external_bart_cvd_sum"
saepath <- "external_bart_sae_sum"
cvdpath <- "external_sf_cvd_sum"
saepath <- "external_sf_sae_sum"
cvdpath <- "deepsurv_external_cvd_sum"
saepath <- "deepsurv_external_sae_sum"
cvdpath <- "deepsurv_external_cvd"
saepath <- "deepsurv_external_sae"

# Deepsurv load estimates
sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)
#sum_cvd$index <- sum_cvd$index + 1
#sum_cvd <- sum_cvd[order(sum_cvd$index),];head(sum_cvd)

sum_sae <- read.csv(paste0(saepath,".csv"));head(sum_sae)
#sum_sae$index <- sum_sae$index + 1
#sum_sae <- sum_sae[order(sum_sae$index),];head(sum_sae)

# sum_cvd <- NULL
# for (i in 1:10){
#   tmpres <- read.csv(paste0(cvdpath,i,".csv"))
#   sum_cvd <- rbind(sum_cvd, tmpres)
# }
# sum_cvd <- sum_cvd[order(sum_cvd$testIndexes),];head(sum_cvd)
# 
# sum_sae <- NULL
# for (i in 1:10){
#   tmpres <- read.csv(paste0(saepath,i,".csv"))
#   sum_sae <- rbind(sum_sae, tmpres)
# }
# sum_sae <- sum_sae[order(sum_sae$testIndexes),];head(sum_sae)


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
  bh <- bh_dat[bh_dat$time==1191,]$hazard
  est_r <- predict(cvdwfit, newdata = testdat[,vars_all], type="risk")
  cvd_ipcw <- c(cvd_ipcw, 1/exp(-bh)^est_r)

  Wdata <- traindat[,colnames(traindat) %in% c(vars_all, "t_saes","sae")]
  saewfit <- coxph(Surv(t_saes, 1-sae) ~., data = Wdata)
  bh_dat <- basehaz(saewfit)
  bh <- bh_dat[bh_dat$time==1217,]$hazard
  est_r <- predict(saewfit, newdata = testdat[,vars_all], type="risk")
  sae_ipcw <- c(sae_ipcw, 1/exp(-bh)^est_r)
}

# Trim extreme weights
cvd_ipcw[cvd_ipcw > 5] <- 5
sae_ipcw[sae_ipcw > 5] <- 5
tmpdata$cvd_ipcw <- cvd_ipcw
tmpdata$sae_ipcw <- sae_ipcw

### Using an RF + IPCW to estimate the CATE function
# X-learner out-of-bag estimates
vars <- vars_all[-which(vars_all=="INTENSIVE")]
traindatUC <- tmpdata[tmpdata$uncensored_cvd==1, ]
traindatUC$cvd <- ifelse(traindatUC$t_cvds>t50, 0, traindatUC$cvd)
fit1_cvd <- ranger(Y ~ ., data = data.frame(Y = traindatUC$cvd[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_cvd[traindatUC$INTENSIVE==1],
                                        traindatUC[traindatUC$INTENSIVE==1,vars]),
               keep.inbag = TRUE, num.trees = 1000, mtry = 6, min.node.size = 1,
               case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
tau1 <- predict(fit1_cvd, tmpdata[,vars])$predictions
tau1[tmpdata$uncensored_cvd==1 & tmpdata$INTENSIVE==1] <- fit1_cvd$predictions
cvd1_imp <- fit1_cvd$variable.importance

fit0_cvd <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_cvd[traindatUC$INTENSIVE==0] - traindatUC$cvd[traindatUC$INTENSIVE==0],
                                        traindatUC[traindatUC$INTENSIVE==0,vars]),
               keep.inbag = TRUE, num.trees = 1000, mtry = 6, min.node.size = 1,
               case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
tau0 <- predict(fit0_cvd, tmpdata[,vars])$predictions
tau0[tmpdata$uncensored_cvd==1 & tmpdata$INTENSIVE==0] <- fit0_cvd$predictions
cvd0_imp <- fit0_cvd$variable.importance
XL_RD_cvd <- (tau1+tau0)/2

# SAE
traindatUC <- tmpdata[tmpdata$uncensored_sae==1, ]
traindatUC$sae <- ifelse(traindatUC$t_saes>t50, 0, traindatUC$sae)
fit1_sae <- ranger(Y ~ ., data = data.frame(Y = traindatUC$sae[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_sae[traindatUC$INTENSIVE==1],
                                            traindatUC[traindatUC$INTENSIVE==1,vars]),
                   keep.inbag = TRUE, num.trees = 1000, mtry = 6, min.node.size = 1,
                   case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
tau1 <- predict(fit1_sae, tmpdata[,vars])$predictions
tau1[tmpdata$uncensored_sae==1 & tmpdata$INTENSIVE==1] <- fit1_sae$predictions
sae1_imp <- fit1_sae$variable.importance

fit0_sae <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_sae[traindatUC$INTENSIVE==0] - traindatUC$sae[traindatUC$INTENSIVE==0],
                                            traindatUC[traindatUC$INTENSIVE==0,vars]),
                   keep.inbag = TRUE, num.trees = 1000, mtry = 6, min.node.size = 1,
                   case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
tau0 <- predict(fit0_sae, tmpdata[,vars])$predictions
tau0[tmpdata$uncensored_sae==1 & tmpdata$INTENSIVE==0] <- fit0_sae$predictions
sae0_imp <- fit0_sae$variable.importance
XL_RD_sae <- (tau1+tau0)/2

# Evaluate 
T_RD_cvd <- sum_cvd[,4] - sum_cvd[,3]
cvd_S_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = XL_RD_cvd)
CB_Scvd <- cforbenefit(cvd_S_RD);CB_Scvd
calib_Scvd <- HTEcalib(data = cvd_S_RD, times = t50, name = "");calib_Scvd

sae_S_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_saes, D = tmpdata$sae, pred_benefit = XL_RD_sae)
CB_Ssae <- cforbenefit(sae_S_RD);CB_Ssae
calib_Ssae <- HTEcalib(data = sae_S_RD, times = t50, name = "");calib_Ssae

# Save results
sum_cvd$XL_RD_cvd <- XL_RD_cvd
sum_sae$XL_RD_sae <- XL_RD_sae
write.csv(sum_cvd,paste0(cvdpath,"_sum.csv"), row.names = F)
write.csv(sum_sae,paste0(saepath,"_sum.csv"), row.names = F)

# Summarize variable importance
cvd_imp <- (cvd1_imp + cvd0_imp)/2
cvd_imp_dat <- data.frame(vars, cvd_imp)
cvd_imp_dat <- cvd_imp_dat[order(-cvd_imp_dat$cvd_imp),]

sae_imp <- (sae1_imp + sae0_imp)/2
sae_imp_dat <- data.frame(vars, sae_imp)
sae_imp_dat <- sae_imp_dat[order(-sae_imp_dat$sae_imp),]

sum_imp <- data.frame(cvd_imp_dat, sae_imp_dat)
colnames(sum_imp) <- c("cvd_vars", "cvd_imp","sae_vars", "sae_imp")
write.csv(sum_imp,"ex_imp_deepsurv.csv", row.names = F)

# Save all the models for building HTE calculator
saveRDS(fit1_cvd, "./fit1_cvd.rds")
saveRDS(fit0_cvd, "./fit0_cvd.rds")
saveRDS(fit1_sae, "./fit1_sae.rds")
saveRDS(fit0_sae, "./fit0_sae.rds")


### External Validation (HTE Calculator)
tmpdataUC <- tmpdata[tmpdata$uncensored==1, ]

# Tune ranger for CATE cvd outcome
nodesizes <- c(1,35,100,200)
ntrees <- c(1000)
mtrys <- 6
pet1 <- pet0 <- NULL
for (j in 1:length(ntrees)){
  num.trees <- ntrees[j]
  for (k in 1:length(mtrys)){
    mtry <- mtrys[k]
    for (l in 1:length(nodesizes)){
      min.node.size <- nodesizes[l]
      fit1 <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$cvd[tmpdataUC$INTENSIVE==1] - tmpdataUC$T_risk_control_cvd[tmpdataUC$INTENSIVE==1],
                                              tmpdataUC[tmpdataUC$INTENSIVE==1,vars]),
                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                     case.weights = tmpdataUC$cvd_ipcw[tmpdataUC$INTENSIVE==1], importance = "impurity", seed=1)
      pe1 <- fit1$prediction.error
      pet1 <- rbind(pet1, c(num.trees, mtry, min.node.size, pe1))

      fit0 <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$T_risk_treated_cvd[tmpdataUC$INTENSIVE==0] - tmpdataUC$cvd[tmpdataUC$INTENSIVE==0],
                                              tmpdataUC[tmpdataUC$INTENSIVE==0,vars]),
                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                     case.weights = tmpdataUC$cvd_ipcw[tmpdataUC$INTENSIVE==0], importance = "impurity", seed=1)
      pe0 <- fit0$prediction.error
      pet0 <- rbind(pet0, c(num.trees, mtry, min.node.size, pe0))
    }
  }
}
opt_params1 <- as.list(pet1[which.min(pet1[,4]),]);print(paste0("cvd1: ",opt_params1))
fit1_cvd <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$cvd[tmpdataUC$INTENSIVE==1] - tmpdataUC$T_risk_control_cvd[tmpdataUC$INTENSIVE==1],
                                        tmpdataUC[tmpdataUC$INTENSIVE==1,vars]),
               keep.inbag = TRUE, num.trees = opt_params1[[1]], mtry = opt_params1[[2]], min.node.size = opt_params1[[3]],
               case.weights = tmpdataUC$cvd_ipcw[tmpdataUC$INTENSIVE==1], importance = "impurity", seed=1)
tau1 <- predict(fit1_cvd, testdat[,vars])$predictions
cvd1_imp <- fit1_cvd$variable.importance

opt_params0 <- as.list(pet0[which.min(pet0[,4]),]);print(paste0("cvd0: ",opt_params0))
fit0_cvd <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$T_risk_treated_cvd[tmpdataUC$INTENSIVE==0] - tmpdataUC$cvd[tmpdataUC$INTENSIVE==0],
                                        tmpdataUC[tmpdataUC$INTENSIVE==0,vars]),
               keep.inbag = TRUE, num.trees = opt_params0[[1]], mtry = opt_params0[[2]], min.node.size = opt_params0[[3]],
               case.weights = tmpdataUC$cvd_ipcw[tmpdataUC$INTENSIVE==0], importance = "impurity", seed=1)
tau0 <- predict(fit0_cvd, testdat[,vars])$predictions
cvd0_imp <- fit0_cvd$variable.importance

XL_RD_cvd <- (tau1+tau0)/2


# SAE
nodesizes <- c(1,150)
pet1 <- pet0 <- NULL
for (j in 1:length(ntrees)){
  num.trees <- ntrees[j]
  for (k in 1:length(mtrys)){
    mtry <- mtrys[k]
    for (l in 1:length(nodesizes)){
      min.node.size <- nodesizes[l]
      fit1 <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$sae[tmpdataUC$INTENSIVE==1] - tmpdataUC$T_risk_control_sae[tmpdataUC$INTENSIVE==1],
                                              tmpdataUC[tmpdataUC$INTENSIVE==1,vars]),
                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                     case.weights = tmpdataUC$sae_ipcw[tmpdataUC$INTENSIVE==1], importance = "impurity", seed=1)
      pe1 <- fit1$prediction.error
      pet1 <- rbind(pet1, c(num.trees, mtry, min.node.size, pe1))

      fit0 <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$T_risk_treated_sae[tmpdataUC$INTENSIVE==0] - tmpdataUC$sae[tmpdataUC$INTENSIVE==0],
                                              tmpdataUC[tmpdataUC$INTENSIVE==0,vars]),
                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                     case.weights = tmpdataUC$sae_ipcw[tmpdataUC$INTENSIVE==0], importance = "impurity", seed=1)
      pe0 <- fit0$prediction.error
      pet0 <- rbind(pet0, c(num.trees, mtry, min.node.size, pe0))
    }
  }
}
opt_params1 <- as.list(pet1[which.min(pet1[,4]),]);print(paste0("sae1: ",opt_params1))
fit1_sae <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$sae[tmpdataUC$INTENSIVE==1] - tmpdataUC$T_risk_control_sae[tmpdataUC$INTENSIVE==1],
                                        tmpdataUC[tmpdataUC$INTENSIVE==1,vars]),
               keep.inbag = TRUE, num.trees = opt_params1[[1]], mtry = opt_params1[[2]], min.node.size = opt_params1[[3]],
               case.weights = tmpdataUC$sae_ipcw[tmpdataUC$INTENSIVE==1], importance = "impurity", seed=1)
tau1 <- predict(fit1_sae, testdat[,vars])$predictions
sae1_imp <- fit1_sae$variable.importance

opt_params0 <- as.list(pet0[which.min(pet0[,4]),]);print(paste0("sae0: ",opt_params0))
fit0_sae <- ranger(Y ~ ., data = data.frame(Y = tmpdataUC$T_risk_treated_sae[tmpdataUC$INTENSIVE==0] - tmpdataUC$sae[tmpdataUC$INTENSIVE==0],
                                        tmpdataUC[tmpdataUC$INTENSIVE==0,vars]),
               keep.inbag = TRUE, num.trees = opt_params0[[1]], mtry = opt_params0[[2]], min.node.size = opt_params0[[3]],
               case.weights = tmpdataUC$sae_ipcw[tmpdataUC$INTENSIVE==0], importance = "impurity", seed=1)
tau0 <- predict(fit0_sae, testdat[,vars])$predictions
sae0_imp <- fit0_sae$variable.importance

XL_RD_sae <- (tau1+tau0)/2

# Save all the models for building HTE calculator
saveRDS(fit1_cvd, "./fit1_cvd.rds")
saveRDS(fit0_cvd, "./fit0_cvd.rds")
saveRDS(fit1_sae, "./fit1_sae.rds")
saveRDS(fit0_sae, "./fit0_sae.rds")



# # save results
# setwd("../External Results")
# cvdpath <- "deepsurv_external_cvd"
# saepath <- "deepsurv_external_sae"
# sum_cvd <- read.csv(paste0(cvdpath,".csv")); head(sum_cvd)
# sum_sae <- read.csv(paste0(saepath,".csv")); head(sum_sae)
# sum_cvd <- data.frame(sum_cvd, XL_RD_cvd)
# sum_sae <- data.frame(sum_sae, XL_RD_sae)
# write.csv(sum_cvd, paste0(cvdpath, ".csv"), row.names = F)
# write.csv(sum_sae, paste0(saepath, ".csv"), row.names = F)


# logistic model for weights
# vars_all <- colnames(tmpdata)[2:39]
# nfold <- 10
# n <- round(dim(tmpdata)[1]/nfold,0)
# cvd_ipcw <- sae_ipcw <- NULL
# for (i in 1:nfold){
#   if(i==nfold){
#     traindat <- tmpdata[-((1+(i-1)*n):(dim(tmpdata)[1])),]
#     testdat  <- tmpdata[((1+(i-1)*n):(dim(tmpdata)[1])),]
#   }else{
#     traindat <- tmpdata[-((1+(i-1)*n):(i*n)),]
#     testdat  <- tmpdata[((1+(i-1)*n):(i*n)),]
#   }
#   
#   Wdata <- traindat[,colnames(traindat) %in% c(vars_all, "uncensored_cvd")]
#   cvdwfit <- glm(uncensored_cvd ~., family=binomial, data = Wdata)
#   est_r <- predict(cvdwfit, newdata = testdat[,vars_all], type="response")
#   cvd_ipcw <- c(cvd_ipcw, 1/est_r)
#   
#   Wdata <- traindat[,colnames(traindat) %in% c(vars_all, "uncensored_sae")]
#   saewfit <- glm(uncensored_sae ~., family=binomial, data = Wdata)
#   est_r <- predict(saewfit, newdata = testdat[,vars_all], type="response")
#   sae_ipcw <- c(sae_ipcw, 1/est_r)
# }
# 
# # Trim extreme weights
# tmpdata$cvd_ipcw <- cvd_ipcw
# tmpdata$sae_ipcw <- sae_ipcw

nfold <- 10
n <- round(dim(tmpdata)[1]/nfold,0)
XL_RD_sae <- XL_RD_cvd <- NULL
cvd1_imp <- cvd0_imp <- sae1_imp <- sae0_imp <- matrix(NA, length(vars), nfold)
for (i in 1:nfold){
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

  # Tune ranger for CATE cvd outcome
  # nodesizes <- c(1)
  #nodesizes <- c(1,35,100,200)
  # ntrees <- c(1000)
  # mtrys <- 6
  # pet1 <- pet0 <- NULL
  # for (j in 1:length(ntrees)){
  #   num.trees <- ntrees[j]
  #   for (k in 1:length(mtrys)){
  #     mtry <- mtrys[k]
  #     for (l in 1:length(nodesizes)){
  #       min.node.size <- nodesizes[l]
  #       fit1 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$cvd[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_cvd[traindatUC$INTENSIVE==1],
  #                                              traindatUC[traindatUC$INTENSIVE==1,vars]),
  #                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
  #                     case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
  #       pe1 <- fit1$prediction.error
  #       pet1 <- rbind(pet1, c(num.trees, mtry, min.node.size, pe1))
  # 
  #       fit0 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_cvd[traindatUC$INTENSIVE==0] - traindatUC$cvd[traindatUC$INTENSIVE==0],
  #                                              traindatUC[traindatUC$INTENSIVE==0,vars]),
  #                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
  #                     case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
  #       pe0 <- fit0$prediction.error
  #       pet0 <- rbind(pet0, c(num.trees, mtry, min.node.size, pe0))
  #     }
  #   }
  # }
  opt_params1 <- as.list(pet1[which.min(pet1[,4]),]);print(paste0("cvd1: ",opt_params1))
  fit1 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$cvd[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_cvd[traindatUC$INTENSIVE==1],
                                          traindatUC[traindatUC$INTENSIVE==1,vars]),
                 keep.inbag = TRUE, num.trees = 1000, mtry = 6, min.node.size = 1,
                 case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
  tau1 <- predict(fit1, testdat[,vars])$predictions
  cvd1_imp[,i] <- fit1$variable.importance

  opt_params0 <- as.list(pet0[which.min(pet0[,4]),]);print(paste0("cvd0: ",opt_params0))
  fit0 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_cvd[traindatUC$INTENSIVE==0] - traindatUC$cvd[traindatUC$INTENSIVE==0],
                                          traindatUC[traindatUC$INTENSIVE==0,vars]),
                 keep.inbag = TRUE, num.trees = 1000, mtry = 6, min.node.size = 1,
                 case.weights = traindatUC$cvd_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
  tau0 <- predict(fit0, testdat[,vars])$predictions
  cvd0_imp[,i] <- fit0$variable.importance

  XL_RD_cvd <- c(XL_RD_cvd, (tau1+tau0)/2)

  # # SAE
  # traindatUC <- traindat[traindat$uncensored_sae==1, ]
  # traindatUC$sae <-  ifelse(traindatUC$t_saes>t50, 0, traindatUC$sae)
  #
  # nodesizes <- c(1,150)
  # pet1 <- pet0 <- NULL
  # for (j in 1:length(ntrees)){
  #   num.trees <- ntrees[j]
  #   for (k in 1:length(mtrys)){
  #     mtry <- mtrys[k]
  #     for (l in 1:length(nodesizes)){
  #       min.node.size <- nodesizes[l]
  #       fit1 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$sae[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_sae[traindatUC$INTENSIVE==1],
  #                                              traindatUC[traindatUC$INTENSIVE==1,vars]),
  #                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
  #                     case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
  #       pe1 <- fit1$prediction.error
  #       pet1 <- rbind(pet1, c(num.trees, mtry, min.node.size, pe1))
  #
  #       fit0 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_sae[traindatUC$INTENSIVE==0] - traindatUC$sae[traindatUC$INTENSIVE==0],
  #                                              traindatUC[traindatUC$INTENSIVE==0,vars]),
  #                     keep.inbag = TRUE, num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
  #                     case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
  #       pe0 <- fit0$prediction.error
  #       pet0 <- rbind(pet0, c(num.trees, mtry, min.node.size, pe0))
  #     }
  #   }
  # }
  # opt_params1 <- as.list(pet1[which.min(pet1[,4]),]);print(paste0("sae1: ",opt_params1))
  # fit1 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$sae[traindatUC$INTENSIVE==1] - traindatUC$T_risk_control_sae[traindatUC$INTENSIVE==1],
  #                                         traindatUC[traindatUC$INTENSIVE==1,vars]),
  #                keep.inbag = TRUE, num.trees = opt_params1[[1]], mtry = 6, min.node.size = opt_params1[[3]],
  #                case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==1], importance = "impurity", seed=1)
  # tau1 <- predict(fit1, testdat[,vars])$predictions
  # sae1_imp[,i] <- fit1$variable.importance
  #
  # opt_params0 <- as.list(pet0[which.min(pet0[,4]),]);print(paste0("sae0: ",opt_params0))
  # fit0 <- ranger(Y ~ ., data = data.frame(Y = traindatUC$T_risk_treated_sae[traindatUC$INTENSIVE==0] - traindatUC$sae[traindatUC$INTENSIVE==0],
  #                                         traindatUC[traindatUC$INTENSIVE==0,vars]),
  #                keep.inbag = TRUE, num.trees = opt_params0[[1]], mtry = 6, min.node.size = opt_params0[[3]],
  #                case.weights = traindatUC$sae_ipcw[traindatUC$INTENSIVE==0], importance = "impurity", seed=1)
  # tau0 <- predict(fit0, testdat[,vars])$predictions
  # sae0_imp[,i] <- fit0$variable.importance
  #
  # XL_RD_sae <- c(XL_RD_sae, (tau1+tau0)/2)
}
# cvd_imp <- (rowMeans(cvd1_imp) + rowMeans(cvd0_imp))/2
# cvd_imp_dat <- data.frame(vars, cvd_imp)
# cvd_imp_dat <- cvd_imp_dat[order(-cvd_imp_dat$cvd_imp),]
# 
# sae_imp <- (rowMeans(sae1_imp) + rowMeans(sae0_imp))/2
# sae_imp_dat <- data.frame(vars, sae_imp)
# sae_imp_dat <- sae_imp_dat[order(-sae_imp_dat$sae_imp),]
# 
# sum_imp <- data.frame(cvd_imp_dat, sae_imp_dat)
# colnames(sum_imp) <- c("cvd_vars", "cvd_imp","sae_vars", "sae_imp")
# write.csv(sum_imp,"combined_sum_imp_bart.csv", row.names = F)
# cvdpath <- "bart_cvd_tuned_results/combined_risk_T1bart_cvd"
# saepath <- "bart_sae_tuned_results/combined_risk_T0bart_sae"
# # Others load estimates
# sum_cvd <- NULL
# for (i in 1:10){
#   tmpres <- read.csv(paste0(cvdpath,i,".csv"))
#   sum_cvd <- rbind(sum_cvd, tmpres)
# }
# sum_cvd_S <- sum_cvd
# sum_cvd_T1 <- sum_cvd
# sum_cvd_T0 <- sum_cvd
# sum_cvd <- data.frame(sum_cvd_S, sum_cvd_T0, sum_cvd_T1)
# write.csv(sum_cvd, "combined_risk_bart_cvd.csv", row.names = F)
# 
# sum_sae <- NULL
# for (i in 1:10){
#   tmpres <- read.csv(paste0(saepath,i,".csv"))
#   sum_sae <- rbind(sum_sae, tmpres)
# }
# sum_sae_S <- sum_sae
# sum_sae_T1 <- sum_sae
# sum_sae_T0 <- sum_sae
# sum_sae <- data.frame(sum_sae_S, sum_sae_T0, sum_sae_T1)
# write.csv(sum_sae, "combined_risk_bart_sae.csv", row.names = F)
