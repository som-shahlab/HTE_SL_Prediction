
##---------------------- Apply Cox PH method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021     -------------------##
library(survival); library(glmnet)
# compute the baseline survival at all unique event times using a penalized cox model
base_surv <- function(fit, Y, D, x){
  data <- data.frame(t_event=Y, event=D, x)
  tab <- data.frame(table(data[data$event == 1, "t_event"])) 
  y <- as.numeric(as.character(sort(unique(tab[,1]))))
  d <- tab[,2]  # number of events at each unique time                               
  
  betaHat <- as.vector((fit$glmnet.fit$beta)[,fit$lambda==fit$lambda.1se])
  h0 <- rep(NA, length(y))
  for(l in 1:length(y)){
    h0[l] <- d[l] / sum(exp(x[data$t_event >= y[l], rownames(fit$glmnet.fit$beta)] %*% betaHat))
  }    
  
  S0 <- exp(-cumsum(h0))
  outcome <- data.frame(time=y,survival=S0)
  outcome
}
# predict survival at time t
pred_surv <- function(fit, S0, x, times){
  link <- predict(fit$glmnet.fit,x,type = "link")[,fit$lambda==fit$lambda.1se] # link = x*beta
  colnames(link) <- NULL
  S0_t <- S0$survival[S0$time>=times][1]
  surv <- S0_t^exp(link)
  surv
}

t50 <- floor(365.25*3.26)

setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
sprint_blvars <- read.csv("blvars.csv")
sprint_outcome <- read.csv("outcome.csv")

# Combine baseline variables and outcomes into one data frame 
sprint_train <- merge(sprint_blvars, sprint_outcome, by="MASKID",all=T)
sprint_train <- sprint_train[complete.cases(sprint_train),]; dim(sprint_train) # 8969 43

# Recode treatment variable according to Tian (2015)
# "A simple method for estimating interactions between a treatment and a large number of covariates"
sprint_train$INTENSIVE[sprint_train$INTENSIVE==0] <- -1  
sprint_train$INTENSIVE <- sprint_train$INTENSIVE/2  
Xvar <- sprint_train[,2:dim(sprint_blvars)[2]]
Xvar_int <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar)
Xvar_int[,2:dim(Xvar_int)[2]] <- scale(Xvar_int[,2:dim(Xvar_int)[2]]) # standardize the covariates

## ACCORD data as testing 
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_outcome$t_cvds <- ceiling(accord_outcome$t_cvds)
accord_outcome$t_saes <- ceiling(accord_outcome$t_saes)
accord_test <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_test <- accord_test[complete.cases(accord_test),]; dim(accord_test) # 4448 43

# Set counterfacutal treatment levels 
testXvar <- accord_test[,2:dim(accord_blvars)[2]]
testXvar1 <- testXvar; testXvar1$INTENSIVE <- 1/2
testXvar0 <- testXvar; testXvar0$INTENSIVE <- -1/2

# Derive treatment*covariates terms 
testXvar_int <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, testXvar)
testXvar_int1 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, testXvar1)
testXvar_int0 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, testXvar0)
testXvar_int[,2:dim(testXvar_int)[2]] <- scale(testXvar_int[,2:dim(testXvar_int)[2]])
testXvar_int1[,2:dim(testXvar_int1)[2]] <- scale(testXvar_int1[,2:dim(testXvar_int1)[2]]) # Standardize the covariates
testXvar_int0[,2:dim(testXvar_int0)[2]] <- scale(testXvar_int0[,2:dim(testXvar_int0)[2]])


#-------------------- Benefit model ----------------------------#
sprint_cvd <- data.frame(sprint_train$t_cvds, sprint_train$cvd, Xvar_int)
names(sprint_cvd)[1:2] <- c("t_cvds","cvd")
accord_cvd <- data.frame(accord_test$t_cvds, accord_test$cvd, testXvar_int)
names(accord_cvd)[1:2] <- c("t_cvds","cvd")
accord_cvd1 <- data.frame(accord_test$t_cvds, accord_test$cvd, testXvar_int1)
names(accord_cvd1)[1:2] <- c("t_cvds","cvd")
accord_cvd0 <- data.frame(accord_test$t_cvds, accord_test$cvd, testXvar_int0)
names(accord_cvd0)[1:2] <- c("t_cvds","cvd")

# S-Learner
# Tuning parameter: shrinkage level 
fpf <- c(0,rep(1, (dim(sprint_cvd)[2]-3)))
fpf[regexpr("INTENSIVE.", colnames(sprint_cvd[,3:dim(sprint_cvd)[2]]))!=-1] <- 1

# Training model using SPRINT
nFolds <- 10 
foldid <- rep(seq(nFolds), length.out = nrow(sprint_cvd))
cvd_fit <- cv.glmnet(x=as.matrix(Xvar_int),
                     y=Surv(sprint_cvd$t_cvds,sprint_cvd$cvd),
                     family="cox",
                     alpha=0.75,
                     penalty.factor=fpf,
                     foldid = foldid, 
                     standardize = F)

# Model on accord data to get the S0 in accord
foldid <- rep(seq(nFolds), length.out = nrow(accord_cvd))
cvd_fit_accord <- cv.glmnet(x=as.matrix(testXvar_int),
                            y=Surv(accord_cvd$t_cvds, accord_cvd$cvd),
                            family="cox",
                            alpha=0.75,
                            penalty.factor=fpf,
                            foldid = foldid, 
                            standardize = F)
S0_accord <- base_surv(fit = cvd_fit_accord, 
                         Y = accord_cvd$t_cvds,
                         D = accord_cvd$cvd,
                         x = as.matrix(testXvar_int))

pred_S1_cvd <- pred_surv(fit = cvd_fit, S0 = S0_accord, x = as.matrix(testXvar_int1), times = t50)
pred_S0_cvd <- pred_surv(fit = cvd_fit, S0 = S0_accord, x = as.matrix(testXvar_int0), times = t50)
Scoxph_treated_cvd <- 1 - pred_S1_cvd
Scoxph_control_cvd <- 1 - pred_S0_cvd

# save results 
risk_coxph_EN_cvd <- data.frame(Scoxph_treated_cvd, Scoxph_control_cvd)
write.csv(risk_coxph_EN_cvd, "../Analyses Stanford Team/Analysis Results/External Results/coxphEN_cvd_external.csv", row.names = F)
# S_predinc_cvd <- risk_coxph_EN_cvd[,1]*accord_cvd$INTENSIVE + risk_coxph_EN_cvd[,2]*(1-accord_cvd$INTENSIVE)
# Spredgrp <- as.numeric(cut2(S_predinc_cvd,g=10))
# GND_Scvd <- GND.calib(pred=S_predinc_cvd, tvar=accord_cvd$t_cvds, out=accord_cvd$cvd,
#                       cens.t=t50, groups=Spredgrp, adm.cens=t50, name = "")




#------------------------------ Risk model ---------------------------------#
sprint_sae <- data.frame(sprint_train$t_saes, sprint_train$sae, Xvar_int)
names(sprint_sae)[1:2] <- c("t_saes","sae")
accord_sae <- data.frame(accord_test$t_saes, accord_test$sae, testXvar_int)
names(accord_sae)[1:2] <- c("t_saes","sae")
accord_sae1 <- data.frame(accord_test$t_saes, accord_test$sae, testXvar_int1)
names(accord_sae1)[1:2] <- c("t_saes","sae")
accord_sae0 <- data.frame(accord_test$t_saes, accord_test$sae, testXvar_int0)
names(accord_sae0)[1:2] <- c("t_saes","sae")

# S-Learner
# Tuning parameter: shrinkage level 
fpf <- c(0,rep(1, (dim(sprint_sae)[2]-3)))
fpf[regexpr("INTENSIVE.", colnames(sprint_sae[,3:dim(sprint_sae)[2]]))!=-1] <- 1

# Training model using SPRINT
nFolds <- 10
foldid <- rep(seq(nFolds), length.out = nrow(sprint_sae))
sae_fit <- cv.glmnet(x=as.matrix(Xvar_int),
                     y=Surv(sprint_sae$t_saes,sprint_sae$sae),
                     family="cox",
                     alpha=0.75,
                     penalty.factor=fpf,
                     foldid = foldid, 
                     standardize = F)

# Model on accord data to get the S0 in accord
foldid <- rep(seq(nFolds), length.out = nrow(accord_sae))
sae_fit_accord <- cv.glmnet(x=as.matrix(testXvar_int),
                            y=Surv(accord_sae$t_saes, accord_sae$sae),
                            family="cox",
                            alpha=0.75,
                            penalty.factor=fpf,
                            foldid = foldid, 
                            standardize = F)

S0_accord <- base_surv(fit = sae_fit_accord, 
                       Y = accord_sae$t_saes,
                       D = accord_sae$sae,
                       x = as.matrix(testXvar_int))

pred_S1_sae <- pred_surv(fit = sae_fit, S0 = S0_accord, x = as.matrix(testXvar_int1), times = t50)
pred_S0_sae <- pred_surv(fit = sae_fit, S0 = S0_accord, x = as.matrix(testXvar_int0), times = t50)
Scoxph_treated_sae <- 1 - pred_S1_sae
Scoxph_control_sae <- 1 - pred_S0_sae

# save results 
risk_coxph_EN_sae <- data.frame(Scoxph_treated_sae, Scoxph_control_sae)
write.csv(risk_coxph_EN_sae, "../Analyses Stanford Team/Analysis Results/External Results/coxphEN_sae_external.csv", row.names = F)
# S_predinc_sae <- risk_coxph_EN_sae[,1]*accord_sae$INTENSIVE + risk_coxph_EN_sae[,2]*(1-accord_sae$INTENSIVE)
# Spredgrp <- as.numeric(cut2(S_predinc_sae,g=10))
# GND_Ssae <- GND.calib(pred=S_predinc_sae, tvar=accord_sae$t_saes, out=accord_sae$sae, 
#                       cens.t=t50, groups=Spredgrp, adm.cens=t50, name="")


##----------------------- Define treatment policy  -----------------------------------##
##-- A subject is treated if the treatment reduces 20-25% of her/his baseline risk ---##
##-- g=1 if delta_R >= R0*0.2; g=0 o/w                                             ---##
##------------------------------------------------------------------------------------##
# g_est = ifelse(-coxph_p_RD_cvd >= risk_control_cvd*0.25,1,0)
# table(g_est)
# 
# # Compare to estimated treatment rule based on PCE
# g_ASCVD = ifelse(risk_control_cvd > 0.075,1,0)
# table(g_est, g_ASCVD)
# 
# # define the utility function using mean survival time 
# # compute censoring weight
# censor <- 1-tmpdata_cvd$cvd
# KMfit <- survfit(Surv(t_cvds, censor)~1, data = tmpdata_cvd)
# 
# censor_weight <- rep(NA, dim(tmpdata_cvd)[1])
# for (i in 1:dim(tmpdata_cvd)[1]){
#   censor_weight[i] <- summary(KMfit,time=tmpdata_cvd$t_cvds[i])$surv
# }
# 
# tmpdata_cvd$censor_weight <- censor_weight
# tmpdata_cvd$g_est <- g_est
# tmpdata_cvd_sub <- tmpdata_cvd[tmpdata_cvd$INTENSIVE==tmpdata_cvd$g_est,]
# 
# 
# # compute utility 
# utility <- tmpdata_cvd_sub$t_cvds*tmpdata_cvd_sub$cvd/tmpdata_cvd_sub$censor_weight
# mean(utility, na.rm=T)  # 49.71