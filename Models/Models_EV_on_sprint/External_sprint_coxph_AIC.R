
##---------------------- External validation for Cox PH method using ACCORD -------------------##
##----------------------       Crystal Xu               06/3/2021     -------------------##
library(survival); library(MASS)
#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
t50 <- floor(365.25*3.26)

## ACCORD data as training 
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_train <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_train <- accord_train[complete.cases(accord_train),]; dim(accord_train)
accord_train$t_cvds <- ceiling(accord_train$t_cvds)
accord_train$t_saes <- ceiling(accord_train$t_saes)
Xvar <- accord_train[,2:dim(accord_blvars)[2]]
Xvar_int <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar)

## SPRINT data as testing 
sprint_blvars <- read.csv("blvars.csv")
sprint_outcome <- read.csv("outcome.csv")
sprint_test <- merge(sprint_blvars, sprint_outcome, by="MASKID",all=T)
sprint_test <- sprint_test[complete.cases(sprint_test),]

# Set counterfacutal treatment levels 
testXvar <- sprint_test[,2:dim(sprint_blvars)[2]]
testXvar1 <- testXvar; testXvar1$INTENSIVE <- 1
testXvar0 <- testXvar; testXvar0$INTENSIVE <- 0

# Derive treatment*covariates terms 
testXvar_int <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, testXvar)
testXvar_int1 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, testXvar1)
testXvar_int0 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, testXvar0)

#-------------------- Benefit model ----------------------------#
accord_cvd <- data.frame(accord_train$t_cvds, accord_train$cvd, Xvar_int)
names(accord_cvd)[1:2] <- c("t_cvds","cvd")
sprint_cvd <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar_int)
names(sprint_cvd)[1:2] <- c("t_cvds","cvd")
sprint_cvd1 <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar_int1)
names(sprint_cvd1)[1:2] <- c("t_cvds","cvd")
sprint_cvd0 <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar_int0)
names(sprint_cvd0)[1:2] <- c("t_cvds","cvd")

# S-Learner 
cvd_coxph_init <- coxph(formula = Surv(t_cvds, cvd) ~ .,data = accord_cvd)

# Find the final training model with covariate selection 
set.seed(199)
cvd_coxph <- stepAIC(cvd_coxph_init,
                     scope = list(upper = ~.,
                                 lower = ~ INTENSIVE))   # standardize all covariates by default

# Final list of predictors
vars <- colnames(model.matrix(cvd_coxph))  
formula <- as.formula(paste0("Surv(t_cvds, cvd) ~ ", paste(vars, sep=" ", collapse = "+")))

# SPRINT S0
cvd_coxph_sprint <- coxph(formula, data = sprint_cvd)
cumhaz <- basehaz(cvd_coxph_sprint)
S0_sprint <- exp(-cumhaz[,1])
S0dat <- data.frame(S0_sprint, cumhaz$time); colnames(S0dat) <- c("S0", "time")
S0_accord_t50 <- S0dat[S0dat$time==t50,]$S0

# Predictions 
pred_S1_cvd <- predict(cvd_coxph, newdata=sprint_cvd1[,vars], type="risk")
pred_S0_cvd <- predict(cvd_coxph, newdata=sprint_cvd0[,vars], type="risk")
Scoxph_treated_cvd <- 1 - S0_accord_t50^pred_S1_cvd
Scoxph_control_cvd <- 1 - S0_accord_t50^pred_S0_cvd

# save results
risk_coxph_AIC_cvd <- data.frame(Scoxph_control_cvd, Scoxph_treated_cvd)
write.csv(risk_coxph_AIC_cvd, "../Analyses Stanford Team/Analysis Results/SPRINT_external_results/coxph_AIC_cvd_sprint_external.csv", row.names = F)

# S_predinc_cvd <- risk_coxph_AIC_cvd[,2]*accord_cvd$INTENSIVE + risk_coxph_AIC_cvd[,1]*(1-accord_cvd$INTENSIVE)
# Spredgrp <- as.numeric(cut2(S_predinc_cvd,g=10))
# GND_Scvd <- GND.calib(pred=S_predinc_cvd, tvar=accord_cvd$t_cvds, out=accord_cvd$cvd,
#                       cens.t=t50, groups=Spredgrp, adm.cens=t50, name="")


# #--------------------------------- Risk model --------------------------------#
# accord_sae <- data.frame(accord_train$t_saes, accord_train$sae, Xvar_int)
# names(accord_sae)[1:2] <- c("t_saes","sae")
# sprint_sae <- data.frame(sprint_test$t_saes, sprint_test$sae, testXvar_int)
# names(sprint_sae)[1:2] <- c("t_saes","sae")
# sprint_sae1 <- data.frame(sprint_test$t_saes, sprint_test$sae, testXvar_int1)
# names(sprint_sae1)[1:2] <- c("t_saes","sae")
# sprint_sae0 <- data.frame(sprint_test$t_saes, sprint_test$sae, testXvar_int0)
# names(sprint_sae0)[1:2] <- c("t_saes","sae")
# 
# # S-Learner 
# sae_coxph_init <- coxph(formula = Surv(t_saes, sae) ~ .,data = accord_sae)
# 
# # Find the final training model with covariate selection 
# set.seed(199)
# sae_coxph <- stepAIC(sae_coxph_init,
#                      scope = list(upper = ~.,
#                                   lower = ~ INTENSIVE))   # standardize all covariates by default
# 
# # Final list of predictors
# vars <- colnames(model.matrix(sae_coxph))  
# formula <- as.formula(paste0("Surv(t_saes, sae) ~ ", paste(vars, sep=" ", collapse = "+")))
# 
# # SPRINT S0
# sae_coxph_sprint <- coxph(formula, data = sprint_sae)
# cumhaz <- basehaz(sae_coxph_sprint)
# S0_sprint <- exp(-cumhaz[,1])
# S0dat <- data.frame(S0_sprint, cumhaz$time); colnames(S0dat) <- c("S0", "time")
# S0_accord_t50 <- S0dat[S0dat$time==t50,]$S0
# 
# # Predictions 
# pred_S1_sae <- predict(sae_coxph, newdata=sprint_sae1[,vars], type="risk")
# pred_S0_sae <- predict(sae_coxph, newdata=sprint_sae0[,vars], type="risk")
# Scoxph_treated_sae <- 1 - S0_accord_t50^pred_S1_sae
# Scoxph_control_sae <- 1 - S0_accord_t50^pred_S0_sae
# 
# # save results
# risk_coxph_AIC_sae <- data.frame(Scoxph_control_sae, Scoxph_treated_sae)
# write.csv(risk_coxph_AIC_sae, "../Analyses Stanford Team/Analysis Results/SPRINT_external_results/coxph_AIC_sae_sprint_external.csv", row.names = F)

# S_predinc_sae <- risk_coxph_AIC_sae[,2]*accord_sae$INTENSIVE + risk_coxph_AIC_sae[,1]*(1-accord_sae$INTENSIVE)
# Spredgrp <- as.numeric(cut2(S_predinc_sae,g=10))
# GND_Ssae <- GND.calib(pred=S_predinc_sae, tvar=accord_sae$t_saes, out=accord_sae$sae,
#                       cens.t=t50, groups=Spredgrp, adm.cens=t50, name="")
# 



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
