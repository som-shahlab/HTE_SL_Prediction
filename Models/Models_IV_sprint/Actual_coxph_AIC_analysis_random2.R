
##---------------------- Apply Cox PH method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021     -------------------##

library(survival); library(MASS)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969 43
t50 <- floor(365.25*3.26)

# set counterfacutal treatment levels 
Xvar <- tmpdata[,2:dim(blvars)[2]]
Xvar1 <- Xvar
Xvar1$INTENSIVE <- 1
Xvar0 <- Xvar
Xvar0$INTENSIVE <- 0

# Derive treatment*covariates terms 
Xvar_int  <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar)
Xvar_int1 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar1) 
Xvar_int0 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar0)

# Random splitting seed
if(length(args <- commandArgs(T))>0){
  caseid <- as.integer(args[[1]])
  message("running for parameter set ", caseid)
}
seed_num <- caseid+1


# #-------------------- Benefit model ----------------------------#
# tmpdata_cvd <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar_int)
# names(tmpdata_cvd)[1:2] <- c("t_cvds","cvd")
# tmpdata_cvd1 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar_int1)
# names(tmpdata_cvd1)[1:2] <- c("t_cvds","cvd")
# tmpdata_cvd0 <- data.frame(tmpdata$t_cvds, tmpdata$cvd, Xvar_int0)
# names(tmpdata_cvd0)[1:2] <- c("t_cvds","cvd")
# 
# # 10-fold CV to get the out-of-bag estimates
# # Randomly shuffle the data
# set.seed(seed_num)
# index <- sample(nrow(tmpdata_cvd))
# tmpdata_cvd <- tmpdata_cvd[index,]
# tmpdata_cvd1 <- tmpdata_cvd1[index,]
# tmpdata_cvd0 <- tmpdata_cvd0[index,]
# 
# folds <- cut(seq(1,nrow(tmpdata_cvd)),breaks=10,labels=FALSE)
# Scoxph_treated_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# Scoxph_control_cvd <- rep(NA,dim(tmpdata_cvd)[1])
# for (i in 1:10){
#   
#   testIndexes <- which(folds==i,arr.ind=TRUE)
#   traindat <- tmpdata_cvd[-testIndexes, ]
#   testdat <- tmpdata_cvd[testIndexes, ]
#   testdat1 <- tmpdata_cvd1[testIndexes, ]
#   testdat0 <- tmpdata_cvd0[testIndexes, ]
#   
#   # S-Learner -- ALL the data is used in the covariate selection step
#   set.seed(99)
#   cvd_coxph_init <- coxph(formula = Surv(t_cvds, cvd) ~ .,data = traindat)
#   CS <- stepAIC(cvd_coxph_init,
#                 scope = list(upper = ~.,
#                              lower = ~ INTENSIVE))   # standardize all covariates by default
#   vars <- colnames(model.matrix(CS))
#   formula <- as.formula(paste0("Surv(t_cvds, cvd) ~ ", paste(vars, sep=" ", collapse = "+")))
#   
#   # Fit the final model
#   cvd_coxph <- coxph(formula, data = traindat)
#   
#   # baseline cumulative hazard at t50
#   bh_dat <- basehaz(cvd_coxph)
#   bh <- bh_dat[bh_dat$time==t50,]$hazard
#   
#   # risk score under each treatment arm on exp(lp)
#   est_r1 <- predict(cvd_coxph, newdata = testdat1, type="risk")
#   est_r0 <- predict(cvd_coxph, newdata = testdat0, type="risk")
#   
#   # survival probabilities
#   est_S1_cvd <- exp(-bh)^est_r1
#   est_S0_cvd <- exp(-bh)^est_r0
#   
#   # risk estimates
#   Scoxph_treated_cvd[testIndexes] <- 1 - est_S1_cvd
#   Scoxph_control_cvd[testIndexes] <- 1 - est_S0_cvd
# }
# # save results
# risk_coxph_AIC_cvd <- data.frame(index, Scoxph_treated_cvd,Scoxph_control_cvd)
# write.csv(risk_coxph_AIC_cvd, paste0("./cvd_AIC_results/random_split_coxph_AIC_cvd",seed_num,".csv"), row.names = F)




#-------------------- Risk model ----------------------------#
tmpdata_sae <- data.frame(tmpdata$t_saes, tmpdata$sae, Xvar_int)
names(tmpdata_sae)[1:2] <- c("t_saes","sae")
tmpdata_sae1 <- data.frame(tmpdata$t_saes, tmpdata$sae, Xvar_int1)
names(tmpdata_sae1)[1:2] <- c("t_saes","sae")
tmpdata_sae0 <- data.frame(tmpdata$t_saes, tmpdata$sae, Xvar_int0)
names(tmpdata_sae0)[1:2] <- c("t_saes","sae")

# 10-fold CV to get the out-of-bag estimates
# Randomly shuffle the data
set.seed(seed_num)
index <- sample(nrow(tmpdata_sae))
tmpdata_sae <- tmpdata_sae[index,]
tmpdata_sae1 <- tmpdata_sae1[index,]
tmpdata_sae0 <- tmpdata_sae0[index,]

folds <- cut(seq(1,nrow(tmpdata_sae)),breaks=10,labels=FALSE)
Scoxph_treated_sae <- rep(NA,dim(tmpdata_sae)[1])
Scoxph_control_sae <- rep(NA,dim(tmpdata_sae)[1])
for (i in 1:10){
  
  testIndexes <- which(folds==i,arr.ind=TRUE)
  traindat <- tmpdata_sae[-testIndexes, ]
  testdat <- tmpdata_sae[testIndexes, ]
  testdat1 <- tmpdata_sae1[testIndexes, ]
  testdat0 <- tmpdata_sae0[testIndexes, ]
  
  # S-Learner -- ALL the data is used in the covariate selection step
  set.seed(99)
  sae_coxph_init <- coxph(formula = Surv(t_saes, sae) ~ .,data = traindat)
  CS <- stepAIC(sae_coxph_init,
                scope = list(upper = ~.,
                             lower = ~ INTENSIVE))   # standardize all covariates by default
  vars <- colnames(model.matrix(CS))
  formula <- as.formula(paste0("Surv(t_saes, sae) ~ ", paste(vars, sep=" ", collapse = "+")))
  
  # Fit the final model
  sae_coxph <- coxph(formula, data = traindat)
  
  # baseline cumulative hazard at t50
  bh_dat <- basehaz(sae_coxph)
  bh <- bh_dat[bh_dat$time==t50,]$hazard
  
  # risk score under each treatment arm on exp(lp)
  est_r1 <- predict(sae_coxph, newdata = testdat1, type="risk")
  est_r0 <- predict(sae_coxph, newdata = testdat0, type="risk")
  
  # survival probabilities
  est_S1_sae <- exp(-bh)^est_r1
  est_S0_sae <- exp(-bh)^est_r0
  
  # risk estimates
  Scoxph_treated_sae[testIndexes] <- 1 - est_S1_sae
  Scoxph_control_sae[testIndexes] <- 1 - est_S0_sae
}
# save results
risk_coxph_AIC_sae <- data.frame(index, Scoxph_treated_sae,Scoxph_control_sae)
write.csv(risk_coxph_AIC_sae, paste0("./sae_AIC_results/random_split_coxph_AIC_sae",seed_num,".csv"), row.names = F)



