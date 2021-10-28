
##---------------------- Apply Cox PH method to the actual SPRINT data -------------------##
##----------------------       Crystal Xu               02/23/2021     -------------------##

library(survival); library(glmnet)
#source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
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

#setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 

# Combine baseline variables and outcomes into one data frame 
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969 43

# Recode treatment variable according to Tian (2015)
# "A simple method for estimating interactions between a treatment and a large number of covariates"
tmpdata$INTENSIVE[tmpdata$INTENSIVE==0] <- -1
tmpdata$INTENSIVE <- tmpdata$INTENSIVE/2

# Set counterfacutal treatment levels & create treatment*covariates terms
Xvar <- tmpdata[,2:dim(blvars)[2]]
Xvar1 <- Xvar
Xvar1$INTENSIVE <- 0.5
Xvar0 <- Xvar
Xvar0$INTENSIVE <- -0.5

Xvar_int <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar)
Xvar_int1 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar1)
Xvar_int0 <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar0)

# Standardize the covariates
Xvar_int[,2:dim(Xvar_int)[2]] <- scale(Xvar_int[,2:dim(Xvar_int)[2]])
Xvar_int1[,2:dim(Xvar_int)[2]] <- scale(Xvar_int1[,2:dim(Xvar_int)[2]])
Xvar_int0[,2:dim(Xvar_int)[2]] <- scale(Xvar_int0[,2:dim(Xvar_int)[2]])

t50 <- floor(365.25*3.26)

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
#   testdat  <- tmpdata_cvd[testIndexes, ]
#   testdat1 <- tmpdata_cvd1[testIndexes, ]
#   testdat0 <- tmpdata_cvd0[testIndexes, ]
#   
#   # covariate selection using elastic net 
#   fpf <- c(0,rep(1, (dim(traindat)[2]-3)))
#   fpf[regexpr("INTENSIVE.", colnames(traindat[,3:dim(traindat)[2]]))!=-1] <- 1
#   
#   nFolds <- 10
#   foldid <- rep(seq(nFolds), length.out = nrow(traindat))
#   cvd_fit <- cv.glmnet(x=as.matrix(traindat[,3:dim(traindat)[2]]),
#                        y=Surv(traindat$t_cvds,traindat$cvd),
#                        family="cox",
#                        alpha=0.75,
#                        penalty.factor=fpf,
#                        foldid = foldid, 
#                        standardize = F)
# 
#   # baseline survival & counterfactuals
#   S0 <- base_surv(fit = cvd_fit, 
#                     Y = traindat$t_cvds,
#                     D = traindat$cvd,
#                     x = as.matrix(traindat[,3:dim(traindat)[2]]))
#   
#   est_S1_cvd <- pred_surv(fit = cvd_fit,
#                            S0 = S0,
#                             x = as.matrix(testdat1[,3:dim(testdat1)[2]]),
#                         times = t50)
#   
#   est_S0_cvd <- pred_surv(fit = cvd_fit,
#                            S0 = S0,
#                             x = as.matrix(testdat0[,3:dim(testdat0)[2]]),
#                         times = t50)
# 
#   Scoxph_treated_cvd[testIndexes] <- 1 - est_S1_cvd
#   Scoxph_control_cvd[testIndexes] <- 1 - est_S0_cvd
# }
# 
# # save results 
# risk_coxph_EN_cvd <- data.frame(index, Scoxph_treated_cvd, Scoxph_control_cvd)
# write.csv(risk_coxph_EN_cvd, paste0("./cvd_EN_results/random_split_coxph_EN_cvd",seed_num, ".csv"), row.names = F)



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
  testdat  <- tmpdata_sae[testIndexes, ]
  testdat1 <- tmpdata_sae1[testIndexes, ]
  testdat0 <- tmpdata_sae0[testIndexes, ]
  
  # covariate selection using elastic net 
  fpf <- c(0,rep(1, (dim(traindat)[2]-3)))
  fpf[regexpr("INTENSIVE.", colnames(traindat[,3:dim(traindat)[2]]))!=-1] <- 1
  
  nFolds <- 10
  foldid <- rep(seq(nFolds), length.out = nrow(traindat))
  sae_fit <- cv.glmnet(x=as.matrix(traindat[,3:dim(traindat)[2]]),
                       y=Surv(traindat$t_saes,traindat$sae),
                       family="cox",
                       alpha=0.75,
                       penalty.factor=fpf,
                       foldid = foldid, 
                       standardize = F)
  
  # baseline survival & counterfactuals
  S0 <- base_surv(fit = sae_fit, 
                  Y = traindat$t_saes,
                  D = traindat$sae,
                  x = as.matrix(traindat[,3:dim(traindat)[2]]))
  
  est_S1_sae <- pred_surv(fit = sae_fit,
                          S0 = S0,
                          x = as.matrix(testdat1[,3:dim(testdat1)[2]]),
                          times = t50)
  
  est_S0_sae <- pred_surv(fit = sae_fit,
                          S0 = S0,
                          x = as.matrix(testdat0[,3:dim(testdat0)[2]]),
                          times = t50)
  
  Scoxph_treated_sae[testIndexes] <- 1 - est_S1_sae
  Scoxph_control_sae[testIndexes] <- 1 - est_S0_sae
}

# save results 
risk_coxph_EN_sae <- data.frame(index, Scoxph_treated_sae, Scoxph_control_sae)
write.csv(risk_coxph_EN_sae, paste0("./sae_EN_results/random_split_coxph_EN_sae",seed_num, ".csv"), row.names = F)





