
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

## ACCORD data as testing 
accord_blvars <- read.csv("ACCORD_blvars.csv")
accord_outcome <- read.csv("ACCORD_outcomes.csv")
accord_outcome$t_cvds <- ceiling(accord_outcome$t_cvds)
accord_outcome$t_saes <- ceiling(accord_outcome$t_saes)
accord_train <- merge(accord_blvars, accord_outcome, by="MASKID",all=T)
accord_train <- accord_train[complete.cases(accord_train),]; dim(accord_train) # 4448 43

# Recode treatment variable according to Tian (2015)
# "A simple method for estimating interactions between a treatment and a large number of covariates"
accord_train$INTENSIVE[accord_train$INTENSIVE==0] <- -1  
accord_train$INTENSIVE <- accord_train$INTENSIVE/2  
Xvar <- accord_train[,2:dim(accord_blvars)[2]]
Xvar_int <- model.matrix(~INTENSIVE + . + INTENSIVE:.-1, Xvar)
Xvar_int[,2:dim(Xvar_int)[2]] <- scale(Xvar_int[,2:dim(Xvar_int)[2]]) # standardize the covariates

# Combine baseline variables and outcomes into one data frame 
sprint_blvars <- read.csv("blvars.csv")
sprint_outcome <- read.csv("outcome.csv")
sprint_test <- merge(sprint_blvars, sprint_outcome, by="MASKID",all=T)
sprint_test <- sprint_test[complete.cases(sprint_test),]; dim(sprint_test) # 8969 43

# Set counterfactual treatment levels 
testXvar <- sprint_test[,2:dim(sprint_blvars)[2]]
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
accord_cvd <- data.frame(accord_train$t_cvds, accord_train$cvd, Xvar_int)
names(accord_cvd)[1:2] <- c("t_cvds","cvd")
sprint_cvd <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar_int)
names(sprint_cvd)[1:2] <- c("t_cvds","cvd")
sprint_cvd1 <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar_int1)
names(sprint_cvd1)[1:2] <- c("t_cvds","cvd")
sprint_cvd0 <- data.frame(sprint_test$t_cvds, sprint_test$cvd, testXvar_int0)
names(sprint_cvd0)[1:2] <- c("t_cvds","cvd")

# S-Learner
# Tuning parameter: shrinkage level 
fpf <- c(0,rep(1, (dim(accord_cvd)[2]-3)))
fpf[regexpr("INTENSIVE.", colnames(accord_cvd[,3:dim(accord_cvd)[2]]))!=-1] <- 1

# Training model using SPRINT
nFolds <- 10 
foldid <- rep(seq(nFolds), length.out = nrow(accord_cvd))
cvd_fit <- cv.glmnet(x=as.matrix(Xvar_int),
                     y=Surv(accord_cvd$t_cvds,accord_cvd$cvd),
                     family="cox",
                     alpha=0.75,
                     penalty.factor=fpf,
                     foldid = foldid, 
                     standardize = F)

# Model on sprint data to get the S0 in sprint
foldid <- rep(seq(nFolds), length.out = nrow(sprint_cvd))
cvd_fit_sprint <- cv.glmnet(x=as.matrix(testXvar_int),
                            y=Surv(sprint_cvd$t_cvds, sprint_cvd$cvd),
                            family="cox",
                            alpha=0.75,
                            penalty.factor=fpf,
                            foldid = foldid, 
                            standardize = F)
S0_sprint <- base_surv(fit = cvd_fit_sprint, 
                         Y = sprint_cvd$t_cvds,
                         D = sprint_cvd$cvd,
                         x = as.matrix(testXvar_int))

pred_S1_cvd <- pred_surv(fit = cvd_fit, S0 = S0_sprint, x = as.matrix(testXvar_int1), times = t50)
pred_S0_cvd <- pred_surv(fit = cvd_fit, S0 = S0_sprint, x = as.matrix(testXvar_int0), times = t50)
Scoxph_treated_cvd <- 1 - pred_S1_cvd
Scoxph_control_cvd <- 1 - pred_S0_cvd

# save results 
risk_coxph_EN_cvd <- data.frame(Scoxph_control_cvd, Scoxph_treated_cvd)
write.csv(risk_coxph_EN_cvd, "../Analyses Stanford Team/Analysis Results/SPRINT_external_results/coxphEN_cvd_sprint_external.csv", row.names = F)


