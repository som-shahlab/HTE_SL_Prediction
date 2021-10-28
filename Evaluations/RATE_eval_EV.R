library(survival);library(ggplot2); library(gtsummary);library(boot);library(survcomp);library(nricens)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
# summarize ACCORD data (external validation)
blvars <- read.csv("ACCORD_blvars.csv")
outcome <- read.csv("ACCORD_outcomes.csv")
outcome$t_cvds <- ceiling(outcome$t_cvds)
outcome$t_saes <- ceiling(outcome$t_saes)
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 4448
t50 <- floor(365.25*3.26)

setwd("../Analyses Stanford Team/Analysis Results/External Results") 
cvdpath <- "external_bart_cvd_sum"
saepath <- "external_bart_sae_sum"
cvdpath <- "external_sf_cvd_sum"
saepath <- "external_sf_sae_sum"
cvdpath <- "deepsurv_external_cvd_sum"
saepath <- "deepsurv_external_sae_sum"
cvdpath <- "external_gbm_cvd_sum"
saepath <- "external_gbm_sae_sum"
cvdpath <- "risk_coxph_AIC_cvd_external"
saepath <- "risk_coxph_AIC_sae_external"
cvdpath <- "coxphEN_cvd_external"
saepath <- "coxphEN_sae_external"

xl_cvd <- read.csv("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/External Results/X-learner/gbm_CVD.csv")
xl_sae <- read.csv("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/External Results/X-learner/gbm_SAE.csv")

# Load estimates 
sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)
sum_sae <- read.csv(paste0(saepath,".csv"));head(sum_sae)
sum_cvd$XL_RD_cvd <- unlist(xl_cvd)
sum_sae$XL_RD_sae <- unlist(xl_sae)
write.csv(sum_cvd, paste0(cvdpath, ".csv"))
write.csv(sum_sae, paste0(saepath, ".csv"))
# sum_cvd <- sum_cvd[order(sum_cvd$index),];head(sum_cvd)
# sum_sae <- sum_sae[order(sum_sae$index),];head(sum_sae)

# ARD estimates 
sum_cvd$S_RD_cvd <- sum_cvd[,1] - sum_cvd[,2]
sum_cvd$T_RD_cvd <- sum_cvd[,4] - sum_cvd[,3]
sum_sae$S_RD_sae <- sum_sae[,1] - sum_sae[,2]
sum_sae$T_RD_sae <- sum_sae[,4] - sum_sae[,3]

# RATE 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/grf/r-package/grf/R")
files.sources = list.files()
sapply(files.sources, source)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/External Results/Fig_RATE")

# RATE for cvd outcomes
set.seed(99)
tmpdata$cvd[tmpdata$t_cvds >= t50] <- 1
tmpdata$t_cvds[tmpdata$t_cvds >= t50] <- t50
cf_cvd <- causal_survival_forest(X = tmpdata[,3:39],
                                 Y = tmpdata$t_cvds, 
                                 D = tmpdata$cvd,
                                 W = tmpdata$INTENSIVE,
                                 failure.times = seq(min(tmpdata$t_cvds), max(tmpdata$t_cvds), length.out = 50))

Stau.hats <- -sum_cvd$S_RD_cvd
Ttau.hats <- -sum_cvd$T_RD_cvd
Xtau.hats <- -sum_cvd$XL_RD_cvd

S_priority_cvd <- cut(Stau.hats, breaks = quantile(Stau.hats), include.lowest = TRUE)
T_priority_cvd <- cut(Ttau.hats, breaks = quantile(Ttau.hats), include.lowest = TRUE)
X_priority_cvd <- cut(Xtau.hats, breaks = quantile(Xtau.hats), include.lowest = TRUE)

set.seed(99)
Srate_cvd <- rank_average_treatment_effect(cf_cvd, S_priority_cvd);Srate_cvd
Trate_cvd <- rank_average_treatment_effect(cf_cvd, T_priority_cvd);Trate_cvd
Xrate_cvd <- rank_average_treatment_effect(cf_cvd, X_priority_cvd);Xrate_cvd

# Plot the Targeting Operator Characteristic curve.
filename <- paste0("./EN_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Srate_cvd, main="CVD Risk Reduction, Elastic Net Regularization")
dev.off()

filename <- paste0("./Tdeepsurv_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Trate_cvd, main="CVD Risk Reduction, T-learner Deep Survival Learning")
dev.off()

filename <- paste0("./Xgbm_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Xrate_cvd, main="CVD Risk Reduction, X-learner Gradient Boosting Machine" )
dev.off()

# for sae outcomes
set.seed(99)
tmpdata$sae[tmpdata$t_saes >= t50] <- 1
tmpdata$t_saes[tmpdata$t_saes >= t50] <- t50
cf_sae <- causal_survival_forest(X = tmpdata[,3:39],
                                 Y = tmpdata$t_saes, 
                                 D = tmpdata$sae,
                                 W = tmpdata$INTENSIVE,
                                 failure.times = seq(min(tmpdata$t_saes), max(tmpdata$t_saes), length.out = 50))

Stau.hats <- sum_sae$S_RD_sae
Ttau.hats <- sum_sae$T_RD_sae
Xtau.hats <- sum_sae$XL_RD_sae

S_priority_sae <- cut(Stau.hats, breaks = quantile(Stau.hats), include.lowest = TRUE)
T_priority_sae <- cut(Ttau.hats, breaks = quantile(Ttau.hats), include.lowest = TRUE)
X_priority_sae <- cut(Xtau.hats, breaks = quantile(Xtau.hats), include.lowest = TRUE)

set.seed(99)
Srate_sae <- rank_average_treatment_effect(cf_sae, S_priority_sae);Srate_sae
Trate_sae <- rank_average_treatment_effect(cf_sae, T_priority_sae);Trate_sae
Xrate_sae <- rank_average_treatment_effect(cf_sae, X_priority_sae);Xrate_sae

# Plot the Targeting Operator Characteristic curve.
filename <- paste0("./EN_rate_sae.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Srate_sae, main="SAE Risk Increase, Elastic Net Regularization")
dev.off()

filename <- paste0("./Tdeepsurv_rate_sae.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Trate_sae, main="SAE Risk Increase, T-learner Deep Survival Learning")
dev.off()

filename <- paste0("./Xgbm_rate_sae.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Xrate_sae, main="SAE Risk Increase, X-learner Gradient Boosting Machine")
dev.off()


RATE <- round(sapply(data.frame(rbind(Xrate_cvd, Xrate_sae)[,1:2]), as.numeric),2)
#RATE <- round(sapply(data.frame(rbind(Srate_cvd, Srate_sae)[,1:2]), as.numeric),2)
RATE <- paste0(RATE[,1], " (", RATE[,2], ")");RATE
# names <- c("SL_cvd","TL_cvd","XL_cvd","SL_sae","TL_sae","XL_sae")
# EN_RATE <- c(RATE[1],NA,NA,RATE[2], NA, NA)
deepsurv_RATE <- RATE

RATE_eval <- data.frame(names, gbm_RATE, bart_RATE,sf_RATE,deepsurv_RATE,AIC_RATE,EN_RATE);RATE_eval
write.csv(RATE_eval, "RATE_evals.csv", row.names = F)






