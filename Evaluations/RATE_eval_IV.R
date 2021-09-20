library(survival);library(ggplot2);library(gtsummary);library(boot);library(survcomp);library(nricens); library(png)
library(gridExtra);library(grid);library(lattice); library(grf)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
t50 <- floor(365.25*3.26)

# SPRINT data
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969

setwd("../Analyses Stanford Team/Analysis Results/Tuned Results") 
cvdpath <- "bart_cvd_tuned_results/stratified_sum_cvd_bart"
saepath <- "bart_sae_tuned_results/stratified_sum_sae_bart"
cvdpath <- "deepsurv_cvd_tuned_results/deepsurv_tuned_cvd_sum"
saepath <- "deepsurv_sae_tuned_results/deepsurv_tuned_sae_sum"
cvdpath <- "gbm_cvd_tuned_results/stratified_split_gbm_cvd_sum"
saepath <- "gbm_sae_tuned_results/stratified_split_gbm_sae_sum"
cvdpath <- "sf_cvd_tuned_results/stratified_split_sf_cvd_sum"
saepath <- "sf_sae_tuned_results/stratified_split_sf_sae_sum"
cvdpath <- "coxph_EN_cvd_results/random_split_coxph_EN_cvd2"
saepath <- "coxph_EN_sae_results/random_split_coxph_EN_sae2"
cvdpath <- "coxph_AIC_cvd_results/random_split_coxph_AIC_cvd2"
saepath <- "coxph_AIC_sae_results/random_split_coxph_AIC_sae2"

# Load estimates 
sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)
sum_sae <- read.csv(paste0(saepath,".csv"));head(sum_sae)
# sum_cvd <- sum_cvd[order(sum_cvd$index),];head(sum_cvd)
# sum_sae <- sum_sae[order(sum_sae$index),];head(sum_sae)

# ARD estimates 
# sum_cvd$S_RD_cvd <- sum_cvd[,2] - sum_cvd[,1]
# sum_cvd$T_RD_cvd <- sum_cvd[,4] - sum_cvd[,3]
# sum_sae$S_RD_sae <- sum_sae[,2] - sum_sae[,1]
# sum_sae$T_RD_sae <- sum_sae[,4] - sum_sae[,3]

# RATE 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/grf/r-package/grf/R")
files.sources = list.files()
sapply(files.sources, source)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/Tuned Results/Fig_RATE")

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


Srate_cvd <- rank_average_treatment_effect(cf_cvd, S_priority_cvd);Srate_cvd
Trate_cvd <- rank_average_treatment_effect(cf_cvd, T_priority_cvd);Trate_cvd
set.seed(99)
Xrate_cvd <- rank_average_treatment_effect(cf_cvd, X_priority_cvd);Xrate_cvd

# Plot the Targeting Operator Characteristic curve.
# filename <- paste0("./Sbart_rate_cvd.png")
# png(filename, width = 6, height = 6, units = 'in', res = 300)
# plot(Srate_cvd, main="CVD Risk Reduction, S-learner BART")
# dev.off()
# 
# filename <- paste0("./Tbart_rate_cvd.png")
# png(filename, width = 6, height = 6, units = 'in', res = 300)
# plot(Trate_cvd, main="CVD Risk Reduction, T-learner BART")
# dev.off()

filename <- paste0("./Xsf_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Xrate_cvd, main="CVD Risk Reduction, X-learner Random Survival Forest" )
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


Srate_sae <- rank_average_treatment_effect(cf_sae, S_priority_sae);Srate_sae
Trate_sae <- rank_average_treatment_effect(cf_sae, T_priority_sae);Trate_sae
set.seed(99)
Xrate_sae <- rank_average_treatment_effect(cf_sae, X_priority_sae);Xrate_sae

# Plot the Targeting Operator Characteristic curve.
# filename <- paste0("./Sbart_rate_sae.png")
# png(filename, width = 6, height = 6, units = 'in', res = 300)
# plot(Srate_sae, main="SAE Risk Increase, S-learner Gradient Boosting Machine")
# dev.off()
# 
# filename <- paste0("./Tbart_rate_sae.png")
# png(filename, width = 6, height = 6, units = 'in', res = 300)
# plot(Trate_sae, main="SAE Risk Increase, T-learner Gradient Boosting Machine")
# dev.off()

filename <- paste0("./Xgbm_rate_sae.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Xrate_sae, main="SAE Risk Increase, X-learner BART")
dev.off()


RATE <- round(sapply(data.frame(rbind(Xrate_cvd,Xrate_sae)[,1:2]), as.numeric),2)
#RATE <- round(sapply(data.frame(rbind(Srate_cvd, Srate_sae)[,1:2]), as.numeric),2)
RATE <- paste0(RATE[,1], " (", RATE[,2], ")");RATE
# names <- c("SL_cvd","TL_cvd","XL_cvd","SL_sae","TL_sae","XL_sae")
# EN_RATE <- c(RATE[1],NA,NA,RATE[2], NA, NA)
sf_RATE <- RATE
#"-2.47 (3.23)" "-4.97 (3.21)" "2.44 (2.88)"  "-3.32 (3.78)" "6.29 (3.72)"  "-6.4 (3.49)" 
RATE_eval <- data.frame(RATE_eval, sf_RATE);RATE_eval
write.csv(RATE_eval, "RATE_evals.csv", row.names = F)






