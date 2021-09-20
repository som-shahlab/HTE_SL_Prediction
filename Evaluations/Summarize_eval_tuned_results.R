
##---------------------- Summarize Raw Results for HTE Estimates -------------------##
##----------------------       Crystal Xu         03/15/2021     -------------------##
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("survcomp")
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

# For combined results 
# tmpdata <- read.csv("sprint_accord.csv")


##------------------------------ Descriptive Table 1 --------------------------------------##
# vars <- tmpdata[,!names(tmpdata) %in% c("MASKID","t_cvds","cvd","t_saes","sae")]
# Table1 <- vars %>%
#   tbl_summary(
#     by = INTENSIVE,
#     type =  list(N_AGENTS ~ "continuous"),
#     statistic = list(all_continuous() ~ "{mean} ({sd})",
#                      all_categorical() ~ "{n} ({p}%)"),
#     digits = all_continuous() ~ 1,
#   )%>%
#   add_p(pvalue_fun = ~style_pvalue(.x, digits = 2))
# 
# write.csv(Table1, "Table1.csv",row.names = F)   

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
# cvdpath <- "coxph_EN_cvd_results/random_split_coxph_EN_cvd2"
# saepath <- "coxph_EN_sae_results/random_split_coxph_EN_sae2"
# cvdpath <- "coxph_AIC_cvd_results/random_split_coxph_AIC_cvd2"
# saepath <- "coxph_AIC_sae_results/random_split_coxph_AIC_sae2"

# Load estimates 
sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)
#sum_cvd <- sum_cvd[order(sum_cvd$index),];head(sum_cvd)
sum_sae <- read.csv(paste0(saepath,".csv"));head(sum_sae)
#sum_sae <- sum_sae[order(sum_sae$index),];head(sum_sae)

# variable importance of Sbart
cvdpath <- "bart_cvd_tuned_results/results_stratified/Sbart_varimp"
varimp_cvd <- matrix(NA, 38, 10)
for (i in 1:10){
  varimp_cvd[,i] <- data.frame(read.csv(paste0(cvdpath,i,".csv")))[,2]
} 
cvd_imp <- data.frame(tmpres[,1], rowMeans(varimp_cvd))
colnames(cvd_imp) <- c("varnames", "varimp")
cvd_imp_Sbart <- cvd_imp[order(-cvd_imp$varimp),]
write.csv(cvd_imp_Sbart, "cvd_imp_Sbart.csv", row.names = F)

# sum_cvd <- sum_cvd[order(sum_cvd$testIndexes),];head(sum_cvd)
# 
# sum_cvd_S <- sum_cvd
# sum_cvd_T0 <- sum_cvd
# sum_cvd_T1 <- sum_cvd
# sum_cvd <- data.frame(sum_cvd_S, sum_cvd_T0[,2],sum_cvd_T1[,2])
# write.csv(sum_cvd, "stratified_sum_cvd_bart.csv", row.names = F)

# sum_sae <- NULL
# for (i in 1:10){
#   tmpres <- read.csv(paste0(saepath,i,".csv"))
#   sum_sae <- rbind(sum_sae, tmpres)
# }
# sum_sae <- sum_sae[order(sum_sae$testIndexes),];head(sum_sae)
# sum_sae_S <- sum_sae
# sum_sae_T0 <- sum_sae
# sum_sae_T1 <- sum_sae
# sum_sae <- data.frame(sum_sae_S, sum_sae_T0[,2],sum_sae_T1[,2])
# write.csv(sum_sae, "stratified_sum_sae_bart.csv", row.names = F)


# ARD estimates 
sum_cvd$S_RD_cvd <- sum_cvd[,2] - sum_cvd[,1]
sum_cvd$T_RD_cvd <- sum_cvd[,4] - sum_cvd[,3]
sum_sae$S_RD_sae <- sum_sae[,2] - sum_sae[,1]
sum_sae$T_RD_sae <- sum_sae[,4] - sum_sae[,3]

cvd_S_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = sum_cvd$S_RD_cvd)
cvd_T_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = sum_cvd$T_RD_cvd)
cvd_X_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = sum_cvd$XL_RD_cvd)
sae_S_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_saes, D = tmpdata$sae, pred_benefit = sum_sae$S_RD_sae)
sae_T_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_saes, D = tmpdata$sae, pred_benefit = sum_sae$T_RD_sae)
sae_X_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_saes, D = tmpdata$sae, pred_benefit = sum_sae$XL_RD_sae)

# Calibration test and plot
S_predinc_cvd <- sum_cvd[,2]*tmpdata$INTENSIVE + sum_cvd[,1]*(1-tmpdata$INTENSIVE)
Spredgrp <- as.numeric(cut2(S_predinc_cvd,g=10))
GND_Scvd <- GND.calib(pred=S_predinc_cvd, tvar=tmpdata$t_cvds, out=tmpdata$cvd,
                      cens.t=t50, groups=Spredgrp, adm.cens=t50, name = "CVD Risk, S-learner Deep Survival Learning", calib_plot = F)

T_predinc_cvd <- sum_cvd[,4]*tmpdata$INTENSIVE + sum_cvd[,3]*(1-tmpdata$INTENSIVE)
Spredgrp <- as.numeric(cut2(T_predinc_cvd,g=10))
GND_Tcvd <- GND.calib(pred=T_predinc_cvd, tvar=tmpdata$t_cvds, out=tmpdata$cvd,
                      cens.t=t50, groups=Spredgrp, adm.cens=t50, name = "CVD Risk, T-learner Deep Survival Learning", calib_plot = F)

S_predinc_sae <- sum_sae[,2]*tmpdata$INTENSIVE + sum_sae[,1]*(1-tmpdata$INTENSIVE)
Spredgrp <- as.numeric(cut2(S_predinc_sae,g=10))
GND_Ssae <- GND.calib(pred=S_predinc_sae, tvar=tmpdata$t_saes, out=tmpdata$sae,
                      cens.t=t50, groups=Spredgrp, adm.cens=t50, name = "SAE Risk, S-learner Deep Survival Learning", calib_plot = F)

T_predinc_sae <- sum_sae[,4]*tmpdata$INTENSIVE + sum_sae[,3]*(1-tmpdata$INTENSIVE)
Spredgrp <- as.numeric(cut2(T_predinc_sae,g=10))
GND_Tsae <- GND.calib(pred=T_predinc_sae, tvar=tmpdata$t_saes, out=tmpdata$sae,
                      cens.t=t50, groups=Spredgrp, adm.cens=t50, name = "SAE Risk, T-learner Deep Survival Learning", calib_plot = F)

GNDall <- round(data.frame(rbind(GND_Scvd, GND_Tcvd, GND_Ssae, GND_Tsae)),3)
GNDall$slope <- paste0(GNDall[,4], " (", GNDall[,5], ", ", GNDall[,6],")")
GNDall$intercept <- paste0(GNDall[,7], " (", GNDall[,8], ", ", GNDall[,9],")")
GND <- GNDall[,c(3,10,11)];GND

# C-statistic (time = t50)
Cindex <- concordance.index(x=S_predinc_cvd, surv.time=tmpdata$t_cvds, surv.event=tmpdata$cvd)
C_Scvd <- round(c(Cindex$c.index, Cindex$lower, Cindex$upper),2)

Cindex <- concordance.index(x=T_predinc_cvd, surv.time=tmpdata$t_cvds, surv.event=tmpdata$cvd)
C_Tcvd <- round(c(Cindex$c.index, Cindex$lower, Cindex$upper),2)

Cindex <- concordance.index(x=S_predinc_sae, surv.time=tmpdata$t_saes, surv.event=tmpdata$sae)
C_Ssae <- round(c(Cindex$c.index, Cindex$lower, Cindex$upper),2)

Cindex <- concordance.index(x=T_predinc_sae, surv.time=tmpdata$t_saes, surv.event=tmpdata$sae)
C_Tsae <- round(c(Cindex$c.index, Cindex$lower, Cindex$upper),2)

Cindexall <- data.frame(rbind(C_Scvd, C_Tcvd, C_Ssae, C_Tsae))
Cindex <- paste0(Cindexall[,1], " (", Cindexall[,2], ", ", Cindexall[,3],")")
risk_eval <- data.frame(Cindex, GND)

# C-for-benefit 
CB_Scvd <- cforbenefit(cvd_S_RD)
CB_Tcvd <- cforbenefit(cvd_T_RD)
CB_Xcvd <- cforbenefit(cvd_X_RD)
CB_Ssae <- cforbenefit(sae_S_RD)
CB_Tsae <- cforbenefit(sae_T_RD)
CB_Xsae <- cforbenefit(sae_X_RD)
C_benefit_sum <- data.frame(rbind(CB_Xcvd,CB_Xsae));C_benefit_sum
Cbenefit <- data.frame(paste0(round(C_benefit_sum[,1],2), " (", round(C_benefit_sum[,2],2), ", ", round(C_benefit_sum[,3],2),")"))
colnames(Cbenefit) <-  "C-for-benefit";Cbenefit

# Calibration slope for HTE
calib_Scvd <- HTEcalib(data = cvd_S_RD, times = t50, name = "CVD Risk Reduction, S-learner Deep Survival Learning", calib_plot = F)
calib_Tcvd <- HTEcalib(data = cvd_T_RD, times = t50, name = "CVD Risk Reduction, T-learner Deep Survival Learning", calib_plot = F)
calib_Xcvd <- HTEcalib(data = cvd_X_RD, times = t50, name = "CVD Risk Reduction, X-learner Random Survival Forest", calib_plot = T)
calib_Ssae <- HTEcalib(data = sae_S_RD, times = t50, name = "SAE Risk Increase, S-learner Deep Survival Learning", calib_plot = F)
calib_Tsae <- HTEcalib(data = sae_T_RD, times = t50, name = "SAE Risk Increase, T-learner Deep Survival Learning", calib_plot = F)
calib_Xsae <- HTEcalib(data = sae_X_RD, times = t50, name = "SAE Risk Increase, X-learner Random Survival Forest", calib_plot = T)
Calib_benefit_sum <- round(data.frame(rbind(calib_Xcvd, calib_Xsae)),3)
slope <- paste0(Calib_benefit_sum[,1], " (", Calib_benefit_sum[,2], ", ", Calib_benefit_sum[,3],")")
intercept <- paste0(Calib_benefit_sum[,4], " (", Calib_benefit_sum[,5], ", ", Calib_benefit_sum[,6],")")
#write.csv(HTE_eval, "SPRINT_HTE_eval_AIC.csv", row.names = F)


# RATE 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/grf/r-package/grf/R")
files.sources = list.files()
sapply(files.sources, source)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/Tuned Results")

# RATE for cvd outcomes
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
Xrate_cvd <- rank_average_treatment_effect(cf_cvd, X_priority_cvd);Xrate_cvd

# Plot the Targeting Operator Characteristic curve.
filename <- paste0("./Xdeepsurv_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Srate_cvd)
text(0.8, -15, paste0("AUTOC = ", round(Srate_cvd$estimate,2)))
dev.off()

png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Trate_cvd)
text(0.8, -10, paste0("AUTOC = ", round(Trate_cvd$estimate,2)))
dev.off()

png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Xrate_cvd)
text(0.8, 15, paste0("AUTOC = ", round(Xrate_cvd$estimate,2)))
dev.off()

# for sae outcomes
tmpdata$sae[tmpdata$t_saes >= t50] <- 1
tmpdata$t_saes[tmpdata$t_saes >= t50] <- t50
cf_sae <- causal_survival_forest(X = tmpdata[,3:39],
                                 Y = tmpdata$t_saes, 
                                 D = tmpdata$sae,
                                 W = tmpdata$INTENSIVE,
                                 failure.times = seq(min(tmpdata$t_saes), max(tmpdata$t_saes), length.out = 100))

Stau.hats <- sum_sae$S_RD_sae
Ttau.hats <- sum_sae$T_RD_sae
Xtau.hats <- sum_sae$XL_RD_sae

S_priority_sae <- cut(Stau.hats, breaks = quantile(Stau.hats), include.lowest = TRUE)
T_priority_sae <- cut(Ttau.hats, breaks = quantile(Ttau.hats), include.lowest = TRUE)
X_priority_sae <- cut(Xtau.hats, breaks = quantile(Xtau.hats), include.lowest = TRUE)

Srate_sae <- rank_average_treatment_effect(cf_sae, S_priority_sae);Srate_sae
Trate_sae <- rank_average_treatment_effect(cf_sae, T_priority_sae);Trate_sae
Xrate_sae <- rank_average_treatment_effect(cf_sae, X_priority_sae);Xrate_sae

# Plot the Targeting Operator Characteristic curve.
filename <- paste0("./Sdeepsurv_rate_sae.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Srate_sae)
text(0.8, -15, paste0("AUTOC = ", round(Srate_sae$estimate,2)))
dev.off()

png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Trate_sae)
text(0.8, -10, paste0("AUTOC = ", round(Trate_sae$estimate,2)))
dev.off()

png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Xrate_sae)
text(0.8, -15, paste0("AUTOC = ", round(Xrate_sae$estimate,2)))
dev.off()

RATE <- round(sapply(data.frame(rbind(Srate_cvd,Trate_cvd,Xrate_cvd, Srate_sae, Trate_sae, Xrate_sae)[,1:2]), as.numeric),2)
RATE <- paste0(RATE[,1], " (", RATE[,2], ")")
HTE_eval <- data.frame(Cbenefit, RATE, slope, intercept)
rownames(HTE_eval) <- c("SL_cvd","TL_cvd","XL_cvd","SL_sae","TL_sae","XL_sae")

# Save all results in one table
eval_sum <- rbind(as.matrix(risk_eval), as.matrix(HTE_eval))
eval_sum <- data.frame(rownames(eval_sum),eval_sum)
write.csv(eval_sum, "deepsurv_iv_eval_sum.csv", row.names = F)


# # Tertiles of ARD
# S_tertile_table <- quan_sum_table(data=cvd_S_RD, times=t50)
# S_cvd_tertiles <- data.frame(rep("SLcvd",3), S_tertile_table)
# colnames(S_cvd_tertiles)[1] <- "Meta-Learner"
# 
# T_tertile_table <- quan_sum_table(data=cvd_T_RD, times=t50)
# T_cvd_tertiles <- data.frame(rep("TLcvd",3), T_tertile_table)
# colnames(T_cvd_tertiles)[1] <- "Meta-Learner"
# 
X_tertile_table <- quan_sum_table(data=cvd_X_RD, times=t50)
X_cvd_tertiles <- data.frame(rep("XLcvd",3), X_tertile_table)
colnames(X_cvd_tertiles)[1] <- "Meta-Learner"
# 
# S_tertile_table <- quan_sum_table(data=sae_S_RD, times=t50)
# S_sae_tertiles <- data.frame(rep("SLsae",3),S_tertile_table)
# colnames(S_sae_tertiles)[1] <- "Meta-Learner"
# 
# T_tertile_table <- quan_sum_table(data=sae_T_RD, times=t50)
# T_sae_tertiles <- data.frame(rep("TLsae",3),T_tertile_table)
# colnames(T_sae_tertiles)[1] <- "Meta-Learner"
# 
X_tertile_table <- quan_sum_table(data=sae_X_RD, times=t50)
X_sae_tertiles <- data.frame(rep("XLsae",3),X_tertile_table)
colnames(X_sae_tertiles)[1] <- "Meta-Learner"
# 
tertiles_sum <- data.frame(rbind(X_cvd_tertiles, X_sae_tertiles))
write.csv(tertiles_sum, "Xbart_tertiles_sum.csv")



# Reclassification table (only for CVD outcome)
# we cannot compare ARD estimates as the PCE and Framingham models did not use treatment as a covariate, so ARD=0
# ASCVDrisks <- read.csv("..//sprint_est_risks.csv")
# pcer <- ASCVDrisks[,1]
# famr <- ASCVDrisks[,2]
# # PCE 
# NRI_Scvd_pcer <- nricens(time = tmpdata$t_cvds, event = tmpdata$cvd, p.std = pcer, p.new = S_predinc_cvd, 
#                          t0 = t50, updown = 'category', cut = 0.07, point.method = 'ipw', niter = 0)
# NRI_Tcvd_pcer <- nricens(time = tmpdata$t_cvds, event = tmpdata$cvd, p.std = pcer, p.new = T_predinc_cvd, 
#                          t0 = t50, updown = 'category', cut = 0.07, point.method = 'ipw', niter = 0)
# 
# # Framingham 
# NRI_Scvd_famr <- nricens(time = tmpdata$t_cvds, event = tmpdata$cvd, p.std = famr, p.new = S_predinc_cvd,
#                          t0 = t50, updown = 'category', cut = 0.07, point.method = 'ipw', niter = 0)
# NRI_Tcvd_famr <- nricens(time = tmpdata$t_cvds, event = tmpdata$cvd, p.std = famr, p.new = T_predinc_cvd,
#                          t0 = t50, updown = 'category', cut = 0.07, point.method = 'ipw', niter = 0)
# 
# nri_risk <- data.frame(cbind(NRI_Scvd_pcer$nri,NRI_Tcvd_pcer$nri,NRI_Scvd_famr$nri,NRI_Tcvd_famr$nri))
# colnames(nri_risk) <- c("Scvd_pce","Tcvd_pce","Scvd_fam","Tcvd_fam")
# nri_risk <- data.frame(nri_risk,rep(NA,7),rep(NA,7))
