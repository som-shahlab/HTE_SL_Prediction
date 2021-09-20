library(survival);library(ggplot2); library(gtsummary);library(boot);library(survcomp);library(nricens)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
# Load SPRINT data (external validation)
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) 
t50 <- floor(365.25*3.26)

setwd("../Analyses Stanford Team/Analysis Results/SPRINT_external_results") 
cvdpath <- "external_sprint_gbm_cvd"

# Load estimates 
sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)

# ARD estimates 
sum_cvd$S_RD_cvd <- sum_cvd[,2] - sum_cvd[,1]
sum_cvd$T_RD_cvd <- sum_cvd[,4] - sum_cvd[,3]

# RATE 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/grf/r-package/grf/R")
files.sources = list.files()
sapply(files.sources, source)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/SPRINT_external_results")

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
filename <- paste0("./Sgbm_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Srate_cvd, main="CVD Risk Reduction, S-learner Gradient Boosting Machine")
dev.off()

filename <- paste0("./Tbart_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Trate_cvd, main="CVD Risk Reduction, T-learner BART")
dev.off()

filename <- paste0("./Xbart_rate_cvd.png")
png(filename, width = 6, height = 6, units = 'in', res = 300)
plot(Xrate_cvd, main="CVD Risk Reduction, X-learner BART")
dev.off()

RATE <- round(sapply(data.frame(rbind(Srate_cvd,Trate_cvd,Xrate_cvd)[,1:2]), as.numeric),2)
rownames(RATE) <- c("SL_cvd","TL_cvd","XL_cvd")
#RATE <- round(sapply(data.frame(Srate_cvd$estimate, Srate_cvd$std.err), as.numeric),2)
bart_RATE <- RATE

RATE_eval <- data.frame(rbind(AIC_RATE,EN_RATE,gbm_RATE,deepsurv_RATE,sf_RATE,bart_RATE))
RATE_eval$Method <- c("AIC","EN", rep("GBM",3),rep("Deepsurv",3),rep("RSF",3),rep("BART",3))
RATE_eval$ML <- c("SL","SL",rep(c("SL","TL","XL"),4))
write.csv(RATE_eval, "EV_sprint_RATE_evals.csv", row.names = F)






