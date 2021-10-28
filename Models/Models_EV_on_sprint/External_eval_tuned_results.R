
##---------------------- Summarize Raw Results for HTE Estimates -------------------##
##----------------------       Crystal Xu         03/15/2021     -------------------##
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("survcomp")
library(survival);library(ggplot2); library(gtsummary);library(boot);library(survcomp);library(nricens)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 

# summarize SPRINT data (external validation)
sprint_blvars <- read.csv("blvars.csv")
sprint_outcome <- read.csv("outcome.csv")
sprint_test <- merge(sprint_blvars, sprint_outcome, by="MASKID",all=T)
tmpdata <- sprint_test[complete.cases(sprint_test),]

t50 <- floor(365.25*3.26)

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
setwd("../Analyses Stanford Team/Analysis Results/SPRINT_external_results") 
cvdpath <- "coxph_AIC_cvd_sprint_external"
cvdpath <- "coxphEN_cvd_sprint_external"
cvdpath <- "external_sprint_bart_cvd"
cvdpath <- "external_sprint_deepsurv_cvd"
cvdpath <- "external_sprint_gbm_cvd"
cvdpath <- "external_sprint_sf_cvd"
#XLpath <- "XL_sf_CVD"

# # Load estimates
# XL_RD_cvd <- read.csv(paste0(XLpath,".csv")); head(XL_RD_cvd)
# sum_cvd_T1 <- read.csv(paste0(cvdpath,".csv")); head(sum_cvd_T1)
# sum_cvd <- data.frame(sum_cvd_S, sum_cvd_T0, sum_cvd_T1,XL_RD_cvd);colnames(sum_cvd)[5]<-"XL_RD_cvd"; head(sum_cvd)
# write.csv(sum_cvd, "external_sprint_bart_cvd.csv", row.names = F)
# 
# sum_cvd <- read.csv(paste0(cvdpath,".csv")); head(sum_cvd)
# XL_RD_cvd <- read.csv(paste0(XLpath,".csv")); head(XL_RD_cvd)
# sum_cvd <- data.frame(sum_cvd, XL_RD_cvd);colnames(sum_cvd)[5]<-"XL_RD_cvd"; head(sum_cvd)
# write.csv(sum_cvd, "external_sprint_sf_cvd.csv", row.names = F)

sum_cvd <- read.csv(paste0(cvdpath,".csv")); head(sum_cvd)

# ARD estimates 
sum_cvd$S_RD_cvd <- sum_cvd[,2] - sum_cvd[,1]
# sum_cvd$T_RD_cvd <- sum_cvd[,4] - sum_cvd[,3]
cvd_S_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = sum_cvd$S_RD_cvd)
# cvd_T_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = sum_cvd$T_RD_cvd)
# cvd_X_RD <- data.frame(W = tmpdata$INTENSIVE, Y = tmpdata$t_cvds, D = tmpdata$cvd, pred_benefit = sum_cvd$XL_RD_cvd)

# Calibration test and plot
S_predinc_cvd <- sum_cvd[,2]*tmpdata$INTENSIVE + sum_cvd[,1]*(1-tmpdata$INTENSIVE)
Spredgrp <- as.numeric(cut2(S_predinc_cvd,g=10))
GND_Scvd <- GND.calib(pred=S_predinc_cvd, tvar=tmpdata$t_cvds, out=tmpdata$cvd,
                      cens.t=t50, groups=Spredgrp, adm.cens=t50, name = "CVD Risk, S-learner Elastic Net Regularization")

# T_predinc_cvd <- sum_cvd[,4]*tmpdata$INTENSIVE + sum_cvd[,3]*(1-tmpdata$INTENSIVE)
# Spredgrp <- as.numeric(cut2(T_predinc_cvd,g=10))
# GND_Tcvd <- GND.calib(pred=T_predinc_cvd, tvar=tmpdata$t_cvds, out=tmpdata$cvd,
#                       cens.t=t50, groups=Spredgrp, adm.cens=t50, name = "CVD Risk, T-learner Random Survival Forest")

GNDall <- round(data.frame(rbind(GND_Scvd)),3)
GND <- GNDall[,c(4,5,6,7,8,9)]
colnames(GND) <- c("slope","slope_lb","slope_ub","intercept","int_lb","int_ub");GND

# C-statistic (time = t50)
Cindex <- concordance.index(x=S_predinc_cvd, surv.time=tmpdata$t_cvds, surv.event=tmpdata$cvd)
C_Scvd <- round(c(Cindex$c.index, Cindex$lower, Cindex$upper),2)
# Cindex <- concordance.index(x=T_predinc_cvd, surv.time=tmpdata$t_cvds, surv.event=tmpdata$cvd)
# C_Tcvd <- round(c(Cindex$c.index, Cindex$lower, Cindex$upper),2)
Cindexall <- data.frame(rbind(C_Scvd))
colnames(Cindexall) <- c("CB","CB_lb","CB_ub")
risk_eval <- data.frame(Cindexall, GND)

# C-for-benefit 
CB_Scvd <- cforbenefit(cvd_S_RD)
# CB_Tcvd <- cforbenefit(cvd_T_RD)
# CB_Xcvd <- cforbenefit(cvd_X_RD)
C_benefit_sum <- data.frame(rbind(CB_Scvd))
colnames(C_benefit_sum) <- c("CB","CB_lb","CB_ub","CB_p")

# Calibration slope for HTE
calib_Scvd <- HTEcalib(data = cvd_S_RD, times = t50, name = "CVD Risk Reduction, S-learner Elastic Net Regularization")
# calib_Tcvd <- HTEcalib(data = cvd_T_RD, times = t50, name = "CVD Risk Reduction, T-learner Random Survival Forest")
# calib_Xcvd <- HTEcalib(data = cvd_X_RD, times = t50, name = "CVD Risk Reduction, X-learner Random Survival Forest")
Calib_benefit_sum <- round(data.frame(rbind(calib_Scvd)),3)
colnames(Calib_benefit_sum) <- c("slope","slope_lb","slope_ub","intercept","int_lb","int_ub")
HTE_eval <- data.frame(C_benefit_sum[,1:3], Calib_benefit_sum)

# Save all results in one table
eval_sum <- rbind(risk_eval, HTE_eval)
eval_sum <- data.frame(rownames(eval_sum),eval_sum)
write.csv(eval_sum, "EN_ev_sprint_sum.csv", row.names = F)

