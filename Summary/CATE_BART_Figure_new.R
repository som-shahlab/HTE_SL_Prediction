library(survival);library(ggplot2);library(gtsummary);library(boot);library(survcomp);library(nricens); library(png)
library(gridExtra);library(grid);library(lattice)
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/eval_funs.R")
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data") 
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969
t50 <- floor(365.25*3.26)

setwd("../Analyses Stanford Team/Analysis Results/Tuned Results") 

# Plot distribution of ARDs 
cvdpath <- "bart_cvd_tuned_results/stratified_sum_cvd_bart"
saepath <- "bart_sae_tuned_results/stratified_sum_sae_bart"
sum_cvd <- read.csv(paste0(cvdpath,".csv"));head(sum_cvd)
sum_sae <- read.csv(paste0(saepath,".csv"));head(sum_sae)

# Scatter plot 
DF <- data.frame(sum_cvd$XL_RD_cvd, sum_sae$XL_RD_sae);colnames(DF) <- c("XL_RD_cvd","XL_RD_sae")
DF$Group <- ifelse(DF$XL_RD_cvd>=0 & DF$XL_RD_sae<=0, "Negative Benefits, Negative Harms", 
                   ifelse(DF$XL_RD_cvd<0 & DF$XL_RD_sae<=0, "Postive Benefits, Negative Harms",
                          ifelse(DF$XL_RD_cvd<0 & DF$XL_RD_sae>0, "Postive Benefits, Postive Harms", "Negative Benefits, Postive Harms")))
DF$Group <- as.factor(DF$Group)
png("./scatter_ARDs.png", width = 4, height = 5, units = 'in', res = 300)
ggplot(DF, aes(x=XL_RD_sae, y=XL_RD_cvd, col=Group)) + 
  theme_classic()+
  geom_point() + 
  xlab("Treatment Harms")+
  ylab("Treatment Benefits")+
  xlim(-0.2,0.2)+
  ylim(-0.2,0.2)+
  geom_hline(yintercept = 0,col="black", linetype="dashed", size=1)+
  geom_vline(xintercept = 0,col="black", linetype="dashed", size=1)+
  theme(legend.position="bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
dev.off()

# DF <- data.frame(c(sum_cvd$XL_RD_cvd, sum_sae$XL_RD_sae), 
#                  c(rep("CVD event",length(sum_cvd$XL_RD_cvd)),rep("Serious adverse event", length(sum_sae$XL_RD_sae))))
# colnames(DF) <- c("ARD", "Outcome")
# png("./dist_ARDs.png", width = 7, height = 5, units = 'in', res = 300)
#  DF %>%
#   ggplot( aes(x=ARD, fill=Outcome)) +
#   
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   xlab("Absolute Risk Difference")+
#   labs(fill="")
# dev.off()


# Top 15 important variables 
#impSbart <- read.csv("./cvd_imp_Sbart.csv")
impvars <- read.csv("./IV_imp_bart.csv")
varimp_cvd <- c(as.character(impvars$cvd_vars[1:15]))
varimp_sae <- c(as.character(impvars$sae_vars[1:15]))

## CVD graph
# UMALCR
tmpdata$UMALCR <-  round(tmpdata$UMALCR,0)
group <- cut(tmpdata$UMALCR, breaks = unique(quantile(tmpdata$UMALCR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_UMALCR_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(group),FUN=mean))
CATE_UMALCR_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(group),FUN=sd))
CATE_UMALCR_ave$LB <- CATE_UMALCR_ave$x-CATE_UMALCR_sd$x
CATE_UMALCR_ave$UB <- CATE_UMALCR_ave$x+CATE_UMALCR_sd$x
CATE_UMALCR_ave$N <- table(group)
colnames(CATE_UMALCR_ave) <- c("UMALCR","CATE","LB","UB","N")
UMALCRp <- ggplot(CATE_UMALCR_ave, aes(x = UMALCR, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Urine Albumin/Creatinine ratio")+
  ylab("Estimated HTEs")

# Age
Age_mod <- cut(tmpdata$AGE, breaks = unique(quantile(tmpdata$AGE, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_age_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(Age_mod),FUN=mean))
CATE_age_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(Age_mod),FUN=sd))
CATE_age_ave$LB <- CATE_age_ave$x - CATE_age_sd$x
CATE_age_ave$UB <- CATE_age_ave$x + CATE_age_sd$x
CATE_age_ave$N <- table(Age_mod)
colnames(CATE_age_ave) <- c("Age","CATE","LB","UB","N")
Agep <- ggplot(CATE_age_ave, aes(x = Age, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Age")+
  ylab("Estimated HTEs")

# eGFR
tmpdata$EGFR <-  round(tmpdata$EGFR,0)
eGFR_mod <- cut(tmpdata$EGFR, breaks = unique(quantile(tmpdata$EGFR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_egfr_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(eGFR_mod),FUN=mean))
CATE_egfr_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(eGFR_mod),FUN=sd))
CATE_egfr_ave$LB <- CATE_egfr_ave$x - CATE_egfr_sd$x
CATE_egfr_ave$UB <- CATE_egfr_ave$x + CATE_egfr_sd$x
CATE_egfr_ave$N <- table(eGFR_mod)
colnames(CATE_egfr_ave) <- c("eGFR","CATE","LB","UB","N")
eGFRp <- ggplot(CATE_egfr_ave, aes(x = eGFR, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Estimated glomerular filtration rate")+
  ylab("Estimated HTEs")

# BMI
tmpdata$BMI <- round(tmpdata$BMI,0)
BMI_mod <- cut(tmpdata$BMI, breaks = unique(quantile(tmpdata$BMI, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_BMI_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(BMI_mod),FUN=mean))
CATE_BMI_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(BMI_mod),FUN=sd))
CATE_BMI_ave$LB <- CATE_BMI_ave$x-CATE_BMI_sd$x
CATE_BMI_ave$UB <- CATE_BMI_ave$x+CATE_BMI_sd$x
CATE_BMI_ave$N <- table(BMI_mod)
colnames(CATE_BMI_ave) <- c("BMI","CATE","LB","UB","N")
BMIp <- ggplot(CATE_BMI_ave, aes(x = BMI, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Boday mass index")+
  ylab("Estimated HTEs")

# CHR
tmpdata$CHR <- round(tmpdata$CHR,0)
CHR_mod <- cut(tmpdata$CHR, breaks = unique(quantile(tmpdata$CHR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_CHR_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(CHR_mod),FUN=mean))
CATE_CHR_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(CHR_mod),FUN=sd))
CATE_CHR_ave$LB <- CATE_CHR_ave$x-CATE_CHR_sd$x
CATE_CHR_ave$UB <- CATE_CHR_ave$x+CATE_CHR_sd$x
CATE_CHR_ave$N <- table(CHR_mod)
colnames(CATE_CHR_ave) <- c("CHR","CATE","LB","UB","N")
CHRp <- ggplot(CATE_CHR_ave, aes(x = CHR, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Total cholesterol")+
  ylab("Estimated HTEs")

# TRR
tmpdata$TRR <- round(tmpdata$TRR,0)
TRR_mod <- cut(tmpdata$TRR, breaks = unique(quantile(tmpdata$TRR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_TRR_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(TRR_mod),FUN=mean))
CATE_TRR_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(TRR_mod),FUN=sd))
CATE_TRR_ave$LB <- CATE_TRR_ave$x-CATE_TRR_sd$x
CATE_TRR_ave$UB <- CATE_TRR_ave$x+CATE_TRR_sd$x
CATE_TRR_ave$N <- table(TRR_mod)
colnames(CATE_TRR_ave) <- c("TRR","CATE","LB","UB","N")
TRRp <- ggplot(CATE_TRR_ave, aes(x = TRR, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Triglycerides")+
  ylab("Estimated HTEs")

# GLUR
tmpdata$GLUR <- round(tmpdata$GLUR,0)
GLUR_mod <- cut(tmpdata$GLUR, breaks = unique(quantile(tmpdata$GLUR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_GLUR_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(GLUR_mod),FUN=mean))
CATE_GLUR_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(GLUR_mod),FUN=sd))
CATE_GLUR_ave$LB <- CATE_GLUR_ave$x - CATE_GLUR_sd$x
CATE_GLUR_ave$UB <- CATE_GLUR_ave$x + CATE_GLUR_sd$x
CATE_GLUR_ave$N <- table(GLUR_mod)
colnames(CATE_GLUR_ave) <- c("GLUR","CATE","LB","UB","N")
GLURp <- ggplot(CATE_GLUR_ave, aes(x = GLUR, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Glucose")+
  ylab("Estimated HTEs")

# DBP
tmpdata$DBP <- round(tmpdata$DBP,0)
DBP_mod <- cut(tmpdata$DBP, breaks = unique(quantile(tmpdata$DBP, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_DBP_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(DBP_mod),FUN=mean))
CATE_DBP_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(DBP_mod),FUN=sd))
CATE_DBP_ave$LB <- CATE_DBP_ave$x - CATE_DBP_sd$x
CATE_DBP_ave$UB <- CATE_DBP_ave$x + CATE_DBP_sd$x
CATE_DBP_ave$N <- table(DBP_mod)
colnames(CATE_DBP_ave) <- c("DBP","CATE","LB","UB","N")
DBPp <- ggplot(CATE_DBP_ave, aes(x = DBP, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Diastolic blood pressure")+
  ylab("Estimated HTEs")

# HDL
tmpdata$HDL <- round(tmpdata$HDL,0)
HDL_mod <- cut(tmpdata$HDL, breaks = unique(quantile(tmpdata$HDL, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_HDL_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(HDL_mod),FUN=mean))
CATE_HDL_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(HDL_mod),FUN=sd))
CATE_HDL_ave$LB <- CATE_HDL_ave$x-CATE_HDL_sd$x
CATE_HDL_ave$UB <- CATE_HDL_ave$x+CATE_HDL_sd$x
CATE_HDL_ave$N <- table(HDL_mod)
colnames(CATE_HDL_ave) <- c("HDL","CATE","LB","UB","N")
HDLp <- ggplot(CATE_HDL_ave, aes(x = HDL, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("High density lipoprotein cholesterol")+
  ylab("Estimated HTEs")

# RESULT_K
RESULT_K_mod <- cut(tmpdata$RESULT_K, breaks = unique(quantile(tmpdata$RESULT_K, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_RESULT_K_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(RESULT_K_mod),FUN=mean))
CATE_RESULT_K_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(RESULT_K_mod),FUN=sd))
CATE_RESULT_K_ave$LB <- CATE_RESULT_K_ave$x-CATE_RESULT_K_sd$x
CATE_RESULT_K_ave$UB <- CATE_RESULT_K_ave$x+CATE_RESULT_K_sd$x
CATE_RESULT_K_ave$N <- table(RESULT_K_mod)
colnames(CATE_RESULT_K_ave) <- c("RESULT_K","CATE","LB","UB","N")
RESULT_Kp <- ggplot(CATE_RESULT_K_ave, aes(x = RESULT_K, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Serum Potassium Level")+
  ylab("Estimated HTEs")

# SCREAT
SCREAT_mod <- cut(tmpdata$SCREAT, breaks = unique(quantile(tmpdata$SCREAT, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_SCREAT_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SCREAT_mod),FUN=mean))
CATE_SCREAT_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SCREAT_mod),FUN=sd))
CATE_SCREAT_ave$LB <- CATE_SCREAT_ave$x - CATE_SCREAT_sd$x
CATE_SCREAT_ave$UB <- CATE_SCREAT_ave$x + CATE_SCREAT_sd$x
CATE_SCREAT_ave$N <- table(SCREAT_mod)
colnames(CATE_SCREAT_ave) <- c("SCREAT","CATE","LB","UB","N")
SCREATp <- ggplot(CATE_SCREAT_ave, aes(x = SCREAT, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Serum creatinine")+
  ylab("Estimated HTEs")

# SEATHEART
tmpdata$SEATHEART <- round(tmpdata$SEATHEART,0)
SEATHEART_mod <- cut(tmpdata$SEATHEART, breaks = unique(quantile(tmpdata$SEATHEART, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_SEATHEART_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SEATHEART_mod),FUN=mean))
CATE_SEATHEART_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SEATHEART_mod),FUN=sd))
CATE_SEATHEART_ave$LB <- CATE_SEATHEART_ave$x - CATE_SEATHEART_sd$x
CATE_SEATHEART_ave$UB <- CATE_SEATHEART_ave$x + CATE_SEATHEART_sd$x
CATE_SEATHEART_ave$N <- table(SEATHEART_mod)
colnames(CATE_SEATHEART_ave) <- c("SEATHEART","CATE","LB","UB","N")
SEATHEARTp <- ggplot(CATE_SEATHEART_ave, aes(x = SEATHEART, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Seat heart rate")+
  ylab("Estimated HTEs")

# SBP
tmpdata$SBP <- round(tmpdata$SBP,0)
SBP_mod <- cut(tmpdata$SBP, breaks = unique(quantile(tmpdata$SBP, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_SBP_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SBP_mod),FUN=mean))
CATE_SBP_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SBP_mod),FUN=sd))
CATE_SBP_ave$LB <- CATE_SBP_ave$x-CATE_SBP_sd$x
CATE_SBP_ave$UB <- CATE_SBP_ave$x+CATE_SBP_sd$x
CATE_SBP_ave$N <- table(SBP_mod)
colnames(CATE_SBP_ave) <- c("SBP","CATE","LB","UB","N")
SBPp <- ggplot(CATE_SBP_ave, aes(x = SBP, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.09,0.02)+
  xlab("Systolic Blood Pressure")+
  ylab("Estimated HTEs")

# N_AGENTS
N_AGENTS_mod <- cut(tmpdata$N_AGENTS, breaks = unique(quantile(tmpdata$N_AGENTS, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_N_AGENTS_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(N_AGENTS_mod),FUN=mean))
CATE_N_AGENTS_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(N_AGENTS_mod),FUN=sd))
CATE_N_AGENTS_ave$LB <- CATE_N_AGENTS_ave$x-CATE_N_AGENTS_sd$x
CATE_N_AGENTS_ave$UB <- CATE_N_AGENTS_ave$x+CATE_N_AGENTS_sd$x
CATE_N_AGENTS_ave$N <- table(N_AGENTS_mod)
colnames(CATE_N_AGENTS_ave) <- c("N_AGENTS","CATE","LB","UB","N")
N_AGENTSp <- ggplot(CATE_N_AGENTS_ave, aes(x = N_AGENTS, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  ylim(-0.09,0.02)+
  xlab("Number of anti-hypertensive medications")+
  ylab("Estimated HTEs")

# SUB_CLINICALCVD
SUB_CLINICALCVD <- as.factor(tmpdata$SUB_CLINICALCVD)
CATE_SUB_CLINICALCVD_ave <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SUB_CLINICALCVD),FUN=mean))
CATE_SUB_CLINICALCVD_sd <- data.frame(aggregate(sum_cvd$XL_RD_cvd,by=list(SUB_CLINICALCVD),FUN=sd))
CATE_SUB_CLINICALCVD_ave$LB <- CATE_SUB_CLINICALCVD_ave$x-CATE_SUB_CLINICALCVD_sd$x
CATE_SUB_CLINICALCVD_ave$UB <- CATE_SUB_CLINICALCVD_ave$x+CATE_SUB_CLINICALCVD_sd$x
CATE_SUB_CLINICALCVD_ave$N <- table(SUB_CLINICALCVD)
colnames(CATE_SUB_CLINICALCVD_ave) <- c("SUB_CLINICALCVD","CATE","LB","UB","N")
SUB_CLINICALCVDp <- ggplot(CATE_SUB_CLINICALCVD_ave, aes(x = SUB_CLINICALCVD, y = CATE, group = 1)) +
  theme_classic()+
  geom_pointrange(aes(ymax = UB, ymin = LB), colour="black")+
  scale_x_discrete(breaks = levels(SUB_CLINICALCVD), labels = c("No","Yes"))+
  geom_line(linetype="dashed")+
  geom_point()+
  ylim(-0.09,0.02)+
  xlab("Subclinical cardiovascular disease")+
  ylab("Estimated HTE")

filename <- paste0("./CATE_BART_CVD.png")
png(filename, width = 12, height = 20, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(UMALCRp, Agep, eGFRp, BMIp, CHRp, TRRp, DBPp, SCREATp,  SBPp,  
                               SEATHEARTp, HDLp, GLURp, RESULT_Kp, N_AGENTSp, SUB_CLINICALCVDp,
                               nrow=5, ncol=3), nrow=5, heights=c(10,1,1,1,1)))
dev.off()

# aggregate by treatment benefit 
bart_deciles_cvd <- cut(sum_cvd$XL_RD_cvd, breaks = unique(quantile(sum_cvd$XL_RD_cvd, probs=seq(0,1,1/3))), include.lowest=TRUE)
cvdlabel <- as.numeric(bart_deciles_cvd)
vars <- tmpdata[,varimp_cvd]
datt <- data.frame(vars, cvdlabel, sum_cvd$XL_RD_cvd)
datt$group <- ifelse(datt$cvdlabel==1,"Most",
                     ifelse(datt$cvdlabel==2, "Moderate", "Least"))
datt2 <- datt[,c(1:15,18)]
datt2$N_AGENTS <- as.numeric(datt2$N_AGENTS)
bart_cvd_subgrps <- datt2 %>%
  tbl_summary(
    by = group,
    type =  list(N_AGENTS ~ "continuous"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
  )%>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Net Benefit Group**") %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  add_overall()%>%       # build gtsummary table
  as_gt() %>%             # convert to gt table
  gt::gtsave(             # save table as image
    filename = "bart_cvd_subgroup.png"
  );bart_cvd_subgrps

# reformat data
ctsvars <- tmpdata[,varimp_cvd]
means <- sds <- matrix(NA, 4, dim(ctsvars)[2])
for (j in 1:dim(ctsvars)[2]){
  means[,j] <- c(mean(ctsvars[,j]), aggregate(ctsvars[,j], by=list(datt2$group), FUN=mean)$x)
  sds[,j] <- c(sd(ctsvars[,j]), aggregate(ctsvars[,j], by=list(datt2$group), FUN=sd)$x)
}
means <- as.vector(means)
sds <- as.vector(sds)
stats <- data.frame(means, sds)
stats$UB <- stats[,1] + stats[,2]
stats$LB <- stats[,1] - stats[,2]
stats$variables <- rep(varimp_cvd, each=4)
stats$Treatment_benefit <- rep(c("Overall","Least","Moderate","Most"),dim(ctsvars)[2])

p2 <- ggplot(stats, aes(x = Treatment_benefit, y=means, colour=Treatment_benefit)) +
  theme_classic()+
  facet_wrap(~variables, ncol = 3, scales = "free") +
  scale_colour_manual(values=c("blue", "green","red","black"))+
  geom_pointrange(aes(ymax = UB, ymin = LB), size=0.5)+
  geom_blank(data=stats) +
  theme(axis.title.x=element_text(size=11, vjust=-2)) +
  theme(axis.title.y=element_text(size=11, angle=90, vjust=2)) +
  theme(plot.title=element_text(size=15, vjust=3,hjust=0.5))+
  xlab("Treatment benefit")+ylab("Mean (+/- sd)")+
  theme(legend.position="bottom")+
  theme(plot.margin = unit(rep(0.5,4), "cm")); print(p2)

filename <- paste0("./bart_cvd.png")
png(filename, width = 8, height =10, units = 'in', res = 300)
print(p2)
dev.off()


## SAE graph 
# UMALCR
tmpdata$UMALCR <-  round(tmpdata$UMALCR,0)
group <- cut(tmpdata$UMALCR, breaks = unique(quantile(tmpdata$UMALCR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_UMALCR_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(group),FUN=mean))
CATE_UMALCR_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(group),FUN=sd))
CATE_UMALCR_ave$LB <- CATE_UMALCR_ave$x-CATE_UMALCR_sd$x
CATE_UMALCR_ave$UB <- CATE_UMALCR_ave$x+CATE_UMALCR_sd$x
CATE_UMALCR_ave$N <- table(group)
colnames(CATE_UMALCR_ave) <- c("UMALCR","CATE","LB","UB","N")
UMALCRp <- ggplot(CATE_UMALCR_ave, aes(x = UMALCR, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Urine Albumin/Creatinine ratio")+
  ylab("Estimated HTEs")

# Age
Age_mod <- cut(tmpdata$AGE, breaks = unique(quantile(tmpdata$AGE, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_age_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(Age_mod),FUN=mean))
CATE_age_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(Age_mod),FUN=sd))
CATE_age_ave$LB <- CATE_age_ave$x - CATE_age_sd$x
CATE_age_ave$UB <- CATE_age_ave$x + CATE_age_sd$x
CATE_age_ave$N <- table(Age_mod)
colnames(CATE_age_ave) <- c("Age","CATE","LB","UB","N")
Agep <- ggplot(CATE_age_ave, aes(x = Age, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Age")+
  ylab("Estimated HTEs")

# eGFR
tmpdata$EGFR <-  round(tmpdata$EGFR,0)
eGFR_mod <- cut(tmpdata$EGFR, breaks = unique(quantile(tmpdata$EGFR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_egfr_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(eGFR_mod),FUN=mean))
CATE_egfr_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(eGFR_mod),FUN=sd))
CATE_egfr_ave$LB <- CATE_egfr_ave$x - CATE_egfr_sd$x
CATE_egfr_ave$UB <- CATE_egfr_ave$x + CATE_egfr_sd$x
CATE_egfr_ave$N <- table(eGFR_mod)
colnames(CATE_egfr_ave) <- c("eGFR","CATE","LB","UB","N")
eGFRp <- ggplot(CATE_egfr_ave, aes(x = eGFR, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Estimated glomerular filtration rate")+
  ylab("Estimated HTEs")

# BMI
tmpdata$BMI <- round(tmpdata$BMI,0)
BMI_mod <- cut(tmpdata$BMI, breaks = unique(quantile(tmpdata$BMI, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_BMI_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(BMI_mod),FUN=mean))
CATE_BMI_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(BMI_mod),FUN=sd))
CATE_BMI_ave$LB <- CATE_BMI_ave$x-CATE_BMI_sd$x
CATE_BMI_ave$UB <- CATE_BMI_ave$x+CATE_BMI_sd$x
CATE_BMI_ave$N <- table(BMI_mod)
colnames(CATE_BMI_ave) <- c("BMI","CATE","LB","UB","N")
BMIp <- ggplot(CATE_BMI_ave, aes(x = BMI, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Boday mass index")+
  ylab("Estimated HTEs")

# CHR
tmpdata$CHR <- round(tmpdata$CHR,0)
CHR_mod <- cut(tmpdata$CHR, breaks = unique(quantile(tmpdata$CHR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_CHR_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(CHR_mod),FUN=mean))
CATE_CHR_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(CHR_mod),FUN=sd))
CATE_CHR_ave$LB <- CATE_CHR_ave$x-CATE_CHR_sd$x
CATE_CHR_ave$UB <- CATE_CHR_ave$x+CATE_CHR_sd$x
CATE_CHR_ave$N <- table(CHR_mod)
colnames(CATE_CHR_ave) <- c("CHR","CATE","LB","UB","N")
CHRp <- ggplot(CATE_CHR_ave, aes(x = CHR, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Total cholesterol")+
  ylab("Estimated HTEs")

# TRR
tmpdata$TRR <- round(tmpdata$TRR,0)
TRR_mod <- cut(tmpdata$TRR, breaks = unique(quantile(tmpdata$TRR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_TRR_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(TRR_mod),FUN=mean))
CATE_TRR_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(TRR_mod),FUN=sd))
CATE_TRR_ave$LB <- CATE_TRR_ave$x-CATE_TRR_sd$x
CATE_TRR_ave$UB <- CATE_TRR_ave$x+CATE_TRR_sd$x
CATE_TRR_ave$N <- table(TRR_mod)
colnames(CATE_TRR_ave) <- c("TRR","CATE","LB","UB","N")
TRRp <- ggplot(CATE_TRR_ave, aes(x = TRR, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Triglycerides")+
  ylab("Estimated HTEs")

# GLUR
tmpdata$GLUR <- round(tmpdata$GLUR,0)
GLUR_mod <- cut(tmpdata$GLUR, breaks = unique(quantile(tmpdata$GLUR, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_GLUR_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(GLUR_mod),FUN=mean))
CATE_GLUR_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(GLUR_mod),FUN=sd))
CATE_GLUR_ave$LB <- CATE_GLUR_ave$x - CATE_GLUR_sd$x
CATE_GLUR_ave$UB <- CATE_GLUR_ave$x + CATE_GLUR_sd$x
CATE_GLUR_ave$N <- table(GLUR_mod)
colnames(CATE_GLUR_ave) <- c("GLUR","CATE","LB","UB","N")
GLURp <- ggplot(CATE_GLUR_ave, aes(x = GLUR, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Glucose")+
  ylab("Estimated HTEs")

# DBP
tmpdata$DBP <- round(tmpdata$DBP,0)
DBP_mod <- cut(tmpdata$DBP, breaks = unique(quantile(tmpdata$DBP, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_DBP_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(DBP_mod),FUN=mean))
CATE_DBP_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(DBP_mod),FUN=sd))
CATE_DBP_ave$LB <- CATE_DBP_ave$x - CATE_DBP_sd$x
CATE_DBP_ave$UB <- CATE_DBP_ave$x + CATE_DBP_sd$x
CATE_DBP_ave$N <- table(DBP_mod)
colnames(CATE_DBP_ave) <- c("DBP","CATE","LB","UB","N")
DBPp <- ggplot(CATE_DBP_ave, aes(x = DBP, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Diastolic blood pressure")+
  ylab("Estimated HTEs")

# HDL
tmpdata$HDL <- round(tmpdata$HDL,0)
HDL_mod <- cut(tmpdata$HDL, breaks = unique(quantile(tmpdata$HDL, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_HDL_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(HDL_mod),FUN=mean))
CATE_HDL_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(HDL_mod),FUN=sd))
CATE_HDL_ave$LB <- CATE_HDL_ave$x-CATE_HDL_sd$x
CATE_HDL_ave$UB <- CATE_HDL_ave$x+CATE_HDL_sd$x
CATE_HDL_ave$N <- table(HDL_mod)
colnames(CATE_HDL_ave) <- c("HDL","CATE","LB","UB","N")
HDLp <- ggplot(CATE_HDL_ave, aes(x = HDL, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("High density lipoprotein cholesterol")+
  ylab("Estimated HTEs")

# RESULT_K
RESULT_K_mod <- cut(tmpdata$RESULT_K, breaks = unique(quantile(tmpdata$RESULT_K, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_RESULT_K_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(RESULT_K_mod),FUN=mean))
CATE_RESULT_K_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(RESULT_K_mod),FUN=sd))
CATE_RESULT_K_ave$LB <- CATE_RESULT_K_ave$x-CATE_RESULT_K_sd$x
CATE_RESULT_K_ave$UB <- CATE_RESULT_K_ave$x+CATE_RESULT_K_sd$x
CATE_RESULT_K_ave$N <- table(RESULT_K_mod)
colnames(CATE_RESULT_K_ave) <- c("RESULT_K","CATE","LB","UB","N")
RESULT_Kp <- ggplot(CATE_RESULT_K_ave, aes(x = RESULT_K, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Serum Potassium Level")+
  ylab("Estimated HTEs")

# SCREAT
SCREAT_mod <- cut(tmpdata$SCREAT, breaks = unique(quantile(tmpdata$SCREAT, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_SCREAT_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(SCREAT_mod),FUN=mean))
CATE_SCREAT_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(SCREAT_mod),FUN=sd))
CATE_SCREAT_ave$LB <- CATE_SCREAT_ave$x - CATE_SCREAT_sd$x
CATE_SCREAT_ave$UB <- CATE_SCREAT_ave$x + CATE_SCREAT_sd$x
CATE_SCREAT_ave$N <- table(SCREAT_mod)
colnames(CATE_SCREAT_ave) <- c("SCREAT","CATE","LB","UB","N")
SCREATp <- ggplot(CATE_SCREAT_ave, aes(x = SCREAT, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Serum creatinine")+
  ylab("Estimated HTEs")

# SEATHEART
tmpdata$SEATHEART <- round(tmpdata$SEATHEART,0)
SEATHEART_mod <- cut(tmpdata$SEATHEART, breaks = unique(quantile(tmpdata$SEATHEART, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_SEATHEART_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(SEATHEART_mod),FUN=mean))
CATE_SEATHEART_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(SEATHEART_mod),FUN=sd))
CATE_SEATHEART_ave$LB <- CATE_SEATHEART_ave$x - CATE_SEATHEART_sd$x
CATE_SEATHEART_ave$UB <- CATE_SEATHEART_ave$x + CATE_SEATHEART_sd$x
CATE_SEATHEART_ave$N <- table(SEATHEART_mod)
colnames(CATE_SEATHEART_ave) <- c("SEATHEART","CATE","LB","UB","N")
SEATHEARTp <- ggplot(CATE_SEATHEART_ave, aes(x = SEATHEART, y = CATE, group=1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Seat heart rate")+
  ylab("Estimated HTEs")

# SBP
tmpdata$SBP <- round(tmpdata$SBP,0)
SBP_mod <- cut(tmpdata$SBP, breaks = unique(quantile(tmpdata$SBP, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_SBP_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(SBP_mod),FUN=mean))
CATE_SBP_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(SBP_mod),FUN=sd))
CATE_SBP_ave$LB <- CATE_SBP_ave$x-CATE_SBP_sd$x
CATE_SBP_ave$UB <- CATE_SBP_ave$x+CATE_SBP_sd$x
CATE_SBP_ave$N <- table(SBP_mod)
colnames(CATE_SBP_ave) <- c("SBP","CATE","LB","UB","N")
SBPp <- ggplot(CATE_SBP_ave, aes(x = SBP, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(-0.03, 0.12)+
  xlab("Systolic Blood Pressure")+
  ylab("Estimated HTEs")

# N_AGENTS
N_AGENTS_mod <- cut(tmpdata$N_AGENTS, breaks = unique(quantile(tmpdata$N_AGENTS, probs=seq(0,1,1/10))), include.lowest=TRUE)
CATE_N_AGENTS_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(N_AGENTS_mod),FUN=mean))
CATE_N_AGENTS_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(N_AGENTS_mod),FUN=sd))
CATE_N_AGENTS_ave$LB <- CATE_N_AGENTS_ave$x-CATE_N_AGENTS_sd$x
CATE_N_AGENTS_ave$UB <- CATE_N_AGENTS_ave$x+CATE_N_AGENTS_sd$x
CATE_N_AGENTS_ave$N <- table(N_AGENTS_mod)
colnames(CATE_N_AGENTS_ave) <- c("N_AGENTS","CATE","LB","UB","N")
N_AGENTSp <- ggplot(CATE_N_AGENTS_ave, aes(x = N_AGENTS, y = CATE, group = 1)) +
  theme_classic()+
  geom_ribbon(aes(ymin = LB, ymax = UB), size=1, linetype="dashed", color="black", fill="grey")+
  geom_line()+
  geom_point()+
  ylim(-0.03, 0.12)+
  xlab("Number of anti-hypertensive medications")+
  ylab("Estimated HTEs")

# LIVEWITHOTHERS
LIVEWITHOTHERS <- as.factor(tmpdata$LIVEWITHOTHERS)
CATE_LIVEWITHOTHERS_ave <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(LIVEWITHOTHERS),FUN=mean))
CATE_LIVEWITHOTHERS_sd <- data.frame(aggregate(sum_sae$XL_RD_sae,by=list(LIVEWITHOTHERS),FUN=sd))
CATE_LIVEWITHOTHERS_ave$LB <- CATE_LIVEWITHOTHERS_ave$x-CATE_LIVEWITHOTHERS_sd$x
CATE_LIVEWITHOTHERS_ave$UB <- CATE_LIVEWITHOTHERS_ave$x+CATE_LIVEWITHOTHERS_sd$x
CATE_LIVEWITHOTHERS_ave$N <- table(LIVEWITHOTHERS)
colnames(CATE_LIVEWITHOTHERS_ave) <- c("LIVEWITHOTHERS","CATE","LB","UB","N")
LIVEWITHOTHERSp <- ggplot(CATE_LIVEWITHOTHERS_ave, aes(x = LIVEWITHOTHERS, y = CATE, group = 1)) +
  theme_classic()+
  geom_pointrange(aes(ymax = UB, ymin = LB), colour="black")+
  scale_x_discrete(breaks = levels(LIVEWITHOTHERS), labels = c("No","Yes"))+
  geom_line(linetype="dashed")+
  geom_point()+
  ylim(-0.03, 0.12)+
  xlab("Whether or not live with others")+
  ylab("Estimated HTEs")

filename <- paste0("./CATE_BART_SAE.png")
png(filename, width = 12, height = 20, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(UMALCRp, eGFRp, CHRp, BMIp, SCREATp, Agep, DBPp, TRRp, HDLp, SEATHEARTp,         
                               SBPp, GLURp, RESULT_Kp, N_AGENTSp, LIVEWITHOTHERSp, nrow=5, ncol=3), nrow=5, heights=c(10,1,1,1,1)))
dev.off()

# aggregate by treatment benefit 
bart_deciles_sae <- cut(sum_sae$XL_RD_sae, breaks = unique(quantile(sum_sae$XL_RD_sae, probs=seq(0,1,1/3))), include.lowest=TRUE)
saelabel <- as.numeric(bart_deciles_sae)
vars <- tmpdata[,varimp_sae]
datt <- data.frame(vars, saelabel, sum_sae$XL_RD_sae)
datt$group <- ifelse(datt$saelabel==1,"Least",
                     ifelse(datt$saelabel==2, "Moderate", "Most"))
datt2 <- datt[,c(1:15,18)]
datt2$N_AGENTS <- as.numeric(datt2$N_AGENTS)
bart_sae_subgrps <- datt2 %>%
  tbl_summary(
    by = group,
    type =  list(N_AGENTS ~ "continuous"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
  )%>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Net Benefit Group**") %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  add_overall()%>%       # build gtsummary table
  as_gt() %>%             # convert to gt table
  gt::gtsave(             # save table as image
    filename = "bart_sae_subgroup.png"
  );bart_sae_subgrps

# reformat data
ctsvars <- tmpdata[,varimp_sae]
means <- sds <- matrix(NA, 4, dim(ctsvars)[2])
for (j in 1:dim(ctsvars)[2]){
  means[,j] <- c(mean(ctsvars[,j]), aggregate(ctsvars[,j], by=list(datt2$group), FUN=mean)$x)
  sds[,j] <- c(sd(ctsvars[,j]), aggregate(ctsvars[,j], by=list(datt2$group), FUN=sd)$x)
}
means <- as.vector(means)
sds <- as.vector(sds)
stats <- data.frame(means, sds)
stats$UB <- stats[,1] + stats[,2]
stats$LB <- stats[,1] - stats[,2]
stats$variables <- rep(varimp_sae, each=4)
stats$Treatment_harms <- rep(c("Overall","Least","Moderate","Most"),dim(ctsvars)[2])

p2 <- ggplot(stats, aes(x = Treatment_harms, y=means, colour=Treatment_harms)) +
  theme_classic()+
  facet_wrap(~variables, ncol = 3, scales = "free") +
  scale_colour_manual(values=c("blue", "green","red","black"))+
  geom_pointrange(aes(ymax = UB, ymin = LB), size=0.5)+
  geom_blank(data=stats) +
  theme(axis.title.x=element_text(size=11, vjust=-2)) +
  theme(axis.title.y=element_text(size=11, angle=90, vjust=2)) +
  theme(plot.title=element_text(size=15, vjust=3,hjust=0.5))+
  xlab("Treatment harms")+ylab("Mean (+/- sd)")+
  theme(plot.margin = unit(rep(0.5,4), "cm"))+
  theme(legend.position="bottom"); print(p2)

filename <- paste0("./bart_sae.png")
png(filename, width = 8, height =10, units = 'in', res = 300)
print(p2)
dev.off()

##-------------------------------------------------- Net benefit ---------------------------------------------------------------##
varimp_NB <- c(intersect(varimp_cvd, varimp_sae),"FEMALE","RACE_BLACK","neversmoker","formersmoker","currentsmoker",
               "SUB_CLINICALCVD", "STATIN", "ASPIRIN")

bart_deciles_cvd <- cut(sum_cvd$XL_RD_cvd, breaks = unique(quantile(sum_cvd$XL_RD_cvd, probs=seq(0,1,1/3))), include.lowest=TRUE)
bart_deciles_sae <- cut(sum_sae$XL_RD_sae, breaks = unique(quantile(sum_sae$XL_RD_sae, probs=seq(0,1,1/3))), include.lowest=TRUE)
cvdlabel <- as.numeric(bart_deciles_cvd)
saelabel <- as.numeric(bart_deciles_sae)
table(saelabel, cvdlabel)

vars <- tmpdata[,varimp_NB]
datt <- data.frame(vars, cvdlabel, saelabel, sum_cvd$XL_RD_cvd, sum_sae$XL_RD_sae)
datt$group <- ifelse(datt$cvdlabel==1 & datt$saelabel==1,"Most",
                     ifelse(datt$cvdlabel==3 & datt$saelabel==3, "Least", "Moderate"))
datt2 <- datt[,c(1:22,27)]
datt2$N_AGENTS <- as.numeric(datt2$N_AGENTS)
bart_subgrps <- datt2 %>%
  tbl_summary(
    by = group,
    type =  list(N_AGENTS ~ "continuous"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
  )%>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Net Benefit Group**") %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  add_overall()%>%       # build gtsummary table
  as_gt() %>%             # convert to gt table
  gt::gtsave(             # save table as image
    filename = "bart_subgroup.png"
  );bart_subgrps
   
# reformat data
femaledat <- datt2[,c("FEMALE","group")]; colnames(femaledat)[1] <- "variables"
blackdat <- datt2[,c("RACE_BLACK","group")]; colnames(blackdat)[1] <- "variables"
cat_dat <- rbind(femaledat, blackdat)

catvars <- tmpdata[,c("FEMALE","RACE_BLACK","neversmoker","formersmoker","currentsmoker",
                      "SUB_CLINICALCVD", "STATIN", "ASPIRIN")]
percen <- matrix(NA, 4, dim(catvars)[2])
for (i in 1:dim(catvars)[2]){
  percen[,i] <- c(round(table(catvars[,i])[2]/dim(datt2)[1],2)*100,
                  round(table(catvars[,i], datt2$group)[2,]/table(datt2$group),2)*100)
}

variables <- rep(c("FEMALE","RACE_BLACK","neversmoker","formersmoker","currentsmoker",
                   "SUB_CLINICALCVD", "STATIN", "ASPIRIN"), each=4)
Net_benefit <- rep(c("Overall","Least","Moderate","Most"),dim(catvars)[2])
Percentage <- as.vector(percen)
dataformat <- data.frame(Percentage, variables, Net_benefit)
p1 <- ggplot(dataformat, aes(x = variables, y = Percentage, group = Net_benefit, fill = Net_benefit))+
  theme_classic()+
  scale_fill_manual(values=c("blue", "green","red","black"))+
  geom_bar(stat = "identity", width = 0.5, position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_text(size=11, vjust=-2)) +
  theme(axis.title.y=element_text(size=11, angle=90, vjust=2)) +
  theme(plot.title=element_text(size=15, vjust=3,hjust=0.5))+
  theme(legend.position="none")+
  xlab("Variables")+ylab("Percentage (%)");print(p1)


ctsvarnames <- intersect(varimp_cvd, varimp_sae)
ctsvars <- tmpdata[,ctsvarnames]
means <- sds <- matrix(NA, 4, dim(ctsvars)[2])
for (j in 1:dim(ctsvars)[2]){
  means[,j] <- c(mean(ctsvars[,j]), aggregate(ctsvars[,j], by=list(datt2$group), FUN=mean)$x)
  sds[,j] <- c(sd(ctsvars[,j]), aggregate(ctsvars[,j], by=list(datt2$group), FUN=sd)$x)
}
means <- as.vector(means)
sds <- as.vector(sds)
stats <- data.frame(means, sds)
stats$UB <- stats[,1] + stats[,2]
stats$LB <- stats[,1] - stats[,2]
stats$variables <- rep(ctsvarnames, each=4)
stats$Net_benefit <- rep(c("Overall","Least","Moderate","Most"),dim(ctsvars)[2])

p2 <- ggplot(stats, aes(x = Net_benefit, y=means, colour=Net_benefit)) +
  theme_classic()+
  facet_wrap(~variables, ncol = 4, scales = "free") +
  scale_colour_manual(values=c("blue", "green","red","black"))+
  geom_pointrange(aes(ymax = UB, ymin = LB), size=0.5)+
  geom_blank(data=stats) +
  theme(axis.title.x=element_text(size=11, vjust=-2)) +
  theme(axis.title.y=element_text(size=11, angle=90, vjust=2)) +
  theme(plot.title=element_text(size=15, vjust=3,hjust=0.5))+
  ylab("Mean (+/- sd)")+
  theme(legend.position="bottom")+
  theme(plot.margin = unit(rep(0.5,4), "cm")); print(p2)

filename <- paste0("./bart_NB.png")
png(filename, width = 10, height =13, units = 'in', res = 300)
print(grid.arrange(p1, p2, ncol = 1, heights=c(2,3)))
dev.off()

