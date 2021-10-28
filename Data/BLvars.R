###--------------------------------------------------------------------------###
###-----------  Extract more baseline covariates for SPRINT data ------------###
###-----------             Crystal Xu    02/20/2021              ------------###
###--------------------------------------------------------------------------###
 
library(sas7bdat); library(dataMaid); library(matrixStats)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/SPRINT_2019a/SPRINT/data")

# from labs
labs <- read.sas7bdat("labs.sas7bdat")
labs <- labs[labs$VISITCODE=="RZ1",]
labs_vars <- labs[, c("MASKID", "RESULT_K","RESULT_UMALI","RESULT_CREATR")]
labs_vars$UMALCR <- labs_vars$RESULT_UMALI/labs_vars$RESULT_CREATR  # fewer missing values if we compute it this way 
labs_vars <- labs_vars[,c("MASKID", "RESULT_K","UMALCR")]

# from incl_excl
incl_excl <- read.sas7bdat("incl_excl.sas7bdat"); dim(incl_excl)
incl_excl_vars <- incl_excl[, c("MASKID", "CORONARYREVAS")]

# from bp_manage_base
bp_manage_base <- read.sas7bdat("bp_manage_base.sas7bdat")
bp_manage_base_vars <- bp_manage_base[, c("MASKID", "SEATHEART")]
bp_manage_base_vars$SEATHEART <- as.numeric(as.character(bp_manage_base_vars$SEATHEART))

# from ecg
ecg <- read.sas7bdat("ecg.sas7bdat")
ecg_unique <- ecg[!duplicated(ecg$MASKID),]
ecg_vars <- ecg_unique[, c("MASKID", "LVHANY3")]

# from bl_history
bl_history <- read.sas7bdat("bl_history.sas7bdat"); dim(bl_history)
bl_history_vars <- bl_history[, c("MASKID","ULCER","LIVEWITHOTHERS","DEPRESS","ALCOHOL","UNINSURED")]


##---------------------------------- Medication -------------------------------------------------##
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data")
bl_meds_vars_summary <- read.csv("./bl_meds_vars_summary.csv")  # no missing value here
bl_meds_vars_summary <- bl_meds_vars_summary[, !names(bl_meds_vars_summary) %in% c("NSAID")]

##----------------------------------------- Existing Baseline Vars -------------------------------------------------##
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/SPRINT_2019a/SPRINT-POP/data")
baseline <- read.sas7bdat("baseline.sas7bdat")

# Derive 3 baseline covariates based on TONY's code below 
baseline$Hispanic <- (baseline$RACE4 == "HISPANIC")
baseline$currentsmoker <- (baseline$SMOKE_3CAT == 3)
baseline$formersmoker <- (baseline$SMOKE_3CAT == 2)
baseline$neversmoker <- (baseline$SMOKE_3CAT == 1)

# remove unnecessary variables
baseline <- subset(baseline, select = -c(NEWSITEID,RISK10YRS,INCLUSIONFRS,SMOKE_3CAT,
                                         RACE4,SUB_CVD,SUB_SENIOR,SBPTERTILE,
                                         NOAGENTS,SUB_SUBCLINICALCVD,UMALCR))

#baseline.comp <- baseline[complete.cases(baseline),]; dim(baseline.comp) # 8746 33

# Merge the additional 22 variables above - extracted by Crystal 
baseline1 <- merge(baseline, bl_meds_vars_summary, by = "MASKID", all.x = T); dim(baseline1) 
baseline2 <- merge(baseline1, labs_vars, by = "MASKID", all.x = T); dim(baseline2)
baseline3 <- merge(baseline2, incl_excl_vars, by = "MASKID", all.x = T); dim(baseline3)
baseline4 <- merge(baseline3, bp_manage_base_vars, by = "MASKID", all.x = T); dim(baseline4)
baseline5 <- merge(baseline4, bl_history_vars, by = "MASKID", all.x = T); dim(baseline5) # 9361 34
#baseline6 <- merge(baseline5, ecg_vars, by = "MASKID", all.x = T); dim(baseline6) 
baseline5$insurance <- ifelse(is.na(baseline5$UNINSURED)==T,1,0)
baseline5 <- baseline5[,!names(baseline5)%in%c("UNINSURED")]

baseline1.comp <- baseline1[complete.cases(baseline1),]; dim(baseline1.comp) # 9171 31  2%
baseline2.comp <- baseline2[complete.cases(baseline2),]; dim(baseline2.comp) # 9079 34  
baseline3.comp <- baseline3[complete.cases(baseline3),]; dim(baseline3.comp) # 9079 35  
baseline4.comp <- baseline4[complete.cases(baseline4),]; dim(baseline4.comp) # 9064 37  
baseline5.comp <- baseline5[complete.cases(baseline5),]; dim(baseline5.comp) # 9067 39  4%  # we take "baseline5" as the final sample
#baseline6.comp <- baseline6[complete.cases(baseline6),]; dim(baseline6.comp) # 8750 55  6%

write.csv(baseline5, "./expanded_baseline_vars.csv", row.names = F)


##---------------------------- Derive CVD and SAE outcomes based on TONY's code ----------------------------##
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/SPRINT_2019a/SPRINT-POP/data") 
outcomes = read.sas7bdat("outcomes.sas7bdat")
safety = read.sas7bdat("safety.sas7bdat")

cvd = (outcomes$EVENT_MI == 1) | (outcomes$EVENT_STROKE == 1) | (outcomes$EVENT_CVDDEATH == 1) | (outcomes$EVENT_HF == 1)
t_censor = rowMaxs(cbind(outcomes$T_MI,outcomes$T_STROKE,outcomes$T_CVDDEATH,outcomes$T_HF))
t_cvds = rowMaxs(cbind(outcomes$T_MI*outcomes$EVENT_MI,
                       outcomes$T_STROKE*outcomes$EVENT_STROKE,
                       outcomes$T_CVDDEATH*outcomes$EVENT_CVDDEATH,
                       outcomes$T_HF*outcomes$EVENT_HF))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
mainY = data.frame(outcomes$MASKID, t_cvds, cvd); names(mainY)[1] <- "MASKID"

sae = (safety$HYP_SAE_EVNT==1)|(safety$SYN_SAE_EVNT==1)|(safety$ELE_SAE_EVNT==1)|(safety$AKI_SAE_EVNT==1)|(safety$BRA_SAE_EVNT==1)
t_censor = rowMaxs(cbind(safety$HYP_SAE_DAYS,safety$SYN_SAE_DAYS,safety$ELE_SAE_DAYS,safety$AKI_SAE_DAYS,safety$BRA_SAE_DAYS))
t_saes = rowMaxs(cbind(safety$HYP_SAE_DAYS*safety$HYP_SAE_EVNT,
                       safety$SYN_SAE_DAYS*safety$SYN_SAE_EVNT,
                       safety$ELE_SAE_DAYS*safety$ELE_SAE_EVNT,
                       safety$AKI_SAE_DAYS*safety$AKI_SAE_EVNT,
                       safety$BRA_SAE_DAYS*safety$BRA_SAE_EVNT))
t_saes[t_saes==0] = t_censor[t_saes==0]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)
saeY = data.frame(safety$MASKID, t_saes, sae); names(saeY)[1] <- "MASKID"

outcome = merge(mainY, saeY, by = "MASKID", all=T); dim(outcome)
outcome$cvd <- as.numeric(outcome$cvd)
outcome$sae <- as.numeric(outcome$sae)
write.csv(outcome,"C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data/outcome.csv", row.names = F)



##----------------------------------- Analysis Data Set  -------------------------------##
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analysis_Data")
blvars <- read.csv("expanded_baseline_vars.csv");dim(blvars) # 8989 54 
blvars$Hispanic <- as.numeric(blvars$Hispanic)
blvars$currentsmoker <- as.numeric(blvars$currentsmoker)
blvars$formersmoker <- as.numeric(blvars$formersmoker)
blvars$neversmoker <- as.numeric(blvars$neversmoker)
#blvars$NSAID <- as.numeric(blvars$NSAID)
blvars$ACE_ARB <- as.numeric(blvars$ACE_ARB)
blvars$CCB <- as.numeric(blvars$CCB)
blvars$ThiazideDiuretic <- as.numeric(blvars$ThiazideDiuretic)
blvars$LoopDiuretic <- as.numeric(blvars$LoopDiuretic)
blvars$BetaBlocker <- as.numeric(blvars$BetaBlocker)
blvars$AlphaOneBlocker <- as.numeric(blvars$AlphaOneBlocker)
blvars$OtherAntiHyperMed <- as.numeric(blvars$OtherAntiHyperMed) 
write.csv(blvars,"blvars.csv", row.names = F)


#---------- create full codebook for expanded_baseline_vars data set -----------------#
attr(expanded_baseline_vars$MASKID, "shortDescription") <- "Masked unique identifier"
attr(expanded_baseline_vars$INTENSIVE, "shortDescription") <- "Assigned to intensive BP arm"
attr(expanded_baseline_vars$NEWSITEID, "shortDescription") <- "Randomization site"
attr(expanded_baseline_vars$RISK10YRS, "shortDescription") <- "Derived: Framingham estimation of 10-year CVD risk"
attr(expanded_baseline_vars$INCLUSIONFRS, "shortDescription") <- "Derived: Framingham 10-year CVD risk > 15%"
attr(expanded_baseline_vars$SBP, "shortDescription") <- "Derived: Seated systolic blood pressure, mm Hg"
attr(expanded_baseline_vars$DBP, "shortDescription") <- "Derived: Seated diastolic blood pressure, mm Hg"
attr(expanded_baseline_vars$N_AGENTS, "shortDescription") <- "Derived: Number of anti-hypertensive medication prescribed"
attr(expanded_baseline_vars$NOAGENTS, "shortDescription") <- "Derived: participants on no anti-hypertensive agents"
attr(expanded_baseline_vars$SMOKE_3CAT, "shortDescription") <- "Derived: Baseline smoking status"
attr(expanded_baseline_vars$ASPIRIN, "shortDescription") <- "BSL Hist: Daily Aspirin Use"
attr(expanded_baseline_vars$EGFR, "shortDescription") <- "Lab: eGFR MDRD (ml/min/1.73m^2)"
attr(expanded_baseline_vars$SCREAT, "shortDescription") <- "Lab: serum creatinine, mh/dL"
attr(expanded_baseline_vars$SUB_CKD, "shortDescription") <- "Derived: Subgroup with CKD (eGFR<60)"
attr(expanded_baseline_vars$RACE_BLACK, "shortDescription") <- "Incl/Excel: Black, African-American"
attr(expanded_baseline_vars$AGE, "shortDescription") <- "Derived: Age at randomization top-coded at 90 years"
attr(expanded_baseline_vars$FEMALE, "shortDescription") <- "Derived: Female gender"
attr(expanded_baseline_vars$SUB_CVD, "shortDescription") <- "Derived: subgroup with history of clinical/subclinical CVD"
attr(expanded_baseline_vars$SUB_CLINICALCVD, "shortDescription") <- "Derived: subgroup with history of clinical CVD"
attr(expanded_baseline_vars$SUB_SUBCLINICALCVD, "shortDescription") <- "Derived: subgroup with history of subclinical CVD"
attr(expanded_baseline_vars$SUB_SENIOR, "shortDescription") <- "Derived: subgroup >=75 years old at randomization"
attr(expanded_baseline_vars$RACE4, "shortDescription") <- "Derived: Four-level race variable (character)"
attr(expanded_baseline_vars$CHR, "shortDescription") <- "Lab: Cholesterol, mg/dL"
attr(expanded_baseline_vars$GLUR, "shortDescription") <- "Lab: Glucose, mg/dL"
attr(expanded_baseline_vars$HDL, "shortDescription") <- "Lab: HDL-cholesterol direct, mg/dL"
attr(expanded_baseline_vars$TRR, "shortDescription") <- "Lab: Triglycerides, mg/dL"
attr(expanded_baseline_vars$UMALCR, "shortDescription") <- "Lab: Urine Albumin/Creatinine ratio - mg Urine Alb / (g Creat * 0.01), mg/g Cr"
attr(expanded_baseline_vars$BMI, "shortDescription") <- "Derived: body mass index (kg/m^2)"
attr(expanded_baseline_vars$STATIN, "shortDescription") <- "Derived: on any statin"
attr(expanded_baseline_vars$SBPTERTILE, "shortDescription") <- "Derived: Systolic BP tertile"

attr(expanded_baseline_vars$Hispanic, "shortDescription") <- "Derived: Whether or not is Hispanic race (RACE4)"
attr(expanded_baseline_vars$currentsmoker, "shortDescription") <- "Derived: Whether or not is current smoker (SMOKE_3CAT)"
attr(expanded_baseline_vars$formersmoker, "shortDescription") <- "Derived: Whether or not is former smoker (SMOKE_3CAT)"

attr(expanded_baseline_vars$NSAID, "shortDescription") <- "Derived: Using NSAID (y/n)"
attr(expanded_baseline_vars$ACE_ARB, "shortDescription") <- "Derived: Using angiotensin converting enzyme inhibitors (ACEI) or angiotensin-receptor blockers (ARB) (y/n)"
attr(expanded_baseline_vars$CCB, "shortDescription") <- "Derived: Using calcium channel blockers (y/n)"
attr(expanded_baseline_vars$ThiazideDiuretic, "shortDescription") <- "Derived: Using Thiazide diuretic (y/n)"
attr(expanded_baseline_vars$LoopDiuretic, "shortDescription") <- "Derived: Using LOOP diuretic (y/n)"
attr(expanded_baseline_vars$BetaBlocker, "shortDescription") <- "Derived: Using Beta-blocker (y/n)"
attr(expanded_baseline_vars$AlphaOneBlocker, "shortDescription") <- "Derived: Using Alpha-blocker (y/n)"
attr(expanded_baseline_vars$OtherAntiHyperMed, "shortDescription") <- "Derived: Using other antihypertensive medication (y/n)"
attr(expanded_baseline_vars$RESULT_K, "shortDescription") <- "Lab: Potassium, mmol/L"
attr(expanded_baseline_vars$RESULT_NA, "shortDescription") <- "Lab: Sodium, mmol/L"
attr(expanded_baseline_vars$CORONARYREVAS, "shortDescription") <- "Baseline Coronary Revascularization(CABG,PCI)"
attr(expanded_baseline_vars$DIZZY, "shortDescription") <- "Baseline Dizziness or light headed feeling when standing"
attr(expanded_baseline_vars$SEATHEART, "shortDescription") <- "Baseline Average of 3 seated heart rate measurements"
attr(expanded_baseline_vars$PVD, "shortDescription") <- "Baseline peripheral vascular disease"
attr(expanded_baseline_vars$PTSD, "shortDescription") <- "Baseline Post-Traumatic Stress Disorder"
attr(expanded_baseline_vars$RHEARTHRITIS, "shortDescription") <- "Baseline Rheumatoid Arthritis"
attr(expanded_baseline_vars$THYROIDDIS, "shortDescription") <- "Baseline Thyroid disease"
attr(expanded_baseline_vars$ULCER, "shortDescription") <- "Baseline gastric or peptic ulcer"
attr(expanded_baseline_vars$GOUT, "shortDescription") <- "Baseline gout"
attr(expanded_baseline_vars$HIPPROB, "shortDescription") <- "Baseline hip problems"
attr(expanded_baseline_vars$OSTEOARTHRITIS, "shortDescription") <- "Baseline Osteoarthritis or degenerative arthritis"
attr(expanded_baseline_vars$BPH, "shortDescription") <- "Derived: Four-level race variable (character)"
attr(expanded_baseline_vars$LOWBKPAIN, "shortDescription") <- "Baseline low back pain"
attr(expanded_baseline_vars$ATRIALFIB, "shortDescription") <- "Baseline Atrial Fibrillation/Atrial Flutter"
attr(expanded_baseline_vars$LIVEWITHOTHERS, "shortDescription") <- "Live with other adults"
attr(expanded_baseline_vars$DEPRESS, "shortDescription") <- "Depression at baseline"
attr(expanded_baseline_vars$ANXIETY, "shortDescription") <- "Anxiety or Panic Disorder at baseline"
attr(expanded_baseline_vars$ANEMIA, "shortDescription") <- "Anemia or low blood count at baseline"
attr(expanded_baseline_vars$ALCOHOL, "shortDescription") <- "Alcohol abuse at baseline"
attr(expanded_baseline_vars$CANCER, "shortDescription") <- "History of cancer (not including skin cancer unless Melanoma)"
makeCodebook(expanded_baseline_vars, replace = T)







####################################
# TONY'S ADDITION TO SANJAY'S CODE #
####################################
# attach(sprint_set)
# 
# hisp = (RACE4 == "HISPANIC")
# currentsmoker = (SMOKE_3CAT == 3)
# formersmoker = (SMOKE_3CAT == 2)
# cvd = (EVENT_MI == 1) | (EVENT_STROKE == 1) | (EVENT_CVDDEATH == 1) | (EVENT_HF == 1)
# t_censor = rowMaxs(cbind(T_MI,T_STROKE,T_CVDDEATH,T_HF))
# t_cvds = rowMaxs(cbind(T_MI*EVENT_MI,T_STROKE*EVENT_STROKE,T_CVDDEATH*EVENT_CVDDEATH,T_HF*EVENT_HF))
# t_cvds[t_cvds==0] = t_censor[t_cvds==0]
# t_cvds[t_cvds==0] = 'NA'
# t_cvds = as.numeric(t_cvds)
# cOutcome = Surv(time=t_cvds, event = cvd)
# sae = (HYP_SAE_EVNT==1)|(SYN_SAE_EVNT==1)|(ELE_SAE_EVNT==1)|(AKI_SAE_EVNT==1)|(BRA_SAE_EVNT==1)
# t_censor = rowMaxs(cbind(HYP_SAE_DAYS,SYN_SAE_DAYS,ELE_SAE_DAYS,AKI_SAE_DAYS,BRA_SAE_DAYS))
# t_saes = rowMaxs(cbind(HYP_SAE_DAYS*HYP_SAE_EVNT,SYN_SAE_DAYS*SYN_SAE_EVNT,ELE_SAE_DAYS*ELE_SAE_EVNT,AKI_SAE_DAYS*AKI_SAE_EVNT,BRA_SAE_DAYS*BRA_SAE_EVNT))
# t_saes[t_saes==0] = t_censor[t_saes==0]
# t_saes[t_saes==0] = 'NA'
# t_saes = as.numeric(t_saes)
# dOutcome = Surv(time=t_saes, event = sae)
