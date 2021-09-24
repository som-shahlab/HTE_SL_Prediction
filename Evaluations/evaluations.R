##-------------------------------- Functions for evaluation metrics -----------------------------##
library(tidyverse)
library(dplyr)
library(Hmisc)
library(survival)
library(grf)
library(pec)
library(gbm) 
library(BART)
library(nricens)
library(MatchIt)
library(glmnet)

# ARR summary table in tertiles
quan_sum_table <- function(data, times, grp=3){
  
  # Compute the difference b/w expected ARR and obs ARR in tertiles
  probs <- seq(0,1,1/grp)
  data$groups <- cut(data$pred_benefit, breaks = unique(quantile(data$pred_benefit, probs=seq(0,1,1/grp))), include.lowest=TRUE)

  # Frequencies
  data1 <- data[data$D==1,]
  freq <- cbind(table(data$groups, data$W), table(data1$groups), round(table(data1$groups)/table(data$groups)*100,1), 
                table(data1$groups, data1$W),round(table(data1$groups, data1$W)/table(data$groups,data$W)*100,1))

  diff <- function(data, indicies){
    d <- data[indicies,]
    exp.delta <- aggregate(d$pred_benefit, by=list(d$groups), FUN="mean")$x

    sfit.t0.delta <- summary(survfit(Surv(d$Y[d$W==0], d$D[d$W==0]) ~ d$groups[d$W==0], se.fit = TRUE), times = times)
    sfit.t1.delta <- summary(survfit(Surv(d$Y[d$W==1], d$D[d$W==1]) ~ d$groups[d$W==1], se.fit = TRUE), times = times)
    obs.risk.t0.tmp <- 1- sfit.t0.delta$surv
    obs.risk.t1.tmp <- 1- sfit.t1.delta$surv

    # observed ARD
    obs.delta <-  obs.risk.t1.tmp - obs.risk.t0.tmp

    diff <- exp.delta-obs.delta

    return(diff)
  }
  boot_obj <- boot(data = data, statistic = diff, R=1000)
  diff <- round(boot_obj$t0,4)
  diff_ci <- matrix(NA,3,2)
  for (z in 1:3){
    diff_ci[z,] <- round(boot.ci(boot_obj, conf=0.95, type="perc", index=z)$percent[4:5],4)
  }
  
  # Expected ARR and 95% CI 
  exp.delta <- function(data, indicies){
    d <- data[indicies,]
    exp.delta <- aggregate(d$pred_benefit, by=list(d$groups), FUN="mean")$x
    return(exp.delta)
  }
  boot_obj <- boot(data = data, statistic = exp.delta, R=1000)
  exp.delta <- round(boot_obj$t0,4)
  exp.delta_ci <- matrix(NA,3,2)
  for (z in 1:3){
    exp.delta_ci[z,] <- round(boot.ci(boot_obj, conf=0.95, type="perc", index=z)$percent[4:5],4)
  }
  
  # Observed ARR and 95% CI 
  obs.delta <- function(data, indicies){
    d <- data[indicies,]
    sfit.t0.delta <- summary(survfit(Surv(d$Y[d$W==0], d$D[d$W==0]) ~ d$groups[d$W==0], se.fit = TRUE), times = times)
    sfit.t1.delta <- summary(survfit(Surv(d$Y[d$W==1], d$D[d$W==1]) ~ d$groups[d$W==1], se.fit = TRUE), times = times)
    obs.risk.t0.tmp <- 1- sfit.t0.delta$surv
    obs.risk.t1.tmp <- 1- sfit.t1.delta$surv
    obs.delta <-  obs.risk.t1.tmp - obs.risk.t0.tmp
    return(obs.delta)
  }
  boot_obj <- boot(data = data, statistic = obs.delta, R=1000)
  obs.delta <- round(boot_obj$t0,4)
  obs.delta_ci <- matrix(NA,3,2)
  for (z in 1:3){
    obs.delta_ci[z,] <- round(boot.ci(boot_obj, conf=0.95, type="perc", index=z)$percent[4:5],4)
  }
  
  labels <- c(paste0("[0-",round(probs[1],2)*100,"%)"),
              paste0("[",round(probs[1],2)*100,"-",round(probs[2],2)*100,"%)"),
              paste0("[",round(probs[2],2)*100,"-100%]"))
  tertile_table <- data.frame(labels, freq, exp.delta, exp.delta_ci,obs.delta, obs.delta_ci,diff, diff_ci)
  colnames(tertile_table) <- c("group", "n0","n1","nevent","nevent_perc",
                               "nevent0","nevent1","nevent0_perc","nevent1_perc",
                                        "exp_ARD","exp_ARD_LB","exp_ARD_UB",
                                        "obs_ARD","obs_ARD_LB","obs_ARD_UB",
                                        "diff_ARD","diff_ARD_LB","diff_ARD_UB")
  tertile_table
}


### HTE evaluation metrics (survival outcome): 1) C-for-benefit, 2) calibration 
cforbenefit <- function(data){

  set.seed(1)
  match.it <- matchit(W ~ pred_benefit, data = data, method="nearest", ratio=1)
  df.match <- match.data(match.it)
  data1 <- df.match[df.match$W==1,]
  names(data1) <- paste(names(data1),1,sep="")
  data1 <- data1[order(data1$subclass1),]
  data0 <- df.match[df.match$W==0,]
  names(data0) <- paste(names(data0),0,sep="")
  data0 <- data0[order(data0$subclass0),]
  tmp <- data.frame(data0,data1)
  
  # Observed treatment effect (survival outcome)
  # Reference:  van Klaveren D, Steyerberg EW, Serruys PW, Kent DM. 
  #             The proposed 'concordance-statistic for benefit' 
  #             provided a useful metric when modeling heterogeneous
  #             treatment effects. J Clin Epidemiol. 2018 Feb;94:59-68. 
  tmp$obs_benefit <- ifelse(tmp$D1==1 & tmp$Y1<tmp$Y0, 1,
                            ifelse(tmp$D0==1 & tmp$Y1>tmp$Y0, -1,
                                   ifelse(tmp$D0==1 & tmp$D1==1 & abs(tmp$Y1-tmp$Y0)==0,0,0)))
  
  tmp$avg_benefit <- (tmp$pred_benefit0 + tmp$pred_benefit1)/2
  tmp <- tmp%>%filter(!is.na(obs_benefit))
  out <- concordance.index(x = tmp$avg_benefit, cl = tmp$obs_benefit, method="noether")
  
  c_for_benefit <- out$c.index
  se <- out$se
  LB <- out$lower
  UB <- out$upper
  pvalue <- out$p.value
  c(c_for_benefit, LB, UB, pvalue)
}

# HTE calibration 
HTEcalib <- function(data, times, g=5, name="", calib_plot=T){
  
  # bootstrap-sampling calibration slope and intercept
  bootS <- function(data,times,indicies){
    d <- data[indicies,]
    d$cut.delta <- cut(d$pred_benefit, breaks = unique(quantile(d$pred_benefit, probs=seq(0,1,1/g))), include.lowest=TRUE)
    
    # observed ARD
    sfit.t0.delta <- summary(survfit(Surv(d$Y[d$W==0], d$D[d$W==0]) ~ d$cut.delta[d$W==0], se.fit = TRUE), extend=T,times = times)
    sfit.t1.delta <- summary(survfit(Surv(d$Y[d$W==1], d$D[d$W==1]) ~ d$cut.delta[d$W==1], se.fit = TRUE), extend=T,times = times)
    obs.risk.t0.tmp <- 1- sfit.t0.delta$surv
    obs.risk.t1.tmp <- 1- sfit.t1.delta$surv
    obs.delta <-  obs.risk.t1.tmp - obs.risk.t0.tmp
    
    # expected ARD 
    exp.delta <- aggregate(d$pred_benefit, by = list(d$cut.delta), FUN = "mean")$x
    
    # Fitted line for expected and observed ARRs
    calline <- lm(obs.delta ~ exp.delta)   
    intercept <- round(calline$coefficients[1],4)
    slope <- round(calline$coefficients[2],4)
    c(obs.delta, exp.delta, slope, intercept)
  } 
  
  bootout <- boot(data = data, statistic = bootS, R=150, times=times)
  means <- bootout$t0
  CIs <- matrix(NA,length(means),2)
  for (k in 1:length(means)){
    CIs[k,] <- boot.ci(bootout, conf=0.95, type="perc", index=k)$percent[4:5]
  }
  
  HTEests <- data.frame(obs_delta = means[1:5],
                        LB_obs = CIs[1:5,1],
                        UB_obs = CIs[1:5,2],
                        exp_delta = means[6:10],
                        LB_exp = CIs[6:10,1],
                        UB_exp = CIs[6:10,2])
  
  calibs <- data.frame(slope = means[11],
                       LB_slope = CIs[11,1],
                       UB_slope = CIs[11,2],
                       intercept = means[12],
                       LB_intercept = CIs[12,1],
                       UB_intercept = CIs[12,2])
  
  if(calib_plot==T){
    # HTE calibration plot 
    p <- ggplot(HTEests, aes(x = exp_delta, y = obs_delta)) +
      theme_classic()+
      #geom_smooth(method=lm , color="darkgrey", linetype=0)+
      geom_point() +
      xlim(-0.25,0.25)+ylim(-0.25,0.25)+
      geom_pointrange(aes(ymax = UB_obs, ymin = LB_obs), size=0.3, colour="black")+
      geom_errorbarh(aes(xmax = UB_exp, xmin = LB_exp, height = 0), size=1, colour="blue")+
      geom_abline(intercept = 0, slope = 1, col="red", size=0.5)+
      geom_abline(intercept = calibs$intercept, slope = calibs$slope, col="black", size=0.5)+
      #geom_ribbon(aes(ymin = LB_exp, ymax = UB_exp), alpha = 0.7, size=1, linetype="dashed", color="black")+
      ggtitle(name)+
      theme(axis.title.x=element_text(size=11, vjust=-2)) +
      theme(axis.title.y=element_text(size=11, angle=90, vjust=2)) +
      theme(plot.title=element_text(size=15, vjust=3,hjust=0.5)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))+
      xlab("Predicted Absolute Risk Reduction")+ylab("Observed Absolute Risk Reduction");print(p)
    filename <- paste0("./",name, " .png")
    png(filename, width = 6, height = 6, units = 'in', res = 600)
    print(p)
    dev.off()
  }
  
  return(c(calib_slope = slope, LB_slope = LB_slope, UB_slope = UB_slope, 
           calib_intercept = intercept, LB_int = LB_int, UB_int = UB_int))
}


### GND test (using existing code from this paper "https://pubmed.ncbi.nlm.nih.gov/25684707/")
#######################################################################
# R FUNCTION TO CALCULATE GREENWOOD-NAM-D'AGOSTINO CALIBRATION TEST FOR SURVIVAL MODEL
# Most up-to date version of this code is available at http://ncook.bwh.harvard.edu/r-code.html
# Author: Olga Demler, BWH, HMS
# Version 2 - Updated 8/4/2015
# FOR MORE DETAILS SEE Demler, Paynter, Cook "Tests of Calibration and Goodness of Fit 
# in the Survival Setting" Stat Med 2015; 34(10):1659-80. PMID: 25684707
# TO RUN:
# GND.calib(pred,tvar,out,cens.t, groups, adm.cens)
# PARAMETERS:
# pred - PREDICTED PROBABILITIES OF AN EVENT CALCULATED FOR THE FIXED TIME WHICH IS THE SAME FOR ALL OBSERVATIONS (=adm.cens)
# out  - OUTCOME 0/1 1=EVENT
# cens.t - CENSORED/NOT CENSORED INDICATOR 1=CENSORED
# groups - GROUPING ASSIGNMENT FOR EACH OBSERVATION
# adm.cens - END OF STUDY TIME 
# REQUIRES AT LEAST 2 EVENTS PER GROUP, AT LEAST 5 EVENTS PER GROUP IS RECOMMENDED
# IF <2 EVENTS PER GROUP THEN QUITS
#######################################################################
kmdec <- function(dec.num,dec.name, datain, adm.cens){
  stopped=0
  data.sub=datain[datain[,dec.name]==dec.num,]
  if (sum(data.sub$out)>1){
    avsurv=survfit(Surv(tvar,out) ~ 1, data=datain[datain[,dec.name]==dec.num,], error="g")
    avsurv.est=ifelse(min(avsurv$time)<=adm.cens,avsurv$surv[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],1)
    avsurv.stderr=ifelse(min(avsurv$time)<=adm.cens,avsurv$std.err[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
    avsurv.stderr=avsurv.stderr*avsurv.est
    avsurv.num=ifelse(min(avsurv$time)<=adm.cens,avsurv$n.risk[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
  } else {
    return(c(0,0,0,0,stopped=-1))
  }
  if (sum(data.sub$out)<5) stopped=1
  c(avsurv.est, avsurv.stderr, avsurv.num, dec.num, stopped) 
}
GNDf_obs  <- function(data, adm.cens, kmdec, indicies){
  data   <- data[indicies,]
  groups <- sort(unique(data$dec))
  kmtab  <- matrix(unlist(lapply(groups, kmdec,"dec", datain=data, adm.cens)),ncol=5, byrow=TRUE)
  obs_risk <- 1-kmtab[,1]
  return(obs_risk)
}

GND.calib <- function(pred, tvar, out, cens.t, groups, adm.cens, name="", calib_plot=T){
  tvar.t  <- ifelse(tvar>adm.cens, adm.cens, tvar)
  out.t   <- ifelse(tvar>adm.cens, 0, out)
  datause <- data.frame(pred=pred, tvar=tvar.t, out=out.t, count=1, cens.t=cens.t, dec=groups)
  numcat  <- length(unique(datause$dec))
  groups  <- sort(unique(datause$dec))
  kmtab   <- matrix(unlist(lapply(groups,kmdec,"dec",datain=datause, adm.cens)),ncol=5, byrow=TRUE)
  if (any(kmtab[,5] == -1)) stop("Stopped because at least one of the groups contains <2 events. Consider collapsing some groups.")
  else if (any(kmtab[,5] == 1)) warning("At least one of the groups contains < 5 events. GND can become unstable.\
                                        (see Demler, Paynter, Cook 'Tests of Calibration and Goodness of Fit in the Survival Setting' DOI: 10.1002/sim.6428) \
                                        Consider collapsing some groups to avoid this problem.")
  
  # bootstrapping 95% CI for observed risk 
  boot_obs <- boot(data = datause, statistic = GNDf_obs, R=1000, kmdec=kmdec, adm.cens=adm.cens)
  obs.risk <- boot_obs$t0
  CI_obs <- matrix(NA,10,2)
  for (q in 1:10){
    CI_obs[q,] <- boot.ci(boot_obs, conf=0.95, type="perc", index=q)$percent[4:5]
  }
  
  # data frame for expected and observed risk estimates 
  hltab <- data.frame(group = kmtab[,4],
                      totaln = tapply(datause$count,datause$dec,sum),
                      censn = tapply(datause$cens.t,datause$dec,sum),
                      numevents = tapply(datause$out,datause$dec,sum),
                      expected = tapply(datause$pred,datause$dec,sum),
                      kmperc = obs.risk,
                      #kmperc = 1-kmtab[,1],
                      kmvar = kmtab[,2]^2, 
                      kmnrisk = kmtab[,3],
                      expectedperc = tapply(datause$pred,datause$dec,mean),
                      LB_obs = CI_obs[,1],
                      UB_obs = CI_obs[,2])
  hltab$kmnum <- hltab$kmperc*hltab$totaln
  hltab$GND_component <- ifelse(hltab$kmvar==0, 0,(hltab$kmperc-hltab$expectedperc)^2/(hltab$kmvar))
  print(hltab[c(1,2,3,4,10,5,6,9,7,13)], digits=4)
  
  # fitted line - observed vs. expected 
  calline <- lm(hltab$kmperc~hltab$expectedperc)
  intercept <- round(calline$coefficients[1],4)
  slope <- round(calline$coefficients[2],4)
  LB_slope <- round(slope-1.96*summary(calline)$coefficients[2,2],4)
  UB_slope <- round(slope+1.96*summary(calline)$coefficients[2,2],4)
  LB_int   <- round(intercept-1.96*summary(calline)$coefficients[1,2],4)
  UB_int   <- round(intercept+1.96*summary(calline)$coefficients[1,2],4)
  
  if(calib_plot==T){
    
    # calibration plots
    p <- ggplot(hltab, aes(x = expectedperc, y = kmperc)) +
      theme_classic()+
      geom_smooth(method=lm, fill="darkgrey", linetype=0)+
      geom_point() +
      xlim(0,0.3)+ylim(0,0.3)+
      geom_pointrange(aes(ymax = UB_obs, ymin = LB_obs), size=0.3, colour="black")+
      geom_abline(intercept = 0, slope = 1, col="red", size=0.5)+
      geom_abline(intercept = intercept, slope = slope, col="black",size=0.5)+
      #geom_ribbon(aes(ymin = LB_exp, ymax = UB_exp), alpha = 0.7, size=1, linetype="dashed", color="black")+
      ggtitle(name)+
      theme(axis.title.x=element_text(size=11, vjust=-2)) +
      theme(axis.title.y=element_text(size=11, angle=90, vjust=2)) +
      theme(plot.title=element_text(size=15, vjust=3,hjust=0.5)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))+
      xlab("Predicted Risk")+ylab("Observed Risk");print(p)
    filename <- paste0("./",name, " .png")
    png(filename, width = 5, height = 5, units = 'in', res = 600)
    print(p)
    dev.off()
  }
  
  return(c(df = numcat-1, chi2gw = round(sum(hltab$GND_component),4),
           pvalgw = round(1-pchisq(sum(hltab$GND_component),numcat-1),4),
           slope = slope, LB_slope = LB_slope, UB_slope = UB_slope, 
           calib_intercept = intercept, LB_int = LB_int, UB_int = UB_int))
}

GND.calib.simp = function(pred, tvar, out, cens.t, groups, adm.cens){
  
  tvar.t=ifelse(tvar>adm.cens, adm.cens, tvar)
  out.t=ifelse(tvar>adm.cens, 0, out)
  
  datause=data.frame(pred=pred, tvar=tvar.t, out=out.t, count=1, cens.t=cens.t, dec=groups)
  numcat=length(unique(datause$dec))
  groups=sort(unique(datause$dec))
  
  kmtab=matrix(unlist(lapply(groups,kmdec,"dec",datain=datause, adm.cens)),ncol=5, byrow=TRUE)
  
  if (any(kmtab[,5] == -1)) stop("Stopped because at least one of the groups contains <2 events. Consider collapsing some groups.")
  else if (any(kmtab[,5] == 1)) warning("At least one of the groups contains < 5 events. GND can become unstable.\ 
(see Demler, Paynter, Cook 'Tests of Calibration and Goodness of Fit in the Survival Setting' DOI: 10.1002/sim.6428) \
Consider collapsing some groups to avoid this problem.")
  
  hltab=data.frame(group=kmtab[,4],
                   totaln=tapply(datause$count,datause$dec,sum),
                   censn=tapply(datause$cens.t,datause$dec,sum),
                   numevents=tapply(datause$out,datause$dec,sum),
                   expected=tapply(datause$pred,datause$dec,sum),
                   kmperc=1-kmtab[,1], 
                   kmvar=kmtab[,2]^2, 
                   kmnrisk=kmtab[,3],
                   expectedperc=tapply(datause$pred,datause$dec,mean))
  
  hltab$kmnum=hltab$kmperc*hltab$totaln
  hltab$GND_component=ifelse(hltab$kmvar==0, 0,(hltab$kmperc-hltab$expectedperc)^2/(hltab$kmvar))
  
  print(hltab[c(1,2,3,4,10,5,6,9,7,11)], digits=4)
  
  c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1))
}#GND.calib

## PCE calculator (based on Syalowsky's revised-pooled-ascvd code and Jenna's PCE code )
# Computing S0 using Cox model 
bsurv <- function(fit, Y, D, x){
  data <- data.frame(t_event=Y, event=D, x)
  tab <- data.frame(table(data[data$event == 1, "t_event"])) 
  y <- as.numeric(as.character(sort(unique(tab[,1]))))
  d <- tab[,2]  # number of events at each unique time                               
  
  betaHat <- as.vector(fit$coefficients)
  h0 <- rep(NA, length(y))
  for(l in 1:length(y)){
    h0[l] <- d[l] / sum(predict(fit, x[data$t_event >= y[l],], type = "risk"))
  }    
  
  S0 <- exp(-cumsum(h0))
  outcome <- data.frame(time=y,survival=S0)
  outcome
}

# PCE: compute X*beta-mean(X*beta)
original.link <- function(x, group) {
  age = x$age
  totchol = x$totchol
  hdl = x$hdl
  sysbp = x$sysbp
  dm = x$dm
  rxbp = x$rxbp
  cursmoke = x$cursmoke
  
  vars = matrix(c(log(age),
                  log(age)^2,
                  log(totchol),
                  log(age)*log(totchol),
                  log(hdl),
                  log(age)*log(hdl),
                  rxbp*log(sysbp),
                  rxbp*log(age)*log(sysbp),
                  (1-rxbp)*log(sysbp),
                  (1-rxbp)*log(age)*log(sysbp),
                  cursmoke,
                  log(age)*cursmoke,
                  dm), ncol=13)
  
  coefs = list(
    c(17.114, 0, 0.94, 0, -18.920, 4.475, 29.291, -6.432, 27.820, -6.087, 0.691, 0, 0.874),      # African American women 
    c(-29.799, 4.884, 13.54, -3.114, -13.578, 3.149, 2.019, 0, 1.957, 0, 7.574, -1.665, 0.661),  # white women 
    c(2.469, 0, 0.302, 0, -0.307, 0, 1.916, 0, 1.809, 0, 0.549, 0, 0.645),                       # African American men
    c(12.344, 0, 11.853, -2.664, -7.990, 1.769, 1.797, 0, 1.7864, 0, 7.837, -1.795, 0.658))      # white men
  
  # Crystal modified this by replacing the mean(X*beta) from the PCE paper with
  # the mean of individual.risk that calculated with the current given data 
  individual.risk = unlist(vars %*% coefs[[group]])
  link = individual.risk - mean(individual.risk)
  return(link)
}

# PCE: compute 1-S_0^exp(X*beta-mean(X*beta))
original.model <- function(x, Y, D, group) {
  
  link <- rep(NA,dim(x)[1])
  for (i in 1:length(unique(group))){
    tmp <- x[group==i,]
    index <- as.numeric(rownames(tmp))
    link[index] <- original.link(tmp,i)
  }
  relative.risk <- exp(link)
  
  # Cox PH model 
  data0 <- data.frame(Y,D,x)
  fit <- coxph(Surv(Y,D) ~ FEMALE + RACE_BLACK, data = data0)

  # compute baseline survival for any given time and given data 
  baseline.survival <- bsurv(fit=fit,
                               Y=Y, 
                               D=D, 
                               x=data0[,c("FEMALE","RACE_BLACK")])
  
  risk <- 1 - baseline.survival ^ relative.risk
  return(risk)
}
 
framingham.link <- function(x, group) {
  age = x$age
  totchol = x$totchol
  hdl = x$hdl
  sysbp = x$sysbp
  dm = x$dm
  rxbp = x$rxbp
  cursmoke = x$cursmoke
  
  vars = matrix(c(log(age),
                  log(totchol),
                  log(hdl),
                  (1-rxbp)*log(sysbp),
                  rxbp*log(sysbp),
                  cursmoke,
                  dm), ncol=7)
  
  coefs = list(
    c(2.32888, 1.20904, -0.70833, 2.76157, 2.82263, 0.52873, 0.69154),
    c(3.06117, 1.12370, -0.93263, 1.93303, 1.99881, 0.65451, 0.57367))
  coefs <- coefs[[group]]
  individual.risk = unlist(vars %*% coefs)
  link = individual.risk - mean(individual.risk)
  return(link)
}

framingham.model <- function(x, Y, D, group) {
  link <- rep(NA,dim(x)[1])
  for (i in 1:length(unique(group))){
    tmp <- x[group==i,]
    index <- as.numeric(rownames(tmp))
    link[index] <- framingham.link(tmp,i)
  }
  relative.risk = exp(link)
  
  # Cox PH model 
  ones <- rep(1, nrow(x))
  data0 <- data.frame(Y, D, ones, group)
  fit <- coxph(Surv(Y,D) ~ ones + group, data = data0)
  
  # compute baseline survival for any given time and given data 
  baseline.survival <- bsurv(fit=fit,
                             Y=Y, 
                             D=D, 
                             x=data0[,c("ones","group")])
  
  risk <- 1 - baseline.survival ^ relative.risk
  return(risk)
}


# ### Add "predictSurvProb" function for ML methods in order to use "pec" package 
# # Function for adding the surv.bart object to pec 
# predictSurvProb.survbart <- function (object, newdata, times, ...) {
#   ptemp <- matrix(BART:::predict.survbart(object, newdata=object$tx.test)$surv.test.mean, ncol = length(object$times))
#   uniquetimes <- object$times
#   pos <- prodlim::sindex(jump.times = uniquetimes, 
#                          eval.times = times)
#   p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
#   if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
#     stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
#                NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
#                NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
#   p
# }
# 
# # Function for adding the gbm object to pec 
# predictSurvProb.gbm <- function (object, newdata, times,...) {
#   
#   # Estimate the cumulative baseline hazard function using training data
#   fx <-  gbm:::predict.gbm(object, newdata = newdata[,3:dim(newdata)[2]], n.trees=object$n.trees) # log hazard constant over time
#   basehaz.cum <- rep(NA,length(sort(unique(newdata[,1]))))
#   timeinterest <- sort(unique(newdata[,1]))
#   for (q in 1:length(timeinterest)){
#     teval <- timeinterest[q]
#     basehaz.cum[q] <- gbm:::basehaz.gbm(newdata[,1], newdata[,2], fx, t.eval = teval, cumulative = TRUE) # cum baseline hazard depends on t
#   }
#   
#   # S(X)
#   ptemp <- matrix(NA,dim(newdata)[1],length(timeinterest))
#   for (m in 1:length(timeinterest)){
#     ptemp[,m] <- exp(-exp(fx)*basehaz.cum[m])
#   }
#   
#   pos <- prodlim::sindex(jump.times = timeinterest, eval.times = times)
#   p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
#   if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
#     stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
#                NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
#                NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
#   p
# }
# 
# # Function for adding the survivalforest object to pec 
# predictSurvProb.survival_forest <- function (object, newdata, times, ...) {
#   ptemp <- grf:::predict.survival_forest(object, newdata = newdata[,3:dim(newdata)[2]])$predictions
#   uniquetimes <- grf:::predict.survival_forest(object, newdata = newdata[,3:dim(newdata)[2]])$failure.times
#   pos <- prodlim::sindex(jump.times = uniquetimes, 
#                          eval.times = times)
#   p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
#   if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
#     stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
#                NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
#                NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
#   p
# }

# inverse weighing for c-for-benefit
# tmp$uncensored <- ifelse(tmp$D1==1 & tmp$Y1<tmp$Y0, 1,
#                          ifelse(tmp$D0==1 & tmp$Y1>tmp$Y0, 1,
#                                 ifelse(tmp$D0==1 & tmp$D1==1 & abs(tmp$Y1-tmp$Y0)==0,1,0)))
# 
# X = tmp[,5:41] + tmp[,49:85]
# fit <- glm(tmp$uncensored~1,family = binomial)
# weights <- 1/predict(fit, X, type="response")
# out <- concordance.index(x = tmp$avg_benefit[tmp$uncensored==1], cl = tmp$obs_benefit[tmp$uncensored==1], weights = weights[tmp$uncensored==1])
# out$c.index