library(survival);library(ggplot2);library(gtsummary);library(boot);library(survcomp);library(nricens); library(png)
library(gridExtra);library(grid);library(lattice)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results") 

# Internal risk forest plot
IV_risks <- read.csv("IV_risks.csv"); colnames(IV_risks)[1:2] <- c("SL", "Outcome")
IV_risks$SL <- c(rep("Stepwise",2),rep("Elastic Net",2),
                 rep("GBM",4),rep("Deepsurv",4),rep("RSF",4),rep("BART",4))
IV_risks$SL <- factor(IV_risks$SL, levels = c("Stepwise","Elastic Net","GBM","Deepsurv","RSF","BART"))
IV_risks$Outcome <- c(rep(c("SL-CVD", "SL-SAE"),2), rep(c("SL-CVD","TL-CVD","SL-SAE","TL-SAE"),4))

# C index
IV_risks_Cindex <- IV_risks[,1:5]
cindex = ggplot(data=IV_risks_Cindex,
           aes(x = SL, y = Cindex, ymin = CLB, ymax = CUB ))+
  geom_pointrange(aes(col=SL))+
  #theme_classic()+
  #geom_hline(aes(fill=SL),yintercept =0.5, linetype=2)+
  xlab('Methods')+ ylab("C index (95% CI)")+
  ylim(0.5,0.8)+
  geom_errorbar(aes(ymin=CLB, ymax=CUB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0.5,vjust = 1,angle=180,face="bold"),
        legend.position = "bottom")+
  coord_flip()
cindex

# Slope
IV_risks_slope <- IV_risks[,c(1:2,6:8)]
slope = ggplot(data=IV_risks_slope,
           aes(x = SL, y = Slope, ymin = Slope.LB, ymax = Slope.UB ))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept =1, linetype=2, size=1)+
  #xlab('Methods')+
  ylab("Slope (95% CI)")+
  ylim(min(IV_risks_slope$Slope.LB), max(IV_risks_slope$Slope.UB))+
  geom_errorbar(aes(ymin=Slope.LB, ymax=Slope.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
slope

# Intercept
IV_risks_intercept <- IV_risks[,c(1:2,9:11)]; colnames(IV_risks_intercept)[4:5] <- c("Int.LB","Int.UB")
intercept = ggplot(data=IV_risks_intercept,
           aes(x = SL, y = Int, ymin = Int.LB, ymax = Int.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2, size=1)+
  #xlab('Methods')+
  ylab("Intercept (95% CI)")+
  ylim(min(IV_risks_intercept$Int.LB), max(IV_risks_intercept$Int.UB))+
  geom_errorbar(aes(ymin=Int.LB, ymax=Int.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=4,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
intercept

# Combine three figures in one row 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(cindex)

filename <- paste0("./IV_risk_forest.png")
png(filename, width = 6.5, height = 10, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(cindex+theme(legend.position="none"), 
                               slope+theme(legend.position="none"),
                               intercept+theme(legend.position="none"),
                               nrow=1, widths = c(4.5,4,4)), mylegend, nrow=3, heights=c(10,1,0)))
dev.off()


# External risk forest plot
EV_risks <- read.csv("EV_risks.csv"); colnames(EV_risks)[1:2] <- c("SL", "Outcome")
EV_risks$SL <- c(rep("Stepwise",1),rep("Elastic Net",1),
                 rep("GBM",2),rep("Deepsurv",2),rep("RSF",2),rep("BART",2))
EV_risks$SL <- factor(EV_risks$SL, levels = c("Stepwise","Elastic Net","GBM","Deepsurv","RSF","BART"))
EV_risks$Outcome <- c(rep("SL-CVD",2), rep(c("SL-CVD","TL-CVD"),4))

# C index
EV_risks_Cindex <- EV_risks[,1:5]
cindex = ggplot(data=EV_risks_Cindex,
                aes(x = SL, y = Cinedx, ymin = CLB, ymax = CUB ))+
  geom_pointrange(aes(col=SL))+
  #theme_classic()+
  #geom_hline(aes(fill=SL),yintercept =0.5, linetype=2)+
  xlab('Methods')+ ylab("C index (95% CI)")+
  ylim(0.5,0.8)+
  geom_errorbar(aes(ymin=CLB, ymax=CUB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0.5,vjust = 1,angle=180,face="bold"),
        legend.position = "bottom")+
  coord_flip()
cindex

# Slope
EV_risks_slope <- EV_risks[,c(1:2,6:8)]
slope = ggplot(data=EV_risks_slope,
               aes(x = SL, y = slope, ymin = SLB, ymax = SUB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept =1, linetype=2, size=1)+
  ylab("Slope (95% CI)")+
  ylim(min(EV_risks_slope$SLB), max(EV_risks_slope$SUB))+
  geom_errorbar(aes(ymin=SLB, ymax=SUB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
slope

# Intercept
EV_risks_intercept <- EV_risks[,c(1:2,9:11)]; colnames(EV_risks_intercept)[4:5] <- c("Int.LB","Int.UB")
intercept = ggplot(data=EV_risks_intercept,
                   aes(x = SL, y = Int, ymin = Int.LB, ymax = Int.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2, size=1)+
  #xlab('Methods')+
  ylab("Intercept (95% CI)")+
  ylim(min(EV_risks_intercept$Int.LB), max(EV_risks_intercept$Int.UB))+
  geom_errorbar(aes(ymin=Int.LB, ymax=Int.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=4,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
intercept

# Combine three figures in one row 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(cindex)

filename <- paste0("./EV_risk_forest.png")
png(filename, width = 10, height = 8, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(cindex+theme(legend.position="none"), 
                               slope+theme(legend.position="none"),
                               intercept+theme(legend.position="none"),
                               nrow=1, widths = c(4.5,4,4)), mylegend, nrow=3, heights=c(10,1,0)))
dev.off()


# External SPRINT risk forest plot
EV_risks <- read.csv("EV_sprint_risks.csv"); colnames(EV_risks)[1:2] <- c("SL", "Outcome")
EV_risks$SL <- c(rep("Stepwise",1),rep("Elastic Net",1),
                 rep("GBM",2),rep("Deepsurv",2),rep("RSF",2),rep("BART",2))
EV_risks$SL <- factor(EV_risks$SL, levels = c("Stepwise","Elastic Net","GBM","Deepsurv","RSF","BART"))

# C index
EV_risks_Cindex <- EV_risks[,1:5]
cindex = ggplot(data=EV_risks_Cindex,
                aes(x = SL, y = Cindex, ymin = C_lb, ymax = C_ub))+
  geom_pointrange(aes(col=SL))+
  #theme_classic()+
  #geom_hline(aes(fill=SL),yintercept =0.5, linetype=2)+
  xlab('Methods')+ ylab("C index (95% CI)")+
  ylim(0.5,0.8)+
  geom_errorbar(aes(ymin=C_lb, ymax=C_ub,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0.5,vjust = 1,angle=180,face="bold"),
        legend.position = "bottom")+
  coord_flip()
cindex

# Slope
EV_risks_slope <- EV_risks[,c(1:2,6:8)]
slope = ggplot(data=EV_risks_slope,
               aes(x = SL, y = slope, ymin = slope_lb, ymax = slope_ub))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept =1, linetype=2, size=1)+
  ylab("Slope (95% CI)")+
  ylim(min(EV_risks_slope$slope_lb), max(EV_risks_slope$slope_ub))+
  geom_errorbar(aes(ymin=slope_lb, ymax=slope_ub,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
slope

# Intercept
EV_risks_intercept <- EV_risks[,c(1:2,9:11)]; colnames(EV_risks_intercept)[3:5] <- c("Int","Int.LB","Int.UB")
intercept = ggplot(data=EV_risks_intercept,
                   aes(x = SL, y = Int, ymin = Int.LB, ymax = Int.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2, size=1)+
  #xlab('Methods')+
  ylab("Intercept (95% CI)")+
  ylim(min(EV_risks_intercept$Int.LB), max(EV_risks_intercept$Int.UB))+
  geom_errorbar(aes(ymin=Int.LB, ymax=Int.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=4,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
intercept

# Combine three figures in one row 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(cindex)

filename <- paste0("./EV_sprint_risk_forest.png")
png(filename, width = 10, height = 8, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(cindex+theme(legend.position="none"), 
                               slope+theme(legend.position="none"),
                               intercept+theme(legend.position="none"),
                               nrow=1, widths = c(4.5,4,4)), mylegend, nrow=3, heights=c(10,1,0)))
dev.off()


# Internal HTE forest plot
IV_HTEs <- read.csv("IV_HTEs.csv"); colnames(IV_HTEs)[1:2] <- c("SL", "Outcome")
IV_HTEs$SL <- c(rep("Stepwise",2),rep("Elastic Net",2),
                 rep("GBM",6),rep("Deepsurv",6),rep("RSF",6),rep("BART",6))
IV_HTEs$SL <- factor(IV_HTEs$SL, levels = c("Stepwise","Elastic Net","GBM","Deepsurv","RSF","BART"))
IV_HTEs$Outcome <- c(rep(c("SL-CVD", "SL-SAE"),2), rep(c("SL-CVD","TL-CVD","XL-CVD","SL-SAE","TL-SAE","XL-SAE"),4))

# AUTOC
IV_HTEs_AUTOC <- IV_HTEs[,c(1:3,5:6)]
AUTOC <- ggplot(data=IV_HTEs_AUTOC,
                aes(x = SL, y = AUTOC, ymin = LB, ymax = UB))+
  geom_pointrange(aes(col=SL))+
  #theme_classic()+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2,size=1)+
  xlab('Methods')+ ylab("AUTOC (95% CI)")+
  ylim(min(IV_HTEs_AUTOC$LB),max(IV_HTEs_AUTOC$UB))+
  geom_errorbar(aes(ymin=LB, ymax=UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0.5,vjust = 1,angle=180,face="bold"),
        legend.position = "bottom")+
  coord_flip()
AUTOC

# Slope
IV_HTEs_slope <- IV_HTEs[,c(1:2,7:9)]; colnames(IV_HTEs_slope)[4:5] <- c("Slope.LB","Slope.UB")
slope = ggplot(data=IV_HTEs_slope,
               aes(x = SL, y = Slope, ymin = Slope.LB, ymax = Slope.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept =1, linetype=2, size=1)+
  ylab("Slope (95% CI)")+
  ylim(min(IV_HTEs_slope$Slope.LB), max(IV_HTEs_slope$Slope.UB))+
  geom_errorbar(aes(ymin=Slope.LB, ymax=Slope.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=6,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
slope

# Intercept
IV_HTEs_intercept <- IV_HTEs[,c(1:2,10:12)]; colnames(IV_HTEs_intercept)[4:5] <- c("Int.LB","Int.UB")
intercept = ggplot(data=IV_HTEs_intercept,
                   aes(x = SL, y = Int, ymin = Int.LB, ymax = Int.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2, size=1)+
  ylab("Intercept (95% CI)")+
  ylim(min(IV_HTEs_intercept$Int.LB), max(IV_HTEs_intercept$Int.UB))+
  geom_errorbar(aes(ymin=Int.LB, ymax=Int.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=6,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
intercept

# Combine three figures in one row 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(AUTOC)

filename <- paste0("./IV_HTEs_forest.png")
png(filename, width = 6.5, height = 10, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(AUTOC+theme(legend.position="none"), 
                               slope+theme(legend.position="none"),
                               intercept+theme(legend.position="none"),
                               nrow=1, widths = c(4.5,4,4)), mylegend, nrow=3, heights=c(10,1,0)))
dev.off()


# External HTE forest plot
EV_HTEs <- read.csv("EV_HTEs.csv"); colnames(EV_HTEs)[1:2] <- c("SL", "Outcome")
EV_HTEs$SL <- c(rep("Stepwise",1),rep("Elastic Net",1),
                rep("GBM",3),rep("Deepsurv",3),rep("RSF",3),rep("BART",3))
EV_HTEs$SL <- factor(EV_HTEs$SL, levels = c("Stepwise","Elastic Net","GBM","Deepsurv","RSF","BART"))
EV_HTEs$Outcome <- c(rep(c("SL-CVD"),2), rep(c("SL-CVD","TL-CVD","XL-CVD"),4))

# AUTOC
EV_HTEs_AUTOC <- EV_HTEs[,c(1:3,5:6)]
AUTOC <- ggplot(data=EV_HTEs_AUTOC,
                aes(x = SL, y = AUTOC, ymin = LB, ymax = UB))+
  geom_pointrange(aes(col=SL))+
  #theme_classic()+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2,size=1)+
  xlab('Methods')+ ylab("AUTOC (95% CI)")+
  ylim(min(EV_HTEs_AUTOC$LB),max(EV_HTEs_AUTOC$UB))+
  geom_errorbar(aes(ymin=LB, ymax=UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0.5,vjust = 1,angle=180,face="bold"),
        legend.position = "bottom")+
  coord_flip()
AUTOC

# Slope
EV_HTEs_slope <- EV_HTEs[,c(1:2,7:9)]; colnames(EV_HTEs_slope)[4:5] <- c("Slope.LB","Slope.UB")
slope = ggplot(data=EV_HTEs_slope,
               aes(x = SL, y = Slope, ymin = Slope.LB, ymax = Slope.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept =1, linetype=2, size=1)+
  ylab("Slope (95% CI)")+
  ylim(min(EV_HTEs_slope$Slope.LB), max(EV_HTEs_slope$Slope.UB))+
  geom_errorbar(aes(ymin=Slope.LB, ymax=Slope.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=6,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
slope

# Intercept
EV_HTEs_intercept <- EV_HTEs[,c(1:2,10:12)]; colnames(EV_HTEs_intercept)[4:5] <- c("Int.LB","Int.UB")
intercept = ggplot(data=EV_HTEs_intercept,
                   aes(x = SL, y = Int, ymin = Int.LB, ymax = Int.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2, size=1)+
  ylab("Intercept (95% CI)")+
  ylim(min(EV_HTEs_intercept$Int.LB), max(EV_HTEs_intercept$Int.UB))+
  geom_errorbar(aes(ymin=Int.LB, ymax=Int.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=6,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
intercept

# Combine three figures in one row 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(AUTOC)

filename <- paste0("./EV_HTEs_forest.png")
png(filename, width = 10, height = 8, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(AUTOC+theme(legend.position="none"), 
                               slope+theme(legend.position="none"),
                               intercept+theme(legend.position="none"),
                               nrow=1, widths = c(4.5,4,4)), mylegend, nrow=3, heights=c(10,1,0)))
dev.off()


# External sprint HTE forest plot
EV_HTEs <- read.csv("EV_sprint_HTEs.csv"); colnames(EV_HTEs)[1:2] <- c("SL", "Outcome")
EV_HTEs$SL <- c(rep("Stepwise",1),rep("Elastic Net",1),
                rep("GBM",3),rep("Deepsurv",3),rep("RSF",3),rep("BART",3))
EV_HTEs$SL <- factor(EV_HTEs$SL, levels = c("Stepwise","Elastic Net","GBM","Deepsurv","RSF","BART"))

# AUTOC
EV_HTEs_AUTOC <- EV_HTEs[,c(1:5)]; colnames(EV_HTEs_AUTOC)[3:5] <- c("AUTOC", "LB","UB")
AUTOC <- ggplot(data=EV_HTEs_AUTOC,
                aes(x = SL, y = AUTOC, ymin = LB, ymax = UB))+
  geom_pointrange(aes(col=SL))+
  #theme_classic()+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2,size=1)+
  xlab('Methods')+ ylab("AUTOC (95% CI)")+
  ylim(min(EV_HTEs_AUTOC$LB),max(EV_HTEs_AUTOC$UB))+
  geom_errorbar(aes(ymin=LB, ymax=UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0.5,vjust = 1,angle=180,face="bold"),
        legend.position = "bottom")+
  coord_flip()
AUTOC

# Slope
EV_HTEs_slope <- EV_HTEs[,c(1:2,6:8)]; colnames(EV_HTEs_slope)[3:5] <- c("Slope","Slope.LB","Slope.UB")
EV_HTEs_slope <- EV_HTEs_slope[-c(3,9),]
slope = ggplot(data=EV_HTEs_slope,
               aes(x = SL, y = Slope, ymin = Slope.LB, ymax = Slope.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept =1, linetype=2, size=1)+
  ylab("Slope (95% CI)")+
  ylim(min(EV_HTEs_slope$Slope.LB), max(EV_HTEs_slope$Slope.UB))+
  geom_errorbar(aes(ymin=Slope.LB, ymax=Slope.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=6,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
slope

# Intercept
EV_HTEs_intercept <- EV_HTEs[,c(1:2,9:11)]; colnames(EV_HTEs_intercept)[3:5] <- c("Int","Int.LB","Int.UB")
EV_HTEs_intercept <- EV_HTEs_intercept[-c(3,9),]
intercept = ggplot(data=EV_HTEs_intercept,
                   aes(x = SL, y = Int, ymin = Int.LB, ymax = Int.UB))+
  geom_pointrange(aes(col=SL))+
  geom_hline(aes(fill=SL),yintercept = 0, linetype=2, size=1)+
  ylab("Intercept (95% CI)")+
  ylim(min(EV_HTEs_intercept$Int.LB), max(EV_HTEs_intercept$Int.UB))+
  geom_errorbar(aes(ymin=Int.LB, ymax=Int.UB,col=SL),width=0.5,cex=1)+ 
  facet_wrap(~Outcome, strip.position="left",nrow=6,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_blank())+
  coord_flip()
intercept

# Combine three figures in one row 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(AUTOC)

filename <- paste0("./EV_sprint_HTEs_forest.png")
png(filename, width = 10, height = 8, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(AUTOC+theme(legend.position="none"), 
                               slope+theme(legend.position="none"),
                               intercept+theme(legend.position="none"),
                               nrow=1, widths = c(4.5,4,4)), mylegend, nrow=3, heights=c(10,1,0)))
dev.off()










