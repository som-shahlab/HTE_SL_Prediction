setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/Tuned Results") 
bart_imp_vars <- read.csv("IV_imp_bart.csv")
deepsurv_imp_vars <- read.csv("IV_imp_deepsurv.csv")
gbm_imp_vars <- read.csv("IV_imp_gbm.csv")
sf_imp_vars <- read.csv("IV_imp_sf.csv")

# CVD
vars <- bart_imp_vars$cvd_vars[1:15]
bart_imp_vars_sub <- bart_imp_vars[bart_imp_vars$cvd_vars %in% vars,1:2]
deepsurv_imp_vars_sub <- deepsurv_imp_vars[deepsurv_imp_vars$cvd_vars %in% vars,1:2]
gbm_imp_vars_sub <- gbm_imp_vars[gbm_imp_vars$cvd_vars %in% vars,1:2]
sf_imp_vars_sub <- sf_imp_vars[sf_imp_vars$cvd_vars %in% vars,1:2]
cvd_data <- rbind(bart_imp_vars_sub, sf_imp_vars_sub, gbm_imp_vars_sub, deepsurv_imp_vars_sub)
cvd_data$Methods <- c(rep("BART",15),rep("SF",15),rep("GBM",15),rep("Deepsurv",15))
cvd_data$cvd_vars <- as.character(cvd_data$cvd_vars)
cvd_data$cvd_vars <- factor(cvd_data$cvd_vars, levels = cvd_data$cvd_vars[1:15])

p1 <- ggplot(cvd_data, aes(x = cvd_vars, y = cvd_imp, group = Methods, fill = Methods))+
  #theme_classic()+
  scale_fill_manual(values=c("bisque4", "darkgoldenrod1","cadetblue3","darkgoldenrod4"))+
  geom_bar(stat = "identity", width = 0.5, position = "dodge")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_text(size=11, vjust=2)) +
  theme(axis.title.y=element_text(size=11, angle=90,vjust=3)) +
  theme(plot.title=element_text(size=15, vjust=3, hjust=0.5))+
  #scale_x_discrete(labels= xlabels)+
  xlab("Variables")+ylab("Variable Importance for CVD Outcome");print(p1)

# SAE
vars <- bart_imp_vars$sae_vars[1:15]
bart_imp_vars_sub <- bart_imp_vars[bart_imp_vars$sae_vars %in% vars,3:4]
deepsurv_imp_vars_sub <- deepsurv_imp_vars[deepsurv_imp_vars$sae_vars %in% vars,3:4]
gbm_imp_vars_sub <- gbm_imp_vars[gbm_imp_vars$sae_vars %in% vars,3:4]
sf_imp_vars_sub <- sf_imp_vars[sf_imp_vars$sae_vars %in% vars,3:4]
sae_data <- rbind(bart_imp_vars_sub, sf_imp_vars_sub, gbm_imp_vars_sub, deepsurv_imp_vars_sub)
sae_data$Methods <- c(rep("BART",15),rep("SF",15),rep("GBM",15),rep("Deepsurv",15))
sae_data$sae_vars <- as.character(sae_data$sae_vars)
sae_data$sae_vars <- factor(sae_data$sae_vars, levels = sae_data$sae_vars[1:15])

p2 <- ggplot(sae_data, aes(x = sae_vars, y = sae_imp, group = Methods, fill = Methods))+
  #theme_classic()+
  scale_fill_manual(values=c("bisque4", "darkgoldenrod1","cadetblue3","darkgoldenrod4"))+
  geom_bar(stat = "identity", width = 0.5, position = "dodge")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_text(size=11, vjust=2)) +
  theme(axis.title.y=element_text(size=11, angle=90,vjust=3)) +
  theme(plot.title=element_text(size=15, vjust=3, hjust=0.5))+
  theme(legend.position="bottom")+
  #scale_x_discrete(labels= xlabels)+
  xlab("Variables")+ylab("Variable Importance for Severe Adverse Events");print(p2)

filename <- paste0("./vars_imp.png")
png(filename, width = 8, height = 10, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2, nrow=2, ncol=1), nrow=2, heights=c(10,1)))
dev.off()
