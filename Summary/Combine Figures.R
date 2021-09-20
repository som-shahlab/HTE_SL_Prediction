# Combine figures for by method and estimand

foldernames <- c("EV_rate_main","EV_rate_sup_S","EV_rate_sup_TX")
foldernames <- c("IV_rate_main","IV_rate_sup_S","IV_rate_sup_T")
foldernames <- c("IV_risk_main","IV_risk_sup_T","IV_hte_main","IV_hte_sup_S","IV_hte_sup_T","IV_rate_main","IV_rate_sup_S","IV_rate_sup_T")
foldernames <- c("Fig_risk_S","Fig_risk_T","Fig_HTE_X","Fig_HTE_T","Fig_HTE_S","EV_rate_main","EV_rate_sup_S","EV_rate_sup_X")
foldernames <- c("EV_sprint_risk","EV_sprint_HTE","EV_sprint_AUTOC")
for (i in 1:2){
  setwd(paste0("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/Analysis Results/SPRINT_external_results/",foldernames[i]))
  files <- list.files(path=".", pattern="*.png", all.files=T, full.names=T)
  filelist <- lapply(files, readPNG)
  names(filelist) <- paste0(basename((files)))
  list2env(filelist, envir=.GlobalEnv)
  
  nameout=paste0("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/",foldernames[i],".png")
  png(nameout, width = 6, height = 4, units = 'in', res=300)
  par(mar=rep(0,4))
  layout(matrix(1:6, ncol=3, byrow=TRUE))
  
  for(j in 1:6) {
    img <- readPNG(names(filelist[j]))
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
  }
  dev.off()
}
