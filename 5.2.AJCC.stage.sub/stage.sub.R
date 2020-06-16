dim(cor.immune.drive.sur)
dim(clin.data.CRC.selected)
length(intersect(rownames(cor.immune.drive.sur), rownames(clin.data.CRC.selected)))
#####
cor.immune.drive.sur.stage=cbind(cor.immune.drive.sur, clin.data.CRC.selected[rownames(cor.immune.drive.sur),])
head(cor.immune.drive.sur.stage)
####
library(survminer)
library(survival)
######
idx.stage=unique(cor.immune.drive.sur.stage$AJCC.stage)[1:4]
for (i in 1:length(idx.stage)) {
  stage.sub=subset(cor.immune.drive.sur.stage, AJCC.stage==idx.stage[i])
  ss=Surv(stage.sub$PFI.time, stage.sub$PFI)
  #
  cox = summary(coxph(ss~factor(stage.sub$TE.cluster.agg,levels = c("TE.low","TE.high"))))
  pvalue=cox$sctest['pvalue']
  hr = round(cox$conf.int[1],2)
  hr_left = round(cox$conf.int[3],2)
  hr_right = round(cox$conf.int[4],2)
  conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
  list(pvalue, hr, hr_left, hr_right)
  #
  fit1<- survfit(Surv(PFI.time, PFI) ~ TE.cluster.agg, data =stage.sub)
  pdf(paste0("5.2.AJCC.stage.sub/",idx.stage[i],".PFI.","pdf"))
  print(ggsurvplot(fit1, data =stage.sub,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = paste0(idx.stage[i],"TE.cluster"),               # Change legend titles
                   legend.labs = c("TE.high","TE.low"),  # Change legend labels
                   palette = c("#ca0020","#0571b0","#00C5CD","tomato2"),  # Use JCO journal color palette
                   #palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
                   risk.table = TRUE,                  # Add No at risk table
                   cumevents = TRUE,                   # Add cumulative No of events table
                   tables.height = 0.15,               # Specify tables height
                   tables.theme = theme_cleantable(),  # Clean theme for tables
                   tables.y.text = FALSE,             # Hide tables y axis text
                   #ylab="Disease specific survival probability",
                   #ylab="Overall survival probability",
                   ylab="progression-free interval probability",
                   #ylab="disease-free interval probability",
                   xlab="Time(days)"
  )
  )
  dev.off()
}
save(cor.immune.drive.sur.stage, file="5.2.AJCC.stage.sub/stage.sub.RData")

