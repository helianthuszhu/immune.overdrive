######validation on KIRC dataset
#####
#####survival
#############
survival.info.kirc.tcga=read.table("8.KIRC/KIRC_survival.txt.gz",header = T,row.names = 1,sep = "\t")
survival.info.kirc.tcga=survival.info.kirc.tcga[,-10]
head(survival.info.kirc.tcga)
############
head(mean.stats.kirc)
dim(mean.stats.kirc)
stat.kirc.kaps.vali=mean.stats.kirc
stat.kirc.kaps.vali$z.of.mean.exp=(stat.kirc.kaps.vali$mean.exp - mean(stat.kirc.kaps.vali$mean.exp))/sd(stat.kirc.kaps.vali$mean.exp)
stat.kirc.kaps.vali=stat.kirc.kaps.vali[order(stat.kirc.kaps.vali$z.of.mean.exp,decreasing = T),]
stat.kirc.kaps.vali$pos=seq(1:nrow(stat.kirc.kaps.vali))
summary(stat.kirc.kaps.vali$z.of.mean.exp)
head(stat.kirc.kaps.vali)
table(stat.kirc.kaps.vali$OS)
#save(stat.kirc.kaps.vali,file="8.KIRC/1.tescore.survival.kaps.validation/for.kaps.classification.RData")
###use the same cutoff value
###
stat.kirc.kaps.vali$kaps.group=cut(stat.kirc.kaps.vali$z.of.mean.exp, 
    breaks = c(min(stat.kirc.kaps.vali$z.of.mean.exp),0.33,0.957,1.347,max(stat.kirc.kaps.vali$z.of.mean.exp)),
    labels = c("set1","set2","set3","set4"), include.lowest = TRUE)

head(stat.kirc.kaps.vali)
table(stat.kirc.kaps.vali$kaps.group)
####load the kaps results
load("8.KIRC/1.tescore.survival.kaps.validation/for.kaps.classification.output.k4.kirc.RData")
fitkaps.kirc
stat.kirc.kaps.vali$kaps.group.kirc=cut(stat.kirc.kaps.vali$z.of.mean.exp, 
                                   breaks = c(min(stat.kirc.kaps.vali$z.of.mean.exp),0.24,0.42,0.74,max(stat.kirc.kaps.vali$z.of.mean.exp)),
                                   labels = c("set1","set2","set3","set4"), include.lowest = TRUE)
table(stat.kirc.kaps.vali$kaps.group.kirc,stat.kirc.kaps.vali$kaps.group)
############
########test the survival
#####OS
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group.kirc, data = stat.kirc.kaps.vali )
kp1=ggsurvplot(fit2, data =stat.kirc.kaps.vali,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
#####
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group.kirc, data = stat.kirc.kaps.vali )
kp2=ggsurvplot(fit2, data =stat.kirc.kaps.vali,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
###
fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group.kirc, data = stat.kirc.kaps.vali )
kp3=ggsurvplot(fit2, data =stat.kirc.kaps.vali,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               ylab="disease-free interval probability",
               xlab="Time(days)"
)
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group.kirc, data = stat.kirc.kaps.vali )
kp4=ggsurvplot(fit2, data =stat.kirc.kaps.vali,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
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
##
generate.PDF <- function(fig) {
  pdf("8.KIRC/1.tescore.survival.kaps.validation/surC.same.cutoff.kirc.test.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
###################

