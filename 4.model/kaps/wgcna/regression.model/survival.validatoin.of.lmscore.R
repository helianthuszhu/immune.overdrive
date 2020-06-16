######
######test on regression model lmscore of kaps classification on CRC
#####################
#####################
load("4.model/kaps/wgcna/regression.model/stat.roctest.clin.CRC.lmscore.for.kaps.test.output.RData")
#######
#######
head(stat.roctest.clin)
summary(stat.roctest.clin$lmscore)
#####
stat.roctest.clin.lmscore.kaps=stat.roctest.clin
stat.roctest.clin.lmscore.kaps$kaps.group.lmscore<- cut(stat.roctest.clin.lmscore.kaps$lmscore, 
                                                      breaks = c(min(stat.roctest.clin.lmscore.kaps$lmscore),0.2,0.3,0.47,max(stat.roctest.clin.lmscore.kaps$lmscore)),
                                                      labels = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"), include.lowest = TRUE)
stat.roctest.clin.lmscore.kaps$kaps.group.lmscore.quantile=cut(as.numeric(as.character(stat.roctest.clin.lmscore.kaps$lmscore)),
                                                               breaks=quantile(as.numeric(as.character(stat.roctest.clin.lmscore.kaps$lmscore)), 
                                                                               c(0, 0.25,0.5,0.75, 1), na.rm=T),labels=c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"), include.lowest = TRUE)
table(stat.roctest.clin.lmscore.kaps$kaps.group.lmscore.quantile)
########
########
#####OS
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group.lmscore.quantile, data = stat.roctest.clin.lmscore.kaps )
kp1=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
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
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group.lmscore.quantile, data = stat.roctest.clin.lmscore.kaps )
kp2=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
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
fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group.lmscore.quantile, data = stat.roctest.clin.lmscore.kaps)
kp3=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
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
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group.lmscore.quantile, data = stat.roctest.clin.lmscore.kaps )
kp4=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
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
  pdf("4.model/kaps/wgcna/regression.model/sur.lmscore.quantile.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
#######################
#################survival test of these 6 markers
######################
head(stat.roctest.clin.lmscore.kaps)
stat.roctest.clin.lmscore.kaps$group.gene=cut(as.numeric(as.character(stat.roctest.clin.lmscore.kaps$STAT2)),
                                              breaks=quantile(as.numeric(as.character(stat.roctest.clin.lmscore.kaps$STAT2)), 
                                                              c(0, 0.25,0.5,0.75, 1), na.rm=T),labels=c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),
                                              include.lowest = TRUE)

#####OS
fit2<- survfit(Surv(OS.time, OS) ~ group.gene, data = stat.roctest.clin.lmscore.kaps )
kp1=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"), # Use JCO journal color palette,
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
fit2<- survfit(Surv(DSS.time, DSS) ~ group.gene, data = stat.roctest.clin.lmscore.kaps )
kp2=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),# Use JCO journal color palette,
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
fit2<- survfit(Surv(DFI.time, DFI) ~ group.gene, data = stat.roctest.clin.lmscore.kaps)
kp3=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"), # Use JCO journal color palette,
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
fit2<- survfit(Surv(PFI.time, PFI) ~ group.gene, data = stat.roctest.clin.lmscore.kaps )
kp4=ggsurvplot(fit2, data =stat.roctest.clin.lmscore.kaps,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),# Use JCO journal color palette,
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
  pdf("4.model/kaps/wgcna/regression.model/sur.gene.JAK3.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
