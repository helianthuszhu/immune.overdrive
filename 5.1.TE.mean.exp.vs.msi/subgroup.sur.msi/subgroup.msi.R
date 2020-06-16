head(cor.immune.drive.sur)
########OS
fit1<- survfit(Surv(OS.time, OS) ~ msi.TE.group.bin, data = cor.immune.drive.sur)
gc1=ggsurvplot(fit1, data =cor.immune.drive.sur,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster",               # Change legend titles
               legend.labs = c("TE.high_MSI", "TE.high_MSS","TE.low_MSI", "TE.low_MSS"),  # Change legend labels
               #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
               palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
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
####DSS
fit1<- survfit(Surv(DSS.time, DSS) ~ msi.TE.group.bin, data = cor.immune.drive.sur)
gc2=ggsurvplot(fit1, data =cor.immune.drive.sur,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster",               # Change legend titles
               legend.labs = c("TE.high_MSI", "TE.high_MSS","TE.low_MSI", "TE.low_MSS"),  # Change legend labels
               #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
               palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
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
###DFI
fit1<- survfit(Surv(DFI.time, DFI) ~ msi.TE.group.bin, data = cor.immune.drive.sur)
gc3=ggsurvplot(fit1, data =cor.immune.drive.sur,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster",               # Change legend titles
               legend.labs = c("TE.high_MSI", "TE.high_MSS","TE.low_MSI", "TE.low_MSS"),  # Change legend labels
               #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
               palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
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
####PFI
fit1<- survfit(Surv(PFI.time, PFI) ~ msi.TE.group.bin, data = cor.immune.drive.sur)
gc4=ggsurvplot(fit1, data =cor.immune.drive.sur,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster",               # Change legend titles
               legend.labs = c("TE.high_MSI", "TE.high_MSS","TE.low_MSI", "TE.low_MSS"),  # Change legend labels
               #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
               palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
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
#####
generate.PDF <- function(fig) {
  pdf("5.1.TE.mean.exp.vs.msi/subgroup.sur.msi/subgroup.msi.pdf",width = 7,height = 8)
  print(gc1)
  print(gc2)
  print(gc3)
  print(gc4)
  dev.off()
}
generate.PDF(fig)
save(cor.immune.drive.sur,file="5.1.TE.mean.exp.vs.msi/subgroup.sur.msi/subgroup.msi.RData")
