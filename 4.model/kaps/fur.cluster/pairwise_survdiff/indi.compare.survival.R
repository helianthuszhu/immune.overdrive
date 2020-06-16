######compare the survival difference between kaps group
######
######
########DSS
setsidx=c("set3","set2","set1")
surstatusidx=c("OS","DSS","PFI","DFI")
surtimeidx=c("OS.time","DSS.time","PFI.time","DFI.time")
datalist1=list()
datalist2=list()
for (i in 1:3) {
  for (j in 1:4) {
    kaps.sub.tmp=subset(kaps.td, kaps.group==setsidx[i]|kaps.group=="set4")
    table(kaps.sub.tmp$kaps.group)
    #
    ss=Surv(kaps.sub.tmp[,surtimeidx[j]], kaps.sub.tmp[, surstatusidx[j]])
    cox = summary(coxph(ss~factor(kaps.sub.tmp$kaps.group,levels = c(setsidx[i],"set4"))))
    pvalue=cox$sctest['pvalue']
    hr = round(cox$conf.int[1],2)
    hr_left = round(cox$conf.int[3],2)
    hr_right = round(cox$conf.int[4],2)
    conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
    #list(pvalue, hr, hr_left, hr_right)
    dfm <- cbind(Hazard_Ratio = hr,CI95 = conf_int,P_Value = round(pvalue, 4))
    rownames(dfm)=paste(setsidx[i], surstatusidx[j],sep = "_")
    datalist1[[j]]=dfm
  }
  aa=do.call(rbind, datalist1)
  datalist2[[i]]=aa
}
ressurcompare=do.call(rbind, datalist2)
ressurcompare=as.data.frame(ressurcompare)
ressurcompare$`HR(95%CI)`=paste0(ressurcompare$Hazard_Ratio, ressurcompare$CI95)
ressurcompare$summary=paste(ressurcompare$Hazard_Ratio, ressurcompare$CI95, ressurcompare$P_Value,sep = "\t")
write.csv(ressurcompare, "4.model/kaps/fur.cluster/pairwise_survdiff/surC.4.cluster.individually.compare.survival.pvalueandHR.set4.as.ref.csv")
########
ress <- pairwise_survdiff(Surv(DFI.time, DFI) ~ kaps.group,
                          data = kaps.td)
ress=as.data.frame(ress$p.value)
ress=round(ress, 4)
write.csv(ress, "4.model/kaps/fur.cluster/pairwise_survdiff/total.DFI.csv")
#
ress <- pairwise_survdiff(Surv(OS.time, OS) ~ kaps.group,
                          data = subset(kaps.td,MSI.status.bin=="MSI" ))
ress=as.data.frame(ress$p.value)
ress=round(ress, 4)
write.csv(ress, "4.model/kaps/fur.cluster/pairwise_survdiff/total.OS.MSI.csv")




kaps.sub.tmp=subset(kaps.td, MSI.status.bin=="MSS" )
ress <- pairwise_survdiff(Surv(DSS.time,DSS) ~ kaps.group,
                          data = kaps.sub.tmp)
ress








kaps.sub.tmp=subset(kaps.td, MSI.status.bin=="MSS" & (kaps.group=="set4" | kaps.group=="set3") )
#
fit1<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data = kaps.sub.tmp)
ggsurvplot(fit1, data =kaps.sub.tmp,
                  # title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster",               # Change legend titles
                   legend.labs = c("set3","set4"),  # Change legend labels
                   palette = c("set4"="#d73027","set3"="#E69F00"),  # Use JCO journal color palette
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
#
kaps.sub.tmp=subset(kaps.td, MSI.status.bin=="MSI" )
fit1<- survfit(Surv(OS.time, OS) ~ kaps.group, data = kaps.sub.tmp)
ggsurvplot(fit1, data =kaps.sub.tmp,
           # title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TE.cluster",               # Change legend titles
           legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
           palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),  # Use JCO journal color palette
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
