######showing example of univariate survival of 7 overlapped
######
seven.overlap.output=cb.data[cb.data$repName %in% teso,]
seven.overlap.output[order(seven.overlap.output$repName),]
dim(seven.overlap.output)
write.csv(seven.overlap.output[order(seven.overlap.output$repName),], "2.survival.TE.screen/survival.seven.example.show/seven.overlap.output.csv")
######
te.tumorexp.sel[1:4,1;4]
head(CRC.survivaldata)
######exp
seven.TE.exp=te.tumorexp.sel[, colnames(te.tumorexp.sel) %in% teso]
head(seven.TE.exp)
######
length(intersect(rownames(seven.TE.exp), rownames(CRC.survivaldata)))
###
seven.ids=intersect(rownames(seven.TE.exp), rownames(CRC.survivaldata))
seven.sur.stat=cbind(seven.TE.exp[seven.ids,], CRC.survivaldata[seven.ids,])
head(seven.sur.stat)
seven.sur.stat$OS.time=as.numeric(paste(seven.sur.stat$OS.time))
seven.sur.stat$DSS.time=as.numeric(paste(seven.sur.stat$DSS.time))
seven.sur.stat$DFI.time=as.numeric(paste(seven.sur.stat$DFI.time))
seven.sur.stat$PFI.time=as.numeric(paste(seven.sur.stat$PFI.time))
##########set group for MSTA-int
seven.sur.stat$group=ifelse(seven.sur.stat$`MSTA-int` > median(seven.sur.stat$`MSTA-int`),"high","low")
##########
###plot
seven.idx.time=c("OS.time", "DSS.time","DFI.time","PFI.time")
seven.idx.status=c("OS", "DSS","DFI","PFI")
seven.idx.ylab=c("Overall survival probability","Disease specific survival probability",
                 "disease-free interval probability","progression-free interval probability"
                 )
for (k in 1:4) {
  aa=data.frame(exp=seven.sur.stat$`MSTA-int`,
           time=seven.sur.stat[,seven.idx.time[k]],
           status=seven.sur.stat[,seven.idx.status[k]])
  aa=subset(aa,time > 0 & status>= 0)
  aa$group=ifelse(aa$exp > median(aa$exp),"high","low")
  head(aa)
  ss=Surv(aa$time, aa$status)
  cox = summary(coxph(ss~factor(aa$group,levels = c("low","high"))))
  pvalue=cox$sctest['pvalue']
  hr = round(cox$conf.int[1],2)
  hr_left = round(cox$conf.int[3],2)
  hr_right = round(cox$conf.int[4],2)
  conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
  #list(pvalue, hr, hr_left, hr_right)
  #
  #
  my.surv <- Surv(aa$time, aa$status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ factor(aa$group,levels = c("low","high")), data =  aa)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- round(exp(beta),3)
  HRse <- HR * se
  p = round(1 - pchisq((beta/se)^2, 1),3)
  HRCILL = round(exp(beta - qnorm(.975, 0, 1) * se),3)
  HRCIUL = round(exp(beta + qnorm(.975, 0, 1) * se),3)
  #summary(m)
  #
  fit1<- survfit(Surv(time, status) ~ group, data = aa)
  #
  
  do=ggsurvplot(fit1, data =aa,
                     title =paste("HR","=", HR," (", HRCILL, " - ", HRCIUL, ")","\n","logrank","P = ", p, sep=""),
                     pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                     surv.median.line = "hv",            # Add median survival lines
                     legend.title = "MSTA-int.group",               # Change legend titles
                     legend.labs = c("high","low"),  # Change legend labels
                     palette = c("tomato2","#00C5CD"),  # Use JCO journal color palette
                     #palette = c("#ca0020","#0571b0","#4daf4a"),
                     risk.table = TRUE,                  # Add No at risk table
                     cumevents = TRUE,                   # Add cumulative No of events table
                     tables.height = 0.15,               # Specify tables height
                     tables.theme = theme_cleantable(),  # Clean theme for tables
                     tables.y.text = FALSE,             # Hide tables y axis text
                     #ylab="Disease specific survival probability",
                     #ylab="Overall survival probability",
                     #ylab="progression-free interval probability",
                     #ylab="disease-free interval probability",
                     ylab=seven.idx.ylab[k],
                     xlab="Time(days)"
  )
  generate.PDF <- function(fig) {
    pdf(paste0("2.survival.TE.screen/survival.seven.example.show/", seven.idx.status[k],".pdf"),width = 4,height = 6)
    print(do)
    dev.off()
  }
  generate.PDF(fig)
}

