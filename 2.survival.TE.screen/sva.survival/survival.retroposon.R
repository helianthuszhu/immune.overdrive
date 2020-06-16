################
head(CRC.survivaldata)
head(kaps.td)
load("2.survival.TE.screen/sva.survival/teclass.tumorexp.filterd.RData")
head(teclass.tumorexp.filterd)
dim(teclass.tumorexp.filterd)
###get the expression of SVFs   only focus on 590 samples
#te.tumorexp.sel[1:4,1:4]
svaexp=te.tumorexp.sel[, grep("SVA",colnames(te.tumorexp.sel))]
#svaexp$z.of.mean.sva=colMeans(svaexp)
#svaexp=as.data.frame(scale(svaexp))
###############
#dim(svaexp)
head(svaexp)
######
length(intersect(rownames(kaps.td), rownames(svaexp)))
length(intersect(rownames(kaps.td), rownames(teclass.tumorexp.filterd)))
surstat.sva.crc=cbind(teclass.tumorexp.filterd[rownames(kaps.td),], kaps.td)
surstat.sva.crc=cbind(svaexp[rownames(kaps.td),c(1,2,4,5)], surstat.sva.crc)
#colnames(surstat.sva.crc)[1]="z.of.mean.sva"
head(surstat.sva.crc)
dim(surstat.sva.crc)
#
head(clin.data.CRC.selected)
dim(clin.data.CRC.selected)
length(intersect(rownames(surstat.sva.crc), rownames(clin.data.CRC.selected)))
surstat.sva.crc=cbind(surstat.sva.crc, clin.data.CRC.selected[rownames(surstat.sva.crc),])
head(surstat.sva.crc)
dim(surstat.sva.crc)
#
################
################
##########
library(survminer)
library(survival)
ss=Surv(surstat.sva.crc$PFI.time, surstat.sva.crc$PFI)
surstat.sva.crc$risk.group.sva.agg=ifelse(surstat.sva.crc$Retroposon>= median(surstat.sva.crc$Retroposon),"2high","1low")
#surstat.sva.crc$risk.group.sva=ifelse(stats.data[, colnames(stats.data) %in% index.sva[i]]> median(stats.data[, colnames(stats.data) %in% index.sva[i]]),"2high","1low")
cox = summary(coxph(ss~surstat.sva.crc$risk.group.sva.agg))
  pvalue=cox$sctest['pvalue']
  hr = round(cox$conf.int[1],2)
  hr_left = round(cox$conf.int[3],2)
  hr_right = round(cox$conf.int[4],2)
  conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep=""); 
  #txt = paste("HR = ", hr, conf_int, "\nlogrank P = ", signif(pvalue, 2), sep="")
  #text(grconvertX(0.98, "npc"), grconvertY(.97, "npc"),labels=txt, adj=c(1, 1))
  list(pvalue, hr, hr_left, hr_right)
  #
fit2 <- survfit(Surv(surstat.sva.crc$PFI.time, surstat.sva.crc$PFI) ~ risk.group.sva.agg, data = surstat.sva.crc)
pdf("2.survival.TE.screen/sva.survival/retroposon.pfi.pdf",width = 5)
ggsurvplot(fit2, data =surstat.sva.crc,
                  title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                  pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                  surv.median.line = "hv",            # Add median survival lines
                  legend.title = "group",               # Change legend titles
                  legend.labs = c( "low-risk","high-risk"),  # Change legend labels
                  #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                  palette = c( "#0571b0","#ca0020"),
                  risk.table = TRUE,                  # Add No at risk table
                  cumevents = TRUE,                   # Add cumulative No of events table
                  tables.height = 0.15,               # Specify tables height
                  tables.theme = theme_cleantable(),  # Clean theme for tables
                  tables.y.text = FALSE,             # Hide tables y axis text
                  #ylab="Disease specific survival probability",         #######################change
                  #ylab="Overall survival probability",
                  ylab="progression-free interval probability",
                  #ylab="disease-free interval probability",
                  xlab="Time(days)"
  )
dev.off()
######
save(surstat.sva.crc,file="2.survival.TE.screen/sva.survival/survival.retroposon.RData")

