####
rcc.RERV.expM=read.csv("6.immnue.treat.data/RCC/ERV.exp.RCC.csv",header = T)
rcc.RERV.expM[1:4,1:4]
#
aa=rep.sequence.index
aa=separate(aa, col="idx",into=c("chr","ni","start","end"),sep = "_")
aa$chr=paste("chr",aa$chr)
aa$chr=gsub(" ","",aa$chr)
head(aa)
aa$ids=paste(aa$chr, aa$start, aa$end,sep = ".")
###
length(intersect(rcc.RERV.expM$id, aa$ids))
table(kaps.td$kaps.group)
mantelhaen.test(Table)
#################
#################mRNA
###
rcc.mRNA=read.csv("6.immnue.treat.data/RCC/mRNA.matrix.RCC.csv",header=T)
dim(rcc.mRNA)
rcc.mRNA= rcc.mRNA  %>% group_by(gene_name) %>% summarise_all(mean)
rcc.mRNA=as.data.frame(rcc.mRNA)
rownames(rcc.mRNA)=rcc.mRNA$gene_name
rcc.mRNA=rcc.mRNA[,-1]
rcc.mRNA[1:3,1:4]
write.table(rcc.mRNA, "6.immnue.treat.data/RCC/mRNA.rcc.txt",quote = F,sep = "\t")
#################clin
rcc.clin=read.csv("6.immnue.treat.data/RCC/RCC.clin.csv",header = T)
rcc.clin=rcc.clin %>% drop_na(RNA_ID)
rcc.clin$RNA_ID=gsub("[-]",".",rcc.clin$RNA_ID)
head(rcc.clin)
#
length(intersect(colnames(rcc.mRNA), rcc.clin$RNA_ID))
#####
#####
rcc.clin.stat=data.frame(RNAid=rcc.clin$RNA_ID,
                         TM_CD8_Density=rcc.clin$TM_CD8_Density,
                         TC_CD8_Density=rcc.clin$TC_CD8_Density,
                         TM_TC_Ratio=rcc.clin$TM_TC_Ratio,
                         ImmunoPhenotype=rcc.clin$ImmunoPhenotype,
                         ORR=rcc.clin$ORR,
                         Benefit=rcc.clin$Benefit,
                         PFS.time=rcc.clin$PFS,
                         PFS=rcc.clin$PFS_CNSR,
                         OS.time=rcc.clin$OS,
                         OS=rcc.clin$OS_CNSR)
rownames(rcc.clin.stat)=rcc.clin.stat$RNAid
head(rcc.clin.stat)
#####GEP calculation
#####
######
icb.gene=c("STK11IP","ZBTB34","TBC1D10B","OAZ1","POLR2A","G6PD","ABCF1","NRDE2","UBB","TBP",
           "SDHA","CD27", "CD274", "CD276", "CD8A",
           "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1",
           "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
           "PSMB10", "STAT1","TIGIT")
icb.expM.rcc=as.data.frame(t(rcc.mRNA[icb.gene,]))
head(icb.expM.rcc)
icb.expM.rcc$housexp=rowMeans(icb.expM.rcc[,1:11],na.rm = T)
#icb.expM.kircfh$housexp=rowSums(icb.expM.kircfh[,1:11])
#icb.expM.kircfh$housexp=icb.expM.kircfh$housexp/11
icb.expM.rcc[1:4,]
dim(icb.expM.rcc)
preexp.rcc=icb.expM.rcc[,-c(1:11)]
head(preexp.rcc)
datalist=list()
for (i in 1:17) {
  aa1=preexp.rcc[,c(i,18)]
  aa2=data.frame(value=aa1[,1]-aa1[,2])
  rownames(aa2)=rownames(preexp.rcc)
  colnames(aa2)=colnames(preexp.rcc)[i]
  datalist[[i]]=aa2
}
preexp.cal.rcc=do.call(cbind, datalist)
preexp.cal.rcc[is.na(preexp.cal.rcc)] <- 0
preexp.cal.rcc[1:4,1:4]
#####
coefs=read.table("4.model/immune/GEP.ICB.predictor.txt",header = T,sep = "\t")
rownames(coefs)=coefs$symbol
coefs=coefs[colnames(preexp.cal.rcc),]
coefs$coef=as.numeric(paste(coefs$coef))
#
preexp.cal.rcc$GEP=preexp.cal.rcc[,1]*coefs[1,2]+preexp.cal.rcc[,2]*coefs[2,2]+preexp.cal.rcc[,3]*coefs[3,2]+
  preexp.cal.rcc[,4]*coefs[4,2]+preexp.cal.rcc[,5]*coefs[5,2]+preexp.cal.rcc[,6]*coefs[6,2]+
  preexp.cal.rcc[,7]*coefs[7,2]+preexp.cal.rcc[,8]*coefs[8,2]+preexp.cal.rcc[,9]*coefs[9,2]+
  preexp.cal.rcc[,10]*coefs[10,2]+preexp.cal.rcc[,11]*coefs[11,2]+preexp.cal.rcc[,12]*coefs[12,2]+
  preexp.cal.rcc[,13]*coefs[13,2]+preexp.cal.rcc[,14]*coefs[14,2]+preexp.cal.rcc[,15]*coefs[15,2]+
  preexp.cal.rcc[,16]*coefs[16,2]+preexp.cal.rcc[,17]*coefs[17,2]
head(preexp.cal.rcc)
########
length(intersect(rownames(rcc.clin.stat), rownames(preexp.cal.rcc)))
rcc.clin.gep.stat=cbind(rcc.clin.stat, preexp.cal.rcc[rownames(rcc.clin.stat),])
head(rcc.clin.gep.stat)
########
pdf("6.immnue.treat.data/RCC/cor.gep.vs.TC.cd6.density.pdf",width = 5,height = 5)
ggscatter(rcc.clin.gep.stat, x = "GEP", y = "TC_CD8_Density", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "GEP", ylab = "TC_CD8_Density")+
  ggtitle("GEP vs TC_CD8_Density in RCC (N=103)")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
#########
#########immunescore
####
library(estimate)
rccmrna <- "6.immnue.treat.data/RCC/mRNA.rcc.txt"
filterCommonGenes(input.f=rccmrna, output.f="6.immnue.treat.data/RCC/RCC.mRNA.gct", id="GeneSymbol")
#
estimateScore("6.immnue.treat.data/RCC/RCC.mRNA.gct", "6.immnue.treat.data/RCC/RCC_estimate_score.illumina.gct", platform="illumina")
#######
######
rcc.immunescore=read.table("6.immnue.treat.data/RCC/RCC_estimate_score.illumina.gct",header = F,sep = "\t")
colnames(rcc.immunescore)=rcc.immunescore[2,]
rcc.immunescore=rcc.immunescore[-c(1:2),]
rownames(rcc.immunescore)=rcc.immunescore$NAME
rcc.immunescore=rcc.immunescore[,-c(1,2)]
rcc.immunescore=as.data.frame(t(rcc.immunescore))
head(rcc.immunescore)
######
length(intersect(rownames(rcc.clin.gep.stat), rownames(rcc.immunescore)))
rcc.clin.gep.stat=cbind(rcc.clin.gep.stat, rcc.immunescore[rownames(rcc.clin.gep.stat),])
rcc.clin.gep.stat$ImmuneScore=as.numeric(paste(rcc.clin.gep.stat$ImmuneScore))
rcc.clin.gep.stat$StromalScore=as.numeric(paste(rcc.clin.gep.stat$StromalScore))
rcc.clin.gep.stat$ESTIMATEScore=as.numeric(paste(rcc.clin.gep.stat$ESTIMATEScore))


table(rcc.clin.gep.stat$Benefit)

head(rcc.clin.gep.stat)


ggboxplot(rcc.clin.gep.stat, x = "Benefit", y = "ImmuneScore",
          color = "Benefit", palette = c("CB"="#d01c8b","ICB"="#4dac26","NCB"="#E69F00"),
          add = "jitter")+stat_compare_means()


######
rcc.clin.gep.stat$GEP.0.7.bin= cut(as.numeric(as.character(rcc.clin.gep.stat$GEP)),
                               breaks=quantile(as.numeric(as.character(rcc.clin.gep.stat$GEP)), c(0, 0.7, 1), na.rm=T),labels=c("low", "high"),include.lowest = T)

rcc.clin.gep.stat$GEP.0.5.bin= cut(as.numeric(as.character(rcc.clin.gep.stat$GEP)),
                               breaks=quantile(as.numeric(as.character(rcc.clin.gep.stat$GEP)), c(0, 0.5, 1), na.rm=T),labels=c("low", "high"),include.lowest = T)

table(rcc.clin.gep.stat$GEP.bin)
#########
#########
pdf("6.immnue.treat.data/RCC/survival.GEP.bin.pdf",width = 5,height = 5)
fit2<- survfit(Surv(OS.time, OS) ~ GEP.0.7.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "GEP.0.7.bin",               # Change legend titles
               legend.labs = c("low","high"),  # Change legend labels
               palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
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
               xlab="Time(months)"
)
#####
fit2<- survfit(Surv(PFS.time, PFS) ~ GEP.0.7.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "GEP.0.7.bin",               # Change legend titles
           legend.labs = c("low","high"),  # Change legend labels
           palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
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
           #ylab="disease-free interval probability",
           ylab="progression-free survival probability",
           xlab="Time(months)"
)
####
fit2<- survfit(Surv(OS.time, OS) ~ GEP.0.5.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "GEP.0.5.bin",               # Change legend titles
           legend.labs = c("low","high"),  # Change legend labels
           palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
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
           xlab="Time(months)"
)
#####
fit2<- survfit(Surv(PFS.time, PFS) ~ GEP.0.5.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "GEP.0.5.bin",               # Change legend titles
           legend.labs = c("low","high"),  # Change legend labels
           palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
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
           #ylab="disease-free interval probability",
           ylab="progression-free survival probability",
           xlab="Time(months)"
)
dev.off()
###########
###########
pdf("6.immnue.treat.data/RCC/cor.gep.vs.immunescore.pdf",width = 5,height = 5)
ggscatter(rcc.clin.gep.stat, x = "GEP", y = "ImmuneScore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "GEP", ylab = "ImmuneScore")+
  ggtitle("GEP vs ImmuneScore in RCC (N=311)")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()

rcc.clin.gep.stat$ImmuneScore.0.5.bin= cut(as.numeric(as.character(rcc.clin.gep.stat$ImmuneScore)),
                                          breaks=quantile(as.numeric(as.character(rcc.clin.gep.stat$ImmuneScore)), c(0, 0.5, 1), na.rm=T),
                                          labels=c("low", "high"),include.lowest = T)

rcc.clin.gep.stat$ImmuneScore.0.7.bin= cut(as.numeric(as.character(rcc.clin.gep.stat$ImmuneScore)),
                                       breaks=quantile(as.numeric(as.character(rcc.clin.gep.stat$ImmuneScore)), c(0, 0.7, 1), na.rm=T),
                                       labels=c("low", "high"),include.lowest = T)

table(rcc.clin.gep.stat$ImmuneScore.bin)
#
pdf("6.immnue.treat.data/RCC/survival.immunescore.bin.pdf",width = 5,height = 5)
fit2<- survfit(Surv(OS.time, OS) ~ ImmuneScore.0.7.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "ImmuneScore.0.7.bin",               # Change legend titles
           legend.labs = c("low","high"),  # Change legend labels
           palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
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
           xlab="Time(months)"
)
#
fit2<- survfit(Surv(PFS.time, PFS) ~ ImmuneScore.0.7.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "ImmuneScore.0.7.bin",               # Change legend titles
           legend.labs = c("low","high"),  # Change legend labels
           palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
           #palette =c("#d73027","#E69F00","#00AFBB"),
           #palette = c("#ca0020","#0571b0","#4daf4a"),
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE,             # Hide tables y axis text
           #ylab="Disease specific survival probability",
           #ylab="Overall survival probability",
           ylab="progression-free survival probability",
           #ylab="disease-free interval probability",
           xlab="Time(months)"
)
######
fit2<- survfit(Surv(OS.time, OS) ~ ImmuneScore.0.5.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "ImmuneScore.0.5.bin",               # Change legend titles
           legend.labs = c("low","high"),  # Change legend labels
           palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
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
           xlab="Time(months)"
)
#
fit2<- survfit(Surv(PFS.time, PFS) ~ ImmuneScore.0.5.bin, data = rcc.clin.gep.stat )
ggsurvplot(fit2, data =rcc.clin.gep.stat,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "ImmuneScore.0.5.bin",               # Change legend titles
           legend.labs = c("low","high"),  # Change legend labels
           palette = c("#00AFBB","#756bb1"),  # Use JCO journal color palette,
           #palette =c("#d73027","#E69F00","#00AFBB"),
           #palette = c("#ca0020","#0571b0","#4daf4a"),
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE,             # Hide tables y axis text
           #ylab="Disease specific survival probability",
           #ylab="Overall survival probability",
           ylab="progression-free survival probability",
           #ylab="disease-free interval probability",
           xlab="Time(months)"
)
dev.off()
####
save(rcc.clin.gep.stat,file="6.immnue.treat.data/RCC/rcc.clin.gep.immunescore.stat.RData")
