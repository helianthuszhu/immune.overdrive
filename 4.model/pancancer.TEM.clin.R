load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.PCA/pca.data.RData")
#########survival screen
surpan=read.csv("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/clinical.data/pan/pan.sur.clinical.full.info.csv",header = T)
rownames(surpan)=surpan$bcr_patient_barcode
head(surpan)
dim(surpan)
######
cc=rclin.agg.tumor
rownames(cc)=cc$patient_ID
head(cc)
ids=intersect(rownames(cc), rownames(surpan))
###
stat.surclin=cbind(surpan[ids,], cc[ids,])
rownames(stat.surclin)=stat.surclin$id
head(stat.surclin)
###
stat.rexp=rexp.agg.T[rownames(stat.surclin),]
dim(stat.rexp)
dim(stat.surclin)
####
save(stat.rexp, stat.surclin,file = "~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/pancancer.TEM.clin.RData")
#########
######
#########
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/model.cor.>=0.4.totoal.surTEs.RData")
pan.exp.8s=stat.rexp[,colnames(stat.rexp) %in%  rownames(cndi.rep)]
pan.exp.8s[1:4,]
###
head(stat.surclin)
table(stat.surclin$type)
stat.surclin.KIRC=subset(stat.surclin, type=="KIRC")
#
#stat.surclin.KIRC=subset(stat.surclin, type=="HNSC")

stat.exp.KIRC=pan.exp.8s[rownames(stat.surclin.KIRC),]
#
pkirc=pheatmap(t(stat.exp.KIRC),show_rownames = T,border_color = NA,
         show_colnames = F,
         #cluster_rows = hc,
         #cluster_cols =hc,
         #annotation_col = ann.col,
         scale="row",
         color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256), 
         annotation_legend = T
         #annotation_colors  = list(response = c("PD" = "#dd1c77","PRCR"="#74c476","SD"="#feb24c"),
         #                         treatment.status=c("On"="#74a9cf", "Pre"="#f768a1"))
)
clu.KIRC=as.data.frame(sort(cutree(pkirc$tree_col, k=4)))   #change
colnames(clu.KIRC)="pt.clu3"
clu.KIRC$pt.clu3=paste("cluster",clu.KIRC$pt.clu3,sep = "_")
clu.KIRC$idss=rownames(clu.KIRC)
table(clu.KIRC$pt.clu3)
head(clu.KIRC)
length(intersect(rownames(clu.KIRC), rownames(stat.surclin.KIRC)))
dim(clu.KIRC)
clu.KIRC=cbind(clu.KIRC, stat.surclin.KIRC[rownames(clu.KIRC),])
table(clu.KIRC$pt.clu3)
ann.col.KIRC=data.frame(type=clu.KIRC$type,
                        cluster=clu.KIRC$pt.clu3)
rownames(ann.col.KIRC)=rownames(clu.KIRC)

head(ann.col.KIRC)
#ann.col.KIRC$cluster.new=ifelse(ann.col.KIRC$cluster=="cluster_4","cluster_2","cluster_other")
ann.col.KIRC$cluster.new=gsub("cluster_4","cluster_2",ann.col.KIRC$cluster)
ann.col.KIRC$cluster.new=gsub("cluster_3","cluster_1",ann.col.KIRC$cluster.new)
ann.col.KIRC$TE.cluster.agg=ifelse(ann.col.KIRC$cluster.new=="cluster_1","TE.low","TE.high")
#
pkirc=pheatmap(t(stat.exp.KIRC),show_rownames = T,border_color = NA,
               show_colnames = F,fontsize = 7,
               #cluster_rows = hc,
               #cluster_cols =hc,
               annotation_col = ann.col.KIRC[,-3],
               annotation_row = cndi.rep[,c(2,3)],
               scale="row",
               color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256), 
               annotation_legend = T,
               annotation_colors  = list(cluster= c("cluster_1" ="#00AFBB" ,"cluster_2"="#d73027","cluster_3"="#E69F00","cluster_4"="#31a354"),
                                         TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                         #,MSI.status=c("MSI-H"="#a65628", "MSI-L"="#ffd92f","MSS"="#8da0cb"),
                                         repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                         repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
                                    )
)

save_pheatmap_pdf(pkirc, "4.model/heatmpa.9.TEs.kirc.agg.pdf",height = 4)
#
head(clu.KIRC)
sur.kirc=cbind(ann.col.KIRC$TE.cluster.agg, clu.KIRC)
colnames(sur.kirc)[1]="TE.cluster.agg"

sur.kirc$DSS.time=as.numeric(paste(sur.kirc$DSS.time))
sur.kirc$DSS=as.numeric(paste(sur.kirc$DSS))
sur.kirc$OS.time=as.numeric(paste(sur.kirc$OS.time))
sur.kirc$OS=as.numeric(paste(sur.kirc$OS))
#
sur.kirc$DFI.time=as.numeric(paste(sur.kirc$DFI.time))
sur.kirc$DFI=as.numeric(paste(sur.kirc$DFI))
sur.kirc$PFI.time=as.numeric(paste(sur.kirc$PFI.time))
sur.kirc$PFI=as.numeric(paste(sur.kirc$PFI))
head(sur.kirc)

#######
surstats.sel.kirc=subset(sur.kirc, pt.clu3=="cluster_1" | pt.clu3=="cluster_2")
########DSS
ss=Surv(surstats.sel.kirc$DSS.time, surstats.sel.kirc$DSS)
cox = summary(coxph(ss~factor(surstats.sel.kirc$pt.clu3,levels = c("cluster_1","cluster_2"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(DSS.time, DSS) ~ pt.clu3, data = surstats.sel.kirc)
dss.two.kirc=ggsurvplot(fit1, data =surstats.sel.kirc,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster",               # Change legend titles
                   legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                   palette = c("#0571b0","#ca0020","#4daf4a"),
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
#########
fit2<- survfit(Surv(DSS.time, DSS) ~ pt.clu3, data = sur.kirc )
dss.three.kirc=ggsurvplot(fit2, data =sur.kirc,
                     #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                     pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                     surv.median.line = "hv",            # Add median survival lines
                     legend.title = "TE.cluster",               # Change legend titles
                     legend.labs = c("cluster_1","cluster_2","cluster_3","cluster_4"),  # Change legend labels
                     #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                     palette = c("#0571b0","#ca0020","#4daf4a","#984ea3"),
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
#
ss=Surv(sur.kirc$DSS.time, sur.kirc$DSS)
cox = summary(coxph(ss~factor(sur.kirc$TE.cluster.agg,levels = c("TE.low","TE.high"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(DSS.time, DSS) ~ TE.cluster.agg, data = sur.kirc )
dss.kirc.agg=ggsurvplot(fit2, data =sur.kirc,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster.agg",               # Change legend titles
                   legend.labs = c("TE.high","TE.low"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                   palette =c("#d73027","#E69F00","#00AFBB"),
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

generate.PDF <- function(fig) {
  pdf("4.model/surC.dss.kirc.agg.pdf",width = 6,height = 8)
  print(dss.two.kirc)
  print(dss.three.kirc)
  print(dss.kirc.agg)
  
  dev.off()
}
generate.PDF(fig)
#################
########OS
ss=Surv(surstats.sel.kirc$OS.time, surstats.sel.kirc$OS)
cox = summary(coxph(ss~factor(surstats.sel.kirc$pt.clu3,levels = c("cluster_1","cluster_2"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(OS.time, OS) ~ pt.clu3, data = surstats.sel.kirc)
os.two.kirc=ggsurvplot(fit1, data =surstats.sel.kirc,
                        title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                        pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                        surv.median.line = "hv",            # Add median survival lines
                        legend.title = "TE.cluster",               # Change legend titles
                        legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                        #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                        palette = c("#0571b0","#ca0020","#4daf4a"),
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
#########
fit2<- survfit(Surv(OS.time, OS) ~ pt.clu3, data = sur.kirc )
os.three.kirc=ggsurvplot(fit2, data =sur.kirc,
                          #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                          pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                          surv.median.line = "hv",            # Add median survival lines
                          legend.title = "TE.cluster",               # Change legend titles
                          legend.labs = c("cluster_1","cluster_2","cluster_3","cluster_4"),  # Change legend labels
                          #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                          palette = c("#0571b0","#ca0020","#4daf4a","#984ea3"),
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
##
ss=Surv(sur.kirc$OS.time, sur.kirc$OS)
cox = summary(coxph(ss~factor(sur.kirc$TE.cluster.agg,levels = c("TE.low","TE.high"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(OS.time, OS) ~ TE.cluster.agg, data = sur.kirc )
os.kirc.agg=ggsurvplot(fit2, data =sur.kirc,
                        title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                        pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                        surv.median.line = "hv",            # Add median survival lines
                        legend.title = "TE.cluster.agg",               # Change legend titles
                        legend.labs = c("TE.high","TE.low"),  # Change legend labels
                        #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                        palette =c("#d73027","#E69F00","#00AFBB"),
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


generate.PDF <- function(fig) {
  pdf("4.model/surC.os.kirc.agg.pdf",width = 6,height = 8)
  print(os.two.kirc)
  print(os.three.kirc)
  print(os.kirc.agg)
  
  dev.off()
}
generate.PDF(fig)
###########
#####DFI
ss=Surv(surstats.sel.kirc$DFI.time, surstats.sel.kirc$DFI)
cox = summary(coxph(ss~factor(surstats.sel.kirc$pt.clu3,levels = c("cluster_1","cluster_2"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(DFI.time, DFI) ~ pt.clu3, data = surstats.sel.kirc)
dfi.two.kirc=ggsurvplot(fit1, data =surstats.sel.kirc,
                       title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                       pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                       surv.median.line = "hv",            # Add median survival lines
                       legend.title = "TE.cluster",               # Change legend titles
                       legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                       #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                       palette = c("#0571b0","#ca0020","#4daf4a"),
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
#########
fit2<- survfit(Surv(DFI.time, DFI) ~ pt.clu3, data = sur.kirc )
dfi.three.kirc=ggsurvplot(fit2, data =sur.kirc,
                         #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                         pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                         surv.median.line = "hv",            # Add median survival lines
                         legend.title = "TE.cluster",               # Change legend titles
                         legend.labs = c("cluster_1","cluster_2","cluster_3","cluster_4"),  # Change legend labels
                         #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                         palette = c("#0571b0","#ca0020","#4daf4a","#984ea3"),
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
##
ss=Surv(sur.kirc$DFI.time, sur.kirc$DFI)
cox = summary(coxph(ss~factor(sur.kirc$TE.cluster.agg,levels = c("TE.low","TE.high"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(DFI.time, DFI) ~ TE.cluster.agg, data = sur.kirc )
dfi.kirc.agg=ggsurvplot(fit2, data =sur.kirc,
                       title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                       pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                       surv.median.line = "hv",            # Add median survival lines
                       legend.title = "TE.cluster.agg",               # Change legend titles
                       legend.labs = c("TE.high","TE.low"),  # Change legend labels
                       #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                       palette =c("#d73027","#E69F00","#00AFBB"),
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

generate.PDF <- function(fig) {
  pdf("4.model/surC.dfi.kirc.agg.pdf",width = 6,height = 8)
  print(dfi.two.kirc)
  print(dfi.three.kirc)
  print(dfi.kirc.agg)
  
  dev.off()
}
generate.PDF(fig)
#########
#PFI
ss=Surv(surstats.sel.kirc$PFI.time, surstats.sel.kirc$PFI)
cox = summary(coxph(ss~factor(surstats.sel.kirc$pt.clu3,levels = c("cluster_1","cluster_2"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(PFI.time, PFI) ~ pt.clu3, data = surstats.sel.kirc)
pfi.two.kirc=ggsurvplot(fit1, data =surstats.sel.kirc,
                        title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                        pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                        surv.median.line = "hv",            # Add median survival lines
                        legend.title = "TE.cluster",               # Change legend titles
                        legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                        #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                        palette = c("#0571b0","#ca0020","#4daf4a"),
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
#########
fit2<- survfit(Surv(PFI.time, PFI) ~ pt.clu3, data = sur.kirc )
pfi.three.kirc=ggsurvplot(fit2, data =sur.kirc,
                          #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                          pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                          surv.median.line = "hv",            # Add median survival lines
                          legend.title = "TE.cluster",               # Change legend titles
                          legend.labs = c("cluster_1","cluster_2","cluster_3","cluster_4"),  # Change legend labels
                          #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                          palette = c("#0571b0","#ca0020","#4daf4a","#984ea3"),
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
ss=Surv(sur.kirc$PFI.time, sur.kirc$PFI)
cox = summary(coxph(ss~factor(sur.kirc$TE.cluster.agg,levels = c("TE.low","TE.high"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(PFI.time, PFI) ~ TE.cluster.agg, data = sur.kirc )
pfi.kirc.agg=ggsurvplot(fit2, data =sur.kirc,
                        title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                        pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                        surv.median.line = "hv",            # Add median survival lines
                        legend.title = "TE.cluster.agg",               # Change legend titles
                        legend.labs = c("TE.high","TE.low"),  # Change legend labels
                        #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                        palette =c("#d73027","#E69F00","#00AFBB"),
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
generate.PDF <- function(fig) {
  pdf("4.model/surC.pfi.kirc.agg.pdf",width = 6,height = 8)
  print(pfi.two.kirc)
  print(pfi.three.kirc)
  print(pfi.kirc.agg)
  
  dev.off()
}
generate.PDF(fig)
#########################
##################compare the GEP icb score
#
load("~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/GEP.ICB.score.pancancer.RData")
######
kirc.gep=ann.col.KIRC
kirc.gep$id=rownames(kirc.gep)
kirc.gep$idd=substr(kirc.gep$id, 1,12)
rownames(kirc.gep)=kirc.gep$idd
head(kirc.gep)
id.kirc=intersect(rownames(kirc.gep), rownames(preexp.cal.pancancer))
#####
kirc.gep.stats=cbind(kirc.gep[id.kirc,], preexp.cal.pancancer[id.kirc,])
kirc.gep.stats$TE.cluster.agg=ifelse(kirc.gep.stats$cluster.new=="cluster_1","TE.low","TE.high")
kirc.gep.stats$TE.cluster.agg=factor(kirc.gep.stats$TE.cluster.agg,levels = c("TE.high","TE.low"))
head(kirc.gep.stats)

##
require(ggplot2)
require(reshape2)
library(ggpubr)
kirc.gep.stats$cluster=factor(kirc.gep.stats$cluster,levels = c("cluster_2","cluster_4","cluster_1","cluster_3"))
my_comparisons <- list( c("cluster_1", "cluster_2"), c("cluster_1", "cluster_3"), c("cluster_2", "cluster_3"),
                        c("cluster_1", "cluster_4"),c("cluster_2", "cluster_4",c("cluster_3", "cluster_4")))

gep1=ggviolin(kirc.gep.stats,x = "cluster", y ="GEP" , fill = "cluster",alpha = 1,size = 0.3,
              palette = c("#d73027","#00AFBB", "#E69F00","#9970ab"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+theme_bw()+xlab("TE.cluster")+ylab("GEP") 

gep2=ggviolin(kirc.gep.stats,x = "TE.cluster.agg", y ="GEP" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
              palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("GEP")+ggtitle("KIRC.combined.two.clusters")

gep3=ggviolin(subset(kirc.gep.stats, cluster=="cluster_1" |cluster=="cluster_2" ),x = "cluster", y ="GEP" , fill = "cluster",alpha = 1,size = 0.3,
              palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster")+ylab("GEP") +ggtitle("KIRC.just.two.clusters")

generate.PDF <- function(fig) {
  pdf("4.model/immune/GEP/GEP.clusters.KIRC.agg.pdf",height  = 4,width = 3)
  print(gep1)
  print(gep2)
  print(gep3)
  dev.off()
}
generate.PDF(fig)
####
#################compare the averaage expression of 9 tes
#######
dim(stat.exp.KIRC)
mean.stats.kirc=cbind( stat.exp.KIRC[rownames(sur.kirc),],sur.kirc)
head(mean.stats.kirc)
mean.stats.kirc$mean.exp=rowSums(mean.stats.kirc[,1:9])/9
#####
my_comparisons <- list( c("cluster_1", "cluster_2"), c("cluster_1", "cluster_3"), c("cluster_2", "cluster_3"),
                        c("cluster_1", "cluster_4"),c("cluster_2", "cluster_4",c("cluster_3", "cluster_4")))
mean.stats.kirc$pt.clu3=factor(mean.stats.kirc$pt.clu3,levels = c("cluster_4","cluster_2","cluster_1","cluster_3"))
gepm1=ggviolin(mean.stats.kirc,x = "pt.clu3", y ="mean.exp" , fill = "pt.clu3",alpha = 1,size = 0.3,
              palette = c("#d73027","#00AFBB", "#E69F00","#31a354"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+theme_bw()+xlab("TE.cluster")+ylab("mean.expression.of.9.TEs")

gepm2=ggviolin(mean.stats.kirc,x = "TE.cluster.agg", y ="mean.exp" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
               palette = c("#d73027","#E69F00","#00AFBB","#31a354"),
               add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("mean.expression.of.9.TEs")+ggtitle("KIRC.mean.TE.expression")
########
generate.PDF <- function(fig) {
  pdf("4.model/Mean.expression.of.9.TE.model.agg.pdf",height  = 4,width = 3)
  print(gepm1)
  print(gepm2)
  dev.off()
}
generate.PDF(fig)
#####################
###########compare the immune vairbale from pancancer
###########
tmp.clu.kirc=ann.col.KIRC
rownames(tmp.clu.kirc)=substr(rownames(tmp.clu.kirc),1,12)
####

pan.id.kirc=intersect(rownames(tmp.clu.kirc), rownames(immsubtype))
CB.data.pan.kirc=cbind(tmp.clu.kirc[pan.id.kirc,], immsubtype[pan.id.kirc,])
CB.data.pan.kirc$TE.cluster.agg=factor(CB.data.pan.kirc$TE.cluster.agg,levels = c("TE.high","TE.low"))
#CB.data.pan$TE.cluster.agg=ifelse(CB.data.pan$pt.clu3=="cluster_1","TE.high","TE.low")
colnames(CB.data.pan.kirc)

#
index.pan.kirc=colnames(CB.data.pan.kirc)[c(8:35,40:62,65:67)]
for (i in 1:length(index.pan.kirc)) {
  pdf(file = paste0("4.model/immune/immune.pan/pan.kirc/",index.pan.kirc[i],".pancancer.kirc.pdf"),height  = 4,width = 3)
  print(ggviolin(CB.data.pan.kirc,x = "TE.cluster.agg", y =index.pan[i] , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
          stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab(index.pan.kirc[i]) +
          ggtitle(paste0(index.pan[i],".KIRC.pancancer"))
  )
  dev.off()
}
#############correlation of GEP and mean TE exp in KIRC
head(mean.stats.kirc)
head(kirc.gep.stats)

tmps.kirc=mean.stats.kirc
rownames(tmps.kirc)=substr(rownames(tmps.kirc),1,12)

length(intersect(rownames(kirc.gep.stats),rownames(tmps.kirc)))

cor.GEP.te.mean.kirc=cbind(kirc.gep.stats[,24:25],tmps.kirc[rownames(kirc.gep.stats),])
head(cor.GEP.te.mean.kirc)

gg1=ggscatter(cor.GEP.te.mean.kirc[,-2], x = "mean.exp", y = "GEP", 
              add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
              cor.coef = T, cor.method = "spearman",
              xlab = "mean.TE.score", ylab = "GEP")+
  ggtitle("mean.TE.exp vs GEP in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#annotate(geom="text", x=quantile(cd$score)[1], y=quantile(cd$exp)[4], label=paste0("pcor =",coef, "\n", "pvalue =",pvalue),color="red")
generate.PDF <- function(fig) {
  pdf("4.model/immune/GEP.vs.mean.TE.exp/cor.mean.TE.with.GEP.kirc.pdf",width = 7,height = 7)
  print(gg1)
  dev.off()
}
generate.PDF(fig)
