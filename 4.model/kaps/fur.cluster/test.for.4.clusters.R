kaps.td=cor.immune.drive.sur
kaps.td$kaps.group<- cut(kaps.td$z.of.mean.exp, 
                                     breaks = c(min(kaps.td$z.of.mean.exp),0.33,0.957,1.347,max(kaps.td$z.of.mean.exp)),
                                     labels = c("set1","set2","set3","set4"), include.lowest = TRUE)
kaps.td$kaps.group.three<- cut(kaps.td$z.of.mean.exp, 
                         breaks = c(min(kaps.td$z.of.mean.exp),0.375,0.573,max(kaps.td$z.of.mean.exp)),
                         labels = c("set1","set2","set3"), include.lowest = TRUE)

save(kaps.td,file="4.model/kaps/fur.cluster/kaps.td.RData")
summary(subset(kaps.td, kaps.group=="set4")$z.of.mean.exp)
#####OS
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = kaps.td )
kp1=ggsurvplot(fit2, data =kaps.td,
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
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group, data = kaps.td )
kp2=ggsurvplot(fit2, data =kaps.td,
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
fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group, data = kaps.td )
kp3=ggsurvplot(fit2, data =kaps.td,
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
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data = kaps.td )
kp4=ggsurvplot(fit2, data =kaps.td,
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
  pdf("4.model/kaps/fur.cluster/surC.4.cluster.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
##################heatmap
library(pheatmap)
##
kaps.td=kaps.td[order(kaps.td$z.of.mean.exp,decreasing = T),]
kaps.td$pos=seq(1:nrow(kaps.td))
ppd=as.data.frame(t(kaps.td[,c(3:11)]))
ppd=as.data.frame(t(kaps.td[,c(19:29)]))
#########
kapp=pheatmap(ppd,show_rownames = T,border_color = NA,
         show_colnames = F,fontsize = 4,
         #cluster_rows = hc,
         cluster_cols =F,
         annotation_col = kaps.td[,c(1,13,14,17,38)],
         annotation_row = cndi.rep[rownames(ppd),c(2,3)],
         scale="row",
         #color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8), 
         color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
         annotation_legend = T,
         annotation_colors  = list( TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                    TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                    #MSI.status=c("MSI"="#a65628", "MSI"="#ffd92f","MSS"="#8da0cb"),
                                    MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                    kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027"),
                                    repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                    repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
         )
)
save_pheatmap_pdf(kapp, "4.model/kaps/fur.cluster/heatmap.4.cluster.pdf",height = 3)

######## score distribution
ks=ggplot(kaps.td, aes(x=pos, y=z.of.mean.exp, color=kaps.group)) +
  scale_color_manual(values=c("#00AFBB","#756bb1","#E69F00","#d73027"))+
  geom_vline(xintercept = 51, color = "#d73027", size=0.2)+
  geom_vline(xintercept = 98, color = "#E69F00", size=0.2)+
  geom_vline(xintercept = 211, color = "#756bb1", size=0.2)+
  geom_point(size=0.5) +theme_classic()

generate.PDF <- function(fig) {
  pdf("4.model/kaps/fur.cluster/z.of.mean.exp.dis.pdf",width = 7,height = 2)
  print(ks)
  dev.off()
}
generate.PDF(fig)

table(kaps.td$MSI.status.bin, kaps.td$kaps.group)
table(kaps.td$kaps.group,kaps.td$MSI.status.bin)
chisq.test(kaps.td$MSI.status.bin, kaps.td$kaps.group,correct = T)
###
######subgroup analysis
msisubsur=kaps.td
msisubsur$msi.kaps.sub=paste(msisubsur$MSI.status.bin, msisubsur$kaps.group,sep = "_")
msisubsur$kaps.group.agg=gsub("set1","set_1&2",msisubsur$kaps.group)
msisubsur$kaps.group.agg=gsub("set2","set_1&2",msisubsur$kaps.group.agg)
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = subset(kaps.td, MSI.status.bin=="MSI") )
msi1=ggsurvplot(fit2, data =subset(kaps.td, MSI.status.bin=="MSI"),
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
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group, data = subset(kaps.td, MSI.status.bin=="MSI") )
msi2=ggsurvplot(fit2, data =subset(kaps.td, MSI.status.bin=="MSI"),
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

fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group, data = subset(kaps.td, MSI.status.bin=="MSI") )
msi3=ggsurvplot(fit2, data =subset(kaps.td, MSI.status.bin=="MSI"),
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

fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data = subset(kaps.td, MSI.status.bin=="MSI") )
msi4=ggsurvplot(fit2, data =subset(kaps.td, MSI.status.bin=="MSI"),
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
generate.PDF <- function(fig) {
  pdf("4.model/kaps/fur.cluster/sub.msi.status.4.group//MSI.4.group.pdf",height = 8,width = 5)
  print(msi1)
  print(msi2)
  print(msi3)
  print(msi4)
  dev.off()
}
generate.PDF(fig)
##################################################
####individual compare survival group
######
ress <- pairwise_survdiff(Surv(PFI.time, PFI) ~ kaps.group,
                         data = kaps.td)
ress
#
library(survMisc)
fitss = survfit(Surv(OS.time, OS) ~ kaps.group, data = kaps.td)
comp(ten(fitss))$tests$lrTests

########DSS
kaps.sub.tmp=subset(kaps.td, kaps.group=="set4"|kaps.group=="set3")
table(kaps.sub.tmp$kaps.group)
ss=Surv(kaps.sub.tmp$OS.time, kaps.sub.tmp$OS)

cox = summary(coxph(ss~factor(kaps.sub.tmp$kaps.group,levels = c("set3","set4"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

#
fit1<- survfit(Surv(DSS.time, DSS) ~ pt.clu3, data = surstats.sel)
dss.two=ggsurvplot(fit1, data =surstats.sel,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster",               # Change legend titles
                   legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                   palette = c("#ca0020","#0571b0","#4daf4a"),
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
#################################################
############compare the immune vairables

##############
####pan cancer variables
###
head(immsubtype)
tmp.kaps=msisubsur
rownames(tmp.kaps)=substr(rownames(tmp.kaps),1,12)
pan.id.kaps=intersect(rownames(tmp.kaps), rownames(immsubtype))

CB.data.pan.kaps=cbind(tmp.kaps[pan.id.kaps,], immsubtype[pan.id.kaps,])
#CB.data.pan$TE.cluster.agg=ifelse(CB.data.pan$pt.clu3=="cluster_1","TE.high","TE.low")
colnames(CB.data.pan.kaps)
############
#
index.pan.kaps=colnames(CB.data.pan.kaps)[c(46:73,78:100,103:105)]
my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
CB.data.pan.kaps$kaps.group.agg=factor(CB.data.pan.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.pan.kaps$kaps.group=factor(CB.data.pan.kaps$kaps.group,levels = c("set4","set3","set2","set1"))
for (i in 1:length(index.pan.kaps)) {
  pdf(file = paste0("4.model/kaps/immue.pan.compare/kaps.four/",index.pan[i],".pancancer.kaps.four.pdf"),height  = 4,width = 3)
  print(ggviolin(CB.data.pan.kaps[,-(74:77)],x = "kaps.group", y =index.pan.kaps[i] , fill = "kaps.group",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab(index.pan.kaps[i]) +
          ggtitle(paste0(index.pan.kaps[i],"CRC.pancancer.kaps"))
  )
  dev.off()
}
###############
#######self.genesets.score
###########
head(self.genesets.score)
###
self.id.kaps=intersect(rownames(self.genesets.score), rownames(msisubsur))
#####
CB.data.self.kaps=cbind(msisubsur[self.id.kaps,], self.genesets.score[self.id.kaps,])
colnames(CB.data.self.kaps)
######
index.self.kaps=colnames(CB.data.self.kaps)[c(43:231,309:311)]
my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
CB.data.self.kaps$kaps.group.agg=factor(CB.data.self.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.self.kaps$kaps.group=factor(CB.data.self.kaps$kaps.group,levels = c("set4","set3","set2","set1"))
for (i in 1:length(index.self.kaps)) {
  pdf(file = paste0("4.model/kaps/immune.self.geneset/four.group/",index.self.kaps[i],".self.genesets.kaps.four.pdf"),height  = 4,width = 3)
  print(ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y =index.self.kaps[i] , fill = "kaps.group",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab(index.self.kaps[i]) +
          ggtitle(paste0(index.self.kaps[i],"CRC.self.genesets"))
  )
  dev.off()
}
########
head(HE.stain)
#######
HE.id.kaps=intersect(rownames(HE.stain), rownames(tmp.kaps))
###
CB.data.HE.kaps=cbind(tmp.kaps[HE.id.kaps,], HE.stain[HE.id.kaps,])
head(CB.data.HE.kaps)

table(CB.data.HE.kaps$kaps.group,CB.data.HE.kaps$Immune.Subtype)

index.HE.kaps=colnames(CB.data.HE.kaps)[c(44,48:54,56)]
my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
CB.data.HE.kaps$kaps.group.agg=factor(CB.data.HE.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.HE.kaps$kaps.group=factor(CB.data.HE.kaps$kaps.group,levels = c("set4","set3","set2","set1"))
for (i in 1:length(index.HE.kaps)) {
  pdf(file = paste0("4.model/kaps/HE.staining/four.group/",index.HE.kaps[i],".HE.IHC.kaps.four.pdf"),height  = 4,width = 3)
  print(ggviolin(CB.data.HE.kaps[,-c(46)],x = "kaps.group", y =index.HE.kaps[i] , fill = "kaps.group",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab(index.HE.kaps[i]) +
          ggtitle(paste0("IHC.",index.HE.kaps[i],"CRC"))
  )
  dev.off()
}
#########
tmp.id.kaps=intersect(rownames(tmp.kaps), rownames(ips))
CB.data.ips.kaps=cbind(ips[tmp.id.kaps,51:54], tmp.kaps[tmp.id.kaps,])
head(CB.data.ips.kaps)
##########
index.ips.kaps=colnames(CB.data.ips.kaps)[1:4]
for (i in 1:4) {
  drawdata=CB.data.ips.kaps
  drawdata=drawdata[!(is.na(drawdata[,index.ips.kaps[i]])),]
  drawdata$group=ifelse(drawdata[,index.ips.kaps[i]]>=median(drawdata[,index.ips.kaps[i]]), "ips-high","ips-low")
  ################
  x=as.data.frame.matrix(table(drawdata$kaps.group, drawdata$group))
  x$group=rownames(x)
  head(x)
  write.csv(x, paste0("4.model/kaps/ips/four.group/","stacked.",index.ips.kaps[i],".kaps.four",".csv"))
  library(reshape2)
  library(plyr)
  pval <- chisq.test(drawdata$kaps.group,drawdata$group)$p.value
  #chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
  # long format with column of proportions within each id
  xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
  p1=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
    geom_bar(stat = 'identity')+
    scale_fill_manual(values= c("#41ab5d","#fe9929"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    annotate("text", x=2, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
    ggtitle(paste0(index.ips.kaps[i],"\n","-group as >=median-",median(drawdata[,index.ips.kaps[i]])))
  generate.PDF <- function(fig) {
    pdf(paste0("4.model/kaps/ips/four.group/","stacked.",index.ips.kaps[i],".kaps.four",".pdf"),width = 4,height = 6)
    print(p1)
    dev.off()
  }
  generate.PDF(fig)
}
##########
######
#####ipres
#####
id.ipres.kaps=intersect(rownames(msisubsur), rownames(ipres))
###
CB.data.ipres.kaps=cbind(msisubsur[id.ipres.kaps,], ipres[id.ipres.kaps,])

my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
CB.data.ipres.kaps$kaps.group.agg=factor(CB.data.ipres.kaps$kaps.group.a,levels = c("set4","set3","set_1&2"))
CB.data.ipres.kaps$kaps.group=factor(CB.data.ipres.kaps$kaps.group,levels = c("set4","set3","set2","set1"))
#####
gg1=ggviolin(CB.data.ipres.kaps,x = "kaps.group.agg", y ="z.mean.IPRES" , fill = "kaps.group.agg",alpha = 1,size = 0.3,
             palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
             add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group.agg")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+theme_bw()+xlab("kaps.group.agg")+ylab("z.mean.IPRES") +ggtitle("IPRES.CRC.cohort")
generate.PDF <- function(fig) {
  pdf("4.model/kaps/ipres/ipres.CRC.kaps.agg.pdf",height  = 4,width = 3)
  print(gg1)
  dev.off()
}
generate.PDF(fig)

#############
#######MSI distribution for kapa.four.group
table(kaps.td$MSI.status.bin, kaps.td$TE.cluster.agg)

x=as.data.frame.matrix(table(kaps.td$kaps.group,kaps.td$MSI.status.bin))
x$group=rownames(x)
head(x)
write.csv(x, paste0("4.model/kaps/msi/","MSIvs.kaps.group",".csv"))
library(reshape2)
library(plyr)
pval <- chisq.test(kaps.td$kaps.group,kaps.td$MSI.status.bin ,correct = T)$p.value
#chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
# long format with column of proportions within each id
xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
p1=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
  ggtitle("MSI distribution among four group")

generate.PDF <- function(fig) {
  pdf("4.model/kaps/msi/MSIvs.kaps.group.pdf",width = 4,height = 6)
  print(p1)
  dev.off()
}
generate.PDF(fig)
##########################GEP
head(GEP.stats)
head(kaps.td)
GEP.stats.kaps=cbind(GEP.stats[rownames(kaps.td),]$GEP, kaps.td)
colnames(GEP.stats.kaps)[1]="GEP"
head(GEP.stats.kaps)
#
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
ge1=ggviolin(GEP.stats.kaps,x = "kaps.group", y ="GEP" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.group.kaps")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.group.kaps")+ylab("GEP")+ggtitle("GEP")

generate.PDF <- function(fig) {
  pdf("4.model/kaps/GEP/GEPvs.kaps.group.pdf",width = 4,height = 6)
  print(ge1)
  dev.off()
}
generate.PDF(fig)
#

#
#
#
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group.agg, data = msisubsur )
ggsurvplot(fit2, data =msisubsur,
                #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                surv.median.line = "hv",            # Add median survival lines
                legend.title = "TE.cluster.kaps",               # Change legend titles
                #legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
                #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
                #palette =c("#d73027","#E69F00","#00AFBB"),
                #palette = c("#ca0020","#0571b0","#4daf4a"),
                palette = "lancet",
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


ggviolin(kaps.td,x = "kaps.group", y ="FABP1" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("TE.MSI.group")+ylab("MS4A1")+ggtitle("MS4A1")


ggscatter(cor.immune.drive.sur, x = "z.of.mean.exp", y = "CD274",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black",fill = "lightgray"),#shape = "MSI.status.bin",
          color = "msi.TE.group.bin", palette =  c("#d73027","#00AFBB","#E69F00","#756bb1","#27408B","#FF0000","#2E8B57","#CD00CD","#EE7942","#009ACD")
)+
  stat_cor(method = "spearman", label.x = 2, label.y = 2) +
  geom_vline(xintercept = -0.41, color = "#d73027", size=0.2)+
  geom_hline(yintercept = 6, color = "#d73027", size=0.2)
