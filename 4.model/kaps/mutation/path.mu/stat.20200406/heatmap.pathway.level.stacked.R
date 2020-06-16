####oncogenic pathways gene mutation
########
pmu=read.table("4.model/kaps/mutation/path.mu/ten.paths.pathway.level.alterd.txt",header = T,sep = "\t",row.names = 1)
head(pmu)
sh.id.pmu=intersect(rownames(kaps.td),rownames(pmu))
caldata.pmu=cbind(kaps.td[sh.id.pmu,]$kaps.group,pmu[sh.id.pmu,])
colnames(caldata.pmu)[1]="kaps.group"
head(caldata.pmu)
###calculate p value
###
dfm.pmu=NULL
for (k in 2:ncol(caldata.pmu)){
  p.value=chisq.test(caldata.pmu[,k],caldata.pmu$kaps.group,correct = T)
  val=as.data.frame(p.value$p.value)
  rownames(val)=colnames(caldata.pmu)[k]
  dfm.pmu=rbind(dfm.pmu,val)
}
colnames(dfm.pmu)="pvalue"
dfm.pmu$pvalue=round(dfm.pmu, digits = 4)
dfm.pmu=as.data.frame(dfm.pmu)
head(dfm.pmu)
dim(dfm.pmu)
class(dfm.pmu)
###
caldata.pmu[is.na(caldata.pmu)] <- "not_available"
#####set the color
kaps_col = c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
#snv.col=c("1" = "#E41A1C", "0" = "#377EB8","not_available"="white")

snv.col=c("1" = "#c51b8a", "0" = "#a6d96a","not_available"="white")
####
ha = HeatmapAnnotation(
  kaps.group=caldata.pmu$kaps.group,
  Cell.Cycle=caldata.pmu$Cell.Cycle,
  HIPPO=caldata.pmu$HIPPO,
  MYC  =caldata.pmu$MYC,
  NOTCH  =caldata.pmu$NOTCH,
  NRF2   =caldata.pmu$NRF2,
  PI3K   =caldata.pmu$PI3K,
  RTK.RAS =caldata.pmu$RTK.RAS,
  TP53    =caldata.pmu$TP53,
  TGF.Beta =caldata.pmu$TGF.Beta,
  WNT       =caldata.pmu$WNT,
  #
  col = list(
    kaps.group=kaps_col,
    Cell.Cycle=snv.col,
    HIPPO   =snv.col,
    MYC     =snv.col,
    NOTCH   =snv.col,
    NRF2    =snv.col,
    PI3K    =snv.col,
    RTK.RAS  =snv.col,
    TP53     =snv.col,
    TGF.Beta =snv.col,
    WNT      =snv.col
    
  ),
  na_col = "white", border = TRUE,
  show_legend = c(TRUE, TRUE, FALSE, FALSE, FALSE,FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    kaps.group=list(title="TE.cluster"),
    Cell.Cycle =list(title="Cell.Cycle"),
    HIPPO      =list(title="HIPPO"),
    MYC        =list(title="MYC"),
    NOTCH      =list(title="NOTCH"),
    NRF2       =list(title="NRF2"),
    PI3K       =list(title="PI3K"),
    RTK.RAS    =list(title="RTK.RAS"),
    TP53       =list(title="TP53"),
    TGF.Beta   =list(title="TGF.Beta"),
    WNT        =list(title="WNT")
  )
)
#dim(caldata.pmu)
zero_row_mat = matrix(nrow = 0, ncol = nrow(caldata.pmu))
ht = Heatmap(zero_row_mat, top_annotation = ha,cluster_columns = F,cluster_rows = F, column_title = "pathway level mutation among TE clusters")
pdf("4.model/kaps/mutation/path.mu/stat.20200406/heatmap.pathway.lelvel.pdf",height = 10,width = 15)
draw(ht, padding = unit(c(22, 80, 20,60), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     #,heatmap_legend_side = "right"
)

annotation_titles2 = c(
  kaps.group=paste("TE.cluster"),
  Cell.Cycle =paste("Cell.Cycle"),
  HIPPO      =paste("HIPPO"),
  MYC        =paste("MYC"),
  NOTCH      =paste("NOTCH"),
  NRF2       =paste("NRF2"),
  PI3K       =paste("PI3K"),
  RTK.RAS    =paste("RTK.RAS"),
  TP53       =paste("TP53"),
  TGF.Beta   =paste("TGF.Beta"),
  WNT        =paste("WNT")
)
annotation_titles3 = c(
  kaps.group=paste("TE.cluster"),
  Cell.Cycle =paste0("p=", dfm.pmu$pvalue$pvalue[1]),
  HIPPO      =paste0("p=", dfm.pmu$pvalue$pvalue[2]),
  MYC        =paste0("p=", dfm.pmu$pvalue$pvalue[3]),
  NOTCH      =paste0("p=", dfm.pmu$pvalue$pvalue[4]),
  NRF2       =paste0("p=", dfm.pmu$pvalue$pvalue[5]),
  PI3K       =paste0("p=", dfm.pmu$pvalue$pvalue[6]),
  RTK.RAS    =paste0("p=", dfm.pmu$pvalue$pvalue[7] ),
  TP53       =paste0("p=", dfm.pmu$pvalue$pvalue[8]),
  TGF.Beta   =paste0("p=", dfm.pmu$pvalue$pvalue[9]),
  WNT        =paste0("p=", dfm.pmu$pvalue$pvalue[10])
)
for(an in names(annotation_titles3)) {
  decorate_annotation(an, {
    grid.text(annotation_titles3[an], unit(243, "mm"), just = "left")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
for(an in names(annotation_titles2)) {
  decorate_annotation(an, {
    grid.text(annotation_titles2[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}


dev.off()
save(caldata.pmu, dfm.pmu, file="4.model/kaps/mutation/path.mu/stat.20200406/heatmap.pathway.lelvel.RData")
#########
#########draw using stacked bar plot
#######
table(caldata.pmu$kaps.group, caldata.pmu$Cell.Cycle)
pidx=colnames(caldata.pmu)
datalist=list()
for (i in 2:ncol(caldata.pmu)) {
  x=as.data.frame.matrix(table(caldata.pmu[,1], caldata.pmu[,i]))[,-3]
  x$group=rownames(x)
  #x$pthyName=rep(pidx[i], times=nrow(x))
  xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
  xlong$pathName=rep(pidx[i],times=nrow(xlong))
  datalist[[i]]=xlong
}
stacked.pmu=do.call(rbind, datalist)
#colnames(stacked.pmu)[1:3]=c("WT","Mut","kaps.group")
#stacked.pmu$group2=paste(stacked.pmu$group,stacked.pmu$pathName,sep = "_")
head(stacked.pmu)
###############
library(reshape2)
library(plyr)
#pval <- chisq.test(caldata.pmu[,1], caldata.pmu[,i],correct = T)$p.value
#chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
# long format with column of proportions within each id
pdf("4.model/kaps/mutation/path.mu/stat.20200406/stacked.barplot.pdf",width = 10  ,height = 5)
ggplot(stacked.pmu, aes(x = group, y = prop, fill = variable)) + 
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#a6d96a","#c51b8a","#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotate("text", x=2, y=1, label=paste0("p-value: " ,dfm.pmu$pvalue$pvalue[1]),size=2)+
  annotate("text", x=2, y=0.9, label=paste0("p-value: ", dfm.pmu$pvalue$pvalue[2]),size=2)+
  annotate("text", x=2, y=0.8, label=paste0("p-value: ",dfm.pmu$pvalue$pvalue[3]),size=2)+
  annotate("text", x=2, y=0.7, label=paste0("p-value: ",dfm.pmu$pvalue$pvalue[4]),size=2)+
  annotate("text", x=2, y=0.6, label=paste0("p-value: ",dfm.pmu$pvalue$pvalue[5]),size=2)+
  annotate("text", x=2, y=0.5, label=paste0("p-value: ",dfm.pmu$pvalue$pvalue[6]),size=2)+
  annotate("text", x=2, y=0.4, label=paste0("p-value: ",dfm.pmu$pvalue$pvalue[7]),size=2)+
  annotate("text", x=2, y=0.3, label=paste0("p-value: " ,dfm.pmu$pvalue$pvalue[8]),size=2)+
  annotate("text", x=2, y=0.2, label=paste0("p-value: ",dfm.pmu$pvalue$pvalue[9]),size=2)+
  annotate("text", x=2, y=0.1, label=paste0("p-value: " ,dfm.pmu$pvalue$pvalue[10]),size=2)+
  ggtitle(paste("pathway level mutation"))+xlab("kaps.TE.cluster")+ylab("Proportion")+
  labs(fill = "mutation")+
  facet_grid(.~ pathName)
dev.off()
save(stacked.pmu, dfm.pmu, file="4.model/kaps/mutation/path.mu/stat.20200406/stacked.barplot.RData")
######
