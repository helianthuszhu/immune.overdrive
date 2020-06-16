######
drug.prediction=read.csv("4.model/kaps/drug.prediction/durg_response_prediction.csv",header = T,row.names = 1)
head(drug.prediction)
drug.prediction=drug.prediction[,-c(1:2)]
drug.prediction$id=rownames(drug.prediction)
drug.prediction=drug.prediction[grepl("TCGA",drug.prediction$id),]
drug.prediction$id=gsub("[.]","-",drug.prediction$id)
head(drug.prediction)
dim(drug.prediction)
#
aa=kaps.td
aa$id=substr(rownames(aa),1,12)
head(aa)
#
length(intersect(drug.prediction$id, aa$id))       
#####
stat.drug=merge(drug.prediction, aa, by="id", all.y=TRUE)
stat.drug=stat.drug[order(stat.drug$z.of.mean.exp,decreasing = T),]
dim(stat.drug)
##############
head(stat.drug)
table(stat.drug$FOLFIRI_response)
#######
#######
#######calculate the p value
table(stat.drug$FL_response, stat.drug$kaps.group)
datalist=list()
for (i in 2:14) {
  pvalue=chisq.test(stat.drug[,i], stat.drug$kaps.group,correct = T)$p.value
  pvalue=as.data.frame(pvalue)
  colnames(pvalue)="pvalue"
  rownames(pvalue)=colnames(stat.drug)[i]
  datalist[[i]]=pvalue
}
pvalue.drug=do.call(rbind, datalist)
pvalue.drug$pvalue=round(pvalue.drug$pvalue,digits = 4)
##########
#####set color
kaps_col = c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
drug_col=c("with_response"="#ae017e","with_no_response"="#a6d96a","undetermined"="white")
#####
#
ha = HeatmapAnnotation(
  kaps.group=stat.drug$kaps.group,
  FOLFIRI_response=stat.drug$FOLFIRI_response,
  Vandetanib_response=stat.drug$Vandetanib_response,
  Gefitinib_response=stat.drug$Gefitinib_response,
  AZD8931_response=stat.drug$AZD8931_response,
  afatinib_response=stat.drug$afatinib_response,
  cetuximab_response=stat.drug$cetuximab_response,
  Avastin_response=stat.drug$Avastin_response,
  five_FU_response=stat.drug$five_FU_response,
  Oxaliplatin=stat.drug$Oxaliplatin,
  FL_response=stat.drug$FL_response,
  irinotecan=stat.drug$irinotecan,
  five_fu_cl=stat.drug$five_fu_cl,
  BEV_based=stat.drug$BEV_based,
  #
  col = list(
    kaps.group=kaps_col,
    FOLFIRI_response=drug_col,
    Vandetanib_response=drug_col,
    Gefitinib_response=drug_col,
    AZD8931_response=drug_col,
    afatinib_response=drug_col,
    cetuximab_response=drug_col,
    Avastin_response=drug_col,
    five_FU_response=drug_col,
    Oxaliplatin=drug_col,
    FL_response=drug_col,
    irinotecan=drug_col,
    five_fu_cl=drug_col,
    BEV_based=drug_col
  ),
  na_col = "white", border = TRUE,
  show_legend = c(TRUE, TRUE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE, FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    kaps.group=list(title="TE.cluster"),
    FOLFIRI_response=list(title="FOLFIRI"),
    Vandetanib_response=list(title="Vandetanib"),
    Gefitinib_response=list(title="Gefitinib"),
    AZD8931_response=list(title="AZD8931"),
    afatinib_response=list(title="Afatinib"),
    cetuximab_response=list(title="Cetuximab"),
    Avastin_response=list(title="Avastin"),
    five_FU_response=list(title="five_FU"),
    Oxaliplatin=list(title="xaliplatin"),
    FL_response=list(title="FL"),
    irinotecan=list(title="Irinotecan"),
    five_fu_cl=list(title="five_fu_cl"),
    BEV_based=list(title="BEV")
  )
)
zero_row_mat = matrix(nrow = 0, ncol = 590)
ht.drug = Heatmap(zero_row_mat, top_annotation = ha,cluster_columns = F,cluster_rows = F, column_title = "drug prediction among TE clusters")

pdf("4.model/kaps/drug.prediction/ht.drug.prediction.pdf",height = 10,width = 15)
draw(ht.drug, padding = unit(c(22, 80, 20,80), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     #,heatmap_legend_side = "right"
)

annotation_titles2 = c(
  kaps.group=paste("TE.cluster"),
  kaps.group=paste("TE.cluster"),
  FOLFIRI_response=paste("FOLFIRI"),
  Vandetanib_response=paste("Vandetanib"),
  Gefitinib_response=paste("Gefitinib"),
  AZD8931_response=paste("AZD8931"),
  afatinib_response=paste("Afatinib"),
  cetuximab_response=paste("Cetuximab"),
  Avastin_response=paste("Avastin"),
  five_FU_response=paste("five_FU"),
  Oxaliplatin=paste("xaliplatin"),
  FL_response=paste("FL"),
  irinotecan=paste("Irinotecan"),
  five_fu_cl=paste("five_fu_cl"),
  BEV_based=paste("BEV")
)
annotation_titles3 = c(
  kaps.group=paste("TE.cluster"),
  FOLFIRI_response=paste0("p=", pvalue.drug$pvalue[1]),
  Vandetanib_response=paste0("p=", pvalue.drug$pvalue[2]),
  Gefitinib_response=paste0("p=", pvalue.drug$pvalue[3]),
  AZD8931_response=paste0("p=", pvalue.drug$pvalue[4]),
  afatinib_response=paste0("p=", pvalue.drug$pvalue[5]),
  cetuximab_response=paste0("p=", pvalue.drug$pvalue[6]),
  Avastin_response=paste0("p=", pvalue.drug$pvalue[7]),
  five_FU_response=paste0("p=", pvalue.drug$pvalue[8]),
  Oxaliplatin=paste0("p=", pvalue.drug$pvalue[9]),
  FL_response=paste0("p=", pvalue.drug$pvalue[10]),
  irinotecan=paste0("p=", pvalue.drug$pvalue[11]),
  five_fu_cl=paste0("p=", pvalue.drug$pvalue[12]),
  BEV_based=paste0("p=", pvalue.drug$pvalue[13])
)
for(an in names(annotation_titles3)) {
  decorate_annotation(an, {
    grid.text(annotation_titles3[an], unit(222, "mm"), just = "left")
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
##############
#############
####stacked plot
####
pidx.drug=colnames(stat.drug)[2:14]
datalist=list()
for (i in 2:14) {
  x=as.data.frame.matrix(table(stat.drug[,52], stat.drug[,i]))
  x$group=rownames(x)
  #x$pthyName=rep(pidx[i], times=nrow(x))
  xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
  xlong$drugName=rep(pidx.drug[i-1],times=nrow(xlong))
  datalist[[i]]=xlong
}
stacked.drug=do.call(rbind, datalist)
#colnames(stacked.pmu)[1:3]=c("WT","Mut","kaps.group")
#stacked.pmu$group2=paste(stacked.pmu$group,stacked.pmu$pathName,sep = "_")
head(stacked.drug)
stacked.drug$drugName=factor(stacked.drug$drugName,levels = rownames(pvalue.drug))
###############
library(reshape2)
library(plyr)
#pval <- chisq.test(caldata.pmu[,1], caldata.pmu[,i],correct = T)$p.value
#chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
# long format with column of proportions within each id
pdf("4.model/kaps/drug.prediction/stacked.barplot.pdf",width = 20  ,height = 5)
ggplot(stacked.drug, aes(x = group, y = prop, fill = variable)) + 
  geom_bar(stat = 'identity')+
  #cale_fill_manual(values= c("#a6d96a","#c51b8a","#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  scale_fill_manual(values= c("with_response"="#ae017e","with_no_response"="#a6d96a","undetermined"="gray"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotate("text", x=2, y=1, label=paste0("p-value: " ,pvalue.drug$pvalue[1]),size=3)+
  annotate("text", x=2, y=0.95, label=paste0("p-value: ", pvalue.drug$pvalue[2]),size=3)+
  annotate("text", x=2, y=0.9, label=paste0("p-value: ",pvalue.drug$pvalue[3]),size=3)+
  annotate("text", x=2, y=0.85, label=paste0("p-value: ",pvalue.drug$pvalue[4]),size=3)+
  annotate("text", x=2, y=0.8, label=paste0("p-value: ",pvalue.drug$pvalue[5]),size=3)+
  annotate("text", x=2, y=0.75, label=paste0("p-value: ",pvalue.drug$pvalue[6]),size=3)+
  annotate("text", x=2, y=0.7, label=paste0("p-value: ",pvalue.drug$pvalue[7]),size=3)+
  annotate("text", x=2, y=0.65, label=paste0("p-value: ",pvalue.drug$pvalue[8]),size=3)+
  annotate("text", x=2, y=0.6, label=paste0("p-value: ",pvalue.drug$pvalue[9]),size=3)+
  annotate("text", x=2, y=0.55, label=paste0("p-value: ",pvalue.drug$pvalue[10]),size=3)+
  annotate("text", x=2, y=0.5, label=paste0("p-value: ",pvalue.drug$pvalue[11]),size=3)+
  annotate("text", x=2, y=0.45, label=paste0("p-value: ",pvalue.drug$pvalue[12]),size=3)+
  annotate("text", x=2, y=0.4, label=paste0("p-value: ",pvalue.drug$pvalue[13]),size=3)+
  
  ggtitle(paste("drug response prediction"))+xlab("kaps.TE.cluster")+ylab("Proportion")+
  labs(fill = "response")+
  facet_grid(.~ drugName)
dev.off()
save(stacked.drug,pvalue.drug,file="4.model/kaps/drug.prediction/stacked.barplot.RData")
