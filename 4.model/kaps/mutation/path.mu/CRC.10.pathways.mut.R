####oncogenic pathways gene mutation
########
pmu=read.table("4.model/kaps/mutation/path.mu/ten.paths.pathway.level.alterd.txt",header = T,sep = "\t",row.names = 1)
head(pmu)
#
sh.id.pmu=intersect(rownames(kaps.td),rownames(pmu))
caldata.pmu=cbind(kaps.td[sh.id.pmu,]$kaps.group,pmu[sh.id.pmu,])
colnames(caldata.pmu)[1]="kaps.group"
head(caldata.pmu)
#
dfm.pmu=NULL
for (k in 2:ncol(caldata.pmu)){
  p.value=chisq.test(caldata.pmu[,k],caldata.pmu$kaps.group,correct = T)
  val=as.data.frame(p.value$p.value)
  rownames(val)=colnames(caldata.pmu)[k]
  dfm.pmu=rbind(dfm.pmu,val)
}
head(dfm.pmu)
write.csv(dfm.pmu,"4.model/kaps/mutation/path.mu/chisq.results.for.10.pathways.kaps.four.csv")
########
table(caldata.pmu$kaps.group, caldata.pmu$Cell.Cycle)
pidx=colnames(caldata.pmu)
for (i in 2:ncol(caldata.pmu)) {
  x=as.data.frame.matrix(table(caldata.pmu[,1], caldata.pmu[,i]))
  x$group=rownames(x)
  head(x)
  write.csv(x, paste0("4.model/kaps/mutation/path.mu/","10.pathways_",pidx[i],".csv"))
  library(reshape2)
  library(plyr)
  pval <- chisq.test(caldata.pmu[,1], caldata.pmu[,i],correct = T)$p.value
  #chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
  # long format with column of proportions within each id
  xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
  p1=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
    geom_bar(stat = 'identity')+
    scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
    ggtitle(paste("pathway level_",pidx[i]))
  
  generate.PDF <- function(fig) {
    pdf(paste0("4.model/kaps/mutation/path.mu/","10.pathways_",pidx[i],".pdf"),width = 4,height = 6)
    print(p1)
    dev.off()
  }
  generate.PDF(fig)
  
}
###################
###########at pathway gene level
################
pmu=read.table("4.model/kaps/mutation/path.mu/ten.paths.gene.level.alterd.txt",header = T,sep = "\t",row.names = 1)
head(pmu)
dim(pmu)
#
sh.id.pmu=intersect(rownames(kaps.td),rownames(pmu))
caldata.pmu=cbind(kaps.td[sh.id.pmu,]$kaps.group,pmu[sh.id.pmu,])
colnames(caldata.pmu)[1]="kaps.group"
head(caldata.pmu)
###select out genes at least one mutation
###
gc=as.data.frame(colSums(caldata.pmu[,-1]))
colnames(gc)="count"
gcsub=subset(gc, count>0)
head(gcsub)
###
caldata.pmu.sel=caldata.pmu[,colnames(caldata.pmu) %in% rownames(gcsub)]
caldata.pmu.sel=cbind(caldata.pmu$kaps.group, caldata.pmu.sel)
colnames(caldata.pmu.sel)[1]="kaps.group"
head(caldata.pmu.sel)
#
dfm.pmu=NULL
for (k in 2:ncol(caldata.pmu.sel)){
  p.value=chisq.test(caldata.pmu.sel[,k],caldata.pmu.sel$kaps.group,correct = T)
  val=as.data.frame(p.value$p.value)
  rownames(val)=colnames(caldata.pmu.sel)[k]
  dfm.pmu=rbind(dfm.pmu,val)
}
head(dfm.pmu)
write.csv(dfm.pmu,"4.model/kaps/mutation/path.mu/chisq.results.for.10.pathways.genelevel.kaps.four.csv")
length(dfm.pmu[which(dfm.pmu$`p.value$p.value`<0.05),])
sig.dfm=subset(dfm.pmu,dfm.pmu$`p.value$p.value`<0.05)
head(sig.dfm)
dim(sig.dfm)
######
##########################################################final hetamp
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(circlize)
######################
head(caldata.pmu.sel)
#caldata.draw=caldata.mu[order(caldata.mu$kaps.group),]
caldata.draw=caldata.pmu.sel
######pvalue anova
pvalueanova = dfm.pmu$`p.value$p.value`
is_siganova = pvalueanova < 0.05
pchanova = rep("*", length(pvalueanova))
pchanova[!is_siganova] = NA
###
col = list(kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027"),
           MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e")
)
top.anno <- columnAnnotation(kaps.group=caldata.draw$kaps.group,
                             MSI.status.bin=kaps.td[rownames(caldata.draw),]$MSI.status.bin,
                             col=col)
####
mycol=colorRampPalette(c("white", "white", "black"))(50)
mu.matrix.total=as.data.frame(t(caldata.draw[,-c(1)]))
mu.matrix.total[1:4,1:4]
dim(mu.matrix.total)
#
####choose significant mutations to show
dim(sig.dfm)
showinggs=rownames(sig.dfm)
idx=match(showinggs,rownames(mu.matrix.total))
#
# color mapping
pvalue_col_funanova = colorRamp2(c(10e-5, 0.05, 0.1), c("red", "white","green")) 
row_ha_right = rowAnnotation(
  pvalueanova = anno_simple(pvalueanova, col = pvalue_col_funanova, pch = pchanova),
  sig.genes = anno_mark(at = idx, labels = showinggs,labels_gp = gpar(fontsize = 6),link_width = unit(5, "mm")),
  #width = unit(2, "mm"),
  annotation_name_side = "bottom")
#####
genetopath=read.table("4.model/kaps/mutation/path.mu/genetopathway.txt",header = T,sep = "\t")
genetopath=genetopath[!(duplicated(genetopath$symbol)),]
rownames(genetopath)=genetopath$symbol
dim(genetopath)
head(genetopath)
#####
length(intersect(rownames(genetopath), rownames(sig.dfm)))
genetopath.sel=genetopath[rownames(sig.dfm),]
#row_ha_right = rowAnnotation(sig.genes = anno_mark(at = idx, labels = showinggs))
####

#
ht=Heatmap(as.matrix(mu.matrix.total), name = "CRC.10.pathways.gene.mutations",
           col = mycol,border = F,
           cluster_columns = F,
           show_column_names = F,show_row_names = F,
           row_names_gp = gpar(fontsize = 6),
           column_names_gp = gpar(fontsize = 10),
           #left_annotation = row_ha.left,
           right_annotation = row_ha_right,
           top_annotation = top.anno
           
)
# now we generate two legends, one for the p-value
#
lgd_pvalue.anova = Legend(title = "anova.p-value", col = pvalue_col_funanova, at = c(10e-5, 0.05, 0.1), 
                          labels = c("10e-5", "0.05", "0.1"))
lgd_sig.anova = Legend(pch = "*", type = "points", labels = "< 0.05")


ht1=draw(ht, annotation_legend_list = list(lgd_pvalue.anova, lgd_sig.anova))

generate.PDF <- function(fig) {
  pdf("4.model/kaps/mutation/path.mu//ht.153.genes.with.anova.p.sig.annoted.3.pdf",
      height = 5,width = 10)
  print(ht1)
  dev.off()
}
generate.PDF(fig)
####################################
table(caldata.pmu.sel$FAT1,caldata.pmu.sel$kaps.group)

pidx=colnames(caldata.pmu.sel)
for (i in 2:ncol(caldata.pmu.sel)) {
  x=as.data.frame.matrix(table(caldata.pmu.sel[,1], caldata.pmu.sel[,i]))
  x$group=rownames(x)
  head(x)
  write.csv(x, paste0("4.model/kaps/mutation/path.mu/indi.gene.chisq.table/","10.pathways_",pidx[i],".csv"))
  library(reshape2)
  library(plyr)
  pval <- chisq.test(caldata.pmu.sel[,1], caldata.pmu.sel[,i],correct = T)$p.value
  #chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
  # long format with column of proportions within each id
  xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
  p1=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
    geom_bar(stat = 'identity')+
    scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
    ggtitle(paste("pathway level_",pidx[i]))
  
  generate.PDF <- function(fig) {
    pdf(paste0("4.model/kaps/mutation/path.mu/indi.gene.chisq.table/","10.pathways_",pidx[i],".pdf"),width = 4,height = 6)
    print(p1)
    dev.off()
  }
  generate.PDF(fig)
  
}
##########################
###############
###oncogenic pathway genes
#my_comparisons <- list( c("set1", "set2"), c("set3", "cluster_3"), c("cluster_2", "cluster_3") )
gidx.TGFB=c("TGFBR2", "ACTA2", "COL4A1", "TAGLN", "SH3PXD2A",rownames(genetopath))

marid.tgfb=cbind(col_ha.top, genestats[rownames(col_ha.top), colnames(genestats) %in% gidx.TGFB])
marid.tgfb$kaps.group=factor(marid.tgfb$kaps.group,levels = c("set4","set3","set2","set1"))
head(marid.tgfb)
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))

for (i in 6:ncol(marid.tgfb)) {
  pdf(paste0("4.model/kaps/mutation/path.mu/indi.exp/",colnames(marid.tgfb)[i],".kaps.group.pdf"))
  #cateN=markerpanel[which(markerpanel$symbol==colnames(marid)[i]),]$category[1]
  print(ggviolin(marid.tgfb,x = "kaps.group", y = colnames(marid.tgfb)[i], fill = "kaps.group",alpha = 1,size = 0.3,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group.agg")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+
          theme_bw()+xlab("kaps.group")+ylab(colnames(marid.tgfb)[i])+ggtitle(paste0(colnames(marid.tgfb)[i]))
  )
  dev.off()
}
