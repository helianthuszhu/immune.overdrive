########validation on GSE107422
###
rdsdata=readRDS("vali.GSE107422/RE.output/RE_intergenic_2_counts_normalized.RDS")
rdsdata
rdsexp=as.data.frame(exprs(rdsdata))
dim(rdsexp)
rdsexp.sel=as.data.frame(t(rdsexp[rownames(cndi.rep),]))
rdsexp.sel$mean.exp=rowMeans(rdsexp.sel)
rdsexp.sel$z.of.mean.exp=(rdsexp.sel$mean.exp - mean(rdsexp.sel$mean.exp))/sd(rdsexp.sel$mean.exp)
#rdsexp.sel$mean.exp=rowSums(rdsexp.sel)/9
#rdsexp.sel=as.data.frame(scale(rdsexp.sel))
head(rdsexp.sel)
summary(rdsexp.sel$z.of.mean.exp)
###
rdsclin=read.table("vali.GSE107422/clin.GSE107422.txt",header = T,row.names = 1)
rdsclin=rdsclin[rownames(rdsexp.sel),]
length(intersect(rownames(rdsclin), rownames(rdsexp.sel)))
rdsclin$sampleID=gsub("-",".",rdsclin$sampleID)
head(rdsclin)
#######
datadraw.gse170422=cbind(rdsexp.sel,rdsclin )
datadraw.gse170422=datadraw.gse170422[order(datadraw.gse170422$z.of.mean.exp,decreasing = T),]
head(datadraw.gse170422)
#
datadraw.gse170422$group=cut(datadraw.gse170422$z.of.mean.exp, 
                             breaks = c(min(datadraw.gse170422$z.of.mean.exp),0.33,0.957,1.347,max(datadraw.gse170422$z.of.mean.exp)),
                             labels = c("set1","set2","set3","set4"), include.lowest = TRUE)


table(datadraw.gse170422$DFS.status, datadraw.gse170422$group)
chisq.test(datadraw.gse170422$DFS.status, datadraw.gse170422$group,correct = T)
####
######TE annotation
TE.ann.data=cndi.rep[rownames(cndi.rep),c(2,3)]

TE.an.ha=rowAnnotation(repClass=TE.ann.data$repClass,
                       repFamily=TE.ann.data$repFamily,
                       col=list(repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")),show_annotation_name = FALSE
)

min_cor = min(as.vector(datadraw.gse170422[,c(1:9)]))
max_cor = max(as.vector(datadraw.gse170422[,c(1:9)]))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#00F5FF", "white","#FF3E96"))(50))
#

anno.gse107422=columnAnnotation(z.of.mean.exp.score=anno_lines(datadraw.gse170422$z.of.mean.exp),
                                z.of.mean.exp=datadraw.gse170422$z.of.mean.exp,
                                DFS.status=datadraw.gse170422$DFS.status,
                                group=datadraw.gse170422$group,
                                col=list(z.of.mean.exp=colorRamp2(c(-3,-2,0,2,3), 
                                                                  c("#ffffcc","#d9f0a3","#addd8e","#78c679","#006837")),
                                         DFS.status=c("Yes"="#d01c8b","No"="#4dac26"),
                                         group= c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                                ),show_annotation_name = TRUE)



ht.gse17=Heatmap(t(scale(datadraw.gse170422[,c(1:9)])),#col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256)
                       col=col.pal_cor,
                       top_annotation = anno.gse107422,
                       show_column_names = F,left_annotation = TE.an.ha,
                       cluster_columns = F,cluster_rows = T, column_title = "clinical.comparison among TE clusters.GSE107422")
pdf("vali.GSE107422/ht.TEs.gse107422.pdf",width = 10,height = 7)
draw(ht.gse17, padding = unit(c(40, 20, 40,20), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     #,heatmap_legend_side = "right"
)

dev.off()
#####

#########
rds.stats.gse170422=datadraw.gse170422
rds.stats.gse170422$group=factor(rds.stats.gse170422$group,levels = c("set4", "set3","set2","set1"))


pdf("vali.GSE107422/test.boxplot.pdf",width = 4,height = 5)
ggboxplot(subset(rds.stats.gse170422, group=="set4"), x = "group", y = "z.of.mean.exp",
               color = "DFS.status", palette = c("Yes"="#d01c8b","No"="#4dac26"),
               add = "jitter")+stat_compare_means(aes(group = DFS.status))
ggboxplot(subset(rds.stats.gse170422), x = "DFS.status", y = "z.of.mean.exp",
          color = "DFS.status", palette = c("Yes"="#d01c8b","No"="#4dac26"),
          add = "jitter")+stat_compare_means(aes(group = DFS.status))
dev.off()





ggviolin(subset(rds.stats.gse170422),x = "group", y ="z.of.mean.exp" , fill = "DFS.status",alpha = 1,size = 0.01,width = 1,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )#,
         #add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2)
         )+
  labs(fill = "DFS.status")+
  stat_compare_means(label = "p.signif")+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+
  xlab("kaps.group")+ylab("z.of.mean.exp")+ggtitle(paste0("z.of.mean.exp.comparison","\n","in.DFS.Yes.samples (n=38)"))




###########
save(datadraw.gse170422, rds.stats.gse170422,file="vali.GSE107422/ht.TEs.box.gse107422.RData")
################
#########
##########
##########deal with TPM data
####
fi<-list.files(path = "~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/vali.GSE107422/TPM/file",full.names=T)
fi
#datalist = list()
df=read.table("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/vali.GSE107422/TPM/file/GSM2866686_AMC.R1_log2_TPM.txt.gz",header = T)
head(df)
#
for (k in 2:length(fi)){
  a <- read.table(fi[k],header = T)
  #head(a)
  #
  #datalist[[k]] <- a
  df=merge(a,df,by.x = "unique_id",all = T)
}
dim(df)
rownames(df)=df.anno$refseq
head(df)
df.anno=as.data.frame(df$unique_id)
colnames(df.anno)="id"
df.anno$id2=df.anno$id
df.anno=separate(df.anno, "id2",into = c("giN","refs","ref","NMid"),sep = "[|]")
df.anno$NMid2=df.anno$NMid
df.anno=separate(df.anno, "NMid2",into = c("refseq","point"),sep = "[.]")
head(df.anno)
dim(df.anno)
length(unique(df.anno$refseq))
#
eg=bitr(df.anno$refseq, fromType="REFSEQ", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
colnames(eg)[1]="refseq"
head(eg)
##########
df.anno.symbol=merge(eg,df.anno, by="refseq",all.x=TRUE)
rownames(df.anno.symbol)=df.anno.symbol$refseq
head(df.anno.symbol)
dim(df.anno.symbol)
#
df.sel=df[rownames(df.anno.symbol),]
df.sel=cbind(df.anno.symbol$SYMBOL, df.sel)
colnames(df.sel)[1]="symbol"
df.sel=df.sel[,-2]
dim(df.sel)
head(df.sel)
########agg
df.sel.agg= df.sel %>% group_by(symbol) %>% summarise_all(mean)
df.sel.agg=as.data.frame(df.sel.agg)
rownames(df.sel.agg)=df.sel.agg$symbol
df.sel.agg=df.sel.agg[,-1]
head(df.sel.agg)
dim(df.sel.agg)
########
save(df.sel.agg,file = "vali.GSE107422/TPM/TPM.expMatrix.GSE107422.RData")
############
length(intersect(colnames(df.sel.agg), datadraw.gse170422$sampleID))
###################
###################
