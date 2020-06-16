########## 5-AZA teatment data
##########
azaexpm=read.table("9.AZAdata.GBM/GSE80137_normalized_counts-source.txt",header = T,sep = "\t")
head(azaexpm)
rownames(azaexpm)=azaexpm$Gene.ID
dim(azaexpm)
####
azaexpm.sel=azaexpm
rownames(azaexpm.sel)=azaexpm$Gene.ID
azaexpm.sel=azaexpm.sel[,-c(1:3)]
azaexpm.sel=log(azaexpm.sel+1)
head(azaexpm.sel)
####filter out
azanacount=as.data.frame(rowSums(azaexpm.sel == 0))
colnames(azanacount)="count"
head(azanacount)
azanacount.sel=subset(azanacount, count<=8)
dim(azanacount.sel)
#########
azaexpm.sel.filter=azaexpm.sel[rownames(azanacount.sel),]
dim(azaexpm.sel.filter)
#########
azapdata=as.data.frame(colnames(azaexpm.sel.filter))
colnames(azapdata)="id"
azapdata$treament=rep(c("untreatment","untreatment","aza","aza"),times=3)
#write.csv(azapdata,"9.AZAdata.GBM/aza.GBM.pdata.csv")
head(azapdata)
##################
############
azapdata.modify=read.csv("9.AZAdata.GBM/aza.GBM.pdata.csv",header = T,row.names = 1)
head(azapdata.modify)
#pca
library("FactoMineR")
library("factoextra")
#
aza.pca <- PCA(as.matrix(t(azaexpm.sel.filter)), graph = FALSE)
#
pdf("9.AZAdata.GBM/pca.before.batch.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca,pointshape = 19,geom = "point",invisible="quali",
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = azapdata.modify$treatment, # color by groups
               palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups",title="PCA.GBM.cell.line")+theme_bw()
fviz_pca_ind(aza.pca,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = azapdata.modify$cell.line, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GBM.cell.line")+theme_bw()
dev.off()


##########remove batch effect
library(sva)
library(pamr)
library(limma)
aza.batch = azapdata.modify$cell.line
aza.modcombat = model.matrix(~1, data=azapdata.modify)
combat_edata.azaexp = ComBat(dat=as.matrix(azaexpm.sel.filter), batch=aza.batch, mod=aza.modcombat, par.prior=TRUE, prior.plots=FALSE)
head(combat_edata.azaexp)
##########
aza.pca.after <- PCA(as.matrix(t(combat_edata.azaexp)), graph = FALSE)

pdf("9.AZAdata.GBM/pca.after.batch.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = azapdata.modify$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GBM.cell.line")+theme_bw()
fviz_pca_ind(aza.pca.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = azapdata.modify$cell.line, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GBM.cell.line")+theme_bw()
dev.off()
##################degs
SibShip <- factor(azapdata.modify$sample)
Treat <- factor(azapdata.modify$treatment, levels=c("untreatment","aza"))
design.aza <- model.matrix(~SibShip+Treat)
fit.aza <- lmFit(combat_edata.azaexp, design.aza)
fit.aza <- eBayes(fit.aza)
##################
res.deg.aza<- topTable(fit.aza, adjust.method="BH", coef="Treataza",number = Inf,resort.by="p")
res.deg.aza$FC=2^(res.deg.aza$logFC)
all(rownames(res.deg.aza) %in% rownames(azaexpm))
res.deg.aza=cbind(azaexpm[rownames(res.deg.aza),1:2], res.deg.aza, combat_edata.azaexp[rownames(res.deg.aza),] )
head(res.deg.aza)
write.csv(res.deg.aza, "9.AZAdata.GBM/degs.genes.AZA.GBM.csv")
dim(res.deg.aza)
##################
save(res.deg.aza, azapdata.modify, azaexpm,azaexpm.sel.filter,file = "9.AZAdata.GBM/degs.genes.AZA.GBM.RData")
##################
##################map the te inhibitor 
#
head(inhibition.markers)
aa=inhibition.markers
aa$Gene.name=aa$symbol
res.deg.aza.inhibitors=merge(aa, res.deg.aza, by="Gene.name", all.x=TRUE)
dim(res.deg.aza.inhibitors)
write.csv(res.deg.aza.inhibitors,"9.AZAdata.GBM/degs.genes.AZA.GBM.res.deg.aza.inhibitors.csv")
#################
#################heatmap of significant te inhibitors
#################
res.deg.aza.inhibitors.draw=subset(res.deg.aza.inhibitors, adj.P.Val < 0.05)
head(res.deg.aza.inhibitors.draw)
colnames(res.deg.aza.inhibitors.draw)
#####
#########draw heatmap
library(circlize)
head(azapdata.modify)
col_ha.aza=columnAnnotation(treatment=azapdata.modify$treatment,
                            cell.line=azapdata.modify$cell.line,
                            col=list(treatment=c("aza"="#E7B800", "untreatment"="#FC4E07"),
                                     cell.line=c("LNT.229"="#00AFBB","T98G"="#e41a1c","U.87"="#377eb8")
                            )
                            ,show_annotation_name = TRUE)
library(ComplexHeatmap)
#
aa=res.deg.aza.inhibitors.draw[,c(12:23)]
rownames(aa)=res.deg.aza.inhibitors.draw$Gene.name
aaz=as.data.frame(t(scale(t(aa))))
#
htaza=Heatmap(aaz, column_title = "TE.inhibitors.incrased.after.aza.treatment", name = "exp",
               #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               col = viridis(10),border = F,
               #col=col.pal_cor,
               show_row_dend = F,
               show_column_names = F,show_row_names = T,
               cluster_columns = T,
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               top_annotation  = col_ha.aza#, 
               #right_annotation = row_ha.right,
               #show_row_dend = T,show_column_dend = T,
               #row_names_side = "left",
               #left_annotation = row_ha.module.gene
)
pdf("9.AZAdata.GBM/heatmap.sig.increased.TE.inhibitors.after.AZA.treatment.pdf",width = 5,height = 6)
draw(htaza, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     ,heatmap_legend_side = "bottom"
)
dev.off()
###############
###############
###############enrichment analysis
###############
## 转换
library(clusterProfiler)
gene.aza = bitr(res.deg.aza$Gene.name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene.aza <- dplyr::distinct(gene.aza,SYMBOL,.keep_all=TRUE)
gene_df.aza <- data.frame(logFC=res.deg.aza$logFC,
                          SYMBOL = res.deg.aza$Gene.name)
gene_df.aza=gene_df.aza[!(duplicated(gene_df.aza$SYMBOL)),]
gene_df.aza <- merge(gene_df.aza,gene.aza,by="SYMBOL")
## geneList 三部曲
## 1.获取基因logFC
geneList.aza <- gene_df.aza$logFC
## 2.命名
names(geneList.aza) = gene_df.aza$ENTREZID
## 3.排序很重要
geneList.aza = sort(geneList.aza, decreasing = TRUE)
#
hallmarks <- read.gmt("~/nas/Xiaoqiang/opti.data/TCGA.data/liver/stat.20200418/c2.all.v7.1.entrez.gmt")
hallmarks <- read.gmt("~/nas/Xiaoqiang/opti.data/TCGA.data/liver/stat.20200418/output.correlation.GSEA/h.all.v7.1.entrez.gmt")
y.aza <- GSEA(geneList.aza,TERM2GENE =hallmarks)
yd.aza <- data.frame(y.aza)
write.csv(yd.aza,"9.AZAdata.GBM/gsea.based.on.ranklogFC.hallmarker.csv")
#############
#############
#############deal with expMatrix for ssgsea of hallmarker
head(res.deg.aza)
aza.expmatrix.normd=res.deg.aza[,-c(1,3:9)]
head(aza.expmatrix.normd)
library(dplyr)
aza.expmatrix.normd.symbol= aza.expmatrix.normd %>% group_by(Gene.name) %>% summarise_all(mean)
aza.expmatrix.normd.symbol=as.data.frame(aza.expmatrix.normd.symbol)
rownames(aza.expmatrix.normd.symbol)=aza.expmatrix.normd.symbol$Gene.name
aza.expmatrix.normd.symbol=aza.expmatrix.normd.symbol[,-1]
#
tail(aza.expmatrix.normd.symbol)
dim(aza.expmatrix.normd.symbol)
save(aza.expmatrix.normd.symbol,file="9.AZAdata.GBM/aza.expmatrix.normd.symbol.RData")
############compre hallmarker geneset score
load("9.AZAdata.GBM/aza.hallmarker.ssgsea/genesets.hallmarker.score.aza.RData")
dim(genesets.hall.marker.score.aza)
#####
stat.hallmarkerscore.aza=cbind(azapdata.modify, genesets.hall.marker.score.aza[rownames(azapdata.modify),])
head(stat.hallmarkerscore.aza)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza[,4:ncol(stat.hallmarkerscore.aza)])){  
  datalist[[i]] <- wilcox.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza)
}
pvalue.hm.aza=do.call(rbind, datalist)
pvalue.hm.aza=as.data.frame(pvalue.hm.aza)

pvalue.hm.sig.aza=subset(pvalue.hm.aza, p.value< 0.05)
########
########draw heatmap
stat.hallmarkerscore.aza[1:4,1:4]
#ht.hallM.sig.aza=as.data.frame(t(stat.hallmarkerscore.aza[,-c(1:3)]))[rownames(pvalue.hm.sig.aza),]
ht.hallM.sig.aza=as.data.frame(t(stat.hallmarkerscore.aza[,-c(1:3)]))
#zscore
ht.hallM.sig.aza.z=as.data.frame(t(scale(t(ht.hallM.sig.aza))))
ht.hallM.sig.aza.z[ht.hallM.sig.aza.z< -2] <- -2
ht.hallM.sig.aza.z[ht.hallM.sig.aza.z>2] <- 2
ht.hallM.sig.aza.z[1:4,1:4]

min_cor = min(as.vector(ht.hallM.sig.aza))
max_cor = max(as.vector(ht.hallM.sig.aza))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
###
row_ha.left.hm.aza = rowAnnotation(wilcox.pvalue=as.numeric(pvalue.hm.aza$p.value),  #change
                               col=list(
                                 #wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001,0.0001), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354")),
                                 wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001), c("#f1b6da","#de77ae","#c51b7d","#8e0152"))
                               ),show_annotation_name = FALSE)

col_ha_top.hm.aza = columnAnnotation(
  treatment=stat.hallmarkerscore.aza$treatment,
  cell.line=stat.hallmarkerscore.aza$cell.line,
  col=list(treatment=c("aza"="#E7B800", "untreatment"="#FC4E07"),
           cell.line=c("LNT.229"="#00AFBB","T98G"="#e41a1c","U.87"="#377eb8")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
hmht.aza=Heatmap(ht.hallM.sig.aza.z, name = "hall.marker", 
             #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             #col = rev(viridis(10)),
             col = viridis(10),
             width = unit(2, "cm"),
             height = unit(12, "cm"),
             border = F,
             #col=col.pal_cor,
             show_column_names = T,show_row_names = T,
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             top_annotation = col_ha_top.hm.aza,
             #right_annotation = row_ha.right,
             show_row_dend = F,show_column_dend = F,
             #row_names_side = "left",
             left_annotation = row_ha.left.hm.aza,
             column_title=paste0("hall.marker.genesets"," GSE80137 GMB"),
             column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/aza.hallmarker.ssgsea/ht.hallM.ssgsea.aza.ALL.pdf",width = 8,height = 7)
draw(hmht.aza, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#
save(stat.hallmarkerscore.aza, pvalue.hm.aza,ht.hallM.sig.aza,file="9.AZAdata.GBM/aza.hallmarker.ssgsea/ht.hallM.ssgsea.aza.RData")
#############
#############