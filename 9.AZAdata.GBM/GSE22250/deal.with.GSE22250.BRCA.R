##########GSE22250
library(GEOquery)
library(limma)
GSEdata <- getGEO('GSE22250', destdir=".",getGPL = F)
exprSet=exprs(GSEdata[[1]])
GSEdata[[1]]
pdata.array=pData(GSEdata[[1]])
exprSet=exprs(GSEdata[[1]])
exprSet=as.data.frame(exprSet)
exprSet[1:4,1:4]
dim(exprSet)
###############
library(hgu133plus2.db)
gene.symbol = unlist(as.list(hgu133plus2SYMBOL))
length(gene.symbol)
GEO.data_exp<- cbind(gene.symbol, exprSet)
GEO.data_exp[1:4,1:4]
GEO.data_exp=na.omit(GEO.data_exp, cols=c("gene.symbol"))
GEO.data_exp.agg=GEO.data_exp %>% group_by(gene.symbol) %>% summarise_all(mean)
GEO.data_exp.agg=as.data.frame(GEO.data_exp.agg)
rownames(GEO.data_exp.agg)=GEO.data_exp.agg$gene.symbol
GEO.data_exp.agg=GEO.data_exp.agg[,-1]
GEO.data_exp.agg[1:4,1:4]
#
head(pdata.array)
pdata.array.sel=data.frame(treatment=pdata.array$`conditions:ch1`,
                           cell.line=pdata.array$`cell line:ch1`)
rownames(pdata.array.sel)=rownames(pdata.array)
pdata.array.sel.ok=pdata.array.sel[c(2:15),]
pdata.array.sel.ok$cell.line=gsub("-",".",pdata.array.sel.ok$cell.line)
head(pdata.array.sel)
pdata.array.sel.ok
#
GEO.data_exp.agg=GEO.data_exp.agg[, colnames(GEO.data_exp.agg) %in% rownames(pdata.array.sel.ok)]
dim(GEO.data_exp.agg)
##########pca
library("FactoMineR")
library("factoextra")
#
aza.pca.GEOdata <- PCA(as.matrix(t(GEO.data_exp.agg)), graph = FALSE)
#
pdf("9.AZAdata.GBM/GSE22250/pca.before.batch.gse22250.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.GEOdata,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE22250.BRCA.cell.line")+theme_bw()
fviz_pca_ind(aza.pca.GEOdata,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$cell.line, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE22250.BRCA.cell.line")+theme_bw()
dev.off()
#############
#############
##########remove batch effect
library(sva)
library(pamr)
library(limma)
aza.batch.GEOdata = pdata.array.sel.ok$cell.line
aza.modcombat.GEOdata = model.matrix(~1, data=pdata.array.sel.ok)
combat_edata.azaexp.GEOdata = ComBat(dat=as.matrix(GEO.data_exp.agg), batch=aza.batch.GEOdata, mod=aza.modcombat.GEOdata, par.prior=TRUE, prior.plots=FALSE)
head(combat_edata.azaexp.GEOdata)
#########PCA again
#
aza.pca.GEOdata.after <- PCA(as.matrix(t(combat_edata.azaexp.GEOdata)), graph = FALSE)
pdf("9.AZAdata.GBM/GSE22250/pca.after.batch.gse22250.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.GEOdata.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE22250.BRCA.cell.line")+theme_bw()
fviz_pca_ind(aza.pca.GEOdata.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$cell.line, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE22250.BRCA.cell.line")+theme_bw()
dev.off()
##############
##############
save(combat_edata.azaexp.GEOdata, file="9.AZAdata.GBM/GSE22250/combat_edata.azaexp.GEOdata.GSE22250.RData")
#############
#ssgsea
#
load("9.AZAdata.GBM/GSE22250/genesets.hallmarker.score.aza.gse22250.RData")
#
dim(pdata.array.sel.ok)
#####
stat.hallmarkerscore.aza.GEOdata=cbind(pdata.array.sel.ok, genesets.hall.marker.score.aza.gse22250[rownames(pdata.array.sel.ok),])
head(stat.hallmarkerscore.aza.GEOdata)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza.GEOdata[,3:ncol(stat.hallmarkerscore.aza.GEOdata)])){  
  datalist[[i]] <- wilcox.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza.GEOdata)
}
pvalue.hm.aza.GEOdata=do.call(rbind, datalist)
pvalue.hm.aza.GEOdata=as.data.frame(pvalue.hm.aza.GEOdata)

#pvalue.hm.sig.aza.GEOdata=subset(pvalue.hm.aza.GEOdata, p.value< 0.05)
########
########draw heatmap
stat.hallmarkerscore.aza.GEOdata[1:4,1:4]
stat.hallmarkerscore.aza.GEOdata=stat.hallmarkerscore.aza.GEOdata[order(stat.hallmarkerscore.aza.GEOdata$treatment,decreasing = T),]
#ht.hallM.sig.aza=as.data.frame(t(stat.hallmarkerscore.aza[,-c(1:3)]))[rownames(pvalue.hm.sig.aza),]
ht.hallM.sig.aza.GEOdata=as.data.frame(t(stat.hallmarkerscore.aza.GEOdata[,-c(1:2)]))
#zscore
ht.hallM.sig.aza.z.GEOdata=as.data.frame(t(scale(t(ht.hallM.sig.aza.GEOdata))))
ht.hallM.sig.aza.z.GEOdata[ht.hallM.sig.aza.z.GEOdata< -2] <- -2
ht.hallM.sig.aza.z.GEOdata[ht.hallM.sig.aza.z.GEOdata>2] <- 2
ht.hallM.sig.aza.z.GEOdata[1:4,1:4]

min_cor = min(as.vector(ht.hallM.sig.aza.z.GEOdata))
max_cor = max(as.vector(ht.hallM.sig.aza.z.GEOdata))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
###
row_ha.left.hm.aza.GEOdata = rowAnnotation(wilcox.pvalue=as.numeric(pvalue.hm.aza.GEOdata$p.value),  #change
                                   col=list(
                                     #wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001,0.0001), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354")),
                                     wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001), c("#f1b6da","#de77ae","#c51b7d","#8e0152"))
                                   ),show_annotation_name = FALSE)

col_ha_top.hm.aza.GEOdata = columnAnnotation(
  treatment=stat.hallmarkerscore.aza.GEOdata$treatment,
  cell.line=stat.hallmarkerscore.aza.GEOdata$cell.line,
  col=list(treatment=c("5-AZA"="#E7B800", "WT"="#FC4E07"),
           cell.line=c("MCF.7"="#1b9e77","T47D"="#d95f02","SKBR3"="#7570b3","BT20"="#e7298a","MDA.MB.231"="#66a61e","MDA.MB.361"="#e6ab02","ZR.75.1"="#a6761d")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
hmht.aza.GEOdata=Heatmap(ht.hallM.sig.aza.z.GEOdata, name = "hall.marker", 
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
                 top_annotation = col_ha_top.hm.aza.GEOdata,
                 #right_annotation = row_ha.right,
                 show_row_dend = F,show_column_dend = F,
                 #row_names_side = "left",
                 left_annotation = row_ha.left.hm.aza.GEOdata,
                 column_title=paste0("hall.marker.genesets"," GSE22250 BRCA"),
                 column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/GSE22250/ht.hallM.ssgsea.aza.GSE22250.BRCAs.pdf",width = 8,height = 7)
draw(hmht.aza.GEOdata, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#####################
##########
#degs
##################degs
head(combat_edata.azaexp.GEOdata)
head(pdata.array.sel.ok)
pdata.array.sel.ok$treatment=gsub("-","",pdata.array.sel.ok$treatment)

SibShip <- factor(pdata.array.sel.ok$cell.line)
Treat <- factor(pdata.array.sel.ok$treatment, levels=c("WT","5AZA"))
design.aza <- model.matrix(~SibShip+Treat)
fit.aza.GEOdata <- lmFit(combat_edata.azaexp.GEOdata, design.aza)
fit.aza.GEOdata <- eBayes(fit.aza.GEOdata)
##################
res.deg.aza.GEOdata<- topTable(fit.aza.GEOdata, adjust.method="BH", coef="Treat5AZA",number = Inf,resort.by="p")
res.deg.aza.GEOdata$FC=2^(res.deg.aza.GEOdata$logFC)
head(res.deg.aza.GEOdata)
write.csv(res.deg.aza.GEOdata, "9.AZAdata.GBM/GSE22250/degs.genes.AZA.gse22250.BRCA.csv")
dim(res.deg.aza.GEOdata)
##################
##################
