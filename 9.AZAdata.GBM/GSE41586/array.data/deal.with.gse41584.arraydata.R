##########GSE41364
library(GEOquery)
library(limma)
GSEdata <- getGEO('GSE41364', destdir="9.AZAdata.GBM/GSE41586/array.data",getGPL = F)
exprSet=exprs(GSEdata[[1]])
GSEdata[[1]]
pdata.array=pData(GSEdata[[1]])
exprSet=exprs(GSEdata[[1]])
exprSet=as.data.frame(exprSet)
exprSet[1:4,1:4]
dim(exprSet)
#
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
dim(GEO.data_exp.agg)
#
head(pdata.array)
pdata.array.sel=data.frame(id1=pdata.array$description,
                           id2=pdata.array$`treatment:ch1`)
rownames(pdata.array.sel)=rownames(pdata.array)
head(pdata.array.sel)
length(intersect(colnames(GEO.data_exp.agg), rownames(pdata.array.sel)))

write.csv(pdata.array.sel, "9.AZAdata.GBM/GSE41586/array.data/pdata.array.sel.csv")
#######
pdata.array.sel.ok=read.csv("9.AZAdata.GBM/GSE41586/array.data/pdata.array.sel.csv",header = T,row.names = 1)
head(pdata.array.sel.ok)
pdata.array.sel.ok$dose=as.factor(pdata.array.sel.ok$dose)
table(pdata.array.sel.ok$cell.line,pdata.array.sel.ok$treatment)
length(intersect(colnames(GEO.data_exp.agg), rownames(pdata.array.sel.ok)))
#
#GEO.data_exp.agg=GEO.data_exp.agg[, colnames(GEO.data_exp.agg) %in% rownames(pdata.array.sel.ok)]
dim(GEO.data_exp.agg)
##########pca
library("FactoMineR")
library("factoextra")
#
aza.pca.GEOdata <- PCA(as.matrix(t(GEO.data_exp.agg)), graph = FALSE)
#
pdf("9.AZAdata.GBM/GSE41586/array.data/pca.before.batch.gse41586.arraydata.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.GEOdata,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE41586.CRC.HT29")+theme_bw()
fviz_pca_ind(aza.pca.GEOdata,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$dose, # color by groups
             palette = c( "#ED2891","#B2509E","#D49DC7","#C1A72F","#E8C51D","#F9ED32",
                          "#104A7F","#9EDDF9","#007EB5","#CACCDB","#6E7BA2","#DAF1FC","#00AEEF",
                          "#F6B667","#D97D25","#FBE3C7"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE41586.CRC.HT29")+theme_bw()
dev.off()
########
########
save(GEO.data_exp.agg,pdata.array.sel.ok, file="9.AZAdata.GBM/GSE41586/array.data/combat_edata.azaexp.GEOdata.GSE41586.arraydata.RData")
#######
#######ssgsea
#
load("9.AZAdata.GBM/GSE41586/array.data/genesets.hallmarker.score.aza.gse41586.arraydata.RData")
##############
stat.hallmarkerscore.aza.gse41586.arraydata=cbind(pdata.array.sel.ok, genesets.hall.marker.score.aza.gse41586.arraydata[rownames(pdata.array.sel.ok),])
stat.hallmarkerscore.aza.gse41586.arraydata$treatment=as.factor(stat.hallmarkerscore.aza.gse41586.arraydata$treatment)
head(stat.hallmarkerscore.aza.gse41586.arraydata)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza.gse41586.arraydata[,4:ncol(stat.hallmarkerscore.aza.gse41586.arraydata)])){    ######change
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza.gse41586.arraydata)
}
pvalue.hm.aza.gse41586.arraydata=do.call(rbind, datalist)
pvalue.hm.aza.gse41586.arraydata=as.data.frame(pvalue.hm.aza.gse41586.arraydata)

#pvalue.hm.sig.aza.GEOdata=subset(pvalue.hm.aza.GEOdata, p.value< 0.05)
########
########draw heatmap
stat.hallmarkerscore.aza.gse41586.arraydata[1:4,1:4]
#ht.hallM.sig.aza=as.data.frame(t(stat.hallmarkerscore.aza[,-c(1:3)]))[rownames(pvalue.hm.sig.aza),]
ht.hallM.sig.aza.gse41586.arraydata=as.data.frame(t(stat.hallmarkerscore.aza.gse41586.arraydata[,-c(1:3)]))    ####change
#zscore
ht.hallM.sig.aza.z.gse41586.arraydata=as.data.frame(t(scale(t(ht.hallM.sig.aza.gse41586.arraydata))))
ht.hallM.sig.aza.z.gse41586.arraydata[ht.hallM.sig.aza.z.gse41586.arraydata< -2] <- -2
ht.hallM.sig.aza.z.gse41586.arraydata[ht.hallM.sig.aza.z.gse41586.arraydata>2] <- 2
ht.hallM.sig.aza.z.gse41586.arraydata[1:4,1:4]

min_cor = min(as.vector(ht.hallM.sig.aza.z.gse41586.arraydata))
max_cor = max(as.vector(ht.hallM.sig.aza.z.gse41586.arraydata))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
###
row_ha.left.hm.aza.gse41586.arraydata = rowAnnotation(wilcox.pvalue=as.numeric(pvalue.hm.aza.gse41586.arraydata$p.value),  #change
                                            col=list(
                                              #wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001,0.0001), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354")),
                                              wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001), c("#f1b6da","#de77ae","#c51b7d","#8e0152"))
                                            ),show_annotation_name = FALSE)

col_ha_top.hm.aza.gse41586.arraydata = columnAnnotation(
  treatment=pdata.array.sel.ok$treatment,
  cell.line=pdata.array.sel.ok$cell.line,
  col=list(treatment=c("High"="#cc4c02","Low"="#E7B800", "Control"="#FC4E07"),
           cell.line=c("HT29"="#8dd3c7")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))


####
hmht.aza.gse41586.arraydata=Heatmap(ht.hallM.sig.aza.z.gse41586.arraydata, name = "hall.marker", 
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
                          top_annotation = col_ha_top.hm.aza.gse41586.arraydata,
                          #right_annotation = row_ha.right,
                          show_row_dend = F,show_column_dend = F,
                          #row_names_side = "left",
                          left_annotation = row_ha.left.hm.aza.gse41586.arraydata,
                          column_title=paste0("hall.marker.genesets"," GSE41586 CRC HT29"),
                          column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/GSE41586/array.data/ht.hallM.ssgsea.aza.GSE41586.CRC.HT29.arraydata.pdf",width = 8,height = 7)
draw(hmht.aza.gse41586.arraydata, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
############
############
save(pvalue.hm.aza.gse41586.arraydata,stat.hallmarkerscore.aza.gse41586.arraydata,ht.hallM.sig.aza.gse41586.arraydata,
     file="9.AZAdata.GBM/GSE41586/array.data/ht.hallM.ssgsea.aza.GSE41586.CRC.ht29.arraydata.RData")
##############
##############
##############degs
##############
#####simply compare
library(limma)
pdata.array.sel.ok.6samples=subset(pdata.array.sel.ok, treatment=="Control"|treatment=="Low")
degdata.gse41586.array=GEO.data_exp.agg[, colnames(GEO.data_exp.agg) %in% rownames(pdata.array.sel.ok.6samples)]
#
design2=model.matrix(~0+ factor(pdata.array.sel.ok.6samples$treatment))
colnames(design2)=c('Control','Low')
head(design2)
fit=lmFit(degdata.gse41586.array,design2)
cont.matrix=makeContrasts('Low-Control',levels = design2)
fit=contrasts.fit(fit,cont.matrix)
fit=eBayes(fit)
#
res.deg.aza.gse41586.array<- topTable(fit, adjust.method="BH", coef=1,number = Inf,resort.by="p")
res.deg.aza.gse41586.array$FC=2^(res.deg.aza.gse41586.array$logFC)
head(res.deg.aza.gse41586.array)
#############
############
write.csv(res.deg.aza.gse41586.array, "9.AZAdata.GBM/GSE41586/array.data/degs.genes.AZA.gse41586.arraydata.lowVScontrol.csv")
