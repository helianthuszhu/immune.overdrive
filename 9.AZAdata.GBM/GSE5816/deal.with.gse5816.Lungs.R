##########GSE5816
library(GEOquery)
library(limma)
GSEdata <- getGEO('GSE5816', destdir="9.AZAdata.GBM/GSE5816/",getGPL = F)
exprSet=exprs(GSEdata[[1]])
GSEdata[[1]]
pdata.array=pData(GSEdata[[1]])
exprSet=exprs(GSEdata[[1]])
exprSet=as.data.frame(exprSet)
exprSet[1:4,1:4]
dim(exprSet)
#
library(affy)
GSEdata2 <- ReadAffy(celfile.path = "~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/9.AZAdata.GBM/GSE5816/raw.file/")
data_rma.GSEdata <- rma(GSEdata2)
save(data_rma.GSEdata, file="9.AZAdata.GBM/GSE5816/GSE5816_eset.Rdata")
#
exprSet.GSEdata <- exprs(data_rma.GSEdata)
exprSet.GSEdata <- as.data.frame(exprSet.GSEdata)
colnames(exprSet.GSEdata)=substr(colnames(exprSet.GSEdata),1,nchar(as.character(colnames(exprSet.GSEdata)))-7)
exprSet.GSEdata[1:4,1:4]
dim(exprSet.GSEdata)
###############
library(hgu133plus2.db)
gene.symbol = unlist(as.list(hgu133plus2SYMBOL))
length(gene.symbol)
GEO.data_exp<- cbind(gene.symbol, exprSet.GSEdata)
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
                           id2=pdata.array$treatment_protocol_ch1)
rownames(pdata.array.sel)=rownames(pdata.array)
head(pdata.array.sel)
length(intersect(colnames(GEO.data_exp.agg), rownames(pdata.array.sel)))

write.csv(pdata.array.sel, "9.AZAdata.GBM/GSE5816/pdata.array.sel.csv")
#######
pdata.array.sel.ok=read.csv("9.AZAdata.GBM/GSE5816/pdata.array.sel.csv",header = T,row.names = 1)
head(pdata.array.sel.ok)
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
pdf("9.AZAdata.GBM/GSE5816/pca.before.batch.gse5816.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.GEOdata,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE5816.Lung.cell.line")+theme_bw()
fviz_pca_ind(aza.pca.GEOdata,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$cell.line, # color by groups
             palette = c( "#ED2891","#B2509E","#D49DC7","#C1A72F","#E8C51D","#F9ED32",
                          "#104A7F","#9EDDF9","#007EB5","#CACCDB","#6E7BA2","#DAF1FC","#00AEEF",
                          "#F6B667","#D97D25","#FBE3C7"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE5816.Lung.cell.line")+theme_bw()

fviz_pca_ind(aza.pca.GEOdata,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$cell.type, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE5816.Lung.cell.line")+theme_bw()
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

pdf("9.AZAdata.GBM/GSE5816/pca.after.batch.gse5816.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.GEOdata.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE5816.Lung.cell.line")+theme_bw()

fviz_pca_ind(aza.pca.GEOdata.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$cell.line, # color by groups
             palette = c( "#ED2891","#B2509E","#D49DC7","#C1A72F","#E8C51D","#F9ED32",
                          "#104A7F","#9EDDF9","#007EB5","#CACCDB","#6E7BA2","#DAF1FC","#00AEEF",
                          "#F6B667","#D97D25","#FBE3C7"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE5816.Lung.cell.line")+theme_bw()

fviz_pca_ind(aza.pca.GEOdata.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata.array.sel.ok$cell.type, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE5816.Lung.cell.line")+theme_bw()

dev.off()
##############
##############
save(combat_edata.azaexp.GEOdata,pdata.array.sel.ok, file="9.AZAdata.GBM/GSE5816/combat_edata.azaexp.GEOdata.GSE5816.RData")
#############
#ssgsea
#
load("9.AZAdata.GBM/GSE5816/ssgsea/genesets.hallmarker.score.aza.gse5816.RData")
#
#####
stat.hallmarkerscore.aza.GEOdata=cbind(pdata.array.sel.ok, genesets.hall.marker.score.aza.gse5816[rownames(pdata.array.sel.ok),])
stat.hallmarkerscore.aza.GEOdata$treatment=as.factor(stat.hallmarkerscore.aza.GEOdata$treatment)
stat.hallmarkerscore.aza.GEOdata$cell.type=gsub(" ","",stat.hallmarkerscore.aza.GEOdata$cell.type)
head(stat.hallmarkerscore.aza.GEOdata)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza.GEOdata[,4:ncol(stat.hallmarkerscore.aza.GEOdata)])){    ######change
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza.GEOdata)
}
pvalue.hm.aza.GEOdata=do.call(rbind, datalist)
pvalue.hm.aza.GEOdata=as.data.frame(pvalue.hm.aza.GEOdata)

#pvalue.hm.sig.aza.GEOdata=subset(pvalue.hm.aza.GEOdata, p.value< 0.05)
########
########draw heatmap
stat.hallmarkerscore.aza.GEOdata[1:4,1:4]
stat.hallmarkerscore.aza.GEOdata$treatment=gsub("High","3High",stat.hallmarkerscore.aza.GEOdata$treatment)
stat.hallmarkerscore.aza.GEOdata$treatment=gsub("Low","2Low",stat.hallmarkerscore.aza.GEOdata$treatment)
stat.hallmarkerscore.aza.GEOdata$treatment=gsub("Control","1Control",stat.hallmarkerscore.aza.GEOdata$treatment)
stat.hallmarkerscore.aza.GEOdata=stat.hallmarkerscore.aza.GEOdata[order(stat.hallmarkerscore.aza.GEOdata$treatment,stat.hallmarkerscore.aza.GEOdata$cell.type),]
#ht.hallM.sig.aza=as.data.frame(t(stat.hallmarkerscore.aza[,-c(1:3)]))[rownames(pvalue.hm.sig.aza),]
ht.hallM.sig.aza.GEOdata=as.data.frame(t(stat.hallmarkerscore.aza.GEOdata[,-c(1:3)]))    ####change
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
  cell.type=stat.hallmarkerscore.aza.GEOdata$cell.type,
  col=list(treatment=c("3High"="#cc4c02","2Low"="#E7B800", "1Control"="#FC4E07"),
           cell.type=c("Breat.cancer"="#8dd3c7","Bronchial.epithelial"="#ffffb3","Colon.cancer"="#bebada","Lung.cancer"="#fb8072"),
           cell.line=c("A549"="#8dd3c7","H1299"="#ffffb3","H157"="#bebada","H1819"="#80b1d3","H1993"="#fdb462","H2347"="#b3de69",
                       "H460"="#fccde5","H526"="#d9d9d9","HBEC2"="#fb8072", "HBEC2.Rep2"="#e31a1c", "HBEC3"="#bc80bd","HBEC3.Rep2"="#6a3d9a", 
                       "HBEC4"="#f1b6da", "HBEC4.Rep2"="#c51b7d", "HCT116"="#ccebc5","MCF7"="#ffed6f")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))


####
hmht.aza.GEOdata=Heatmap(ht.hallM.sig.aza.z.GEOdata, name = "hall.marker", 
                         #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                         #col = rev(viridis(10)),
                         col = viridis(10),
                         width = unit(5, "cm"),
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
                         column_title=paste0("hall.marker.genesets"," GSE5816 lung.colon.brca"),
                         column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/GSE5816/ht.hallM.ssgsea.aza.GSE5816.lungs.pdf",width = 8,height = 7)
draw(hmht.aza.GEOdata, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#save data
pvalue.hm.aza.GEOdata.GSE5816=pvalue.hm.aza.GEOdata
stat.hallmarkerscore.aza.GEOdata.GSE5816=stat.hallmarkerscore.aza.GEOdata
ht.hallM.sig.aza.GEOdata.GSE5816=ht.hallM.sig.aza.GEOdata
#
save(pvalue.hm.aza.GEOdata.GSE5816,stat.hallmarkerscore.aza.GEOdata.GSE5816,ht.hallM.sig.aza.GEOdata.GSE5816,file="9.AZAdata.GBM/GSE5816/ht.hallM.ssgsea.aza.GSE5816.lungs.RData")
#####################
##########only compare between Control and high dose group
#degs
##################degs
head(combat_edata.azaexp.GEOdata)
head(pdata.array.sel.ok)
dim(combat_edata.azaexp.GEOdata)
dim(pdata.array.sel.ok)
table(pdata.array.sel.ok$cell.line)
#
targets.GEOdata.cell=subset(pdata.array.sel.ok, treatment=="Control" | treatment=="High")
targets.GEOdata.cell=targets.GEOdata.cell[-c(7,18),]
targets.GEOdata=combat_edata.azaexp.GEOdata[,colnames(combat_edata.azaexp.GEOdata) %in% rownames(targets.GEOdata.cell)]

rownames(targets.GEOdata.cell)
colnames(targets.GEOdata)
table(targets.GEOdata.cell$cell.line,targets.GEOdata.cell$treatment)
####deg
SibShip <- factor(targets.GEOdata.cell$cell.line)
Treat <- factor(targets.GEOdata.cell$treatment, levels=c("Control","High"))
design.aza <- model.matrix(~SibShip+Treat)
fit.aza.GEOdata <- lmFit(targets.GEOdata, design.aza)
fit.aza.GEOdata <- eBayes(fit.aza.GEOdata)
##################
res.deg.aza.GEOdata<- topTable(fit.aza.GEOdata, adjust.method="BH", coef="TreatHigh",number = Inf,resort.by="p")
res.deg.aza.GEOdata$FC=2^(res.deg.aza.GEOdata$logFC)
head(res.deg.aza.GEOdata)
write.csv(res.deg.aza.GEOdata, "9.AZAdata.GBM/GSE5816/degs.genes.AZA.gse5816.Lungs.csv")
dim(res.deg.aza.GEOdata)
#
res.deg.aza.GEOdata.GSE5816=res.deg.aza.GEOdata
targets.GEOdata.cell.GSE5816=targets.GEOdata.cell
targets.GEOdata.GSE5816=targets.GEOdata
save(res.deg.aza.GEOdata.GSE5816, targets.GEOdata.cell.GSE5816, targets.GEOdata.GSE5816,file = "9.AZAdata.GBM/GSE5816/degs.genes.AZA.gse5816.Lungs.RData")
##################
##################
