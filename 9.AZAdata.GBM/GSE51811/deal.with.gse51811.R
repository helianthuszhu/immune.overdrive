###########
azaexpm.gse51811=read.table("9.AZAdata.GBM/GSE51811/GSE51811_series_matrix.probe.clean.txt",sep = "\t",header = T)
#colnames(azaexpm.gse41676)[1]="ID"
azaexpm.gse51811[1:4,1:5]
dim(azaexpm.gse51811)
###########
probe.gse51811=read.table("9.AZAdata.GBM/GSE41676/GPL10558-50081.txt",header = T,sep = "\t")
probe.gse51811[1:4,1:2]
##########
length(intersect(azaexpm.gse51811$ID_REF, probe.gse51811$ID_REF))
####
azaexpm.gse51811=merge(probe.gse51811, azaexpm.gse51811, by="ID_REF",all=FALSE)
azaexpm.gse51811[1:4,1:4]
dim(azaexpm.gse51811)
azaexpm.gse51811.sel=azaexpm.gse51811[!(grepl("phage_lambda",azaexpm.gse51811$Symbol)),]
azaexpm.gse51811.sel=azaexpm.gse51811.sel[!(is.na(azaexpm.gse51811.sel$Symbol)),]
azaexpm.gse51811.sel[1:4,1:7]
dim(azaexpm.gse51811.sel)
#####
azaexpm.gse51811.sel.agg=azaexpm.gse51811.sel[,-1] %>% group_by(Symbol) %>% summarise_all(mean)
azaexpm.gse51811.sel.agg=as.data.frame(azaexpm.gse51811.sel.agg)[-1,]
azaexpm.gse51811.sel.agg.logd=log(azaexpm.gse51811.sel.agg)
azaexpm.gse51811.sel.agg.logd=na.omit(azaexpm.gse51811.sel.agg.logd)
azaexpm.gse51811.sel.agg.logd[1:4,1:4]
save(azaexpm.gse51811.sel.agg.logd,file="9.AZAdata.GBM/GSE51811/azaexpm.gse51811.sel.agg.logd.RData")
#

dim(azaexpm.gse51811.sel.agg.logd)
#############change the wrong id
#############
syidx=c("Mar","Sep","Dec")
syidx.new=c("MARCH","SEPT","DEC")
for (i in 15:1) {
  for (j in 1:3) {
    azaexpm.gse51811.sel.agg$Symbol=gsub(paste(i,syidx[j],sep = "-"),
                                         paste0(syidx.new[j],i),
                                         azaexpm.gse51811.sel.agg$Symbol)
  }
}
###########
rownames(azaexpm.gse51811.sel.agg)=azaexpm.gse51811.sel.agg$Symbol
azaexpm.gse51811.sel.agg=azaexpm.gse51811.sel.agg[,-1]
save(azaexpm.gse51811.sel.agg,file="9.AZAdata.GBM/GSE51811/azaexpm.gse51811.agg.corrected.Matrix.RData")
write.csv(azaexpm.gse51811.sel.agg,"9.AZAdata.GBM/GSE51811/azaexpm.gse51811.agg.corrected.Matrix.csv")
#############
#############
####### cell information
#
cell.line.info.gse51811=read.table("9.AZAdata.GBM/GSE51811/GSE51811.cell.info.txt",header = T,row.names = 1,sep = "\t")
rownames(cell.line.info.gse51811)=cell.line.info.gse51811$ID_REF
head(cell.line.info.gse51811)
dim(cell.line.info.gse51811)
#
length(intersect(rownames(cell.line.info.gse51811), colnames(azaexpm.gse51811.sel.agg)))
#############
######PCA
######
#pca
library("FactoMineR")
library("factoextra")
aza.pca.gse51811 <- PCA(as.matrix(t(na.omit(azaexpm.gse51811.sel.agg.logd))), graph = FALSE)
#
pdf("9.AZAdata.GBM/GSE51811/pca.before.batch.gse51811.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.gse51811,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse51811$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE51811.CRC.HCT116")+theme_bw()

fviz_pca_ind(aza.pca.gse51811,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse51811$time.point, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","black"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE51811.CRC.HCT116")+theme_bw()
dev.off()
########################ssgsea 
load("9.AZAdata.GBM/GSE51811/ssgsea/genesets.hallmarker.score.aza.GSE51811.RData")

head(genesets.hall.marker.score.aza.GSE51811)
########################
#####
stat.hallmarkerscore.aza.gse51811=cbind(cell.line.info.gse51811, genesets.hall.marker.score.aza.GSE51811[rownames(cell.line.info.gse51811),])
head(stat.hallmarkerscore.aza.gse51811)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza.gse51811[,5:ncol(stat.hallmarkerscore.aza.gse51811)])){  
  datalist[[i]] <- wilcox.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza.gse51811)
}
pvalue.hm.aza.gse51811=do.call(rbind, datalist)
pvalue.hm.aza.gse51811=as.data.frame(pvalue.hm.aza.gse51811)

#pvalue.hm.sig.aza.gse51811=subset(pvalue.hm.aza.gse51811, p.value< 0.05)

########draw heatmap
stat.hallmarkerscore.aza.gse51811[1:4,1:5]
ht.hallM.sig.aza.gse51811=as.data.frame(t(stat.hallmarkerscore.aza.gse51811[,-c(1:4)]))
#zscore
ht.hallM.sig.aza.z.gse51811=as.data.frame(t(scale(t(ht.hallM.sig.aza.gse51811))))
ht.hallM.sig.aza.z.gse51811[ht.hallM.sig.aza.z.gse51811< -2] <- -2
ht.hallM.sig.aza.z.gse51811[ht.hallM.sig.aza.z.gse51811>2] <- 2
ht.hallM.sig.aza.z.gse51811[1:4,1:4]

min_cor = min(as.vector(ht.hallM.sig.aza.z.gse51811))
max_cor = max(as.vector(ht.hallM.sig.aza.z.gse51811))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
###
row_ha.left.hm.aza.gse51811 = rowAnnotation(wilcox.pvalue=as.numeric(pvalue.hm.aza.gse51811$p.value),
                                            col=list(
                                              #wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001,0.0001), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354")),
                                              wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001), c("#f1b6da","#de77ae","#c51b7d","#8e0152"))
                                            ),show_annotation_name = FALSE)

col_ha_top.hm.aza = columnAnnotation(
  treatment=stat.hallmarkerscore.aza.gse51811$treatment,
  time.point=stat.hallmarkerscore.aza.gse51811$time.point,
  col=list(treatment=c("aza"="#E7B800", "untreatment"="#FC4E07"),
           time.point=c("baseline"="#d6604d","D5"="#92c5de","D14"="#4393c3","D24"="#2166ac","D42"="#053061")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
hmht.aza.gse51811=Heatmap(ht.hallM.sig.aza.z.gse51811, name = "hall.marker", 
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
                          left_annotation = row_ha.left.hm.aza.gse51811,
                          column_title=paste0("hall.marker.genesets"," GSE51811.HCT116"),
                          column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/GSE51811/ssgsea/ht.hallM.ssgsea.aza.gse51811.pdf",width = 8,height = 6)
draw(hmht.aza.gse51811, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#

#load("9.AZAdata.GBM/GSE51811/neqc.mored/expdata.gse51811.neqcnormed.RData")
#############
#############
library(GEOquery)
library(limma)
GSE51811 <- getGEO('GSE51811', destdir=".",getGPL = F)
exprSet=exprs(GSE51811[[1]])
GSE51811[[1]]
pdata.array=pData(GSE51811[[1]])
exprSet=exprs(GSE51811[[1]])
exprSet[1:4,1:4]
dim(exprSet)
###############
library(limma)
x <- read.ilmn("GSE51811_non-normalized_samples_1-10.txt.gz",expr="SAMPLE ",probeid="ID_REF")
y <- neqc(x)
########
#BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)
#####
#illuminaHumanv4SYMBOL
#probe <- illuminaHumanv4SYMBOL
# Get the probe identifiers that are mapped to a gene symbol
#mapped_probes <- mappedkeys(probe)
# Convert to a list
#xx <- as.list(probe[mapped_probes])
#########################
ids <- rownames(y$E)
qual <- unlist(mget(ids,illuminaHumanv4SYMBOL, ifnotfound = NA))
qual=as.data.frame(qual)
head(qual)
dim(qual)
#########################
dataexp=as.data.frame(y$E)
dataexp[1:4,1:4]
dim(dataexp)
colnames(dataexp)=colnames(exprSet)[1:10]
#
dataexp=cbind(qual,dataexp)
dataexp=na.omit(dataexp, cols="qual")
library(dplyr)
dataexp.agg= dataexp %>% group_by(qual) %>% summarise_all(mean)
dataexp.agg=as.data.frame(dataexp.agg)
rownames(dataexp.agg)=dataexp.agg$qual
dataexp.agg=dataexp.agg[,-1]
dataexp.agg[1:4,1:4]
######
pdata.sel=pdata.array[colnames(dataexp.agg),]
head(pdata.sel)
#########
expdata.gse51811.neqcnormed=dataexp.agg
save(expdata.gse51811.neqcnormed,file="9.AZAdata.GBM/GSE51811/neqc.mored/expdata.gse51811.neqcnormed.RData")
###########
###########load expdata
load("9.AZAdata.GBM/GSE51811/neqc.mored/expdata.gse51811.neqcnormed.RData")
head(expdata.gse51811.neqcnormed)
###########PCA
aza.pca.gse51811.neqcnormed <- PCA(as.matrix(t(expdata.gse51811.neqcnormed)), graph = FALSE)

pdf("9.AZAdata.GBM/GSE51811/neqc.mored/pca.before.batch.gse51811.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.gse51811.neqcnormed,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse51811$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE51811.CRC.HCT116")+theme_bw()

fviz_pca_ind(aza.pca.gse51811.neqcnormed,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse51811$time.point, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","black"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE51811.CRC.HCT116")+theme_bw()
dev.off()
################################
#################remove batch effect of time point
#
library(sva)
library(pamr)
library(limma)
aza.batch.gse51881 = cell.line.info.gse51811$time.point
aza.modcombat.gse51881 = model.matrix(~1, data=cell.line.info.gse51811)
combat_edata.azaexp.gse51881 = ComBat(dat=as.matrix(expdata.gse51811.neqcnormed), batch=aza.batch.gse51881, mod=aza.modcombat.gse51881, par.prior=TRUE, prior.plots=FALSE)
head(combat_edata.azaexp.gse51881)
save(combat_edata.azaexp.gse51881, file="9.AZAdata.GBM/GSE51811/neqc.mored/batch.corrected.necq/combat_edata.azaexp.gse51881.normnecq.RData")
######PCA again

aza.pca.gse51811.neqcnormed.after <- PCA(as.matrix(t(combat_edata.azaexp.gse51881)), graph = FALSE)

pdf("9.AZAdata.GBM/GSE51811/neqc.mored/pca.after.batch.gse51811.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.gse51811.neqcnormed.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse51811$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE51811.CRC.HCT116")+theme_bw()

fviz_pca_ind(aza.pca.gse51811.neqcnormed.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse51811$time.point, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","black"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE51811.CRC.HCT116")+theme_bw()
dev.off()

###############run ssgsea
#
load("9.AZAdata.GBM/GSE51811/neqc.mored/genesets.hallmarker.score.aza.GSE51811.necqnormed.RData")
head(genesets.hall.marker.score.aza.GSE51811.necqnormed)

#load("9.AZAdata.GBM/GSE51811/neqc.mored/batch.corrected.necq/genesets.hallmarker.score.aza.gse51881.normnecq.RData")
#head(genesets.hall.marker.score.aza.GSE51811.necqnormed.batchcorrected)
########################
#####
stat.hallmarkerscore.aza.gse51811.necqnormed=cbind(cell.line.info.gse51811, genesets.hall.marker.score.aza.GSE51811.necqnormed[rownames(cell.line.info.gse51811),])
#stat.hallmarkerscore.aza.gse51811.necqnormed=cbind(cell.line.info.gse51811, genesets.hall.marker.score.aza.GSE51811.necqnormed.batchcorrected[rownames(cell.line.info.gse51811),])
head(stat.hallmarkerscore.aza.gse51811.necqnormed)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza.gse51811.necqnormed[,5:ncol(stat.hallmarkerscore.aza.gse51811.necqnormed)])){  
  datalist[[i]] <- wilcox.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza.gse51811.necqnormed)
}
pvalue.hm.aza.gse51811.necqnormed=do.call(rbind, datalist)
pvalue.hm.aza.gse51811.necqnormed=as.data.frame(pvalue.hm.aza.gse51811.necqnormed)

#pvalue.hm.sig.aza.gse51811=subset(pvalue.hm.aza.gse51811, p.value< 0.05)

########draw heatmap
stat.hallmarkerscore.aza.gse51811.necqnormed[1:4,1:5]
stat.hallmarkerscore.aza.gse51811.necqnormed=cbind(time=stat.hallmarkerscore.aza.gse51811.necqnormed$time.point,stat.hallmarkerscore.aza.gse51811.necqnormed)
stat.hallmarkerscore.aza.gse51811.necqnormed$time=gsub("baseline","D1",stat.hallmarkerscore.aza.gse51811.necqnormed$time)

stat.hallmarkerscore.aza.gse51811.necqnormed$time=substr(stat.hallmarkerscore.aza.gse51811.necqnormed$time, 2, nchar(as.character(stat.hallmarkerscore.aza.gse51811.necqnormed$time)))

stat.hallmarkerscore.aza.gse51811.necqnormed$time=as.numeric(paste(stat.hallmarkerscore.aza.gse51811.necqnormed$time))
stat.hallmarkerscore.aza.gse51811.necqnormed=stat.hallmarkerscore.aza.gse51811.necqnormed[order(
                                                                                                stat.hallmarkerscore.aza.gse51811.necqnormed$time,decreasing = F),]
#save(stat.hallmarkerscore.aza.gse51811.necqnormed,file="9.AZAdata.GBM/GSE51811/neqc.mored/batch.corrected.necq/stat.hallmarkerscore.aza.gse51811.necqnormed.RData")
save(stat.hallmarkerscore.aza.gse51811.necqnormed,pvalue.hm.aza.gse51811.necqnormed,file="9.AZAdata.GBM/GSE51811/neqc.mored/stat.hallmarkerscore.aza.gse51811.necqnormed.RData")

ht.hallM.sig.aza.gse51811.necqnormed=as.data.frame(t(stat.hallmarkerscore.aza.gse51811.necqnormed[,-c(1:5)]))
#zscore
ht.hallM.sig.aza.z.gse51811.necqnormed=as.data.frame(t(scale(t(ht.hallM.sig.aza.gse51811.necqnormed))))
ht.hallM.sig.aza.z.gse51811.necqnormed[ht.hallM.sig.aza.z.gse51811.necqnormed< -2] <- -2
ht.hallM.sig.aza.z.gse51811.necqnormed[ht.hallM.sig.aza.z.gse51811.necqnormed>2] <- 2
ht.hallM.sig.aza.z.gse51811.necqnormed[1:4,1:4]


min_cor = min(as.vector(ht.hallM.sig.aza.z.gse51811.necqnormed))
max_cor = max(as.vector(ht.hallM.sig.aza.z.gse51811.necqnormed))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
###
row_ha.left.hm.aza.gse51811.necqnormed = rowAnnotation(wilcox.pvalue=as.numeric(pvalue.hm.aza.gse51811.necqnormed$p.value),
                                            col=list(
                                              #wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001,0.0001), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354")),
                                              wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001), c("#f1b6da","#de77ae","#c51b7d","#8e0152"))
                                            ),show_annotation_name = FALSE)

col_ha_top.hm.aza.necqnormed = columnAnnotation(
  treatment=stat.hallmarkerscore.aza.gse51811.necqnormed$treatment,
  time.point=stat.hallmarkerscore.aza.gse51811.necqnormed$time.point,
  col=list(treatment=c("aza"="#E7B800", "untreatment"="#FC4E07"),
           time.point=c("baseline"="#d6604d","D5"="#92c5de","D14"="#4393c3","D24"="#2166ac","D42"="#053061")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
hmht.aza.gse51811.necqnormed=Heatmap(ht.hallM.sig.aza.z.gse51811.necqnormed, name = "hall.marker", 
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
                          top_annotation = col_ha_top.hm.aza.necqnormed,
                          #right_annotation = row_ha.right,
                          show_row_dend = F,show_column_dend = F,
                          #row_names_side = "left",
                          left_annotation = row_ha.left.hm.aza.gse51811.necqnormed,
                          column_title=paste0("hall.marker.genesets"," GSE51811.HCT116"),
                          column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/GSE51811/neqc.mored/ht.hallM.ssgsea.aza.gse51811.necqnormed.pdf",width = 8,height = 6)
#pdf("9.AZAdata.GBM/GSE51811/neqc.mored/batch.corrected.necq/ht.hallM.ssgsea.aza.gse51811.necqnormed.batched.pdf",width = 8,height = 6)

draw(hmht.aza.gse51811.necqnormed, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#############################
##############degs using neqc normed data not adjusted by batch effect
#
head(expdata.gse51811.neqcnormed)
head(cell.line.info.gse51811)
#####simply compare
library(limma)
design2=model.matrix(~0+ factor(cell.line.info.gse51811$treatment))
colnames(design2)=c('aza','untreatment')
head(design2)
fit=lmFit(expdata.gse51811.neqcnormed,design2)
cont.matrix=makeContrasts('aza-untreatment',levels = design2)
fit=contrasts.fit(fit,cont.matrix)
fit=eBayes(fit)
#
res.deg.aza.gse51811<- topTable(fit, adjust.method="BH", coef=1,number = Inf,resort.by="p")
res.deg.aza.gse51811$FC=2^(res.deg.aza.gse51811$logFC)
head(res.deg.aza.gse51811)
#############
############
write.csv(res.deg.aza.gse51811, "9.AZAdata.GBM/GSE51811/neqc.mored/degs.genes.AZA.gse51811.hct116.csv")
save(res.deg.aza.gse51811,expdata.gse51811.neqcnormed, cell.line.info.gse51811,file="9.AZAdata.GBM/GSE51811/neqc.mored/degs.genes.AZA.gse51811.hct116.RData" )
