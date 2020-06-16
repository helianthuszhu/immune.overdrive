###########
azaexpm.gse41676=read.table("9.AZAdata.GBM/GSE41676/GSE41676_series_matrix.probe.clean.txt",sep = "\t",header = T)
#colnames(azaexpm.gse41676)[1]="ID"
azaexpm.gse41676[1:4,1:5]
dim(azaexpm.gse41676)
###########
probe.gse41676=read.table("9.AZAdata.GBM/GSE41676/GPL10558-50081.txt",header = T,sep = "\t")
probe.gse41676[1:4,1:2]
##########
length(intersect(azaexpm.gse41676$ID_REF, probe.gse41676$ID_REF))
####
azaexpm.gse41676=merge(probe.gse41676, azaexpm.gse41676, by="ID_REF",all=FALSE)
azaexpm.gse41676[1:4,1:4]
dim(azaexpm.gse41676)
azaexpm.gse41676.sel=azaexpm.gse41676[!(grepl("phage_lambda",azaexpm.gse41676$Symbol)),]
azaexpm.gse41676.sel=azaexpm.gse41676.sel[!(is.na(azaexpm.gse41676.sel$Symbol)),]
azaexpm.gse41676.sel[1:4,1:7]
dim(azaexpm.gse41676.sel)
#####
azaexpm.gse41676.sel.agg=azaexpm.gse41676.sel[,-1] %>% group_by(Symbol) %>% summarise_all(mean)
azaexpm.gse41676.sel.agg=as.data.frame(azaexpm.gse41676.sel.agg)[-1,]
azaexpm.gse41676.sel.agg[1:20,1:4]
dim(azaexpm.gse41676.sel.agg)
#############change the wrong id
#############
syidx=c("Mar","Sep","Dec")
syidx.new=c("MARCH","SEPT","DEC")
for (i in 15:1) {
  for (j in 1:3) {
    azaexpm.gse41676.sel.agg$Symbol=gsub(paste(i,syidx[j],sep = "-"),
                                              paste0(syidx.new[j],i),
                                              azaexpm.gse41676.sel.agg$Symbol)
  }
}
###########
rownames(azaexpm.gse41676.sel.agg)=azaexpm.gse41676.sel.agg$Symbol
azaexpm.gse41676.sel.agg=azaexpm.gse41676.sel.agg[,-1]
save(azaexpm.gse41676.sel.agg,file="9.AZAdata.GBM/GSE41676/azaexpm.gse41676.agg.corrected.Matrix.RData")
write.csv(azaexpm.gse41676.sel.agg,"9.AZAdata.GBM/GSE57341/azaexpm.gse41676.agg.corrected.Matrix.csv")
#############
#############
####### cell information
#
cell.line.info.gse41676=read.table("9.AZAdata.GBM/GSE41676/cell.line.information.GSE41676.txt",header = T,row.names = 1,sep = "\t")
head(cell.line.info.gse41676)
dim(cell.line.info.gse41676)
head(cell.line.info.gse41676)
######select only none and AZA samples
cell.line.info.gse41676.sel=subset(cell.line.info.gse41676, treatment=="none"| treatment=="AZA")
cell.line.info.gse41676.sel$sample=gsub(" ",".",cell.line.info.gse41676.sel$sample)
cell.line.info.gse41676.sel$sample=gsub("-",".",cell.line.info.gse41676.sel$sample)
cell.line.info.gse41676.sel$sample2=substr(cell.line.info.gse41676.sel$sample,1,nchar(as.character(cell.line.info.gse41676.sel$sample))-2)
cell.line.info.gse41676.sel$cell.line=gsub("-","",cell.line.info.gse41676.sel$cell.line)
head(cell.line.info.gse41676.sel)
azaexpm.gse41676.sel.agg.sel=azaexpm.gse41676.sel.agg[, colnames(azaexpm.gse41676.sel.agg) %in% rownames(cell.line.info.gse41676.sel)]
azaexpm.gse41676.sel.agg.sel=na.omit(azaexpm.gse41676.sel.agg.sel)

#
length(intersect(rownames(cell.line.info.gse41676.sel), colnames(azaexpm.gse41676.sel.agg.sel)))
#############
######PCA
######
#pca
library("FactoMineR")
library("factoextra")
aza.pca.gse41676 <- PCA(as.matrix(t(azaexpm.gse41676.sel.agg.sel)), graph = FALSE)
#
pdf("9.AZAdata.GBM/GSE41676/pca.before.batch.gse41676.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.gse41676,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse41676.sel$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE41676")+theme_bw()

fviz_pca_ind(aza.pca.gse41676,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse41676.sel$cell.line, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","black"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE41676")+theme_bw()
dev.off()
###############
################
##########remove batch effect
library(sva)
library(pamr)
library(limma)
aza.batch.gse41676 = cell.line.info.gse41676.sel$cell.line
aza.modcombat.gse41676 = model.matrix(~1, data=cell.line.info.gse41676.sel)
combat_edata.azaexp.gse41676 = ComBat(dat=as.matrix(azaexpm.gse41676.sel.agg.sel), batch=aza.batch.gse41676, mod=aza.modcombat.gse41676, par.prior=TRUE, prior.plots=FALSE)
head(combat_edata.azaexp.gse41676)
save(combat_edata.azaexp.gse41676,file="9.AZAdata.GBM/GSE41676/azaexpm.gse41676.batch.removed.RData")
#############
#############
aza.pca.after.gse41676 <- PCA(as.matrix(t(combat_edata.azaexp.gse41676)), graph = FALSE)
pdf("9.AZAdata.GBM/GSE41676/pca.after.batch.gse41676.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.after.gse41676,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse41676.sel$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE41676")+theme_bw()

fviz_pca_ind(aza.pca.after.gse41676,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse41676.sel$cell.line, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","black"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.cell.line.GSE41676")+theme_bw()
dev.off()
############ssgasea for hallmarker
############
load("9.AZAdata.GBM/GSE41676/ssgsea/genesets.hallmarker.score.aza.gse41676.RData")
head(genesets.hall.marker.score.aza.gse41676)
####
#####
stat.hallmarkerscore.aza.gse41676=cbind(cell.line.info.gse41676.sel, genesets.hall.marker.score.aza.gse41676[rownames(cell.line.info.gse41676.sel),])
head(stat.hallmarkerscore.aza.gse41676)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza.gse41676[,5:ncol(stat.hallmarkerscore.aza.gse41676)])){  
  datalist[[i]] <- wilcox.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza.gse41676)
}
pvalue.hm.aza.gse41676=do.call(rbind, datalist)
pvalue.hm.aza.gse41676=as.data.frame(pvalue.hm.aza.gse41676)

pvalue.hm.sig.aza.gse41676=subset(pvalue.hm.aza.gse41676, p.value< 0.05)
########
########draw heatmap
stat.hallmarkerscore.aza.gse41676[1:4,1:4]
ht.hallM.sig.aza.gse41676=as.data.frame(t(stat.hallmarkerscore.aza.gse41676[,-c(1:4)]))[rownames(pvalue.hm.sig.aza.gse41676),]
#zscore
ht.hallM.sig.aza.z.gse41676=as.data.frame(t(scale(t(ht.hallM.sig.aza.gse41676))))
ht.hallM.sig.aza.z.gse41676[ht.hallM.sig.aza.z.gse41676< -2] <- -2
ht.hallM.sig.aza.z.gse41676[ht.hallM.sig.aza.z.gse41676>2] <- 2
ht.hallM.sig.aza.z.gse41676[1:4,1:4]

min_cor = min(as.vector(ht.hallM.sig.aza.gse41676))
max_cor = max(as.vector(ht.hallM.sig.aza.gse41676))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
###
row_ha.left.hm.aza.gse41676 = rowAnnotation(wilcox.pvalue=as.numeric(pvalue.hm.sig.aza.gse41676$p.value),
                                   col=list(
                                     #wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001,0.0001), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354")),
                                     wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001), c("#f1b6da","#de77ae","#c51b7d","#8e0152"))
                                   ),show_annotation_name = FALSE)

col_ha_top.hm.aza = columnAnnotation(
  treatment=stat.hallmarkerscore.aza.gse41676$treatment,
  cell.line=stat.hallmarkerscore.aza.gse41676$cell.line,
  col=list(treatment=c("AZA"="#E7B800", "none"="#FC4E07"),
           cell.line=c("T24"="#00AFBB","HL60"="#e41a1c","HCT116"="#377eb8")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
hmht.aza.gse41676=Heatmap(ht.hallM.sig.aza.z.gse41676, name = "hall.marker", 
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
                 left_annotation = row_ha.left.hm.aza.gse41676,
                 column_title=paste0("hall.marker.genesets","\n","(significant 12 out of 50)"," GSE41676"),
                 column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/GSE41676/ssgsea//ht.hallM.ssgsea.aza.gse41676.pdf",width = 8,height = 6)
draw(hmht.aza.gse41676, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#
















##############################degs
Treat <- factor(paste(cell.line.info.gse41676.sel$cell.line,cell.line.info.gse41676.sel$treatment,sep="."))
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

corfit <- duplicateCorrelation(combat_edata.azaexp.gse41676,design,block=cell.line.info.gse41676.sel$sample2)
#corfit$consensus
fit <- lmFit(combat_edata.azaexp.gse41676,design,block=cell.line.info.gse41676.sel$sample2,correlation=corfit$consensus)
cm <- makeContrasts(
   DiseasedvsNormalForTissueA = Diseased.A-Normal.A,
   DiseasedvsNormalForTissueB = Diseased.B-Normal.B,
   TissueAvsTissueBForNormal = Normal.B-Normal.A,
   TissueAvsTissueBForDiseased = Diseased.B-Diseased.A,
   levels=design)
> fit2 <- contrasts.fit(fit, cm)
> fit2 <- eBayes(fit2)
topTable(fit2, coef="DiseasedvsNormalForTissueA")



Treat <- factor(paste(targets$Condition,targets$Tissue,sep="."))
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)
> corfit <- duplicateCorrelation(eset,design,block=targets$Subject)
> corfit$consensus
fit <- lmFit(eset,design,block=targets$Subject,correlation=corfit$consensus)
> cm <- makeContrasts(
  + DiseasedvsNormalForTissueA = Diseased.A-Normal.A,
  + DiseasedvsNormalForTissueB = Diseased.B-Normal.B,
  + TissueAvsTissueBForNormal = Normal.B-Normal.A,
  + TissueAvsTissueBForDiseased = Diseased.B-Diseased.A,
  + levels=design)
> fit2 <- contrasts.fit(fit, cm)
> fit2 <- eBayes(fit2)
topTable(fit2, coef="DiseasedvsNormalForTissueA")
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