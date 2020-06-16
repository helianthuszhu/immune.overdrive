#########################
#########################GSE41586
#####
expdata.gse41586=readRDS("9.AZAdata.GBM/GSE41586/RE.output/GENE_3_TPM.RDS")
expdata.gse41586=as.data.frame(expdata.gse41586$counts)
expdata.gse41586=cbind(Gene.stable.ID=rownames(expdata.gse41586),
                       expdata.gse41586)
head(expdata.gse41586)
#
biomartid=read.table("9.AZAdata.GBM/GSE41586/biomart_export.ensg.enst.id.infor.txt",header = T,sep = ",")
head(biomartid)
#########
length(intersect(rownames(expdata.gse41586), biomartid$Gene.stable.ID))
#########
expdata.gse41586=merge(biomartid[,c(1,6)], expdata.gse41586, by="Gene.stable.ID",all=FALSE)
head(expdata.gse41586)
expdata.gse41586.agg=expdata.gse41586[,-1] %>% group_by(Gene.name) %>% summarise_all(mean)
expdata.gse41586.agg=as.data.frame(expdata.gse41586.agg)
rownames(expdata.gse41586.agg)=expdata.gse41586.agg$Gene.name
expdata.gse41586.agg=expdata.gse41586.agg[,-1]
head(expdata.gse41586.agg)
####filter
####filter out
azanacount.gse41586=as.data.frame(rowSums(expdata.gse41586.agg == 0))
colnames(azanacount.gse41586)="count"
head(azanacount.gse41586)
azanacount.gse41586.sel=subset(azanacount.gse41586, count<=6)
dim(azanacount.gse41586.sel)
#########
expdata.gse41586.agg.sel=expdata.gse41586.agg[rownames(azanacount.gse41586.sel),]
dim(expdata.gse41586.agg.sel)
#########cell information
cell.line.info.gse41586=read.table("9.AZAdata.GBM/GSE41586/SRR_meta.GSE41586.txt",header = T,sep = ",")
head(cell.line.info.gse41586)
cell.line.info.gse41586=separate(cell.line.info.gse41586, col="source_name",into = c("cellName","d1","d2","dose","other"),sep = " ")
cell.line.info.gse41586.ok=data.frame(cell.line=cell.line.info.gse41586$cellName,
                                      dose=cell.line.info.gse41586$dose)
rownames(cell.line.info.gse41586.ok)=cell.line.info.gse41586$Run
cell.line.info.gse41586.ok$treatment=gsub("10","high",cell.line.info.gse41586.ok$dose)
cell.line.info.gse41586.ok$treatment=gsub("0","Control",cell.line.info.gse41586.ok$treatment)
cell.line.info.gse41586.ok$treatment=gsub("5","low",cell.line.info.gse41586.ok$treatment)

head(cell.line.info.gse41586.ok)
#########PCA
#########
aza.pca.gse41586 <- PCA(as.matrix(t(expdata.gse41586.agg.sel)), graph = FALSE)
#
pdf("9.AZAdata.GBM/GSE41586/pca.before.batch.gse41586.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.gse41586,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse41586.ok$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE41586.CRC.HT29")+theme_bw()
fviz_pca_ind(aza.pca.gse41586,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse41586.ok$dose, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE41586.CRC.HT29")+theme_bw()

dev.off()
#############
##########remove batch effect
library(sva)
library(pamr)
library(limma)
aza.batch.gse41586 = cell.line.info.gse41586.ok$treatment
aza.modcombat.gse41586 = model.matrix(~1, data=cell.line.info.gse41586.ok)
combat_edata.azaexp.gse41586 = ComBat(dat=as.matrix(expdata.gse41586.agg.sel), batch=aza.batch.gse41586, mod=aza.modcombat.gse41586, par.prior=TRUE, prior.plots=FALSE)
head(combat_edata.azaexp.gse41586)
save(combat_edata.azaexp.gse41586, file="9.AZAdata.GBM/GSE41586/batchcorrected/combat_edata.azaexp.gse41586.RData")
############PCA again
aza.pca.gse41586.after <- PCA(as.matrix(t(combat_edata.azaexp.gse41586)), graph = FALSE)
#
pdf("9.AZAdata.GBM/GSE41586/pca.after.batch.gse41586.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.gse41586.after,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse41586.ok$treatment, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GSE41586.CRC.HT29")+theme_bw()

dev.off()
###########
###########no need to do batch effect for GSE41586
###########
save(expdata.gse41586.agg.sel,cell.line.info.gse41586.ok, file="9.AZAdata.GBM/GSE41586/expdata.gse41586.agg.sel.noneedtodobatcheffect.RData")
#######run ssgsea
#
load("9.AZAdata.GBM/GSE41586/ssgsea/genesets.hallmarker.score.aza.gse41586.RData")
head(genesets.hall.marker.score.aza.gse41586)
#load("9.AZAdata.GBM/GSE41586/batchcorrected/genesets.hallmarker.score.aza.gse41586.batched.RData")
####
#####
stat.hallmarkerscore.aza.gse41586=cbind(cell.line.info.gse41586.ok, genesets.hall.marker.score.aza.gse41586[rownames(cell.line.info.gse41586.ok),])
stat.hallmarkerscore.aza.gse41586$treatment=as.factor(stat.hallmarkerscore.aza.gse41586$treatment)
head(stat.hallmarkerscore.aza.gse41586)
####calculate p value
datalist <- list()
for(i in names(stat.hallmarkerscore.aza.gse41586[,4:ncol(stat.hallmarkerscore.aza.gse41586)])){    ######change
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ treatment")), data = stat.hallmarkerscore.aza.gse41586)
}
pvalue.hm.aza.gse41586=do.call(rbind, datalist)
pvalue.hm.aza.gse41586=as.data.frame(pvalue.hm.aza.gse41586)

#pvalue.hm.sig.aza.GEOdata=subset(pvalue.hm.aza.GEOdata, p.value< 0.05)
########
########draw heatmap
stat.hallmarkerscore.aza.gse41586[1:4,1:4]
#ht.hallM.sig.aza=as.data.frame(t(stat.hallmarkerscore.aza[,-c(1:3)]))[rownames(pvalue.hm.sig.aza),]
ht.hallM.sig.aza.gse41586=as.data.frame(t(stat.hallmarkerscore.aza.gse41586[,-c(1:3)]))    ####change
#zscore
ht.hallM.sig.aza.z.gse41586=as.data.frame(t(scale(t(ht.hallM.sig.aza.gse41586))))
ht.hallM.sig.aza.z.gse41586[ht.hallM.sig.aza.z.gse41586< -2] <- -2
ht.hallM.sig.aza.z.gse41586[ht.hallM.sig.aza.z.gse41586>2] <- 2
ht.hallM.sig.aza.z.gse41586[1:4,1:4]

min_cor = min(as.vector(ht.hallM.sig.aza.z.gse41586))
max_cor = max(as.vector(ht.hallM.sig.aza.z.gse41586))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
###
row_ha.left.hm.aza.gse41586 = rowAnnotation(wilcox.pvalue=as.numeric(pvalue.hm.aza.gse41586$p.value),  #change
                                           col=list(
                                             #wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001,0.0001), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354")),
                                             wilcox.pvalue=colorRamp2(c(0.05,0.03,0.01,0.001), c("#f1b6da","#de77ae","#c51b7d","#8e0152"))
                                           ),show_annotation_name = FALSE)

col_ha_top.hm.aza.gse41586 = columnAnnotation(
  treatment=stat.hallmarkerscore.aza.gse41586$treatment,
  cell.line=stat.hallmarkerscore.aza.gse41586$cell.line,
  col=list(treatment=c("high"="#cc4c02","low"="#E7B800", "Control"="#FC4E07"),
           cell.line=c("HT29"="#8dd3c7")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))


####
hmht.aza.gse41586=Heatmap(ht.hallM.sig.aza.z.gse41586, name = "hall.marker", 
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
                         top_annotation = col_ha_top.hm.aza.gse41586,
                         #right_annotation = row_ha.right,
                         show_row_dend = F,show_column_dend = F,
                         #row_names_side = "left",
                         left_annotation = row_ha.left.hm.aza.gse41586,
                         column_title=paste0("hall.marker.genesets"," GSE41586 CRC HT29"),
                         column_title_gp = gpar(fontsize = 8)
)
pdf("9.AZAdata.GBM/GSE41586/ht.hallM.ssgsea.aza.GSE41586.CRC.HT29.pdf",width = 8,height = 7)
draw(hmht.aza.gse41586, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#save data
#
save(pvalue.hm.aza.gse41586,stat.hallmarkerscore.aza.gse41586,ht.hallM.sig.aza.gse41586,file="9.AZAdata.GBM/GSE41586/ht.hallM.ssgsea.aza.GSE41586.CRC.ht29.RData")
#####################