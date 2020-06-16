###############
###########63 cell lines   
###GSE57341
###########
azaexpm.gse57341=read.table("9.AZAdata.GBM/GSE57341/GSE57341_series_matrix.clean.probe.txt",sep = "\t",header = T)
dim(azaexpm.gse57341)
colnames(azaexpm.gse57341)[1]="ID"
azaexpm.gse57341[1:20,1:5]
###########
probe.gse57341=read.table("9.AZAdata.GBM/GSE57341/GPL4133_GSE57341.clean.txt",header = T,sep = "\t")
probe.gse57341[1:20,1:3]
##########
length(intersect(azaexpm.gse57341$ID, probe.gse57341$ID))
####
azaexpm.gse57341=merge(probe.gse57341, azaexpm.gse57341, by="ID",all=FALSE)
azaexpm.gse57341.sel=azaexpm.gse57341[!(is.na(azaexpm.gse57341$GENE)),]
azaexpm.gse57341.sel[1:4,1:7]
#####
azaexpm.gse57341.sel.agg=azaexpm.gse57341.sel[,-c(1:2)] %>% group_by(GENE_SYMBOL) %>% summarise_all(mean)
azaexpm.gse57341.sel.agg=as.data.frame(azaexpm.gse57341.sel.agg)
azaexpm.gse57341.sel.agg[1:5,1:4]
dim(azaexpm.gse57341.sel.agg)
#############change the wrong id
#############
syidx=c("Mar","Sep","Dec")
syidx.new=c("MARCH","SEPT","DEC")
for (i in 15:1) {
  for (j in 1:3) {
    azaexpm.gse57341.sel.agg$GENE_SYMBOL=gsub(paste(i,syidx[j],sep = "-"),
                                              paste0(syidx.new[j],i),
                                              azaexpm.gse57341.sel.agg$GENE_SYMBOL)
  }
}
###########
rownames(azaexpm.gse57341.sel.agg)=azaexpm.gse57341.sel.agg$GENE_SYMBOL
azaexpm.gse57341.sel.agg=azaexpm.gse57341.sel.agg[,-1]
class(azaexpm.gse57341.sel.agg)
azaexpm.gse57341.sel.agg.nadelete=na.omit(azaexpm.gse57341.sel.agg)
dim(azaexpm.gse57341.sel.agg.nadelete)
#save(azaexpm.gse57341.sel.agg.nadelete,file="9.AZAdata.GBM/GSE57341/azaexpm.gse57341.agg.corrected.Matrix.RData")
#write.csv(azaexpm.gse57341.sel.agg.nadelete,"9.AZAdata.GBM/GSE57341/azaexpm.gse57341.agg.corrected.Matrix.csv")
#############
#############
####### cell information
#
cell.line.info.gse57341=read.table("9.AZAdata.GBM/GSE57341/GSE57341_cell.line.information.txt",header = T,row.names = 1,sep = "\t")
head(cell.line.info.gse57341)
dim(cell.line.info.gse57341)
head(cell.line.info.gse57341)
#
cell.line.info.gse57341$cancer.type=gsub(" ","",cell.line.info.gse57341$cancer.type)
cell.line.info.gse57341$time.point=gsub(" ","",cell.line.info.gse57341$time.point)
#
length(intersect(rownames(cell.line.info.gse57341), colnames(azaexpm.gse57341.sel.agg)))
#########################
#########################
######PCA
######
#pca
library("FactoMineR")
library("factoextra")
####
#cell.line.info.gse57341.colon=subset(cell.line.info.gse57341,cancer.type=="colon cancer cell line")
azaexpm.gse57341.sel.agg.nadelete.matrix=as.numeric(as.matrix(azaexpm.gse57341.sel.agg.nadelete))
dim(azaexpm.gse57341.sel.agg.nadelete.matrix) <- dim(azaexpm.gse57341.sel.agg.nadelete)
#rownames(azaexpm.gse57341.sel.agg.nadelete.matrix)=rownames(azaexpm.gse57341.sel.agg.nadelete)
#colnames(azaexpm.gse57341.sel.agg.nadelete.matrix)=colnames(azaexpm.gse57341.sel.agg.nadelete)
azaexpm.gse57341.sel.agg.nadelete.matrix[1:4,1:4]
dim(azaexpm.gse57341.sel.agg.nadelete.matrix)
#######################
aza.pca.gse57341 <- PCA(as.matrix(t(azaexpm.gse57341.sel.agg.nadelete.matrix)), graph = FALSE)
#
#pdf("9.AZAdata.GBM/pca.before.batch.pdf",width = 5,height = 4)
fviz_pca_ind(aza.pca.gse57341,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse57341$cancer.type, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GBM.cell.line")+theme_bw()

fviz_pca_ind(aza.pca.gse57341,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cell.line.info.gse57341$time.point, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","black"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.GBM.cell.line")+theme_bw()
dev.off()
################################
################################
###load ssgasea score for hall marker genesets
#######
load("9.AZAdata.GBM/GSE57341/ssgsea/genesets.hallmarker.score.aza.gse57341.RData")
dim(genesets.hall.marker.score.aza.gse57341)
genesets.hall.marker.score.aza.gse57341[1:4,1:4]
######draw heatmap
stat.aza.gse57341=cbind(cell.line.info.gse57341,as.data.frame(genesets.hall.marker.score.aza.gse57341)[rownames(cell.line.info.gse57341),])
dim(stat.aza.gse57341)
stat.aza.gse57341[1:4,1:6]
######colon
stat.aza.gse57341.colon=subset(stat.aza.gse57341,cancer.type=="coloncancercellline")
stat.aza.gse57341.colon=cbind(treatment=ifelse(stat.aza.gse57341.colon$time.point=="Day3","mock","treated"), stat.aza.gse57341.colon)
stat.aza.gse57341.colon=cbind(treatment2=paste(stat.aza.gse57341.colon$cell.line,stat.aza.gse57341.colon$treatment,sep = "_"),stat.aza.gse57341.colon)
stat.aza.gse57341.colon.agg=stat.aza.gse57341.colon[,-c(2:7)] %>% group_by(treatment2) %>% summarise_all(mean)
stat.aza.gse57341.colon.agg=as.data.frame(stat.aza.gse57341.colon.agg)
rownames(stat.aza.gse57341.colon.agg)=stat.aza.gse57341.colon.agg$treatment2
stat.aza.gse57341.colon.agg.matrix=stat.aza.gse57341.colon.agg[,-1]
stat.aza.gse57341.colon.agg.matrix=as.data.frame(scale(stat.aza.gse57341.colon.agg.matrix))
stat.aza.gse57341.colon.agg.matrix[1:4,1:4]
stat.aza.gse57341.colon.agg.anno=separate(stat.aza.gse57341.colon.agg,col="treatment2",into = c("cell.line","treat"),sep = "_")
######brca
stat.aza.gse57341.brca=subset(stat.aza.gse57341,cancer.type=="breastcancercellline")
stat.aza.gse57341.brca=cbind(treatment=ifelse(stat.aza.gse57341.brca$time.point=="Day1","mock","treated"), stat.aza.gse57341.brca)
stat.aza.gse57341.brca=cbind(treatment2=paste(stat.aza.gse57341.brca$cell.line,stat.aza.gse57341.brca$treatment,sep = "_"),stat.aza.gse57341.brca)
stat.aza.gse57341.brca.agg=stat.aza.gse57341.brca[,-c(2:7)] %>% group_by(treatment2) %>% summarise_all(mean)
stat.aza.gse57341.brca.agg=as.data.frame(stat.aza.gse57341.brca.agg)
rownames(stat.aza.gse57341.brca.agg)=stat.aza.gse57341.brca.agg$treatment2
stat.aza.gse57341.brca.agg.matrix=stat.aza.gse57341.brca.agg[,-1]
stat.aza.gse57341.brca.agg.matrix=as.data.frame(scale(stat.aza.gse57341.brca.agg.matrix))
stat.aza.gse57341.brca.agg.matrix[1:4,1:4]
stat.aza.gse57341.brca.agg.anno=separate(stat.aza.gse57341.brca.agg,col="treatment2",into = c("cell.line","treat"),sep = "_")
######ovarian
stat.aza.gse57341.ova=subset(stat.aza.gse57341,cancer.type=="ovariancancercellline")
stat.aza.gse57341.ova=cbind(treatment=ifelse(stat.aza.gse57341.ova$time.point=="Day3","mock","treated"), stat.aza.gse57341.ova)
stat.aza.gse57341.ova=cbind(treatment2=paste(stat.aza.gse57341.ova$cell.line,stat.aza.gse57341.ova$treatment,sep = "_"),stat.aza.gse57341.ova)
stat.aza.gse57341.ova.agg= stat.aza.gse57341.ova[,-c(2:7)] %>% group_by(treatment2) %>% summarise_all(mean)
stat.aza.gse57341.ova.agg=as.data.frame(stat.aza.gse57341.ova.agg)
rownames(stat.aza.gse57341.ova.agg)=stat.aza.gse57341.ova.agg$treatment2
stat.aza.gse57341.ova.agg.matrix=stat.aza.gse57341.ova.agg[,-1]
stat.aza.gse57341.ova.agg.matrix=as.data.frame(scale(stat.aza.gse57341.ova.agg.matrix))
stat.aza.gse57341.ova.agg.matrix[1:4,1:4]
stat.aza.gse57341.ova.agg.anno=separate(stat.aza.gse57341.ova.agg,col="treatment2",into = c("cell.line","treat"),sep = "_")
####################
####################
#
col_ha_top.hm.aza.gse57341 = columnAnnotation(
  treatment=stat.aza.gse57341.ova.agg.anno$treat,
  #,
  cell.line=stat.aza.gse57341.ova.agg.anno$cell.line
  #cancer.type=cell.line.info.gse57341$cancer.type,
  #col=list(#cancer.type=c("breastcancercellline"="#E7B800", "coloncancercellline"="#FC4E07","ovariancancercellline"="#00AFBB"),
  #         $treatment=c("Day1"="#a6cee3","Day3"="#1f78b4","Day7"="#b2df8a","Day8"="#33a02c","Day9"="#fb9a99","Day10"="#e31a1c","Day14"="#fdbf6f","Day21"="#ff7f00","Day28"="#cab2d6")
  #),
  #show_annotation_name = FALSE,gp = gpar(col = "black"))
)
#
####
hmht.aza.gse57341=Heatmap(t(stat.aza.gse57341.ova.agg.matrix), name = "hall.marker", 
                 col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                 #col = rev(viridis(10)),
                 #col = viridis(10),
                 width = unit(15, "cm"),
                 height = unit(12, "cm"),
                 border = F,
                 #col=col.pal_cor,
                 show_column_names = T,show_row_names = T,
                 cluster_columns = T,
                 row_names_gp = gpar(fontsize = 5),
                 column_names_gp = gpar(fontsize = 5),
                 top_annotation = col_ha_top.hm.aza.gse57341,
                 #right_annotation = row_ha.right,
                 show_row_dend = F,show_column_dend = F,
                 #row_names_side = "left",
                 #left_annotation = row_ha.left.hm.aza,
                 column_title=paste0("hall.marker.genesets","\n","63 cell lines"," GSE57341"),
                 column_title_gp = gpar(fontsize = 8)
)
#
pdf("9.AZAdata.GBM/GSE57341/ssgsea/ht.hallM.50.ssgsea.aza.gse57341.ova.agg.pdf",width = 10,height = 7)
draw(hmht.aza.gse57341, #padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
##############
##############
save(stat.aza.gse57341.colon.agg.matrix,stat.aza.gse57341.brca.agg.matrix,stat.aza.gse57341.ova.agg.matrix,
     stat.aza.gse57341.colon.agg.anno,stat.aza.gse57341.brca.agg.anno,stat.aza.gse57341.ova.agg.anno,file="9.AZAdata.GBM/GSE57341/ssgsea/ht.hallM.50.ssgsea.aza.gse57341.agg.matrix.anno.RData")

