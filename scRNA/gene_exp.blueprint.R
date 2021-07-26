
pbmc=readRDS("/home/zhuxq/data/scRNA/scPanCompare/blueprint.CRC/blueprint.CRC.count.rds")
pbmc=subset(pbmc, subset= Tissue!="N")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#saveRDS(pbmc, file="/home/zhuxq/data/scRNA/scPanCompare/lineage.CRC/GSE144735.KUL3/KUL3data.cohort.seuratnormalized.rds")
###
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)
saveRDS(pbmc,file = "/home/zhuxq/data/te/scRNA.cd274/blueprint/blueprint_onlyTumor.cohort.seurat.umap.rds")
#########
#########
metadata=pbmc@meta.data
metadata$celltype2=gsub("Cancer","Epithelial",metadata$celltype)
metadata$celltype2=gsub("EC","Fibroblast",metadata$celltype2)
metadata$celltype2=gsub("Enteric_glia","Fibroblast",metadata$celltype2)
table(metadata$celltype2)
pbmc=AddMetaData(pbmc,metadata = metadata$celltype2,col.name = "celltype2")
##
pdf("/home/zhuxq/data/te/scRNA.cd274/blueprint/umap.tsne.blueprintdata.onlyTumor.pdf",width = 8,height = 7)
DimPlot(pbmc, reduction = "umap",pt.size = 0.05#,split.by = "study",
        ,group.by = "celltype2"
        ,cols =rev(brewer.pal(8,"Set3"))
)
DimPlot(pbmc, reduction = "tsne",pt.size = 0.05#,split.by = "study",
        ,group.by = "celltype2"
        ,cols =rev(brewer.pal(8,"Set3"))
)
dev.off()
#####
#####
features.plot=c("CD79A","EPCAM","DCN","KIT", "CD68","CD3D", "JAK2","STAT1","STAT2","IRF1","CD274")
pdf("/home/zhuxq/data/te/scRNA.cd274/blueprint/stacked.violin.blueprintdata.pdf",width = 10,height = 10)
StackedVlnPlot(obj = subset(pbmc, subset= CD274>0), features = features.plot,group.by = "celltype2",cols=rev(brewer.pal(8,"Set3")))
StackedVlnPlot(pbmc, features = features.plot,group.by = "celltype2",cols=rev(brewer.pal(8,"Set3")))
dev.off()
####
a1sce=subset(pbmc, features = as.character(c("CD274","STAT1","STAT2","IRF1","JAK2","CD68","EPCAM","DCN","CD3D","CD79A","KIT")))
a1exp=a1sce[["RNA"]]@data
a1exp=as.data.frame(t(a1exp))
a1stat=cbind(a1exp, as.data.frame(pbmc@meta.data))
head(a1stat)
a1stat$CD274_group=ifelse(a1stat$CD274>0, "positive","negative")
write.csv(a1stat, "/home/zhuxq/data/te/scRNA.cd274/blueprint/percentage.cd274.positive.blueprint.csv")
####
a1stat=read.csv("/home/zhuxq/data/te/scRNA.cd274/blueprint/percentage.cd274.positive.blueprint.csv",header = T)
head(a1stat)
x=as.data.frame.matrix(table(a1stat$CD274_group,a1stat$celltype2))
x
# Transform this data in %
data_percentage <- apply(x, 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage = as.data.frame(data_percentage )
data_percentage$group=rownames(data_percentage)
colnames(data_percentage)=gsub(" ","_",colnames(data_percentage))
#
library(tidyr)
drawdata.per1 <- data_percentage %>% pivot_longer(cols=colnames(data_percentage)[1:6],
                                                  names_to= "celltype", values_to = "percentage")
# Make a stacked barplot--> it will be in %!
pdf("/home/zhuxq/data/te/scRNA.cd274/blueprint/percentage.cd274.positive.blueprint.pdf",width = 5,height = 7)
ggplot(drawdata.per1, aes(x = celltype, y = percentage, fill = group)) +
  scale_y_continuous(labels = scales::percent)+
  geom_col(position = "fill",colour = "white")+scale_fill_manual(values= c("#80cdc1", "#8c510a"))+
  theme_classic2()+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#############
#############
x=as.data.frame.matrix(table(a1stat$celltype2,a1stat$CD274_group))
x
# Transform this data in %
data_percentage <- apply(x, 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage = as.data.frame(data_percentage )
data_percentage$group=rownames(data_percentage)
colnames(data_percentage)=gsub(" ","_",colnames(data_percentage))
library(tidyr)
drawdata.per1 <- data_percentage[,-1] %>% pivot_longer(cols=colnames(data_percentage[,-1])[1],
                                                  names_to= "celltype", values_to = "percentage")
drawdata.per1=drawdata.per1[order(drawdata.per1$percentage,decreasing = F),]
drawdata.per1$group=factor(drawdata.per1$group,levels = drawdata.per1$group)
pdf("/home/zhuxq/data/te/scRNA.cd274/blueprint/percentage.cd274.positive.only.blueprint.coloset.pdf",width = 3,height = 7)
ggplot(drawdata.per1, aes(x = celltype, y = percentage, fill = group)) +
  scale_y_continuous(labels = scales::percent)+
  geom_col(position = "fill",colour = "white")+
  #scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  scale_fill_manual(values= c("B_cell"="#FCCDE5","Mast_cell"="#80B1D3","Fibroblast"="#FDB462",
                              "Epithelial"="#B3DE69","T_cell"="#BEBADA","Myeloid"="#FB8072"))+
  theme_classic2()+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#########
###calculate signature score
pbmc=readRDS("/home/zhuxq/data/te/scRNA.cd274/blueprint/blueprint_onlyTumor.cohort.seurat.umap.rds")
read_gmt <- function(file, tidy = FALSE) {
  con <- file(file, "r")
  gmt_lines <- readLines(file, warn = FALSE)
  close(con)
  rlist <- purrr::map(gmt_lines, parse_gmt_lines)
  rlist_names <- purrr::map_chr(gmt_lines, get_gmt_names)
  names(rlist) <- rlist_names
  if (tidy) rlist <- tidy_gmt(rlist)
  return(rlist)
}

parse_gmt_lines <- function(gl) {
  gl <- unlist(stringr::str_split(gl, "\\\t"), recursive = FALSE)
  gl <- gl[3:length(gl)]
  return(gl)
}

get_gmt_names <- function(gl) {
  gl <- unlist(stringr::str_split(gl, "\\\t"))
  gl_name <- gl[[1]]
  return(gl_name)
}

sel_gmt=read_gmt("/home/zhuxq/data/te/scRNA.cd274/c2.all.v7.2.symbols.sel.gmt")
sets=sel_gmt
datalist=list()
for (i in 1:length(sets)) {
  sce_marker=subset(pbmc, features = as.character(sets[[i]]))
  #scexp=sce_sub[["RNA"]]@counts
  scexp=sce_marker[["RNA"]]@data
  scexp=as.data.frame(scexp)
  #scexp[1:3,1:4]
  #dim(scexp)
  markerscore=as.data.frame(colMeans(scexp))
  colnames(markerscore)=names(sets[i])
  datalist[[i]]=markerscore
}
res.markerscore=do.call(cbind, datalist)
head(res.markerscore)
#
res.markerscore.stat=cbind(res.markerscore,as.data.frame(pbmc@meta.data))
write.csv(res.markerscore.stat, "/home/zhuxq/data/te/scRNA.cd274/blueprint/stacked.violin.all.signature.score.blueprint.csv")
######
res.markerscore.stat$Index= res.markerscore.stat$Cell
#####
pdf("/home/zhuxq/data/te/scRNA.cd274/blueprint/stacked.violin.Jakstat.IFN.score.blueprint.pdf",width = 10,height = 2)
ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, KEGG_JAK_STAT_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, WP_OVERVIEW_OF_INTERFERONSMEDIATED_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, REACTOME_INTERFERON_GAMMA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, REACTOME_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, REACTOME_INTERFERON_ALPHA_BETA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, ST_TYPE_I_INTERFERON_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, WP_TYPE_III_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()
dev.off()
######
pdf("/home/zhuxq/data/te/scRNA.cd274/blueprint/stacked.violin.Jakstat.IFN.score.blueprint.nobox.pdf",width = 10,height = 2)
ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, KEGG_JAK_STAT_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, WP_OVERVIEW_OF_INTERFERONSMEDIATED_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, REACTOME_INTERFERON_GAMMA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, REACTOME_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, REACTOME_INTERFERON_ALPHA_BETA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, ST_TYPE_I_INTERFERON_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype2, WP_TYPE_III_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype2)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype2))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()
dev.off()
###
###
####
####################
#AUCcell
#devtools::install_github("aertslab/AUCell")
library(AUCell)
##1.count matrix
##2.gene sets
##
exprMatrix=pbmc[["RNA"]]@counts
exprMatrix=as.matrix(exprMatrix)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=F)
save(cells_rankings, file="/home/zhuxq/data/te/scRNA.cd274/blueprint/Aucell.analysis/cells_rankings.blueprint.RData")
##
library(GSEABase)
sets= getGmt("/home/zhuxq/data/te/scRNA.cd274/c2.all.v7.2.symbols.gmt")
c1=getGmt("/home/zhuxq/data/te/scRNA.cd274/c2.all.v7.2.symbols.sel.gmt")
#geneSets <- list(geneSet1=sel_gmt$BIOCARTA_FEEDER_PATHWAY[1:5])
sets_sel=sets[names(c1)]
cells_AUC <- AUCell_calcAUC(sets_sel, cells_rankings,nCores = 5)
save(cells_AUC, file="/home/zhuxq/data/te/scRNA.cd274/blueprint/Aucell.analysis/cells_AUC.blueprint.RData")
###
load("/home/zhuxq/data/te/scRNA.cd274/blueprint/Aucell.analysis/cells_AUC.blueprint.RData")
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
######
aucscore.blueprint=as.data.frame(cells_AUC@assays@data$AUC)
aucscore.blueprint=as.data.frame(t(aucscore.blueprint))
aucscore.blueprint[1:4,1:4]
dim(aucscore.blueprint)
metacell=cbind(pbmc@meta.data, pbmc[["umap"]]@cell.embeddings, pbmc[["tsne"]]@cell.embeddings)
head(metacell)
stat.auc.blueprint=cbind(metacell, aucscore.blueprint[rownames(metacell),])
head(stat.auc.blueprint)
save(stat.auc.blueprint, file="/home/zhuxq/data/te/scRNA.cd274/blueprint/Aucell.analysis/stat.aucscore.blueprint.RData")
###
###
pdf("/home/zhuxq/data/te/scRNA.cd274/blueprint/Aucell.analysis/WP_TYPE_II_INTERFERON_SIGNALING_IFNG.pdf",width = 10,height = 4)
ggplot(subset(stat.auc.blueprint,Cell %in%  colnames(subset(pbmc, subset= CD274>0)) ), 
       aes(celltype, WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype2),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()
dev.off()
