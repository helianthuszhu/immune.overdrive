
########
########
fi<-list.files(path = 'zhu/5.TE/scRNA/E-MTAB-8410/scTE',full.names=T,pattern = 'csv')
fi
fisample= sapply(stringr::str_split(fi, "_"), `[`, 2)
fisample= sapply(stringr::str_split(fisample, "[.]"), `[`, 1)
fisample
########
clin= read.table("zhu/5.TE/scRNA/E-MTAB-8410/clean_sampleinfor.txt", header = F, sep='\t')
colnames(clin)=c('tissue','sample')
rownames(clin)= sapply(stringr::str_split(clin$sample, "_"), `[`, 1)
clin$pt= sapply(stringr::str_split(clin$tissue, "[.]"), `[`, 1)
clin$hist= sapply(stringr::str_split(clin$tissue, "[.]"), `[`, 2)
clin$hist= gsub('tumour_core','T',clin$hist)
clin$hist= gsub('tumour_border','B',clin$hist)
clin$hist= gsub('normal_tissue','N',clin$hist)
clin$cdx= paste(clin$pt, clin$hist, sep = '-')
head(clin)

clin= clin[fisample,]

datalist = list()
metalist=list()
#KUL01-T_AAACCTGGTCTTTCAT

for (k in 1:length(fi)){
  require(data.table)
  a <- fread(fi[k],sep = ',',header = T, check.names = F)
  a$barcodes= sapply(stringr::str_split(a$barcodes, "[-]"), `[`, 1)
  a$barcodes= paste(rep(clin$cdx[k], times= nrow(a)),a$barcodes, sep = '_')
  datalist[[k]] <- a
  #
  b= data.frame(Index= a$barcodes, ID= rep(fisample[k],times= nrow(a)), tissue= rep(clin$tissue[k],times= nrow(a)))
  metalist[[k]]= b
  #
}
######
scTEgene_matrix= do.call(rbind, datalist)
scTEgene_matrix=as.data.frame(scTEgene_matrix)
rownames(scTEgene_matrix)= scTEgene_matrix$barcodes
scTEgene_matrix= scTEgene_matrix[,-1]
scTEgene_matrix= t(scTEgene_matrix)
scTEgene_matrix= as.data.frame(scTEgene_matrix)
scTEgene_matrix[1:4,1:4]
#
scTEgene_meta= do.call(rbind, metalist)
row.names(scTEgene_meta)= scTEgene_meta$Index
head(scTEgene_meta)
save(scTEgene_matrix, scTEgene_meta, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/merged.countsANDmeta.RData')
####
####
library(Seurat)
sce= CreateSeuratObject(scTEgene_matrix, project = "scTE_KUL_cohort", assay = "RNA",
                   min.cells = 0, min.features = 0,  meta.data = scTEgene_meta)
sce
saveRDS(sce, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/sceurat_merged_scTEgene.rds' )
############
#############
#####***
#####*
#####
########
########
hg38TE= read.table('zhu/5.TE/scRNA/ref/hg38_rmsk.txt.gz', header=F, sep='\t')
head(hg38TE)
length(unique(hg38TE$V11))
#####run in conda
library('rtracklayer')
gtf <- rtracklayer::import('/home/setup/zhu/5.TE/scRNA/ref/gencode.v30.annotation.gtf')
gtf_df=as.data.frame(gtf)
save(gtf_df, file='/home/setup/zhu/5.TE/scRNA/ref/gencode.v30.annotation.gtf.data.frame.RData' )
######
load('/home/setup/zhu/5.TE/scRNA/ref/gencode.v30.annotation.gtf.data.frame.RData')
head(gtf_df)
gtf_df_sel= gtf_df[!(duplicated(gtf_df$gene_name)),]
head(gtf_df_sel)
dim(gtf_df_sel)
#
length(intersect(rownames(sce), gtf_df_sel$gene_name))
idxgene= gtf_df_sel$gene_name
idxTE= rownames(sce)[!(rownames(sce) %in% idxgene)]
length(idxgene)
length(idxTE)
#
######
######
sce_gene= sce[idxgene, ]
sce_te= sce[idxTE,]

sce_te
head(sce_gene@meta.data)
head(sce_te@meta.data)
head(sce@meta.data)
##
##
saveRDS(sce_gene, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/sceurat_merged_sce_gene.rds')
saveRDS(sce_te, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/sceurat_merged_sce_te.rds')
###
###########normalize the data
#
sce
head(sce@meta.data)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce_filter <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.mt < 20)
sce_filter
#
sce_filter <- NormalizeData(sce_filter)
######
######
nineTE= c('AluSq4',	'HERV1_LTRd',	'LTR21B',	'MER57F',	'MER65C',	'MER92-int',	'SVA_C',	'SVA_F',	'Tigger12A')
oktelist= intersect(nineTE, rownames(sce_filter))
oktelist
#
#TE_ave= AverageExpression(object = sce_filter, features = oktelist)
#head(TE_ave)
expr= as.matrix(sce_filter[['RNA']]@data)
TE_ave= as.data.frame(t(expr[oktelist,]))
TE_ave$z_score= rowMeans(TE_ave)
head(TE_ave)
TE_ave=cbind(TE_ave, sce_filter@meta.data)

######
####cell type information
ctype= read.table('zhu/5.TE/scRNA/E-MTAB-8410/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz', header = T, sep='\t')
head(ctype)
######
length(intersect(TE_ave$Index, ctype$Index))
TE_ave_cellann= merge(TE_ave, ctype, by='Index', all=F )
head(TE_ave_cellann)
dim(TE_ave_cellann)
#
library(ggplot2)
library(ggpubr)
pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_seurat_pipeline.pdf',width = 5,height = 7)
ggpubr::ggboxplot(TE_ave_cellann, x = "Cell_type", y = "z_score",
          color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
          #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
          add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_allcells'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")


ggpubr::ggboxplot(subset(TE_ave_cellann,Class!='Normal' ), x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_only_tumor_tissue'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")

#################

ggpubr::ggboxplot(subset(TE_ave_cellann,Cell_type =='Epithelial cells'), x = "Class", y = "z_score",
                  color = "Class", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_EPI cell only'))+ggplot2::xlab("Class")+ggplot2::ylab("TE score")

dev.off()
#################CPM
####
expr= as.matrix(sce_filter[['RNA']]@counts)
expr <- t(t(expr)/colSums(expr))*1000000
expr= log2(expr+1)
expr[1:4,1:4]

TE_ave= as.data.frame(t(expr[oktelist,]))
TE_ave$z_score= rowMeans(TE_ave)
head(TE_ave)
TE_ave=cbind(TE_ave, sce_filter@meta.data)

######
######
length(intersect(TE_ave$Index, ctype$Index))
TE_ave_cellann= merge(TE_ave, ctype, by='Index', all=F )
TE_ave_cellann$Cell_type= gsub(' ', '_',TE_ave_cellann$Cell_type)
head(TE_ave_cellann)
dim(TE_ave_cellann)
save(TE_ave_cellann, file= 'zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_pipeline.RData')
#


#
pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_pipeline.pdf',width = 5,height = 7)
ggpubr::ggboxplot(TE_ave_cellann, x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  palette=rev(brewer.pal(8,"Set3")),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  #add = "jitter"
                  )+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_allcells'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")


ggpubr::ggboxplot(subset(TE_ave_cellann,Class!='Normal' ), x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  palette=rev(brewer.pal(8,"Set3")),
                  #add = "jitter"
                  )+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_only_tumor_tissue'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")

#################

ggpubr::ggboxplot(subset(TE_ave_cellann,Cell_type =='Epithelial_cells'), x = "Class", y = "z_score",
                  color = "Class", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  palette=rev(brewer.pal(8,"Set3")),
                  #add = "jitter"
                  )+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_EPI cell only'))+ggplot2::xlab("Class")+ggplot2::ylab("TE score")

dev.off()
#
#################
#################filter score > 0
library(RColorBrewer)
library(dplyr)
mevalue = subset(TE_ave_cellann, z_score> 0) %>% group_by(Cell_type) %>%  summarise(Median=median(z_score))
mevalue= mevalue[order(mevalue$Median, decreasing = T),]
head(mevalue)
TE_ave_cellann$Cell_type= factor(TE_ave_cellann$Cell_type, levels = mevalue$Cell_type)
TE_ave_cellann$Class= factor(TE_ave_cellann$Class,  levels = c('Tumor', 'Border','Normal'))

colorset1=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf")

pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_pipeline_filtered.pdf',width = 5,height = 7)

ggpubr::ggboxplot(subset(TE_ave_cellann, z_score >0 ), x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #palette=rev(brewer.pal(8,"Set3")),
                  palette= colorset1, 
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  #add = "jitter"
                  )+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_allcells'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")


ggpubr::ggboxplot(subset(TE_ave_cellann,Class!='Normal' & z_score >0), x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  palette=colorset1,
                  #add = "jitter"
                  )+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_only_tumor_tissue'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")

#################

ggpubr::ggboxplot(subset(TE_ave_cellann,Cell_type =='Epithelial_cells' & z_score > 0), x = "Class", y = "z_score",
                  color = "Class", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  palette=colorset1,
                  #add = "jitter"
                  )+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_EPI cell only'))+ggplot2::xlab("Class")+ggplot2::ylab("TE score")

dev.off()
######dropout rate
######
dropdata= TE_ave_cellann
dropdata$drops= ifelse(dropdata$z_score > 0, 'detected','not.detected')
table(dropdata$drops)
df <- data.frame(
  group = c('detected','non.detected'),
  value = c(19348,5808 )
  )
df$perc= c(df$value[1]/sum(df$value), df$value[2]/sum(df$value))
library(ggplot2)
library(scales)
# Barplot
pdf('piechart_dropout_rate.pdf')
ggplot(df, aes(x="", y=perc, fill=group))+
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+
  scale_fill_brewer("Blues") + scale_fill_brewer(palette="Dark2")+#blank_theme +
    theme(axis.text.x=element_blank())+
    geom_text(aes(y = perc/3 + c(0, cumsum(perc)[-length(perc)]), 
                  label = percent(perc/100)), size=5)
dev.off()



####################
####################seurat pipeline
sce_te
head(ctype)
aa=ctype
rownames(aa)= aa$Index  
head(aa)
iddx= intersect(rownames(aa),  rownames(sce_te@meta.data))
length(iddx)
aa_sel= aa[iddx,]
##
sce_te_sel= subset(sce_te, subset= Index %in% iddx)
sce_te_sel= AddMetaData(sce_te_sel, metadata=aa_sel[,-1])
sce_te_sel
head(sce_te_sel@meta.data)
##
pbmc= sce_te_sel
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#saveRDS(pbmc, file="/home/zhuxq/data/scRNA/scPanCompare/lineage.CRC/GSE144735.KUL3/KUL3data.cohort.seuratnormalized.rds")
###
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
saveRDS(pbmc, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/sce_te_sel.run_seurat.rds')
#
colorset1=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf")
pbmc$Cell_type2= factor(pbmc$Cell_type, levels = c('Epithelial cells', 'Mast cells', 'Stromal cells', 'T cells', 'Myeloids', 'B cells'), )

pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/seurat_umap.pdf',width = 7,height = 6)

DimPlot(pbmc, reduction = "umap",pt.size = 0.1#,split.by = "study",
        ,group.by = "Cell_type2"
        ,cols =colorset1
)

DimPlot(pbmc, reduction = "umap",pt.size = 0.1#,split.by = "study",
        ,group.by = "Class"
        ,cols =colorset1
)

DimPlot(pbmc, reduction = "tsne",pt.size = 0.1#,split.by = "study",
        ,group.by = "Cell_type"
        ,cols =colorset1
)

DimPlot(pbmc, reduction = "tsne",pt.size = 0.1#,split.by = "study",
        ,group.by = "Class"
        ,cols =colorset1
)
dev.off()
#
#
pp= pbmc
Idents(pp)= 'Cell_type'
pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/seurat_violin.pdf',width = 6,height = 6)

VlnPlot(pp, features = oktelist)

dev.off()
##############
##############
#TE fraction proportion
#####
read_total= sce@meta.data[, c('nCount_RNA', 'nFeature_RNA')]
colnames(read_total)= paste('total', colnames(read_total), sep = '.')
head(read_total)

read_gene= sce_gene@meta.data[, c('nCount_RNA', 'nFeature_RNA')]
colnames(read_gene)= paste('gene', colnames(read_gene), sep = '.')
head(read_gene)
#
read_te= sce_te@meta.data[, c('nCount_RNA', 'nFeature_RNA')]
colnames(read_te)= paste('te', colnames(read_te), sep = '.')
head(read_te)
##
readFraction= cbind(read_total, read_gene, read_te)
readFraction$percentage=readFraction$te.nCount_RNA/readFraction$total.nCount_RNA
readFraction$Index= rownames(readFraction)
head(readFraction)
##
head(TE_ave_cellann)
statFraction= merge(readFraction, TE_ave_cellann, my= 'Index', all=F)
dim(statFraction)
#
head(statFraction)
write.csv(statFraction, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_pipeline_filtered_readFraction.csv')



pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_pipeline_filtered_readFraction.pdf',width = 5,height = 7)

ggpubr::ggboxplot(subset(statFraction, z_score >0 ), x = "Cell_type", y = "percentage",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #palette=rev(brewer.pal(8,"Set3")),
                  palette= colorset1, 
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  #add = "jitter"
)+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_allcells'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("percentage")


ggpubr::ggboxplot(subset(statFraction,Class!='Normal' & z_score >0), x = "Cell_type", y = "percentage",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  palette=colorset1,
                  #add = "jitter"
)+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_only_tumor_tissue'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("percentage")

#################

ggpubr::ggboxplot(subset(statFraction,Cell_type =='Epithelial_cells' & z_score > 0), x = "Class", y = "percentage",
                  color = "Class", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  palette=colorset1,
                  #add = "jitter"
)+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_EPI cell only'))+ggplot2::xlab("Class")+ggplot2::ylab("percentage")

dev.off()

####
pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_pipeline_filtered_readFraction_cor.pdf',width = 5,height = 5)

ggpubr::ggscatter(subset(statFraction, z_score > 0), x = "z_score", y = "percentage",
          add = "reg.line", size = 1,                                # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "#253494")
          #add.params = list(color = "black",fill = "lightgray")#,
          #color = "kaps.group", palette =  c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086")
)+ggplot2::xlab("TE score")+ggplot2::ylab("the proportion of reads mapping to TEs")+
  ggplot2::geom_point(fill="#dd1c77",color="#dd1c77")+
  #scale_y_continuous(breaks=seq(0, 0.02, 0.005),limits = c(0, 0.02))+
  ggpubr::stat_cor(method = "spearman", label.x = 0, label.y = 0.6)+ggplot2::theme(legend.position = "right")+ggplot2::ggtitle("scRNA CRC KUL3")

ggpubr::ggscatter(statFraction, x = "z_score", y = "percentage",
                  add = "reg.line", size = 1,                                # Add regression line
                  conf.int = TRUE,                                  # Add confidence interval
                  add.params = list(color = "#253494")
                  #add.params = list(color = "black",fill = "lightgray")#,
                  #color = "kaps.group", palette =  c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086")
)+ggplot2::xlab("TE score")+ggplot2::ylab("the proportion of reads mapping to TEs")+
  ggplot2::geom_point(fill="#dd1c77",color="#dd1c77")+
  #scale_y_continuous(breaks=seq(0, 0.02, 0.005),limits = c(0, 0.02))+
  ggpubr::stat_cor(method = "spearman", label.x = 0, label.y = 0.6)+ggplot2::theme(legend.position = "right")+ggplot2::ggtitle("scRNA CRC KUL3")


dev.off()

#
saveRDS(sce_gene, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/sce_gene.seurat.rds')
saveRDS(sce_te, file='zhu/5.TE/scRNA/E-MTAB-8410/scTE/sce_te.seurat.rds')
#
####
#END
#transfer data



########################get the lib size of combined data, TE, and GENE using CGS cloud server
#######counts ratio not running
####not running
###############################
### using edgeR analysis
########################################
library(edgeR)
y=DGEList(as.matrix(sce_filter[['RNA']]@counts), genes=rownames(sce_filter))
y <- calcNormFactors( y ) #CALCULATE NORMALIZATION FACTOR
####
normalized.log2.cpm=cpm(y,normalized.lib.sizes = TRUE,log = TRUE,prior.count = 5) # following the REdiscoverTE paper
normalized.log2.cpm[1:4,1:4]
dim(normalized.log2.cpm)
###

TE_ave= as.data.frame(t(normalized.log2.cpm[oktelist,]))
TE_ave$z_score= rowMeans(TE_ave)
head(TE_ave)
TE_ave=cbind(TE_ave, sce_filter@meta.data)

######
######
length(intersect(TE_ave$Index, ctype$Index))
TE_ave_cellann= merge(TE_ave, ctype, by='Index', all=F )
head(TE_ave_cellann)
dim(TE_ave_cellann)
#
pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_edgeR_pipeline.pdf',width = 5,height = 7)
ggpubr::ggboxplot(TE_ave_cellann, x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_allcells'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")


ggpubr::ggboxplot(subset(TE_ave_cellann,Class!='Normal' ), x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_only_tumor_tissue'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")

#################

ggpubr::ggboxplot(subset(TE_ave_cellann,Cell_type =='Epithelial cells'), x = "Class", y = "z_score",
                  color = "Class", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_EPI cell only'))+ggplot2::xlab("Class")+ggplot2::ylab("TE score")

dev.off()

#########filtered score >0 
pdf('zhu/5.TE/scRNA/E-MTAB-8410/scTE/output_logCPM_edgeR_pipeline_filtered.pdf',width = 5,height = 7)
ggpubr::ggboxplot(subset(TE_ave_cellann, z_score >9.640414), x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_allcells'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")


ggpubr::ggboxplot(subset(TE_ave_cellann,Class!='Normal' & z_score >9.640414), x = "Cell_type", y = "z_score",
                  color = "Cell_type", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_only_tumor_tissue'))+ggplot2::xlab("Cell_type")+ggplot2::ylab("TE score")

#################

ggpubr::ggboxplot(subset(TE_ave_cellann,Cell_type =='Epithelial cells' & z_score >9.640414), x = "Class", y = "z_score",
                  color = "Class", #palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                  #c("set1" ="#00AFBB" ,"set2"="#d73027","set3"="#E69F00","set4"="#31a354"),
                  add = "jitter")+ggpubr::stat_compare_means(label = "p.signif")+
  ggplot2::theme(legend.position = "right",axis.text.x = ggplot2::element_text(size=12,colour="black",angle=45,hjust=1,vjust=1))+
  #stat_compare_means(comparisons = my_com,label = "p.signif")+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggplot2::ggtitle(paste0('scRNA_KUL3_EPI cell only'))+ggplot2::xlab("Class")+ggplot2::ylab("TE score")

dev.off()













#
#
#

##sce_gene
sce_gene[["percent.mt"]] <- PercentageFeatureSet(sce_gene, pattern = "^MT-")
sce_gene_filter <- subset(sce_gene, subset = nFeature_RNA > 200 & nFeature_RNA < 6000  & nCount_RNA > 1000 & percent.mt < 20)
sce_gene_filter
##



###
###################
####cell type information
ctype= read.table('zhu/5.TE/scRNA/E-MTAB-8410/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz', header = T, sep='\t', row.names=1)
head(ctype)
###
dim(ctype)
table(ctype$Cell_type)
#
head(sce_te@meta.data)
length(intersect(rownames(ctype), rownames(sce_gene_filter@meta.data)))
###
# transform UMI counta into CPM 
expr= as.matrix(sce_te[['RNA']]@counts)
expr <- t(t(expr)/colSums(expr))*1000000
expr[1:4,1:4]
#
nineTE= c('AluSq4',	'HERV1_LTRd',	'LTR21B',	'MER57F',	'MER65C',	'MER92-int',	'SVA_C',	'SVA_F',	'Tigger12A')
#
intersect(nineTE, rownames(expr))
intersect(nineTE, unique(hg38TE$V11))

intersect(nineTE,rownames(sce))

expr_sel= as.data.frame(t(expr[nineTE,]))
head(expr_sel)
