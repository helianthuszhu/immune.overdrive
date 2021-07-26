pbmc=readRDS("/home/zhuxq/data/scRNA/scPanCompare/lineage.CRC/GSE144735.KUL3/KUL3data.cohort.count.rds")
pbmc=readRDS("/home/zhuxq/data/scRNA/scPanCompare/lineage.CRC/GSE132465.SMC/SMCdata.cohort.count.rds")
pbmc=subset(pbmc, subset= Tissue!="Normal")
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
saveRDS(pbmc,file = "/home/zhuxq/data/te/scRNA.cd274/SMCdata_onlyTumor.cohort.seurat.umap.rds")
pbmc=readRDS("/home/zhuxq/data/te/scRNA.cd274/SMCdata_onlyTumor.cohort.seurat.umap.rds")
pbmc_sel=subset(pbmc, subset= celltype!="Mast cells")
#
colorset1=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf")
pdf("/home/zhuxq/data/te/scRNA.cd274/umap.tsne.smcdata.onlyTumor.pdf",width = 8,height = 7)
DimPlot(pbmc_sel, reduction = "umap",pt.size = 0.05#,split.by = "study",
        ,group.by = "celltype"
        ,cols =rev(brewer.pal(8,"Set3"))
)
DimPlot(pbmc_sel, reduction = "tsne",pt.size = 0.05#,split.by = "study",
        ,group.by = "celltype"
        ,cols =rev(brewer.pal(8,"Set3"))
)
dev.off()
#
#FeaturePlot(pbmc, reduction = "tsne",features = c("CD274","CD68","IRF1","STAT1","STAT2","JAK2"),pt.size = 0.01)
####
features.plot=c("CD79A","EPCAM", "CD68","DCN", "CD3D", "JAK2","STAT1","STAT2","IRF1","CD274")

#DoHeatmap(pbmc_sel, features = features.plot, group.by = "celltype", slot="scale.data")

pdf("/home/zhuxq/data/te/scRNA.cd274/ht.SMCdata.marker.pdf",width = 10,height = 10)
#png("/home/zhuxq/data/te/scRNA.cd274/ht.SMCdata.marker.png",width = 10,height = 5)
DoHeatmap(object = pbmc_sel,features =features.plot, group.by = "celltype",
          size = 3, disp.max = 6, disp.min = -1
          )+ 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")


DoHeatmap(object =subset(pbmc_sel, subset= CD274>0),features =features.plot, 
          group.by = "celltype",slot="counts",
          size = 3, disp.max = 6, disp.min = -1,draw.lines=T,lines.width=3
)+ scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")
dev.off()


#StackedVlnPlot(subset(pbmc_sel, subset= CD274>0), features = features.plot,group.by = "celltype",slot="counts")

#VlnPlot(object =subset(pbmc_sel, subset= CD274>0), features  = "CD274",group.by = "celltype" )+
 # geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")

#StackedVlnPlot(obj = subset(pbmc_sel, subset= CD274>0), features = features.plot,group.by = "celltype", cols=rev(brewer.pal(8,"Set3")))
  

pdf("/home/zhuxq/data/te/scRNA.cd274/stacked.violin.pdf",width = 10,height = 10)
StackedVlnPlot(obj = subset(pbmc_sel, subset= CD274>0), features = features.plot,group.by = "celltype",cols=rev(brewer.pal(8,"Set3")))
  
StackedVlnPlot(pbmc_sel, features = features.plot,group.by = "celltype",cols=rev(brewer.pal(8,"Set3")))
  
dev.off()

##########3
#sce_sub=subset(sce, subset = rownames(sce) %in% c("CD274","STAT1","STAT2"))
sce_sub=subset(pbmc, features = as.character(c("CD274","STAT1","STAT2","IRF1","JAK2","CD68","EPCAM","DCN","CD3D","CD79A","KIT")))
#scexp=sce_sub[["RNA"]]@counts
scexp=sce_sub[["RNA"]]@data
scexp=as.data.frame(scexp)
scexp[1:3,1:4]
dim(scexp)
#scexp_sel=scexp[rownames(scexp) %in% c("CD274","STAT1","STAT2"),]
#####
scemeta=pbmc@meta.data
head(scemeta)
#####
ht.hallM.sig.aza=scexp
ht.hallM.sig.aza=t(ht.hallM.sig.aza)
ht.hallM.sig.aza=scale(ht.hallM.sig.aza)
ht.hallM.sig.aza=t(ht.hallM.sig.aza)

#ht.hallM.sig.aza[ht.hallM.sig.aza< -4] <- -4
#ht.hallM.sig.aza[ht.hallM.sig.aza>4] <- 4

min_cor = min(as.vector(ht.hallM.sig.aza))
max_cor = max(as.vector(ht.hallM.sig.aza))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
col_ha_top.hm.aza = columnAnnotation(
  #treatment=stat.hallmarkerscore.aza$treatment,
  celltype=scemeta$celltype#,
  #col=list(#treatment=c("aza"="#E7B800", "untreatment"="#FC4E07"),
   #        cell.line=c("LNT.229"="#00AFBB","T98G"="#e41a1c","U.87"="#377eb8")
  #),
  #show_annotation_name = FALSE#,gp = gpar(col = "black")
  )

####
Heatmap(ht.hallM.sig.aza, name = "marker", use_raster = F,
                 #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                 #col = rev(viridis(10)),
                 #width = unit(2, "cm"),
                 #height = unit(12, "cm"),
                 border = F,
                 col=col.pal_cor,
                 show_column_names = F,show_row_names = T,
                 cluster_columns = F,
                 row_names_gp = gpar(fontsize = 5),
                 column_names_gp = gpar(fontsize = 5),
                 top_annotation = col_ha_top.hm.aza,
                 #right_annotation = row_ha.right,
                 show_row_dend = F,show_column_dend = F,
                 #row_names_side = "left",
                 #left_annotation = row_ha.left.hm.aza,
                 column_title=paste0("hall.marker.genesets","\n","(significant 42 out of 50)"),
                 column_title_gp = gpar(fontsize = 8)
)

######
#####
#####
sel_gmt=read_gmt("/home/zhuxq/data/te/scRNA.cd274/c2.all.v7.2.symbols.sel.gmt")
sets=sel_gmt
sets$KEGG_JAK_STAT_SIGNALING_PATHWAY
sce_marker=subset(pbmc_sel, features = as.character(sets$KEGG_JAK_STAT_SIGNALING_PATHWAY))
#scexp=sce_sub[["RNA"]]@counts
scexp=sce_marker[["RNA"]]@data
scexp=as.data.frame(scexp)
scexp[1:3,1:4]
dim(scexp)
#####
markerscore=as.data.frame(colMeans(scexp))
head(markerscore)
stat.markerscore=data.frame(Jakstatscore=markerscore$`colMeans(scexp)`,as.data.frame(pbmc_sel@meta.data))
head(stat.markerscore)
write.csv(stat.markerscore, "/home/zhuxq/data/te/scRNA.cd274/stacked.violin.Jakstat.score.csv")
#####
pdf("/home/zhuxq/data/te/scRNA.cd274/stacked.violin.Jakstat.score.pdf",width = 10,height = 4)
ggplot(subset(stat.markerscore, Index %in%  colnames(subset(pbmc_sel, subset= CD274>0))), aes(celltype, Jakstatscore)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  #stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(stat.markerscore, aes(celltype, Jakstatscore)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()
dev.off()
########
########calculate sig score
datalist=list()
for (i in 1:length(sets)) {
  sce_marker=subset(pbmc_sel, features = as.character(sets[[i]]))
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
res.markerscore.stat=cbind(res.markerscore,as.data.frame(pbmc_sel@meta.data))
write.csv(res.markerscore.stat, "/home/zhuxq/data/te/scRNA.cd274/stacked.violin.all.signature.score.csv")
res.markerscore.stat=read.csv("/home/zhuxq/data/te/scRNA.cd274/stacked.violin.all.signature.score.csv",header = T)
#####
pdf("/home/zhuxq/data/te/scRNA.cd274/stacked.violin.Jakstat.IFN.score.pdf",width = 10,height = 2)
ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, KEGG_JAK_STAT_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, WP_OVERVIEW_OF_INTERFERONSMEDIATED_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, REACTOME_INTERFERON_GAMMA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, REACTOME_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, REACTOME_INTERFERON_ALPHA_BETA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, ST_TYPE_I_INTERFERON_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, WP_TYPE_III_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  geom_boxplot(aes(fill=celltype),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()
dev.off()
#######
#######
pdf("/home/zhuxq/data/te/scRNA.cd274/stacked.violin.Jakstat.IFN.score.noboxplot.pdf",width = 10,height = 2)
ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, KEGG_JAK_STAT_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, WP_OVERVIEW_OF_INTERFERONSMEDIATED_SIGNALING_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, REACTOME_INTERFERON_GAMMA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, REACTOME_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, REACTOME_INTERFERON_ALPHA_BETA_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, ST_TYPE_I_INTERFERON_PATHWAY)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()


ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()

ggplot(subset(res.markerscore.stat,Index %in%  colnames(subset(pbmc_sel, subset= CD274>0)) ), 
       aes(celltype, WP_TYPE_III_INTERFERON_SIGNALING)) +
  geom_violin(aes(fill = factor(celltype)), scale = "width")+#facet_grid(. ~ celltype)+
  #geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.2,col="black", fill="white")+
  stat_compare_means(label = "p.signif",aes(group = celltype))+
  #geom_jitter(width = 0.1,size=1)+
  scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  theme_classic2()
dev.off()
####
save.image(file = "/home/zhuxq/data/te/scRNA.cd274/cd274.marker.SMCdata.RData")
####
a1sce=subset(pbmc_sel, features = as.character(c("CD274","STAT1","STAT2","IRF1","JAK2","CD68","EPCAM","DCN","CD3D","CD79A","KIT")))
a1exp=a1sce[["RNA"]]@data
a1exp=as.data.frame(t(a1exp))
a1stat=cbind(a1exp, as.data.frame(pbmc_sel@meta.data))
head(a1stat)
a1stat$CD274_group=ifelse(a1stat$CD274>0, "positive","negative")
write.csv(a1stat, "/home/zhuxq/data/te/scRNA.cd274/percentage.cd274.positive.smcdata.csv")
####
a1stat=read.csv("/home/zhuxq/data/te/scRNA.cd274/percentage.cd274.positive.smcdata.csv",header = T)
head(a1stat)
x=as.data.frame.matrix(table(a1stat$CD274_group,a1stat$celltype))
x
# Transform this data in %
data_percentage <- apply(x, 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage = as.data.frame(data_percentage )
data_percentage$group=rownames(data_percentage)
colnames(data_percentage)=gsub(" ","_",colnames(data_percentage))
#
library(tidyr)
drawdata.per1 <- data_percentage %>% pivot_longer(cols=colnames(data_percentage)[1:5],
                    names_to= "celltype", values_to = "percentage")
# Make a stacked barplot--> it will be in %!
pdf("/home/zhuxq/data/te/scRNA.cd274/percentage.cd274.positive.smcdata.pdf",width = 5,height = 7)
ggplot(drawdata.per1, aes(x = celltype, y = percentage, fill = group)) +
  scale_y_continuous(labels = scales::percent)+
  geom_col(position = "fill",colour = "white")+scale_fill_manual(values= c("#80cdc1", "#8c510a"))+
  theme_classic2()+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#####
#####
x=as.data.frame.matrix(table(a1stat$celltype,a1stat$CD274_group))
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
pdf("/home/zhuxq/data/te/scRNA.cd274/percentage.cd274.positive.only.smcdata.coloset.pdf",width = 3,height = 7)
ggplot(drawdata.per1, aes(x = celltype, y = percentage, fill = group)) +
  scale_y_continuous(labels = scales::percent)+
  geom_col(position = "fill",colour = "white")+
  #scale_fill_manual(values= rev(brewer.pal(8,"Set3")))+
  scale_fill_manual(values= c("B cells"="#FCCDE5","Stromal cells"="#80B1D3",
                              "Epithelial cells"="#B3DE69","T cells"="#FB8072","Myeloids"="#FDB462"))+
  theme_classic2()+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
