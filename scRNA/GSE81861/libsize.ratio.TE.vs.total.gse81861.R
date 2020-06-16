#####gse81861
#
gse81861.nmall.cells=read.csv("~/nas/Xiaoqiang/opti.data/scRNA/GSE81861/GSE81861_CRC_NM_all_cells_FPKM.csv",header = T,row.names = 1)
gse81861.nmall.cells[1:4,1:4]
gse81861.nmall.cells.label=data.frame(id=colnames(gse81861.nmall.cells))
gse81861.nmall.cells.label$batch=rep("normal.cells",times=nrow(gse81861.nmall.cells.label))
head(gse81861.nmall.cells.label)
#
gse81861.tuall.cells=read.csv("~/nas/Xiaoqiang/opti.data/scRNA/GSE81861/GSE81861_CRC_tumor_all_cells_FPKM.csv",header = T,row.names = 1)
gse81861.tuall.cells[1:4,1:4]
gse81861.tuall.cells.label=data.frame(id=colnames(gse81861.tuall.cells))
gse81861.tuall.cells.label$batch=rep("tumor.cells",times=nrow(gse81861.tuall.cells.label))
head(gse81861.tuall.cells.label)
####
gse81861.all.cell.label=rbind(gse81861.nmall.cells.label, gse81861.tuall.cells.label)
table(gse81861.all.cell.label$batch)
####
gse81861.all.cell.label$id1=gse81861.all.cell.label$id
gse81861.all.cell.label=separate(gse81861.all.cell.label, col="id1",into = c("RHCid","celltype","color"),sep = "__")
head(gse81861.all.cell.label)
rownames(gse81861.all.cell.label)=gse81861.all.cell.label$RHCid
#######
#######
write.csv(gse81861.all.cell.label, "scRNA/GSE81861/gse81861.all.cell.label.csv")
#######
#######TE counts
#######
gse81861.TEcounts=readRDS("scRNA/GSE81861/RE.output.ver2/RE.output/RE_intergenic_1_raw_counts.RDS")
gse81861.TEcounts.libsize=gse81861.TEcounts$samples[,c(2,4)]
colnames(gse81861.TEcounts.libsize)[1]="lib.size.TE"
head(gse81861.TEcounts.libsize)
#
gse81861.GENEcounts=readRDS("scRNA/GSE81861/RE.output.ver2/RE.output/GENE_1_raw_counts.RDS")
gse81861.GENEcounts.libsize=gse81861.GENEcounts$samples[,c(2,4)]
colnames(gse81861.GENEcounts.libsize)[1]="lib.size.GENE"
head(gse81861.GENEcounts.libsize)
######
gse81861.libsize=merge(gse81861.TEcounts.libsize, gse81861.GENEcounts.libsize, by ="sample",all=FALSE)
rownames(gse81861.libsize)=gse81861.libsize$sample
head(gse81861.libsize)
######
######combine cell information
length(intersect(rownames(gse81861.libsize), rownames(gse81861.all.cell.label)))
stat.libsize.gse81861=cbind(gse81861.libsize, gse81861.all.cell.label[rownames(gse81861.libsize),])
#
stat.libsize.gse81861$lib.size.ratio.total=stat.libsize.gse81861$lib.size.TE+stat.libsize.gse81861$lib.size.GENE
stat.libsize.gse81861$lib.size.ratio.TEvsGENE= stat.libsize.gse81861$lib.size.TE/stat.libsize.gse81861$lib.size.GENE
stat.libsize.gse81861$lib.size.ratio.TEvstotal=stat.libsize.gse81861$lib.size.TE/stat.libsize.gse81861$lib.size.ratio.total
head(stat.libsize.gse81861)
stat.libsize.gse81861.sel=subset(stat.libsize.gse81861, !(celltype=="NA"))
table(stat.libsize.gse81861.sel$batch)
write.csv(stat.libsize.gse81861.sel,"scRNA/GSE81861/libsize.ratio.TE.vs.total.gse81861.csv")
stat.libsize.gse81861.sel.tumor=subset(stat.libsize.gse81861.sel, batch=="tumor.cells")
stat.libsize.gse81861.sel.normal=subset(stat.libsize.gse81861.sel, batch=="normal.cells")
####
####
############
pdf("scRNA/GSE81861/libsize.ratio.TE.vs.total.gse81861.pdf",width = 5)

ggboxplot(stat.libsize.gse81861.sel.normal, x = "celltype", y = "lib.size.ratio.TEvstotal",
          color = "celltype", palette = c("Epithelial"="#d01c8b","Fibroblast"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Macrophage"="#fdc086", 
                                          "Endothelial"="#a6761d","MastCell"="#e41a1c"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggtitle(paste0("library.size.ratio.TE vs GENE","\n","from CRC scRNA GSE81861.normal"))+xlab("cell.type")+ylab("lib.size.ratio.TEvsTotal")

ggboxplot(stat.libsize.gse81861.sel.tumor, x = "celltype", y = "lib.size.ratio.TEvstotal",
          color = "celltype", palette = c("Epithelial"="#d01c8b","Fibroblast"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Macrophage"="#fdc086", 
                                        "Endothelial"="#a6761d","MastCell"="#e41a1c"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggtitle(paste0("library.size.ratio.TE vs GENE","\n","from CRC scRNA GSE81861.tumor"))+xlab("cell.type")+ylab("lib.size.ratio.TEvsTotal")
dev.off()
