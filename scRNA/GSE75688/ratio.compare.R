##############
##############compare the ratio counts in single cell BRCA
#GSE75688
##############
scBRCA.TEcounts=readRDS("scRNA/GSE75688/RE.output/RE_intergenic_1_raw_counts.RDS")
scBRCA.TEcounts
scBRCA.TEcounts.libsize=scBRCA.TEcounts$samples[,c(2,4)]
colnames(scBRCA.TEcounts.libsize)[1]="lib.size.TE"
#
scBRCA.GENEcounts=readRDS("scRNA/GSE75688/RE.output/GENE_1_raw_counts.RDS")
scBRCA.GENEcounts
scBRCA.GENEcounts.libsize=scBRCA.GENEcounts$samples[,c(2,4)]
colnames(scBRCA.GENEcounts.libsize)[1]="lib.size.GENE"
#
scBRCA.libsizes=merge(scBRCA.TEcounts.libsize, scBRCA.GENEcounts.libsize, by ="sample",all=FALSE)

scBRCA.libsizes1=scBRCA.libsizes[grepl("clean",scBRCA.libsizes$sample),]
scBRCA.libsizes1$sample=substr(scBRCA.libsizes1$sample, 1, nchar(as.character(scBRCA.libsizes1$sample))-6)

scBRCA.libsizes2=scBRCA.libsizes[!grepl("clean",scBRCA.libsizes$sample),]

scBRCA.libsizes.clean=rbind(scBRCA.libsizes1, scBRCA.libsizes2)

colnames(scBRCA.libsizes.clean)[1]="srrid.gse75688"
dim(scBRCA.libsizes.clean)
head(scBRCA.libsizes.clean)
##############
scbrca.cell1=read.csv("scRNA/GSE75688/single.gsm.515.id.csv",header = T,row.names = 1)
rownames(scbrca.cell1)=scbrca.cell1$da2
head(scbrca.cell1)
scbrca.cell2=read.table("scRNA/GSE75688/SraRunTable.GSE75688.gsm.srr.id.txt",header = T,sep = "\t")
colnames(scbrca.cell2)[2]="da2"
head(scbrca.cell2)
#
scbrca.cellinfo=merge(scbrca.cell2, scbrca.cell1, by="da2",all=FALSE)
head(scbrca.cellinfo)
dim(scbrca.cellinfo)
#############
length(intersect(scBRCA.libsizes.clean$srrid.gse75688, scbrca.cellinfo$srrid.gse75688))
stat.libsize.brca=merge(scBRCA.libsizes.clean, scbrca.cellinfo, by="srrid.gse75688",all=FALSE)
stat.libsize.brca$lib.size.ratio.total=stat.libsize.brca$lib.size.TE+stat.libsize.brca$lib.size.GENE
stat.libsize.brca$lib.size.ratio.TEvsGENE=stat.libsize.brca$lib.size.TE/stat.libsize.brca$lib.size.GENE
stat.libsize.brca$lib.size.ratio.TEvstotal=stat.libsize.brca$lib.size.TE/stat.libsize.brca$lib.size.ratio.total
stat.libsize.brca$index2=factor(stat.libsize.brca$index2,levels = c("Tumor","Immune","Stromal"))
stat.libsize.brca$index3=factor(stat.libsize.brca$index3,levels = c("Tumor","Myeloid","Bcell","Tcell","Stromal"))
head(stat.libsize.brca)
table(stat.libsize.brca$index3)
############
pdf("scRNA/GSE75688/libsize.ratio.TE.vs.GENE.pdf",width = 5)
ggboxplot(stat.libsize.brca, x = "index2", y = "lib.size.ratio.TEvsGENE",
          color = "index2", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Immune"="#386cb0"),
          add = "jitter")+stat_compare_means()+xlab("cell.type")+ylab("lib.size.ratio.TEvsGENE")+
          ggtitle(paste0("library.size.ratio.TE vs GENE","\n","from BRCA scRNA GSE75688"))+theme(legend.position = "right")

ggboxplot(stat.libsize.brca, x = "index3", y = "lib.size.ratio.TEvsGENE",
          color = "index3", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
          ggtitle(paste0("library.size.ratio.TE vs GENE","\n","from BRCA scRNA GSE75688"))+xlab("cell.type")+ylab("lib.size.ratio.TEvsGENE")
dev.off()
############
pdf("scRNA/GSE75688/libsize.ratio.TE.vs.total.pdf",width = 5)
ggboxplot(stat.libsize.brca, x = "index2", y = "lib.size.ratio.TEvstotal",
          color = "index2", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Immune"="#386cb0"),
          add = "jitter")+stat_compare_means()+xlab("cell.type")+ylab("lib.size.ratio.TEvsTotal")+
          ggtitle(paste0("library.size.ratio.TE vs GENE","\n","from BRCA scRNA GSE75688"))+theme(legend.position = "right")

ggboxplot(stat.libsize.brca, x = "index3", y = "lib.size.ratio.TEvstotal",
          color = "index3", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
  ggtitle(paste0("library.size.ratio.TE vs GENE","\n","from BRCA scRNA GSE75688"))+xlab("cell.type")+ylab("lib.size.ratio.TEvsTotal")
dev.off()
###########
######compare across sample tisses
pdf("scRNA/GSE75688/libsize.ratio.TE.vs.total.across.samples.only.tumorcells.pdf",width = 5)
ggboxplot(subset(stat.libsize.brca, index=="Tumor"), x = "da4", y = "lib.size.ratio.TEvstotal",
          color = "da4", #palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1))+
          ggtitle(paste0("library.size.ratio.TE vs GENE","\n","from BRCA scRNA GSE75688"))+xlab("cell.type")+ylab("lib.size.ratio.TEvsTotal")
dev.off()
#####
pdf("scRNA/GSE75688/libsize.ratio.TE.vs.total.across.samples.Turmor.vs.nonTumor.pdf",width = 10)
ggplot(stat.libsize.brca, aes(x=index, y=lib.size.ratio.TEvstotal, group=index)) + 
  geom_boxplot(aes(fill=index),outlier.colour = "black",outlier.size = 0.5)+
  stat_compare_means(label = "p.signif")+
  geom_jitter(width = 0.2,size=1.5)+
  facet_grid(. ~ da4)+theme(strip.text.x = element_text(size=4))+
  scale_fill_manual(values= c("Tumor"="#d01c8b","nonTumor"="#4dac26"))+
  theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1))
dev.off()
###########################
############## counts only included 1204 TEs
##############
scBRCA.TEcounts.matrix=as.data.frame(scBRCA.TEcounts$counts)
scBRCA.GENEcounts.matrix=as.data.frame(scBRCA.GENEcounts$counts)
scBRCA.TEcounts.matrix.sel=scBRCA.TEcounts.matrix[rownames(cdrep),]
dim(scBRCA.TEcounts.matrix.sel)
scBRCA.TEcounts.matrix.sel[1:4,1:4]
#
stat.genecounts.scBRCA=as.data.frame(colSums(scBRCA.GENEcounts.matrix))
colnames(stat.genecounts.scBRCA)="geneCount"
stat.genecounts.scBRCA$id=rownames(stat.genecounts.scBRCA)
head(stat.genecounts.scBRCA)

stat.tecounts.scBRCA=as.data.frame(colSums(scBRCA.TEcounts.matrix.sel))
colnames(stat.tecounts.scBRCA)="teCount"
stat.tecounts.scBRCA$id=rownames(stat.tecounts.scBRCA)
head(stat.tecounts.scBRCA)
#
stat.counts.scBRCA.cb=cbind(stat.genecounts.scBRCA, stat.tecounts.scBRCA)
stat.counts.scBRCA.cb=stat.counts.scBRCA.cb[,-2]
stat.counts.scBRCA.cb$ratio.TEvsGENE=stat.counts.scBRCA.cb$teCount/stat.counts.scBRCA.cb$geneCount
stat.counts.scBRCA.cb$totalCount=stat.counts.scBRCA.cb$geneCount+stat.counts.scBRCA.cb$teCount
stat.counts.scBRCA.cb$ratio.TEvsTotal=stat.counts.scBRCA.cb$teCount/stat.counts.scBRCA.cb$totalCount
head(stat.counts.scBRCA.cb)
#
stat.counts.scBRCA.cb1=stat.counts.scBRCA.cb[grepl("clean",stat.counts.scBRCA.cb$id),]
stat.counts.scBRCA.cb1$id=substr(stat.counts.scBRCA.cb1$id, 1, nchar(as.character(stat.counts.scBRCA.cb1$id))-6)
stat.counts.scBRCA.cb2=stat.counts.scBRCA.cb[!grepl("clean",stat.counts.scBRCA.cb$id),]
head(stat.counts.scBRCA.cb1)
stat.counts.scBRCA.cb.clean=rbind(stat.counts.scBRCA.cb1,stat.counts.scBRCA.cb2)
colnames(stat.counts.scBRCA.cb.clean)[3]="srrid.gse75688"
head(stat.counts.scBRCA.cb.clean)
###############
stat.ratiocounts.brca=merge(stat.counts.scBRCA.cb.clean, scbrca.cellinfo, by="srrid.gse75688",all=FALSE)
head(stat.ratiocounts.brca)
dim(stat.ratiocounts.brca)
##############
stat.ratiocounts.brca$index2=factor(stat.ratiocounts.brca$index2,levels = c("Tumor","Immune","Stromal"))
stat.ratiocounts.brca$index3=factor(stat.ratiocounts.brca$index3,levels = c("Tumor","Myeloid","Bcell","Tcell","Stromal"))
#############
############# draw again
pdf("scRNA/GSE75688/ratiocounts.TE.vs.GENE.pdf",width = 5)
ggboxplot(stat.ratiocounts.brca, x = "index2", y = "ratio.TEvsGENE",
          color = "index2", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Immune"="#386cb0"),
          add = "jitter")+stat_compare_means()+xlab("cell.type")+ylab("ratiocounts.TEvsGENE")
ggboxplot(stat.ratiocounts.brca, x = "index3", y = "ratio.TEvsGENE",
          color = "index3", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
          ggtitle("CountsRatio.TE vs GENE")
dev.off()
############
pdf("scRNA/GSE75688/ratiocounts.TE.vs.total.pdf",width = 5)
ggboxplot(stat.ratiocounts.brca, x = "index2", y = "ratio.TEvsTotal",
          color = "index2", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Immune"="#386cb0"),
          add = "jitter")+stat_compare_means()+xlab("cell.type")+ylab("ratiocounts.TEvsTotal")

ggboxplot(stat.ratiocounts.brca, x = "index3", y = "ratio.TEvsTotal",
          color = "index3", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1,vjust=0.5))+
          ggtitle("CountsRatio.TE vs TE+GENE")
dev.off()
#################
#################compare the 9 te expression
#################
scBRCA.TEexpression=readRDS("scRNA/GSE75688/RE.output/RE_intergenic_2_counts_normalized.RDS")
scBRCA.TEexpression=as.data.frame(exprs(scBRCA.TEexpression))
scBRCA.TEexpression.sel=as.data.frame(t(scBRCA.TEexpression[rownames(cndi.rep),]))
scBRCA.TEexpression.sel$id=rownames(scBRCA.TEexpression.sel)
scBRCA.TEexpression.sel$mean.exp=rowMeans(scBRCA.TEexpression.sel[,c(1:9)])
scBRCA.TEexpression.sel$z.of.mean.exp=(scBRCA.TEexpression.sel$mean.exp-mean(scBRCA.TEexpression.sel$mean.exp))/sd(scBRCA.TEexpression.sel$mean.exp)
head(scBRCA.TEexpression.sel)
scBRCA.TEexpression.sel1=scBRCA.TEexpression.sel[grepl("clean",scBRCA.TEexpression.sel$id),]
scBRCA.TEexpression.sel1$id=substr(scBRCA.TEexpression.sel1$id, 1, nchar(as.character(scBRCA.TEexpression.sel1$id))-6)

scBRCA.TEexpression.sel2=scBRCA.TEexpression.sel[!grepl("clean",scBRCA.TEexpression.sel$id),]
scBRCA.TEexpression.clean=rbind(scBRCA.TEexpression.sel1,scBRCA.TEexpression.sel2)
colnames(scBRCA.TEexpression.clean)[10]="srrid.gse75688"
head(scBRCA.TEexpression.clean)
#########
stat.scbrca.9tes=merge(scBRCA.TEexpression.clean, scbrca.cellinfo, by="srrid.gse75688",all=FALSE)
head(stat.scbrca.9tes)
########################
pdf("scRNA/GSE75688/TEscore.pdf",width = 5)
ggboxplot(stat.scbrca.9tes, x = "index2", y = "z.of.mean.exp",
          color = "index2", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Immune"="#386cb0"),
          add = "jitter")+stat_compare_means()+xlab("cell.type")+ylab("TE score")+
  ggtitle(paste0("nine TE score","\n","from BRCA scRNA GSE75688"))+theme(legend.position = "right")

ggboxplot(stat.scbrca.9tes, x = "index3", y = "z.of.mean.exp",
          color = "index3", palette = c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086"),
          add = "jitter")+stat_compare_means()+theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1))+
  ggtitle(paste0("nine TE score","\n","from BRCA scRNA GSE75688"))+xlab("cell.type")+ylab("TE score")
dev.off()
#
pdf("scRNA/GSE75688/TE.score.across.samples.pdf",width = 10)
ggplot(stat.scbrca.9tes, aes(x=index, y=z.of.mean.exp, group=index)) + 
  geom_boxplot(aes(fill=index),outlier.colour = "black",outlier.size = 0.5)+
  stat_compare_means(label = "p.signif")+
  geom_jitter(width = 0.2,size=1.5)+
  facet_grid(. ~ da4)+theme(strip.text.x = element_text(size=4))+
  scale_fill_manual(values= c("Tumor"="#d01c8b","nonTumor"="#4dac26"))+
  theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1))
dev.off()

save(stat.libsize.brca,stat.ratiocounts.brca,stat.scbrca.9tes,file = "scRNA/GSE75688/ratio.compare.RData" )
#########################
########################
########################tSNE plot based on 1204 TEs
########################
scBRCA.TEexpression.total=readRDS("scRNA/GSE75688/RE.output/RE_intergenic_2_counts_normalized.RDS")
scBRCA.TEexpression.total=as.data.frame(exprs(scBRCA.TEexpression.total))
scBRCA.TEexpression.1204s=as.data.frame(t(scBRCA.TEexpression.total[rownames(repinfo.sel),]))
rownames(scBRCA.TEexpression.1204s)=gsub(".clean","",rownames(scBRCA.TEexpression.1204s))
############
length(intersect(rownames(scBRCA.TEexpression.1204s), scbrca.cellinfo$srrid.gse75688)) #515
scBRCA.TEexpression.1204s.515s=scBRCA.TEexpression.1204s[scbrca.cellinfo$srrid.gse75688,]
dim(scBRCA.TEexpression.1204s.515s)
scBRCA.TEexpression.1204s.515s[1:4,1:4]
dim(scbrca.cellinfo)
############
############tSNE plot
library(Rtsne)
set.seed(12345)
tsne_out2 <- Rtsne(as.matrix(scBRCA.TEexpression.1204s.515s),initial_dims=5,theta=0.0,perplexity = 10)
tmp2 <- data.frame(x=tsne_out2$Y[,1],y=tsne_out2$Y[,2],scbrca.cellinfo)
head(tmp2)
save(tmp2, file="scRNA/GSE75688/Rtsne.GSE75688.BRCA.1204.TEs.output.RData")
#####
tsnetes=ggplot(tmp2, aes(x, y, label=index3)) + 
  #geom_point(aes(colour=factor(K8.detailed),shape=factor(K8.detailed)), size=2) +
  geom_point(aes(colour=factor(index3)), size=0.5) +
  #scale_shape_manual(values=c(21,22,23,25,24))+
  #scale_shape_manual(values=c(10,4,8,3,6))+
  #scale_color_manual(values=wes_palette(n=5, name="GrandBudapest"))+
  #geom_point(aes(shape=factor(site)))+
  #geom_text(aes(colour=factor(control)), check_overlap = TRUE, size=2.2, hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.005) + 
  #scale_colour_discrete(name = "Patient.ID") +
  scale_color_manual(values=c("Tumor"="#e31a1c", "Tcell"="#377eb8", 
                              "Bcell"="#4daf4a", "Myeloid"="#984ea3", 
                              "Stromal"="#f781bf"))+
  labs(x="tSNE1", y="tSNE2", title="Rtsne.GSE75688.BRCA.1204.TEs") + theme_bw()+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(0.3, "cm"))

ggsave(file.path("scRNA/GSE75688/Rtsne.GSE75688.BRCA.1204.TEs.pdf"),tsnetes,width=4,height=3)
#
tsnetes2=ggplot(tmp2, aes(x, y, label=da4)) + 
  #geom_point(aes(colour=factor(K8.detailed),shape=factor(K8.detailed)), size=2) +
  geom_point(aes(colour=factor(da4)), size=0.5)+
  labs(x="tSNE1", y="tSNE2", title="Rtsne.GSE75688.BRCA.1204.TEs") + theme_bw()+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(0.3, "cm"))

ggsave(file.path("scRNA/GSE75688/Rtsne.GSE75688.BRCA.1204.TEs.samples.pdf"),tsnetes2,width=4,height=3)
#######
tsnetes3=ggplot(tmp2, aes(x, y, label=index)) + 
  #geom_point(aes(colour=factor(K8.detailed),shape=factor(K8.detailed)), size=2) +
  geom_point(aes(colour=factor(index)), size=0.5)+
  scale_color_manual(values=c("Tumor"="#e31a1c",  
                              "nonTumor"="#377eb8"))+
  labs(x="tSNE1", y="tSNE2", title="Rtsne.GSE75688.BRCA.1204.TEs") + theme_bw()+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(0.3, "cm"))

ggsave(file.path("scRNA/GSE75688/Rtsne.GSE75688.BRCA.1204.TEs.Tumor.vs.nonTumor.pdf"),tsnetes3,width=4,height=3)
#############################
###########################
##########correlation TEscore and ratioCounts
head(stat.scbrca.9tes)
head(stat.libsize.brca)
length(intersect(stat.scbrca.9tes$srrid.gse75688, stat.libsize.brca$srrid.gse75688))
#####
stat.scbrca.TEscore.cor.ratiocounts=merge(stat.scbrca.9tes, stat.libsize.brca[,-c(4:14)], by="srrid.gse75688",all=TRUE)
head(stat.scbrca.TEscore.cor.ratiocounts)
summary(stat.scbrca.TEscore.cor.ratiocounts$z.of.mean.exp)
pdf("scRNA/GSE75688/scBRCA.correlation.TEscore.vs.ratioCounts.pdf",width = 5,height = 5)
ggscatter(stat.scbrca.TEscore.cor.ratiocounts, x = "z.of.mean.exp", y = "lib.size.ratio.TEvstotal",
          add = "reg.line", size = 1,                                # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black",fill = "lightgray"),
          color = "index3", palette =  c("Tumor"="#d01c8b","Stromal"="#4dac26","Bcell"="#386cb0","Tcell"="#beaed4","Myeloid"="#fdc086")
)+xlab("TE score")+ylab("Count ratio")+
  stat_cor(method = "spearman", label.x = -1, label.y = 0.2)+theme(legend.position = "right")+ggtitle("TEscore.vs.CountRatio.in.scRNA.BRCA")
dev.off()
