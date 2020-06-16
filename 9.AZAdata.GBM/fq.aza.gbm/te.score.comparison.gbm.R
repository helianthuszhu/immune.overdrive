######TE in GSE80137
###
teexp.gse80137=readRDS("9.AZAdata.GBM/fq.aza.gbm/RE.output/RE_intergenic_2_counts_normalized.RDS")
teexp.gse80137
####
head(cndi.rep)
##############
teexp.gse80137.9.tes=as.data.frame(exprs(teexp.gse80137))[rownames(cndi.rep),]
teexp.gse80137.9.tes=as.data.frame(t(teexp.gse80137.9.tes))
teexp.gse80137.9.tes$mean.exp=rowMeans(teexp.gse80137.9.tes,na.rm = T)
teexp.gse80137.9.tes$z.of.mean.exp=(teexp.gse80137.9.tes$mean.exp - mean(teexp.gse80137.9.tes$mean.exp))/sd(teexp.gse80137.9.tes$mean.exp)
######cell information
#head(ascore.gse80137)
teexp.gse80137.9.tes.cellinfo=read.table("9.AZAdata.GBM/fq.aza.gbm/SRR_Acc_List.GBM.cell.line.meta.txt",header = T,sep = ",")
teexp.gse80137.9.tes.cellinfo=teexp.gse80137.9.tes.cellinfo[!(duplicated(teexp.gse80137.9.tes.cellinfo$GEO_Accession..exp.)),]
teexp.gse80137.9.tes.cellinfo=teexp.gse80137.9.tes.cellinfo[,c("Cell_Line","GEO_Accession..exp.", "treatment")]
rownames(teexp.gse80137.9.tes.cellinfo)=teexp.gse80137.9.tes.cellinfo$GEO_Accession..exp.
head(teexp.gse80137.9.tes.cellinfo)
######
teexp.gse80137.9.tes.stat=cbind(teexp.gse80137.9.tes.cellinfo, teexp.gse80137.9.tes[rownames(teexp.gse80137.9.tes.cellinfo),])
teexp.gse80137.9.tes.stat$treatment=ifelse(teexp.gse80137.9.tes.stat$treatment=="Not treated", "control","low.dose")
save(teexp.gse80137.9.tes.stat, file="9.AZAdata.GBM/fq.aza.gbm/te.score.comparison.gbm.RData")
######
######
pdf("9.AZAdata.GBM/fq.aza.gbm/te.score.comparison.gbm.pdf",width = 4,height = 4)
ggboxplot(teexp.gse80137.9.tes.stat, x = "treatment", y = "z.of.mean.exp",
          color = "treatment", palette = c("low.dose"="#E7B800","control"="#FC4E07"),
          add = "jitter")+stat_compare_means()+ggtitle("gbm.5aza.cells")+ylab("normalized TE score") + xlab("treatment")+
  theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=0,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()
