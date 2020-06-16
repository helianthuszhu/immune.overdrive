############
############
########run GSEA based on correlation z.of.mean.exp variable
head(TE.mean.corM.with.genes)
#
gene.corTEscore <- TE.mean.corM.with.genes$symbol
## 转换
library(clusterProfiler)
gene.corTEscore = bitr(gene.corTEscore, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene.corTEscore,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=TE.mean.corM.with.genes$cor,
                      SYMBOL = TE.mean.corM.with.genes$symbol)
gene_df <- merge(gene_df,gene,by="SYMBOL")
## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
########
####C2
## 读入hallmarks gene set，从哪来？
c2.all.gmt <- read.gmt("~/nas/Xiaoqiang/opti.data/TCGA.data/liver/stat.20200418/c2.all.v7.1.entrez.gmt")
# 需要网络
y.GSEA.TE.score <- GSEA(geneList,TERM2GENE =c2.all.gmt)
yd.GSEA.TE.score <- data.frame(y.GSEA.TE.score)
write.csv(yd.GSEA.TE.score,"5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/cor.gsea/cor.GSEA.C2all.TE.score.based.csv")
####
c5.all.gmt <- read.gmt("5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/cor.gsea/c5.all.v7.1.entrez.gmt")
# 需要网络
y.C5.GSEA.TE.score <- GSEA(geneList,TERM2GENE =c5.all.gmt)
yd.C5.GSEA.TE.score <- data.frame(y.C5.GSEA.TE.score)
write.csv(yd.C5.GSEA.TE.score,"5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/cor.gsea/cor.GSEA.C5.GO.TE.score.based.csv")
####
save(TE.mean.corM.with.genes, y.GSEA.TE.score, y.C5.GSEA.TE.score, file = "5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/cor.gsea/cor.GSEA.TE.score.based.RData")
####