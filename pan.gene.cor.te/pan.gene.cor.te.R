####pan cancer TE correlation with genes
#######
library(Biobase)
library(dplyr)
library(tidyr)
aa=readRDS("~/nas/Xiaoqiang/R.pacakages/TE/REdiscoverTEdata/inst/Fig4_data/eset_TCGA_TE_intergenic_logCPM.RDS")
rexp=as.data.frame(t(exprs(aa)))
rexp[1:4,1:4]
rclin=pData(aa)
rclin$control=substr(rownames(rclin),14,15)
rclin$control=ifelse(rclin$control=="01","tumor","normal")
head(rclin)
#write.csv(rclin, "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/repCLIN.REdis.csv")
rclin.tumor=subset(rclin, control=="tumor")
head(rclin.tumor)
#
rexp.tumor=rexp[rownames(rclin.tumor),]
dim(rexp.tumor)
#
bb=cbind(rclin.tumor$patient_ID, rexp.tumor[rownames(rclin.tumor),])
colnames(bb)[1]="id"
rexp.agg= bb %>% group_by(id) %>% summarise_all(mean)
rexp.agg=as.data.frame(rexp.agg)
rownames(rexp.agg)=rexp.agg$id
rexp.agg=rexp.agg[,-1]
rexp.agg[1:4,1:4]
dim(rexp.agg)
#
rclin.agg=rclin.tumor[!(duplicated(rclin.tumor$patient_ID)),]
rownames(rclin.agg)=rclin.agg$patient_ID
head(rclin.agg)
dim(rclin.agg)
#
length(intersect(rownames(rexp.agg), rownames(rclin.agg)))
#
#
#
gexp=read.table("~/nas/Xiaoqiang/opti.data/TCGA.data/TCGA.pancancer/mRNA.pan/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.txt",header = T,sep = "\t")
gexp=gexp[!(duplicated(gexp$sample)),]
rownames(gexp)=gexp$sample
gexp=gexp[,-1]
gexp=as.data.frame(t(gexp))
gexp[1:4,1:4]
##########
gexp.2=cbind(substr(rownames(gexp),14,15), gexp)
colnames(gexp.2)[1]="control"
gexp.2[1:4,1:4]
gexp.tumor=subset(gexp.2, control=="01")
dim(gexp.tumor)
gexp.tumor=cbind(substr(rownames(gexp.tumor),1,12), gexp.tumor)
colnames(gexp.tumor)[1]="id"
gexp.tumor=gexp.tumor[,-2]
gexp.tumor=gexp.tumor[!(duplicated(gexp.tumor$id)),]
rownames(gexp.tumor)=gexp.tumor$id
gexp.tumor=gexp.tumor[,-1]
rownames(gexp.tumor)=gsub("[.]","-",rownames(gexp.tumor))
gexp.tumor[1:4,1:4]
#
idpan=intersect(rownames(gexp.tumor),rownames(rexp.agg) )
####
reM=rexp.agg[idpan,]
geM=gexp.tumor[idpan,]
#geM=matrix(as.numeric(unlist(gexp.tumor[idpan,])),nrow=nrow(gexp.tumor[idpan,]))
reM[1:4,1:4]
geM[1:4,1:4]
#
#####get the R and pvalue of spearman correlation
datalist=list()
oklist=list()
cor.exp=reM
cor.score=geM
for (i in 1:ncol(cor.exp)) {
  for (j in 1:ncol(cor.score)) {
    res <- cor.test(cor.exp[,i], cor.score[,j],  method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    rownames(res)=colnames(cor.score)[j]
    rownames(res)=paste(colnames(cor.exp)[i], rownames(res),sep = "and")
    datalist[[j]] <- res
  }
  oklist[[i]]=do.call(rbind, datalist)
}
res.cor=do.call(rbind, oklist)
res.cor.gene.with.te=res.cor
save(res.cor.gene.with.te,"/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/pan.gene.cor.te/cor.gene.with.te.pan.cancer.RData")
