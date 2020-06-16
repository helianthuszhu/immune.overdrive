####################################
#############cor of TE matrix and gene expression matrix in CRC
########
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/5.TE.cor.with.gene/TE.gene.CRC.matrix.RData")
cor.exp[1:4,1:4]
gexp=as.data.frame(t(expM))
gexp[1:4,1:4]
idddd=intersect(rownames(cor.exp), rownames(gexp))
#####
gexp.sel=gexp[idddd,]
teexp.sel=cor.exp[idddd,]
dim(gexp.sel)
dim(teexp.sel)
####NA count of genes
gene.na.count=as.data.frame(colSums(is.na(gexp.sel)))
colnames(gene.na.count)="count"
gene.na.count$symbol=rownames(gene.na.count)
gene.na.count.sel=subset(gene.na.count, count<308)
dim(gene.na.count.sel)
####
gexp.sel=gexp.sel[, colnames(gexp.sel) %in% rownames(gene.na.count.sel)]
gexp.sel[1:4,1:4]
######
#####get the R and pvalue of spearman correlation
datalist=list()
oklist=list()
cor.exp=gexp.sel
cor.score=teexp.sel
#######
library(ppcor)
for (i in 1:ncol(cor.exp)) {
  subdata=cbind(cor.exp[,i],cor.score)
  colnames(subdata)[1]=colnames(cor.exp)[i]
  subdata=na.omit(subdata,cols=colnames(subdata)[1])
  subdata[1:3,1:4]
  for (j in 2:ncol(subdata)) {
    aa=cbind(subdata[,j], subdata[,1])
    colnames(aa)=c(colnames(subdata)[2],colnames(subdata)[1])
    aa=na.omit(aa)
    res <- cor.test(aa[,2], aa[,1],method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    rownames(res)=colnames(subdata)[j]
    rownames(res)=paste(colnames(subdata)[1], rownames(res),sep = "&")
    datalist[[j]] <- res
  }
  oklist[[i]]=do.call(rbind, datalist)
}
res.cor.gene.with.te=do.call(rbind, oklist)
save(res.cor.gene.with.te,file = "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/5.TE.cor.with.gene/TE.gene.correlation.output.RData")
