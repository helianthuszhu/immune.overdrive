#######
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/Mean.expression.of.9.TE.model.RData")
#############cor of TE matrix and gene expression matrix in CRC
########
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/5.TE.cor.with.gene/TE.gene.CRC.matrix.RData")
cor.exp[1:4,1:4]
gexp=as.data.frame(t(expM))
gexp[1:4,1:4]
dim(gexp)
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
gexp.sel=gexp.sel[, colnames(gexp.sel) %in% rownames(gene.na.count.sel)]
gexp.sel[1:4,1:4]
######
head(mean.stats)
dim(gexp.sel)
####
length(intersect(rownames(gexp.sel), rownames(mean.stats)))
cordata=cbind(mean.stats[rownames(gexp.sel),]$mean.exp,gexp.sel)
colnames(cordata)[1]="TE.score"
cordata[1:4,1:4]

#####
datalist=list()
for (i in 2:ncol(cordata)) {
    aa=cbind(cordata[,i], cordata[,1])
    colnames(aa)=c(colnames(cordata)[i],colnames(cordata)[1])
    aa=na.omit(aa)
    res <- cor.test(aa[,2], aa[,1],method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    rownames(res)=colnames(cordata)[i]
    rownames(res)=paste(colnames(cordata)[1], rownames(res),sep = "&")
    datalist[[i]] <- res
  }
TE.mean.corM.with.genes=do.call(rbind, datalist)
save(TE.mean.corM.with.genes,"5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/TEmean.cor.with.gene.RData")
