################
gene.set=read.table("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/3.immune.sets.kirc/ipres/IPRES.txt",header = T,sep = "\t",na.strings = (""))
gene.set[1:20,]
sets=as.list(gene.set)
sets=lapply(sets, function(x) x[!is.na(x)])
sets[1]
#
#sypnase.data=read.table("../../../TCGA.data/CRC.577/TCGACRC_expression-merged.txt",header = T,row.names = 1)
#sypnase.data[1:5,1:4]
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/expMatrix.KIRC.RData")
sypnase.data=kirc.expfireh
sypnase.data[1:5,1:4]
#
library(GSVA)
Gsva <- gsva(as.matrix(sypnase.data),sets ,method="ssgsea",
             min.sz=2, max.sz=500, verbose=TRUE)
dim(Gsva)
Gsva[1:4,1:4]
Gsva=(Gsva - rowMeans(Gsva))/apply(Gsva,1,sd)
mean.gsva=as.data.frame(colMeans(Gsva))
head(mean.gsva)
colnames(mean.gsva)="z.mean.IPRES"
mean.gsva$IPRES.status=ifelse(mean.gsva$z.mean.IPRES>0.35,"ipres.enriched","ipres.non.enriched")
ipres.score.kirc=cbind(mean.gsva,as.data.frame(t(Gsva)))
save(ipres.score.kirc, file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/3.immune.sets.kirc/ipres/IPRES.score.KIRC.RData")
