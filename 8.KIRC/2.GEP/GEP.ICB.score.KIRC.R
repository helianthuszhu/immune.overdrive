######calculate GEP for KIRC
######
kirc.expfireh=read.table("~/nas/Xiaoqiang/opti.data/KIRC.firehouse/gdac.broadinstitute.org_KIRC.mRNAseq_Preprocess.Level_3.2016012800.0.0/KIRC.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt",,header=T,sep="\t",check.names=F)
kirc.expfireh=kirc.expfireh[-c(1:29),]
kirc.expfireh=subset(kirc.expfireh, !(gene=="SLC35E2|728661"))
library(tidyr)
geneid=separate(kirc.expfireh[,1:2],gene,into=c("symbol","ensml"),sep="[|]")
#rownames(geneid)=geneid$symbol
rownames(kirc.expfireh)=geneid$symbol
kirc.expfireh=kirc.expfireh[,-1]
save(kirc.expfireh,file='/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/expMatrix.KIRC.RData')
######
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/expMatrix.KIRC.RData")
kirc.expfireh[1:4,1:4]
######
icb.gene=c("STK11IP","ZBTB34","TBC1D10B","OAZ1","POLR2A","G6PD","ABCF1","C14orf102","UBB","TBP",
           "SDHA","CCL5","CD27", "CD274", "CD276", "CD8A",
           "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1",
           "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
           "PSMB10", "STAT1","TIGIT")
icb.expM.kircfh=as.data.frame(t(kirc.expfireh[icb.gene,]))
icb.expM.kircfh$housexp=rowMeans(icb.expM.kircfh[,1:11],na.rm = T)
#icb.expM.kircfh$housexp=rowSums(icb.expM.kircfh[,1:11])
#icb.expM.kircfh$housexp=icb.expM.kircfh$housexp/11
icb.expM.kircfh[1:4,]
dim(icb.expM.kircfh)
preexp.kircfh=icb.expM.kircfh[,-c(1:11)]
head(preexp.kircfh)
datalist=list()
for (i in 1:18) {
  aa1=preexp.kircfh[,c(i,19)]
  aa2=data.frame(value=aa1[,1]-aa1[,2])
  rownames(aa2)=rownames(preexp.kircfh)
  colnames(aa2)=colnames(preexp.kircfh)[i]
  datalist[[i]]=aa2
}
preexp.cal.kircfh=do.call(cbind, datalist)
preexp.cal.kircfh[is.na(preexp.cal.kircfh)] <- 0
preexp.cal.kircfh[1:4,1:4]
#####
coefs=read.table("4.model/immune/GEP.ICB.predictor.txt",header = T,sep = "\t")
rownames(coefs)=coefs$symbol
coefs=coefs[colnames(preexp.cal),]
coefs$coef=as.numeric(paste(coefs$coef))
#
preexp.cal.kircfh$GEP=preexp.cal.kircfh[,1]*coefs[1,2]+preexp.cal.kircfh[,2]*coefs[2,2]+preexp.cal.kircfh[,3]*coefs[3,2]+
  preexp.cal.kircfh[,4]*coefs[4,2]+preexp.cal.kircfh[,5]*coefs[5,2]+preexp.cal.kircfh[,6]*coefs[6,2]+
  preexp.cal.kircfh[,7]*coefs[7,2]+preexp.cal.kircfh[,8]*coefs[8,2]+preexp.cal.kircfh[,9]*coefs[9,2]+
  preexp.cal.kircfh[,10]*coefs[10,2]+preexp.cal.kircfh[,11]*coefs[11,2]+preexp.cal.kircfh[,12]*coefs[12,2]+
  preexp.cal.kircfh[,13]*coefs[13,2]+preexp.cal.kircfh[,14]*coefs[14,2]+preexp.cal.kircfh[,15]*coefs[15,2]+
  preexp.cal.kircfh[,16]*coefs[16,2]+preexp.cal.kircfh[,17]*coefs[17,2]+preexp.cal.kircfh[,18]*coefs[18,2]
save(preexp.cal.kircfh, file="8.KIRC/2.GEP/GEP.ICB.score.KIRC.RData")
###############