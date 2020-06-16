#GEP calculation of pan cancer
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
##################
icb.gene=c("STK11IP","ZBTB34","TBC1D10B","OAZ1","POLR2A","G6PD","ABCF1","C14orf102","UBB","TBP",
           "SDHA","CCL5","CD27", "CD274", "CD276", "CD8A",
           "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1",
           "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
           "PSMB10", "STAT1","TIGIT")

icb.expM=as.data.frame(t(gexp.tumor))[icb.gene,]
icb.expM=as.data.frame(t(icb.expM))
icb.expM$housexp=rowSums(icb.expM[,1:11])
icb.expM$housexp=icb.expM$housexp/11
icb.expM[1:4,]

dim(icb.expM)
preexp=icb.expM[,-c(1:11)]
head(preexp)
datalist=list()
for (i in 1:18) {
  aa1=preexp[,c(i,19)]
  aa2=data.frame(value=aa1[,1]-aa1[,2])
  rownames(aa2)=rownames(preexp)
  colnames(aa2)=colnames(preexp)[i]
  datalist[[i]]=aa2
}
preexp.cal=do.call(cbind, datalist)
#####
coefs=read.table("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/immune/GEP.ICB.predictor.txt",header = T,sep = "\t")
rownames(coefs)=coefs$symbol
coefs=coefs[colnames(preexp.cal),]
coefs$coef=as.numeric(paste(coefs$coef))
#
preexp.cal$GEP=preexp.cal[,1]*coefs[1,2]+preexp.cal[,2]*coefs[2,2]+preexp.cal[,3]*coefs[3,2]+
  preexp.cal[,4]*coefs[4,2]+preexp.cal[,5]*coefs[5,2]+preexp.cal[,6]*coefs[6,2]+
  preexp.cal[,7]*coefs[7,2]+preexp.cal[,8]*coefs[8,2]+preexp.cal[,9]*coefs[9,2]+
  preexp.cal[,10]*coefs[10,2]+preexp.cal[,11]*coefs[11,2]+preexp.cal[,12]*coefs[12,2]+
  preexp.cal[,13]*coefs[13,2]+preexp.cal[,14]*coefs[14,2]+preexp.cal[,15]*coefs[15,2]+
  preexp.cal[,16]*coefs[16,2]+preexp.cal[,17]*coefs[17,2]+preexp.cal[,18]*coefs[18,2]
preexp.cal.pancancer=preexp.cal
save(preexp.cal.pancancer, file="~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/GEP.ICB.score.pancancer.RData")                                                