###########GEP ICB
load("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/expMatrixCRC.RData")
icb.gene=c("STK11IP","ZBTB34","TBC1D10B","OAZ1","POLR2A","G6PD","ABCF1","C14orf102","UBB","TBP",
           "SDHA","CCL5","CD27", "CD274", "CD276", "CD8A",
           "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1",
           "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
           "PSMB10", "STAT1","TIGIT")
icb.expM=as.data.frame(t(expM[icb.gene,]))
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
preexp.cal[is.na(preexp.cal)] <- 0
#####
coefs=read.table("4.model/immune/GEP.ICB.predictor.txt",header = T,sep = "\t")
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
save(preexp.cal, file="~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/GEP.ICB.score.CRC.RData")
#
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/model.cor.>=0.4.totoal.surTEs.RData")
dim(clu)
GEP.stats=cbind(preexp.cal[rownames(clu),]$GEP,clu)
colnames(GEP.stats)[1]="GEP"
GEP.stats$TE.cluster.agg=ifelse(GEP.stats$pt.clu3=="cluster_1","TE.high","TE.low")
head(GEP.stats)
save(GEP.stats,file="4.model/immune/GEP/GEP.stats.RData")
#####draw plot
farb=c("#d73027","#E69F00","#00AFBB","#ca0020","#0571b0","#4daf4a","#27408B","#FF0000","#2E8B57","#CD00CD","#FC4E07")
gep1=ggplot(GEP.stats, aes(x=TE.cluster.agg, y=GEP, group=TE.cluster.agg)) + 
  geom_boxplot(aes(fill=TE.cluster.agg),alpha = 0.8,outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.4)+
  stat_compare_means(label = "p.signif")+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+xlab("TE.cluster.agg")+ylab("GEP")+labs(fill = "TE.cluster.agg")
require(ggplot2)
require(reshape2)
library(ggpubr)
gep2=ggviolin(GEP.stats,x = "TE.cluster.agg", y ="GEP" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
              palette = c("#d73027", "#E69F00","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("GEP") 

generate.PDF <- function(fig) {
  pdf("4.model/immune/GEP/GEP.three.clusters.agg.pdf",height  = 4,width = 3)
  print(gep1)
  print(gep2)
  dev.off()
}
generate.PDF(fig)
###################