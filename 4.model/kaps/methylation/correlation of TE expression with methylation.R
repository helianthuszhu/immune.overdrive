############
############correlation of TE expression with methylation
############
####BULK CRC
####
colnames(stat.meanbeta.cor)
cordata.crc=stat.meanbeta.cor
#
datalist=list()
for (j in 1:6) {
    res <- cor.test(cordata.crc[,j], cordata.crc[,j+6],  method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    rownames(res)=colnames(cordata.crc)[j]
    datalist[[j]] <- res
  }
res.corEvsM.crc=do.call(rbind, datalist)
colnames(res.corEvsM.crc)=paste("bulk",colnames(res.corEvsM.crc),sep = ".")
res.corEvsM.crc$teName=rownames(res.corEvsM.crc)
res.corEvsM.crc$teName=gsub("[-]",".",res.corEvsM.crc$teName)
res.corEvsM.crc
###########
##########wgbs
head(wgbs.stat)
datalist=list()
for (j in 1:9) {
  res <- cor.test(wgbs.stat[,j], wgbs.stat[,j+9],  method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(wgbs.stat)[j]
  datalist[[j]] <- res
}
res.corEvsM.wgbs=do.call(rbind, datalist)
colnames(res.corEvsM.wgbs)=paste("wgbs",colnames(res.corEvsM.wgbs),sep=".")
rownames(res.corEvsM.wgbs)=substr(rownames(res.corEvsM.wgbs),5,nchar(as.character(rownames(res.corEvsM.wgbs))))
res.corEvsM.wgbs$teName=rownames(res.corEvsM.wgbs)
res.corEvsM.wgbs
##########scRNA
head(stat.single.cor.tumor.pt)
datalist=list()
for (j in 1:9) {
  res <- cor.test(stat.single.cor.tumor.pt[,j], stat.single.cor.tumor.pt[,j+9],  method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(stat.single.cor.tumor.pt)[j]
  datalist[[j]] <- res
}
res.corEvsM.scRNA=do.call(rbind, datalist)
colnames(res.corEvsM.scRNA)=paste("scRNA",colnames(res.corEvsM.scRNA),sep=".")
rownames(res.corEvsM.scRNA)=substr(rownames(res.corEvsM.scRNA),5,nchar(as.character(rownames(res.corEvsM.scRNA))))
res.corEvsM.scRNA$teName=rownames(res.corEvsM.scRNA)
res.corEvsM.scRNA$teName=gsub("[-]",".",res.corEvsM.scRNA$teName)
res.corEvsM.scRNA
############
res.corEvsM.CB=merge(res.corEvsM.wgbs,res.corEvsM.scRNA,by="teName",all=TRUE )
res.corEvsM.CB=merge(res.corEvsM.CB, res.corEvsM.crc,by.x="teName",all=TRUE)
rownames(res.corEvsM.CB)=res.corEvsM.CB$teName
############
#
pdf("4.model/kaps/methylation/correlation of TE expression with methylation.pdf",width = 4,height = 4)
ggballoonplot(res.corEvsM.CB[,c(3,5,7)], fill = "value",size.range = c(1, 5))+
  scale_fill_gradientn(colors = c("#980043", "#dd1c77", "#df65b0","#d7b5d8", "#f1eef6"))+
  ggtitle(paste0("   correlation of TE ","\n","expression with methylation"))+theme_bw()+xlab("cohort")+ylab("TE")


ggballoonplot(res.corEvsM.CB[,c(3,5,7)], fill = "value",size.range = c(1, 5))+
  scale_fill_gradientn(colors = c("#276419", "#4d9221", "#7fbc41","#b8e186", "#e6f5d0"))+
  ggtitle(paste0("   correlation of TE ","\n","expression with methylation"))+theme_bw()+xlab("cohort")+ylab("TE")

ggballoonplot(res.corEvsM.CB[,c(3,5,7)], fill = "value",size.range = c(1, 5))+
  scale_fill_gradientn(colors = c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
                                  "#E16462FF", "#FCA636FF", "#F0F921FF"))+
  ggtitle(paste0("   correlation of TE ","\n","expression with methylation"))+theme_bw()+xlab("cohort")+ylab("TE")

dev.off()

save(res.corEvsM.CB,res.corEvsM.crc,res.corEvsM.scRNA,res.corEvsM.wgbs,stat.meanbeta.cor,wgbs.stat,stat.single.cor.tumor.pt,file="4.model/kaps/methylation/correlation of TE expression with methylation.RData")

#
head(TE.mean.corM.with.genes)
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")
balldata=data.frame(corvalue=TE.mean.corM.with.genes[1:5,]$cor)
rownames(balldata)=TE.mean.corM.with.genes$symbol[1:5]

ggballoonplot(balldata, fill = "corvalue",color = "#0073C2FF")










