#######compare immune genesets in CRC
#######load immune variables
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/immune/immu.associated.signatures.reuslts.combined.RData")
######load CRC clusters
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/model.cor.>=0.4.totoal.surTEs.RData")
######
######
#ips
head(ips)
head(clu)
tmp.clu=clu
tmp.clu$TE.cluster.agg=ifelse(tmp.clu$pt.clu3=="cluster_1","TE.high","TE.low")
rownames(tmp.clu)=substr(rownames(tmp.clu),1,12)
tmp.id=intersect(rownames(tmp.clu), rownames(ips))
CB.data.ips=cbind(ips[tmp.id,51:54], tmp.clu[tmp.id,])
head(CB.data.ips)
##########
index.ips=colnames(CB.data.ips)[1:4]
for (i in 1:4) {
  drawdata=CB.data.ips
  drawdata=drawdata[!(is.na(drawdata[,index.ips[i]])),]
  drawdata$group=ifelse(drawdata[,index.ips[i]]>=median(drawdata[,index.ips[i]]), "ips-high","ips-low")
  ################
  x=as.data.frame.matrix(table(drawdata$TE.cluster.agg, drawdata$group))
  x$group=rownames(x)
  head(x)
  write.csv(x, paste0("4.model/immune/ips/","stacked.",index.ips[i],".agg",".csv"))
  library(reshape2)
  library(plyr)
  pval <- chisq.test(drawdata$pt.clu3,drawdata$group)$p.value
  #chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
  # long format with column of proportions within each id
  xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
  p1=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
    geom_bar(stat = 'identity')+
    scale_fill_manual(values= c("#41ab5d","#fe9929"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
    ggtitle(paste0(index.ips[i],"\n","-group as >=median-",median(drawdata[,index.ips[i]])))
  generate.PDF <- function(fig) {
    pdf(paste0("4.model/immune/ips/","stacked.",index.ips[i],"agg",".pdf"),width = 4,height = 6)
    print(p1)
    dev.off()
  }
  generate.PDF(fig)
}
##################
#####ipres
#####
head(ipres)
head(clu)
id.ipres=intersect(rownames(clu), rownames(ipres))
###
CB.data.ipres=cbind(clu[id.ipres,], ipres[id.ipres,])
CB.data.ipres$TE.cluster.agg=ifelse(CB.data.ipres$pt.clu3=="cluster_1","TE.high","TE.low")
#####
gg1=ggviolin(CB.data.ipres,x = "TE.cluster.agg", y ="z.mean.IPRES" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
         palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster")+ylab("z.mean.IPRES") +ggtitle("IPRES.CRC.cohort")
generate.PDF <- function(fig) {
  pdf("4.model/immune/ipres/ipres.CRC.agg.pdf",height  = 4,width = 3)
  print(gg1)
  dev.off()
}
generate.PDF(fig)
#############
#####HE.stain
#####
head(HE.stain)
#######
HE.id=intersect(rownames(HE.stain), rownames(tmp.clu))
###
CB.data.HE=cbind(tmp.clu[HE.id,], HE.stain[HE.id,])
head(CB.data.HE)
CB.data.HE$TE.cluster.agg=ifelse(CB.data.HE$pt.clu3=="cluster_1","TE.high","TE.low")
table(CB.data.HE$TE.cluster.agg,CB.data.HE$Immune.Subtype)
####
x=as.data.frame.matrix(table(CB.data.HE$TE.cluster.agg, CB.data.HE$Immune.Subtype))
x$group=rownames(x)
head(x)
write.csv(x, paste0("4.model/immune/HE.stain/","immune.subtype.C1.C6.agg",".csv"))
library(reshape2)
library(plyr)
pval <- chisq.test(CB.data.HE$TE.cluster.agg,CB.data.HE$Immune.Subtype,correct = T)$p.value
#chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
# long format with column of proportions within each id
xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
p1=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotate("text", x=2, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
  ggtitle("Immunesubtype.C1-c6")

generate.PDF <- function(fig) {
  pdf(paste0("4.model/immune/HE.stain/","stacked.immune.C1.6.subtype.agg",".pdf"),width = 4,height = 6)
  print(p1)
  dev.off()
}
generate.PDF(fig)
####
index.HE=colnames(CB.data.HE)[c(7,11:16,19)]
for (i in 1:length(index.HE)) {
  pdf(file = paste0("4.model/immune/HE.stain/",index.HE[i],".HE.IHC.agg.pdf"),height  = 4,width = 3)
  print(ggviolin(CB.data.HE,x = "TE.cluster.agg", y =index.HE[i] , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
          stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab(index.HE[i]) +
          ggtitle(paste0("IHC.",index.HE[i],"CRC"))
        )
  dev.off()
}
##############
####pan cancer variables
###
head(immsubtype)
pan.id=intersect(rownames(tmp.clu), rownames(immsubtype))
CB.data.pan=cbind(tmp.clu[pan.id,], immsubtype[pan.id,])
CB.data.pan$TE.cluster.agg=ifelse(CB.data.pan$pt.clu3=="cluster_1","TE.high","TE.low")
colnames(CB.data.pan)

#
index.pan=colnames(CB.data.pan)[c(9:36,41:63,66:68)]
for (i in 1:length(index.pan)) {
  pdf(file = paste0("4.model/immune/immune.pan/pan.agg/",index.pan[i],".pancancer.agg.pdf"),height  = 4,width = 3)
  print(ggviolin(CB.data.pan,x = "TE.cluster.agg", y =index.pan[i] , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
          stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab(index.pan[i]) +
          ggtitle(paste0(index.pan[i],"CRC.pancancer"))
  )
  dev.off()
}
###############
#######self.genesets.score
###########
head(self.genesets.score)
###
self.id=intersect(rownames(self.genesets.score), rownames(clu))
#####
CB.data.self=cbind(clu[self.id,], self.genesets.score[self.id,])
CB.data.self$TE.cluster.agg=ifelse(CB.data.self$pt.clu3=="cluster_1","TE.high","TE.low")
colnames(CB.data.self)
######
index.self=colnames(CB.data.self)[c(7:194,272:274)]
for (i in 1:length(index.self)) {
  pdf(file = paste0("4.model/immune/self.sets/self.agg/",index.self[i],".self.genesets.agg.pdf"),height  = 4,width = 3)
  print(ggviolin(CB.data.self[,-c(2:5)],x = "TE.cluster.agg", y =index.self[i] , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
          stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab(index.self[i]) +
          ggtitle(paste0(index.self[i],"CRC.self.genesets"))
  )
  dev.off()
}
########
#########


