#######TMB
#####
tmb.ratio=read.table("~/nas/Xiaoqiang/opti.data/signatureVSimmune/pan.immune.score/pan.TMB/mutation-load_updated.txt",header = T,sep = '\t')
rownames(tmb.ratio)=tmb.ratio$Tumor_Sample_ID
head(tmb.ratio)
####
totoal.M=read.table("~/nas/Xiaoqiang/opti.data/signatureVSimmune/pan.immune.score/pan.TMB/pancan_indel_mut.txt",header = T,sep = "\t")
rownames(totoal.M)=substr(totoal.M$sample, 1,15)
head(totoal.M)
######
idtmp=intersect(rownames(clus), rownames(tmb.ratio))
####
tmb.ratio.stats=cbind(tmb.ratio[idtmp,], clus[idtmp,])
head(tmb.ratio.stats)
tmb1=ggviolin(tmb.ratio.stats,x = "TE.cluster.agg", y ="Silent.per.Mb" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
         palette = c("#d73027", "#E69F00","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("Silent.per.Mb") 
tmb2=ggviolin(tmb.ratio.stats,x = "TE.cluster.agg", y ="Non.silent.per.Mb" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
              palette = c("#d73027", "#E69F00","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("Non.silent.per.Mb") 
generate.PDF <- function(fig) {
  pdf("4.model/immune/TMB/tmb.ratio.pdf",height  = 4,width = 3)
  print(tmb1)
  print(tmb2)
  dev.off()
}
generate.PDF(fig)
########
idtmp=intersect(rownames(clus), rownames(totoal.M))
####
tmb.total.stats=cbind(totoal.M[idtmp,], clus[idtmp,])
head(tmb.total.stats)
tmb3=ggviolin(tmb.total.stats,x = "TE.cluster.agg", y ="total_mut" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
              palette = c("#d73027", "#E69F00","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("total_mut") 




generate.PDF <- function(fig) {
  pdf("4.model/immune/TMB/tmb.total.pdf",height  = 4,width = 3)
  print(tmb3)
  dev.off()
}
generate.PDF(fig)
###################
cor.tmb.stats=cbind(tmb.total.stats, cor.GEP.te.mean[rownames(tmb.total.stats),])
ggg1=ggscatter(cor.tmb.stats[,-c(5:7)], x = "mean.exp", y = "total_mut", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "mean.TE.score", ylab = "total Mut count")+
  ggtitle("mean.TE.exp vs total Mut count in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggg2=ggscatter(cor.tmb.stats[,-c(5:7)], x = "GEP", y = "total_mut", 
               add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
               cor.coef = T, cor.method = "spearman",
               xlab = "GEP", ylab = "total Mut count")+
  ggtitle("GEP vs total Mut count in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
######
cor.ratio.stats=cbind(tmb.ratio.stats, cor.GEP.te.mean[rownames(tmb.ratio.stats),])
ggg3=ggscatter(cor.ratio.stats[,-c(6:8)], x = "mean.exp", y = "Silent.per.Mb", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "mean.TE.score", ylab = "Silent.per.Mb")+
  ggtitle("mean.TE.exp vs Silent.per.Mb in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#
ggg4=ggscatter(cor.ratio.stats[,-c(6:8)], x = "mean.exp", y = "Non.silent.per.Mb", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "mean.TE.score", ylab = "Non.silent.per.Mb")+
  ggtitle("mean.TE.exp vs Non.silent.per.Mb in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#
generate.PDF <- function(fig) {
  pdf("4.model/immune/TMB/tmb.cor.with.te.MEAN.exp.pdf",height  = 7,width = 7)
  print(ggg1)
  print(ggg2)
  print(ggg3)
  print(ggg4)
  dev.off()
}
generate.PDF(fig)
###################
#############msi
###############
head(cor.tmb.stats)
cor.tmb.stats$TMB=cor.tmb.stats$total_mut/50
summary(cor.tmb.stats$TMB)
ggviolin(cor.tmb.stats[,-c(5:7)],x = "MSI.status", y ="total_mut" , fill = "MSI.status",alpha = 1,size = 0.3,
         palette = c("#d73027", "#E69F00","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("total_mut") 
#
ggviolin(cor.tmb.stats[,-c(5:7)],x = "MSI.status", y ="TMB" , fill = "MSI.status",alpha = 1,size = 0.3,
         palette = c("#d73027", "#E69F00","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster.agg")+ylab("TMB") 
table(cor.tmb.stats$MSI.status, cor.tmb.stats$TE.cluster.agg)
cor.tmb.stats$quart.GEP <- cut(as.numeric(as.character(cor.tmb.stats$GEP)),
                               breaks=quantile(as.numeric(as.character(cor.tmb.stats$GEP)), c(0, 0.5,0.75, 1), na.rm=T),
                               labels=c("GEP-low","GEP-indeterminate","GEP-high"))
cor.tmb.stats$GEP.biGroup=ifelse(cor.tmb.stats$GEP> median(cor.tmb.stats$GEP),"GEP-high","GEP-low")
cor.tmb.stats$quart.TMB=ifelse(cor.tmb.stats$TMB> median(cor.tmb.stats$TMB),"TMB-high","TMB-low")
head(cor.tmb.stats)
table(cor.tmb.stats$GEP.biGroup, cor.tmb.stats$quart.TMB,cor.tmb.stats$TE.cluster.agg)
table(cor.tmb.stats$TE.cluster.agg,)
####
x=as.data.frame.matrix(table(subset(cor.tmb.stats,TE.cluster.agg=="TE.high" )$quart.TMB, subset(cor.tmb.stats,TE.cluster.agg=="TE.high" )$GEP.biGroup))
x$group=rownames(x)
head(x)

y=as.data.frame.matrix(table(subset(cor.tmb.stats,TE.cluster.agg=="TE.low" )$quart.TMB, subset(cor.tmb.stats,TE.cluster.agg=="TE.low" )$GEP.biGroup))
y$group=rownames(y)
head(y)
write.csv(x, paste0("4.model/immune/TMB/","TE.high.group.TMB.GEP",".csv"))
write.csv(y, paste0("4.model/immune/TMB/","TE.low.group.TMB.GEP",".csv"))
#
library(reshape2)
library(plyr)
pval <- chisq.test(subset(cor.tmb.stats,TE.cluster.agg=="TE.high" )$quart.TMB,
                   subset(cor.tmb.stats,TE.cluster.agg=="TE.high" )$GEP.biGroup,correct = T)$p.value
#chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
# long format with column of proportions within each id
xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
p1=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotate("text", x=2, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
  ggtitle("TE.high.group")

pval <- chisq.test(subset(cor.tmb.stats,TE.cluster.agg=="TE.low" )$quart.TMB,
                   subset(cor.tmb.stats,TE.cluster.agg=="TE.low" )$GEP.biGroup,correct = T)$p.value
#chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
# long format with column of proportions within each id
xlong <- ddply(melt(y, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
p2=ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotate("text", x=2, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
  ggtitle("TE.low.group")

generate.PDF <- function(fig) {
  pdf(paste0("4.model/immune/TMB/","stacked.count.plot.TMBGEP-TEgroup",".pdf"),width = 4,height = 6)
  print(p1)
  print(p2)
  dev.off()
}
generate.PDF(fig)
########################
datalist1=list()
datalist2=list()
for (i in 1:2) {
  for (j in 1:2) {
    ac1=as.data.frame(x[i,j])
    rownames(ac1)=paste(rownames(x)[i], colnames(x)[j], sep = "-")
    datalist1[[j]]=ac1
  }
  datalist2[[i]]=do.call(rbind, datalist1)
}
ok.high.group=do.call(rbind, datalist2)
colnames(ok.high.group)="TE.high"
ok.high.group
#
datalist1=list()
datalist2=list()
for (i in 1:2) {
  for (j in 1:2) {
    ac1=as.data.frame(y[i,j])
    rownames(ac1)=paste(rownames(y)[i], colnames(y)[j], sep = "-")
    datalist1[[j]]=ac1
  }
  datalist2[[i]]=do.call(rbind, datalist1)
}
ok.low.group=do.call(rbind, datalist2)
colnames(ok.low.group)="TE.low"
ok.low.group
stackdata=cbind(ok.high.group,ok.low.group)
stackdata
#
library(reshape2)
library(ggplot2)
library(scales)
pval1 <- chisq.test(subset(cor.tmb.stats,TE.cluster.agg=="TE.high" )$quart.TMB,
                   subset(cor.tmb.stats,TE.cluster.agg=="TE.high" )$GEP.biGroup,correct = T)$p.value
pval2 <- chisq.test(subset(cor.tmb.stats,TE.cluster.agg=="TE.low" )$quart.TMB,
                   subset(cor.tmb.stats,TE.cluster.agg=="TE.low" )$GEP.biGroup,correct = T)$p.value
datm <- melt(cbind(stackdata, ind = rownames(stackdata)), id.vars = c('ind'))

datm$ind=factor(datm$ind,levels = rev(unique(datm$ind)))

gc1=ggplot(datm,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  ylab("Percentage(%)")+xlab("TE.cluster")+
  annotate("text", x=1, y=0.9, label=paste0("p-value: ","\n" ,signif(pval1,4)),size=3)+
  annotate("text", x=2, y=0.9, label=paste0("p-value: ","\n" ,signif(pval2,4)),size=3)

gc2=ggplot(datm,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right",legend.title = element_text( size = 8))+
  scale_fill_manual(values= c("#feebe2","#fbb4b9","#f768a1","#ae017e"))+ylab("Percentage(%)")+xlab("TE.cluster")+
  annotate("text", x=1, y=0.9, label=paste0("p-value: ","\n" ,signif(pval1,4)),size=3)+
  annotate("text", x=2, y=0.9, label=paste0("p-value: ","\n" ,signif(pval2,4)),size=3)

gc3=ggplot(datm,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right",legend.title = element_text( size = 8))+
  scale_fill_manual(values= c("#fbb4b9","#f768a1","#c51b8a","#7a0177"))+ylab("Percentage(%)")+xlab("TE.cluster")+
  annotate("text", x=1, y=0.9, label=paste0("p-value: ","\n" ,signif(pval1,4)),size=3)+
  annotate("text", x=2, y=0.9, label=paste0("p-value: ","\n" ,signif(pval2,4)),size=3)

generate.PDF <- function(fig) {
  pdf(paste0("4.model/immune/TMB/","stacked.count.plot.TMBGEP-TEgroup.percentage",".pdf"),width = 4,height = 5)
  print(gc1)
  print(gc2)
  print(gc3)
  dev.off()
}
generate.PDF(fig)
save(datm,cor.tmb.stats, file = "4.model/immune/TMB/stacked.count.plot.TMBGEP-TEgroup.percentage.RData")
#########
cor.tmb.stats$four.gorup=paste(cor.tmb.stats$quart.TMB, cor.tmb.stats$GEP.biGroup,sep = "&")
table(cor.tmb.stats$four.gorup, cor.tmb.stats$MSI.status)
