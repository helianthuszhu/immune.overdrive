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
head(clus)
head(kaps.td)
idtmp.kaps=intersect(rownames(kaps.td), rownames(tmb.ratio))
####
tmb.ratio.stats.kaps=cbind(tmb.ratio[idtmp.kaps,], kaps.td[idtmp.kaps,])
tmb.ratio.stats.kaps$kaps.group=factor(tmb.ratio.stats.kaps$kaps.group,levels = c("set4", "set3","set2","set1"))
head(tmb.ratio.stats.kaps)
tmb1=ggviolin(tmb.ratio.stats.kaps,x = "kaps.group", y ="Silent.per.Mb" , fill = "kaps.group",alpha = 1,size = 0.3,
              palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab("Silent.per.Mb") +
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

tmb2=ggviolin(tmb.ratio.stats.kaps,x = "kaps.group", y ="Non.silent.per.Mb" , fill = "kaps.group",alpha = 1,size = 0.3,
              palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab("Non.silent.per.Mb")+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

generate.PDF <- function(fig) {
  pdf("4.model/kaps/TMB/tmb.ratio.kaps.group.new.pdf",height  = 4,width = 3)
  print(tmb1)
  print(tmb2)
  dev.off()
}
generate.PDF(fig)
#######
pdf("4.model/kaps/TMB/Cor.TEscore.vs.non.silent.per.Mb.logged.pdf",width = 5,height = 5)
ggscatter(tmb.ratio.stats.kaps, x = "z.of.mean.exp", y = "Non.silent.per.Mb", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.TE.score", ylab = "Non.silent.per.Mb")+
  ggtitle("mean.TE.exp vs Non.silent.per.Mb in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

ggscatter(tmb.ratio.stats.kaps, x = "z.of.mean.exp", y = "Silent.per.Mb", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.TE.score", ylab = "Silent.per.Mb")+
  ggtitle("mean.TE.exp vs Silent.per.Mb in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))
dev.off()

########
idtmp.kaps=intersect(rownames(kaps.td), rownames(totoal.M))
####
tmb.total.stats.kaps=cbind(totoal.M[idtmp.kaps,], kaps.td[idtmp.kaps,])
tmb.total.stats.kaps$total_mut.loged=log2(tmb.total.stats.kaps$total_mut)
tmb.total.stats.kaps$indel.loged=log2(tmb.total.stats.kaps$indel)
head(tmb.total.stats.kaps)
tmb.total.stats.kaps$kaps.group=factor(tmb.total.stats.kaps$kaps.group,levels = c("set4", "set3","set2","set1"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
tmb3=ggviolin(tmb.total.stats.kaps,x = "kaps.group", y ="total_mut" , fill = "kaps.group",alpha = 1,size = 0.3,
              palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons =my_comparisons ,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab("total_mut") 
tmb4=ggviolin(tmb.total.stats.kaps,x = "kaps.group", y ="total_mut.loged" , fill = "kaps.group",alpha = 1,size = 0.3,
              palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons =my_comparisons ,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab("total_mut.loged") 
generate.PDF <- function(fig) {
  pdf("4.model/kaps/TMB/tmb.total.new.pdf",height  = 4,width = 4)
  print(tmb3)
  print(tmb4)
  dev.off()
}
generate.PDF(fig)
####
tmb3=ggviolin(tmb.total.stats.kaps,x = "kaps.group", y ="indel" , fill = "kaps.group",alpha = 1,size = 0.3,
              palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons =my_comparisons ,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab("indel") 

tmb4=ggviolin(tmb.total.stats.kaps,x = "kaps.group", y ="indel.loged" , fill = "kaps.group",alpha = 1,size = 0.3,
              palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons =my_comparisons ,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab("indel.loged") +
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))
  )
generate.PDF <- function(fig) {
  pdf("4.model/kaps/TMB/tmb.indel.pdf",height  = 4,width = 4)
  print(tmb3)
  print(tmb4)
  dev.off()
}
generate.PDF(fig)
#######
pdf("4.model/kaps/TMB/Cor.TEscore.vs.total.mutation.logged.pdf",width = 5,height = 5)
ggscatter(tmb.total.stats.kaps, x = "z.of.mean.exp", y = "total_mut.loged", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.TE.score", ylab = "total_mut.loged")+
  ggtitle("mean.TE.exp vs total_mut.loged in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
###################
cor.tmb.stats.kaps=cbind(cor.GEP.te.mean[rownames(tmb.total.stats.kaps),]$GEP,tmb.total.stats.kaps)
colnames(cor.tmb.stats.kaps)[1]="GEP"
head(cor.tmb.stats.kaps)
ggg1=ggscatter(cor.tmb.stats.kaps, x = "z.of.mean.exp", y = "total_mut", 
               add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
               cor.coef = T, cor.method = "spearman",
               xlab = "mean.TE.score", ylab = "total Mut count")+
  ggtitle("mean.TE.exp vs total Mut count in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggg2=ggscatter(cor.tmb.stats.kaps, x = "GEP", y = "total_mut", 
               add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
               cor.coef = T, cor.method = "spearman",
               xlab = "GEP", ylab = "total Mut count")+
  ggtitle("GEP vs total Mut count in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
######
#cor.ratio.stats=cbind(tmb.ratio.stats.kaps, cor.GEP.te.mean[rownames(tmb.ratio.stats.kaps),])

ggg3=ggscatter(tmb.ratio.stats.kaps, x = "z.of.mean.exp", y = "Silent.per.Mb", 
               add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
               cor.coef = T, cor.method = "spearman",
               xlab = "mean.TE.score", ylab = "Silent.per.Mb")+
  ggtitle("mean.TE.exp vs Silent.per.Mb in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#
ggg4=ggscatter(tmb.ratio.stats.kaps, x = "mean.exp", y = "Non.silent.per.Mb", 
               add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
               cor.coef = T, cor.method = "spearman",
               xlab = "mean.TE.score", ylab = "Non.silent.per.Mb")+
  ggtitle("mean.TE.exp vs Non.silent.per.Mb in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#
generate.PDF <- function(fig) {
  pdf("4.model/kaps/TMB/tmb.cor.with.te.MEAN.exp.kaps.pdf",height  = 7,width = 7)
  print(ggg1)
  print(ggg2)
  print(ggg3)
  print(ggg4)
  dev.off()
}
generate.PDF(fig)
#######

###################
#############msi
###############
head(cor.tmb.stats.kaps)
cor.tmb.stats.kaps$TMB=cor.tmb.stats.kaps$total_mut/50
summary(cor.tmb.stats.kaps$TMB)
############
#cor.tmb.stats$quart.GEP <- cut(as.numeric(as.character(cor.tmb.stats$GEP)),
 #                              breaks=quantile(as.numeric(as.character(cor.tmb.stats$GEP)), c(0, 0.5,0.75, 1), na.rm=T),
  #                             labels=c("GEP-low","GEP-indeterminate","GEP-high"))
cor.tmb.stats.kaps$GEP.biGroup=ifelse(cor.tmb.stats.kaps$GEP> median(cor.tmb.stats.kaps$GEP),"GEP-high","GEP-low")
cor.tmb.stats.kaps$quart.TMB=ifelse(cor.tmb.stats.kaps$TMB> median(cor.tmb.stats.kaps$TMB),"TMB-high","TMB-low")
head(cor.tmb.stats.kaps)
table(cor.tmb.stats.kaps$GEP.biGroup, cor.tmb.stats.kaps$quart.TMB,cor.tmb.stats.kaps$kaps.group)
cor.tmb.stats.kaps$four.gorup=paste(cor.tmb.stats.kaps$quart.TMB, cor.tmb.stats.kaps$GEP.biGroup,sep = "&")

stackdata.kaps=as.data.frame.matrix(table(cor.tmb.stats.kaps$four.gorup, cor.tmb.stats.kaps$kaps.group))
#####################################
library(reshape2)
library(ggplot2)
library(scales)
x=as.data.frame.matrix(table(cor.tmb.stats.kaps$four.gorup,cor.tmb.stats.kaps$kaps.group))
x$group=rownames(x)
head(x)
write.csv(x, "4.model/kaps/TMB/TE.low.group.TMB.GEP.csv")

pval.kaps <- chisq.test(cor.tmb.stats.kaps$four.gorup,cor.tmb.stats.kaps$kaps.group,correct = T)$p.value


datm.kaps <- melt(cbind(stackdata.kaps, ind = rownames(stackdata.kaps)), id.vars = c('ind'))

datm.kaps$ind=factor(datm.kaps$ind,levels = rev(unique(datm.kaps$ind)))

gc1=ggplot(datm.kaps,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps,4)),size=3)#+
  #annotate("text", x=2, y=0.9, label=paste0("p-value: ","\n" ,signif(pval2,4)),size=3)

gc2=ggplot(datm.kaps,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right",legend.title = element_text( size = 8))+
  scale_fill_manual(values= c("#feebe2","#fbb4b9","#f768a1","#ae017e"))+ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps,4)),size=3)

gc3=ggplot(datm.kaps,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right",legend.title = element_text( size = 8))+
  scale_fill_manual(values= c("#fbb4b9","#f768a1","#c51b8a","#7a0177"))+ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps,4)),size=3)

generate.PDF <- function(fig) {
  pdf(paste0("4.model/kaps/TMB/","stacked.count.plot.TMBGEP-TEgroup.percentage",".pdf"),width = 4,height = 5)
  print(gc1)
  print(gc2)
  print(gc3)
  dev.off()
}
generate.PDF(fig)
save(datm,cor.tmb.stats, file = "4.model/kaps/TMB/stacked.count.plot.TMBGEP-TEgroup.percentage.RData")
#########


