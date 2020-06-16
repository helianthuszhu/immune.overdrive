########class level compare among kaps.group
########
TE.classexp.kirc=read.csv("~/nas/Xiaoqiang/R.pacakages/TE/REdiscoverTEdata/expM.class.csv",header = T)
rownames(TE.classexp.kirc)=TE.classexp.kirc$id
TE.classexp.kirc=TE.classexp.kirc[,-c(1:2)]
TE.classexp.kirc=as.data.frame(t(TE.classexp.kirc))
TE.classexp.kirc$id=substr(rownames(TE.classexp.kirc),1,15)
#TE.classexp.kirc$id2=substr(rownames(TE.classexp.kirc),16,16)
TE.classexp.kirc.agg= TE.classexp.kirc %>% group_by(id) %>% summarise_all(mean)
TE.classexp.kirc.agg=as.data.frame(TE.classexp.kirc.agg)
rownames(TE.classexp.kirc.agg)=TE.classexp.kirc.agg$id
rownames(TE.classexp.kirc.agg)=gsub("[.]","-",rownames(TE.classexp.kirc.agg))
head(TE.classexp.kirc.agg)
save(TE.classexp.kirc.agg,file="8.KIRC/TE.class.compate.kaps.group/TEclass.level.expM.pancancer.RData")
#
head(clin.heat.kirc.var)
colnames(clin.heat.kirc.var)
length(intersect(rownames(clin.heat.kirc.var), rownames(TE.classexp.kirc.agg)))
########
kaps.te.class.compare.kirc=cbind(clin.heat.kirc.var[,c(46,49)], TE.classexp.kirc.agg[rownames(clin.heat.kirc.var),])
######
library(tidyr)
colnames(kaps.te.class.compare.kirc)
drawdata.teclass.kirc <- kaps.te.class.compare.kirc[,c(2,4,6,8,12,16)] %>% pivot_longer(cols=c("DNA","LINE","SINE","LTR","Retroposon"),
                                                                              names_to= "five.class",
                                                                              values_to = "expression")
head(drawdata.teclass.kirc)
#
library(ggplot2)
library(ggpubr)
#farb=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
pdf("8.KIRC/TE.class.compate.kaps.group/class.level.compare.among.kaps.group.KIRC.pdf",width = 5,height = 3)
ggplot(drawdata.teclass.kirc, aes(x=kaps.group.kirc, y=expression, group=kaps.group.kirc)) + 
  geom_boxplot(aes(fill=kaps.group.kirc),outlier.colour = "black",outlier.size = 0.5)+
  stat_compare_means(label = "p.signif")+
  facet_grid(. ~ five.class)+theme(strip.text.x = element_text(size=3))+
  scale_fill_manual(values= farb)+theme_classic()+ylab("TE class level expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#
pdf("8.KIRC/TE.class.compate.kaps.group/class.level.corrlation.with.TE.score.pdf",width = 4,height = 4)
ggscatter(kaps.te.class.compare.kirc, x = "z.of.mean.exp", y = "DNA", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "DNA")+
  ggtitle("z.of.mean.exp vs DNA in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare.kirc, x = "z.of.mean.exp", y = "LINE", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "LINE")+
  ggtitle("z.of.mean.exp vs LINE in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare.kirc, x = "z.of.mean.exp", y = "SINE", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "SINE")+
  ggtitle("z.of.mean.exp vs SINE in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare.kirc, x = "z.of.mean.exp", y = "LTR", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "LTR")+
  ggtitle("z.of.mean.exp vs LTR in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare.kirc, x = "z.of.mean.exp", y = "Retroposon", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "Retroposon")+
  ggtitle("z.of.mean.exp vs Retroposon in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
#####
save(kaps.te.class.compare.kirc,drawdata.teclass.kirc,file="8.KIRC/TE.class.compate.kaps.group/TEclass.level.compare.kaps.groups.KIRC.RData" )
