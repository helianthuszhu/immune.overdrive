########class level compare among kaps.group
########
teclass.tumorexp.filterd[1:4,1:4]
dim(teclass.tumorexp.filterd)
######
colnames(kaps.td)
kaps.te.class.compare=cbind(kaps.td[, c(17,38)], teclass.tumorexp.filterd[rownames(kaps.td),])
head(kaps.te.class.compare)
######
library(tidyr)
drawdata.teclass <- kaps.te.class.compare[,c(2,3,5,7,11,15)] %>% pivot_longer(cols=c("DNA","LINE","SINE","LTR","Retroposon"),
                                                                     names_to= "five.class",
                                                                values_to = "expression")
head(drawdata.teclass)
#
library(ggplot2)
library(ggpubr)
#farb=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/class.level.compare.among.kaps.group.pdf",width = 5,height = 3)
ggplot(drawdata.teclass, aes(x=kaps.group, y=expression, group=kaps.group)) + 
  geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.5)+
  stat_compare_means(label = "p.signif")+
  facet_grid(. ~ five.class)+theme(strip.text.x = element_text(size=3))+
  scale_fill_manual(values= farb)+theme_classic()+ylab("TE class level expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/class.level.corrlation.with.TE.score.pdf",width = 4,height = 4)
ggscatter(kaps.te.class.compare, x = "z.of.mean.exp", y = "DNA", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "DNA")+
  ggtitle("z.of.mean.exp vs DNA in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare, x = "z.of.mean.exp", y = "LINE", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "LINE")+
  ggtitle("z.of.mean.exp vs LINE in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare, x = "z.of.mean.exp", y = "SINE", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "SINE")+
  ggtitle("z.of.mean.exp vs SINE in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare, x = "z.of.mean.exp", y = "LTR", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "LTR")+
  ggtitle("z.of.mean.exp vs LTR in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(kaps.te.class.compare, x = "z.of.mean.exp", y = "Retroposon", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "Retroposon")+
  ggtitle("z.of.mean.exp vs Retroposon in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
#####

################
################
################
######
library(tidyr)
drawdata.teclass.pattern <- kaps.te.class.compare[,c(2,3,5,7,11,13,15)] %>% pivot_longer(cols=c("DNA","LINE","SINE","LTR","Satellite","Retroposon"),
                                                                              names_to= "five.class",
                                                                              values_to = "expression")
head(drawdata.teclass.pattern)
##############
aa=aggregate(expression~five.class,drawdata.teclass.pattern,median)
aa=aa[order(aa$expression,decreasing = T),]
#
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/class.level.compare.across.samples.pdf",width = 4,height = 5)
ggplot(drawdata.teclass.pattern, aes(x=five.class, y=expression, fill=five.class)) + 
  geom_boxplot(outlier.size = 0.5)+scale_x_discrete(limits=aa$five.class)+
  scale_fill_manual(values= c("DNA"="#1f78b4","LINE"="#d95f02","LTR"="#7570b3","Retroposon"="#e7298a","Satellite"="#66a61e","SINE"="#e6ab02"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
  )
dev.off()
####
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/class.level.compare.among.kaps.group.new.pdf",width = 5,height = 3)
ggplot(drawdata.teclass.pattern, aes(x=kaps.group, y=expression, group=kaps.group)) + 
  geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.5)+
  stat_compare_means(label = "p.signif")+
  facet_grid(. ~ five.class)+theme(strip.text.x = element_text(size=3))+
  scale_fill_manual(values= farb)+theme_classic()+ylab("TE class level expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#
save(kaps.te.class.compare,drawdata.teclass,drawdata.teclass.pattern,file="4.model/kaps/TEclass.level.compare.kaps.groups/TEclass.level.compare.kaps.groups.RData" )
