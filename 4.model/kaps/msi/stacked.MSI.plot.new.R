library(reshape2)
library(plyr)
library(ggplot2)
x=read.csv("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/msi/MSIvs.kaps.group.csv",header = T,row.names = 1)
#save(kaps.td, file="4.model/kaps/msi/kaps.td.RData")
load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/msi/kaps.td.RData")
pval <- chisq.test(kaps.td$kaps.group,kaps.td$MSI.status.bin ,correct = T)$p.value
#chisq.test(as.data.frame.matrix(table(drawdata$cluster4, drawdata$group)))
# long format with column of proportions within each id
xlong <- ddply(melt(x, id.vars = 'group'), .(group), mutate, prop = value / sum(value))
xlong$group=factor(xlong$group, levels = c("set4","set3","set2","set1"))
pdf("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/msi/stacked.MSI.plot.new.pdf",width = 4,height = 6)
ggplot(xlong, aes(x = group, y = prop, fill = variable)) + 
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+
  ggtitle("MSI distribution among four group")
dev.off()
