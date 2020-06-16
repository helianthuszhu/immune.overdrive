########GEP vs TMB in KIRC kaps group
head(stat.clin.kirc.TMB.full)
GEP.vs.TMB.kirc=stat.clin.kirc.TMB.full
GEP.vs.TMB.kirc=subset(GEP.vs.TMB.kirc, TMB>0)

GEP.vs.TMB.kirc$GEP.biGroup=ifelse(GEP.vs.TMB.kirc$GEP> median(GEP.vs.TMB.kirc$GEP),"GEP-high","GEP-low")
GEP.vs.TMB.kirc$quart.TMB=ifelse(GEP.vs.TMB.kirc$TMB> median(GEP.vs.TMB.kirc$TMB,na.rm = T),"TMB-high","TMB-low")
head(GEP.vs.TMB.kirc)



table(GEP.vs.TMB.kirc$GEP.biGroup, GEP.vs.TMB.kirc$quart.TMB,GEP.vs.TMB.kirc$kaps.group.kirc)

GEP.vs.TMB.kirc$four.gorup.GEPvsTMB=paste(GEP.vs.TMB.kirc$quart.TMB, GEP.vs.TMB.kirc$GEP.biGroup,sep = "&")

stackdata.kaps.kirc=as.data.frame.matrix(table(GEP.vs.TMB.kirc$four.gorup.GEPvsTMB, GEP.vs.TMB.kirc$kaps.group.kirc))
#####################################
library(reshape2)
library(ggplot2)
library(scales)
x=as.data.frame.matrix(table(GEP.vs.TMB.kirc$four.gorup.GEPvsTMB,GEP.vs.TMB.kirc$kaps.group.kirc))
x$group=rownames(x)
head(x)
write.csv(x, "8.KIRC/2.GEP/GEPvsTMB/kaps.group.TMB.GEP.KIRC.csv")

pval.kaps.kirc <- chisq.test(GEP.vs.TMB.kirc$four.gorup.GEPvsTMB,GEP.vs.TMB.kirc$kaps.group.kirc,correct = T)$p.value


datm.kaps.kirc <- melt(cbind(stackdata.kaps.kirc, ind = rownames(stackdata.kaps.kirc)), id.vars = c('ind'))

datm.kaps.kirc$ind=factor(datm.kaps.kirc$ind,levels = rev(unique(datm.kaps.kirc$ind)))

datm.kaps.kirc$variable=factor(datm.kaps.kirc$variable,levels = c("set4", "set3","set2","set1"))

gc1=ggplot(datm.kaps.kirc,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.kirc,4)),size=3)#+
#annotate("text", x=2, y=0.9, label=paste0("p-value: ","\n" ,signif(pval2,4)),size=3)

gc2=ggplot(datm.kaps.kirc,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right",legend.title = element_text( size = 8))+
  scale_fill_manual(values= c("#feebe2","#fbb4b9","#f768a1","#ae017e"))+ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.kirc,4)),size=3)

gc3=ggplot(datm.kaps.kirc,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right",legend.title = element_text( size = 8))+
  scale_fill_manual(values= c("#fbb4b9","#f768a1","#c51b8a","#7a0177"))+ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.kirc,4)),size=3)

generate.PDF <- function(fig) {
  pdf(paste0("8.KIRC/2.GEP/GEPvsTMB/stacked.count.plot.TMBGEP-TEgroup.percentage.KIRC.pdf"),width = 4,height = 5)
  print(gc1)
  print(gc2)
  print(gc3)
  dev.off()
}
generate.PDF(fig)

save(GEP.vs.TMB.kirc, file="8.KIRC/2.GEP/GEPvsTMB/stacked.count.plot.TMBGEP-TEgroup.percentage.KIRC.RData")
