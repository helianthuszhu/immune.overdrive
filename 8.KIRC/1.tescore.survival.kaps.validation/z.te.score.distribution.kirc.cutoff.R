###z of mean exp distribution
###
stat.kirc.kaps.vali=stat.kirc.kaps.vali[order(stat.kirc.kaps.vali$z.of.mean.exp,decreasing = T),]
stat.kirc.kaps.vali$pos=seq(1:nrow(stat.kirc.kaps.vali))
#
pdf("8.KIRC/1.tescore.survival.kaps.validation/z.te.score.distribution.kirc.cutoff.pdf",width = 7,height = 2)
ggplot(stat.kirc.kaps.vali, aes(x=pos, y=z.of.mean.exp, color=kaps.group.kirc)) +
  scale_color_manual(values=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ))+
  geom_vline(xintercept = 99, color = "#d73027", size=0.2)+
  geom_vline(xintercept = 157, color = "#E69F00", size=0.2)+
  geom_vline(xintercept = 190, color = "#756bb1", size=0.2)+
  geom_point(size=0.5) +theme_classic()
dev.off()

save(stat.kirc.kaps.vali,file="8.KIRC/1.tescore.survival.kaps.validation/z.te.score.distribution.kirc.cutoff.RData")
