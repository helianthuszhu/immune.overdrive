########immune sets
###########
###########
#self genesets
load("8.KIRC/3.immune.sets.kirc/genesets.self.score.kircfh.RData")
head(stat.kirc.kaps.vali)
head(genesets.self.score.kircfh)
length(intersect(rownames(genesets.self.score.kircfh), rownames(stat.kirc.kaps.vali)))

####
self.id.kaps.kirc=intersect(rownames(genesets.self.score.kircfh), rownames(stat.kirc.kaps.vali))
#####
CB.data.self.kaps.kirc=cbind(stat.kirc.kaps.vali[self.id.kaps.kirc,], genesets.self.score.kircfh[self.id.kaps.kirc,])
colnames(CB.data.self.kaps.kirc)
dim(CB.data.self.kaps.kirc)
######
index.self.kaps.kirc=colnames(CB.data.self.kaps.kirc)[c(60:122)]
#my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
#CB.data.self.kaps$kaps.group.agg=factor(CB.data.self.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.self.kaps.kirc$kaps.group.kirc=factor(CB.data.self.kaps.kirc$kaps.group.kirc,levels = c("set4","set3","set2","set1"))

for (i in 1:length(index.self.kaps.kirc)) {
  pdf(file = paste0("8.KIRC/3.immune.sets.kirc/self.genesets/indi/",index.self.kaps.kirc[i],".self.genesets.kaps.four.KIRC.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps.kirc,x = "kaps.group.kirc", y =index.self.kaps.kirc[i] , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
                 palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+xlab("kaps.group")+ylab(index.self.kaps.kirc[i]) +
          ggtitle(paste0(index.self.kaps.kirc[i],"CRC.self"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
#
for (i in 1:length(index.self.kaps.kirc)) {
  pdf(file = paste0("8.KIRC/3.immune.sets.kirc/self.genesets/total/",index.self.kaps.kirc[i],".self.genesets.kaps.four.KIRC.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps.kirc,x = "kaps.group.kirc", y =index.self.kaps.kirc[i] , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
                 palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
          stat_compare_means(label = "p.signif",method = "kruskal.test")+xlab("kaps.group")+ylab(index.self.kaps.kirc[i]) +
          ggtitle(paste0(index.self.kaps.kirc[i],"CRC.self"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
###############
###############
save(CB.data.self.kaps.kirc,file="8.KIRC/3.immune.sets.kirc/self.genesets/CB.data.self.kaps.kirc.RData")
