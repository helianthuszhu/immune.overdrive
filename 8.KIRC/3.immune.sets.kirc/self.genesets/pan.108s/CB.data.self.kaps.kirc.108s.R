########immune sets
###########
###########
#self genesets
load("8.KIRC/3.immune.sets.kirc/gsva.pan.immne.108.geneset.score.KIRC.RData")
head(stat.kirc.kaps.vali)
head(pan108.sets.score)
length(intersect(rownames(pan108.sets.score), rownames(stat.kirc.kaps.vali)))

####
self108s.id.kaps.kirc=intersect(rownames(pan108.sets.score), rownames(stat.kirc.kaps.vali))
#####
CB.data.self.kaps.kirc.108s=cbind(stat.kirc.kaps.vali[self.id.kaps.kirc,], pan108.sets.score[self.id.kaps.kirc,])
colnames(CB.data.self.kaps.kirc.108s)=gsub(" ",".",colnames(CB.data.self.kaps.kirc.108s))
colnames(CB.data.self.kaps.kirc.108s)=gsub("[-]","_",colnames(CB.data.self.kaps.kirc.108s))
colnames(CB.data.self.kaps.kirc.108s)
dim(CB.data.self.kaps.kirc.108s)
######
index.self.kaps.kirc.108s=colnames(CB.data.self.kaps.kirc.108s)[c(60:164)]
#my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
#CB.data.self.kaps$kaps.group.agg=factor(CB.data.self.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.self.kaps.kirc.108s$kaps.group.kirc=factor(CB.data.self.kaps.kirc.108s$kaps.group.kirc,levels = c("set4","set3","set2","set1"))

for (i in 1:length(index.self.kaps.kirc.108s)) {
  pdf(file = paste0("8.KIRC/3.immune.sets.kirc/self.genesets/pan.108s/indi/",index.self.kaps.kirc.108s[i],".self.genesets.kaps.four.KIRC.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps.kirc.108s,x = "kaps.group.kirc", y =index.self.kaps.kirc.108s[i] , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
                 palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+xlab("kaps.group")+ylab(index.self.kaps.kirc.108s[i]) +
          ggtitle(paste0(index.self.kaps.kirc.108s[i],"CRC.self"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
#
for (i in 1:length(index.self.kaps.kirc.108s)) {
  pdf(file = paste0("8.KIRC/3.immune.sets.kirc/self.genesets/pan.108s/total/",index.self.kaps.kirc.108s[i],".self.genesets.kaps.four.KIRC.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps.kirc.108s,x = "kaps.group.kirc", y =index.self.kaps.kirc.108s[i] , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
                 palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
          stat_compare_means(label = "p.signif",method = "kruskal.test")+xlab("kaps.group")+ylab(index.self.kaps.kirc.108s[i]) +
          ggtitle(paste0(index.self.kaps.kirc.108s[i],"CRC.self"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
###############
###############
save(CB.data.self.kaps.kirc.108s,file="8.KIRC/3.immune.sets.kirc/self.genesets/pan.108s/CB.data.self.kaps.kirc.108s.RData")
