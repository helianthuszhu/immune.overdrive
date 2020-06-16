##########
head(rds.stats.gse170422)
rds.stats.gse170422$SRRid=rownames(rds.stats.gse170422)
rownames(rds.stats.gse170422)=rds.stats.gse170422$sampleID
####pan.108s
load("vali.GSE107422/immunesets/gsva.pan.immne.108.geneset.score.GSE107422.RData")
head(pan108.sets.score.GSE107422)
######
####
self108s.id.kaps.gse107422=intersect(rownames(pan108.sets.score.GSE107422), rownames(rds.stats.gse170422))
#####
CB.data.self.kaps.108s.gse107422=cbind(rds.stats.gse170422[self108s.id.kaps.gse107422,], pan108.sets.score.GSE107422[self108s.id.kaps.gse107422,])

colnames(CB.data.self.kaps.108s.gse107422)=gsub(" ",".",colnames(CB.data.self.kaps.108s.gse107422))
colnames(CB.data.self.kaps.108s.gse107422)=gsub("[-]","_",colnames(CB.data.self.kaps.108s.gse107422))
colnames(CB.data.self.kaps.108s.gse107422)
dim(CB.data.self.kaps.108s.gse107422)
######
index.self.kaps.108s.gse107422=colnames(CB.data.self.kaps.108s.gse107422)[c(19:124)]
#my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
#CB.data.self.kaps$kaps.group.agg=factor(CB.data.self.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.self.kaps.108s.gse107422$group=factor(CB.data.self.kaps.108s.gse107422$group,levels = c("set4","set3","set2","set1"))

#
for (i in 1:length(index.self.kaps.108s.gse107422)) {
  pdf(file = paste0("vali.GSE107422/immunesets/pan.108s/total/",index.self.kaps.108s.gse107422[i],".self.genesets.kaps.four.GSE107422.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps.108s.gse107422,x = "group", y =index.self.kaps.108s.gse107422[i] , fill = "group",alpha = 1,size = 0.01,width = 1,
                 palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "group")+
          stat_compare_means(label = "p.signif",method = "kruskal.test")+xlab("kaps.group")+ylab(index.self.kaps.108s.gse107422[i]) +
          ggtitle(paste0(index.self.kaps.108s.gse107422[i],"CRC.self"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
###############
###############
save(CB.data.self.kaps.108s.gse107422,file="vali.GSE107422/immunesets/pan.108s/CB.data.self.kaps.108s.GSE107422.RData")
#########