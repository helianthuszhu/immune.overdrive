###########compare the proliferation score and DNA M
###########
head(globalDNA)
head(stat.kirc.kaps.vali)
length(intersect(rownames(globalDNA), rownames(stat.kirc.kaps.vali)))
globalDNA$id=globalDNA$id1
##########
stat.kirc.kaps.vali.globalDNA=merge(stat.kirc.kaps.vali, globalDNA, by="id",all.x=TRUE)
head(stat.kirc.kaps.vali.globalDNA)
dim(stat.kirc.kaps.vali.globalDNA)
#########
pdf("8.KIRC/4.methylation/DNA.prolifer.score.pdf",height  = 5,width = 4)
ggviolin(stat.kirc.kaps.vali.globalDNA,x = "kaps.group.kirc", y ="Global.methylation.level" , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
         palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("Global.methylation.level") +
  ggtitle(paste0("Global.methylation.level",".KIRC"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

ggviolin(stat.kirc.kaps.vali.globalDNA,x = "kaps.group.kirc", y ="Proliferation.score" , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
         palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("Proliferation.score") +
  ggtitle(paste0("Proliferation.score",".KIRC"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

ggviolin(stat.kirc.kaps.vali.globalDNA,x = "kaps.group.kirc", y ="CD8.T.cells" , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
         palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("CD8.T.cells") +
  ggtitle(paste0("CD8.T.cells",".KIRC"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#########
#########
save(stat.kirc.kaps.vali.globalDNA, file ="8.KIRC/4.methylation/DNA.prolifer.score.RData" )
