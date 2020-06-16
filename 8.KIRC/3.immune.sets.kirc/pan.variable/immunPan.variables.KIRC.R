#################################################
############compare the immune vairables in KIRC

##############
####pan cancer variables
###
head(immsubtype)
tmp.kaps.kirc=stat.kirc.kaps.vali
rownames(tmp.kaps.kirc)=substr(rownames(tmp.kaps.kirc),1,12)
pan.id.kaps.kirc=intersect(rownames(tmp.kaps.kirc), rownames(immsubtype))

CB.data.pan.kaps.kirc=cbind(tmp.kaps.kirc[pan.id.kaps.kirc,], immsubtype[pan.id.kaps.kirc,])
#CB.data.pan$TE.cluster.agg=ifelse(CB.data.pan$pt.clu3=="cluster_1","TE.high","TE.low")
colnames(CB.data.pan.kaps.kirc)
dim(CB.data.pan.kaps.kirc)
############
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))

index.pan.kaps.kirc=colnames(CB.data.pan.kaps.kirc)[c(63:90,95:117,120:122)]
index.pan.kaps.kirc=index.pan.kaps.kirc[-47]
head(CB.data.pan.kaps.kirc)
table(CB.data.pan.kaps.kirc$kaps.group.kirc, CB.data.pan.kaps.kirc$Immune.Subtype)
#
for (i in 1:length(index.pan.kaps.kirc)) {
  pdf(file = paste0("8.KIRC/3.immune.sets.kirc/pan.variable/indi/",index.pan.kaps.kirc[i],".pancancer.kaps.four.kirc.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.pan.kaps.kirc[,-c(91:94)],x = "kaps.group.kirc", y =index.pan.kaps.kirc[i] , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
                 palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps.kirc[i]) +
          ggtitle(paste0(index.pan.kaps.kirc[i],"KIRC.pancancer.kaps"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
####
for (i in 1:length(index.pan.kaps.kirc)) {
  pdf(file = paste0("8.KIRC/3.immune.sets.kirc/pan.variable/total/",index.pan.kaps.kirc[i],".pancancer.kaps.four.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.pan.kaps.kirc[,-c(91:94)],x = "kaps.group.kirc", y =index.pan.kaps[i] , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
                 palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
          stat_compare_means(label = "p.signif",method = "kruskal.test")+xlab("TE.cluster")+ylab(index.pan.kaps.kirc[i]) +
          ggtitle(paste0(index.pan.kaps.kirc[i],"kirc.pancancer.kaps"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
##########################################

########
save(CB.data.pan.kaps.kirc, file="8.KIRC/3.immune.sets.kirc/pan.variable/immunPan.variables.KIRC.RData")
