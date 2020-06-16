#################################################
############compare the immune vairables

##############
####pan cancer variables
###
head(immsubtype)
tmp.kaps=msisubsur
rownames(tmp.kaps)=substr(rownames(tmp.kaps),1,12)
pan.id.kaps=intersect(rownames(tmp.kaps), rownames(immsubtype))

CB.data.pan.kaps=cbind(tmp.kaps[pan.id.kaps,], immsubtype[pan.id.kaps,])
#CB.data.pan$TE.cluster.agg=ifelse(CB.data.pan$pt.clu3=="cluster_1","TE.high","TE.low")
colnames(CB.data.pan.kaps)
dim(CB.data.pan.kaps)
############
#

#my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
#CB.data.pan.kaps$kaps.group.agg=factor(CB.data.pan.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
#CB.data.pan.kaps$kaps.group=factor(CB.data.pan.kaps$kaps.group,levels = c("set4","set3","set2","set1"))
CB.data.pan.kaps=cbind(CB.data.pan.kaps$kaps.group,CB.data.pan.kaps)
colnames(CB.data.pan.kaps)[1]="kaps.group.new"
CB.data.pan.kaps$kaps.group.new=gsub("set4","TE.cluster4",CB.data.pan.kaps$kaps.group.new)
CB.data.pan.kaps$kaps.group.new=gsub("set3","TE.cluster3",CB.data.pan.kaps$kaps.group.new)
CB.data.pan.kaps$kaps.group.new=gsub("set2","TE.cluster2",CB.data.pan.kaps$kaps.group.new)
CB.data.pan.kaps$kaps.group.new=gsub("set1","TE.cluster1",CB.data.pan.kaps$kaps.group.new)
CB.data.pan.kaps$kaps.group.new=factor(CB.data.pan.kaps$kaps.group.new,levels = c("TE.cluster4","TE.cluster3","TE.cluster2","TE.cluster1"))
my_comparisons <- list( c("TE.cluster4", "TE.cluster3"), c("TE.cluster4", "TE.cluster2"), c("TE.cluster4", "TE.cluster1"),
                        c("TE.cluster3", "TE.cluster2"), c("TE.cluster3", "TE.cluster1"), c("TE.cluster2", "TE.cluster1"))
index.pan.kaps=colnames(CB.data.pan.kaps)[c(47:74,79:101,104:106)]
head(CB.data.pan.kaps)
table(CB.data.pan.kaps$kaps.group.new)
for (i in 1:length(index.pan.kaps)) {
  pdf(file = paste0("4.model/kaps/immue.pan.compare/stat.20200407/indi.compare/",index.pan[i],".pancancer.kaps.four.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[i] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps[i]) +
          ggtitle(paste0(index.pan.kaps[i],"CRC.pancancer.kaps"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
####
for (i in 1:length(index.pan.kaps)) {
  pdf(file = paste0("4.model/kaps/immue.pan.compare/stat.20200407/total.compare/",index.pan[i],".pancancer.kaps.four.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[i] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
          stat_compare_means(label = "p.signif",method = "kruskal.test")+xlab("TE.cluster")+ylab(index.pan.kaps[i]) +
          ggtitle(paste0(index.pan.kaps[i],"CRC.pancancer.kaps"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
##########################################
##########################################
#color, shape, size, fill, linetype
b1=ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[19] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps[19]) +
  ggtitle(paste0(index.pan.kaps[19],"CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
b2=ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[20] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps[20]) +
  ggtitle(paste0(index.pan.kaps[20],"CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
b3=ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[21] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps[21]) +
  ggtitle(paste0(index.pan.kaps[21],"CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
b4=ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[22] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps[22]) +
  ggtitle(paste0(index.pan.kaps[22],"CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b5=ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[23] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps[23]) +
  ggtitle(paste0(index.pan.kaps[23],"CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b6=ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group.new", y =index.pan.kaps[24] , fill = "kaps.group.new",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab(index.pan.kaps[24]) +
  ggtitle(paste0(index.pan.kaps[24],"CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

pdf("4.model/kaps/immue.pan.compare/stat.20200407/pan.TCR.TCR.pdf",width = 8,height = 10)
print(
  ggarrange(b2,b3,b1,b4,b5,b6,
            labels = c("A", "B","C","D","E","F"),
            common.legend = T, legend = "bottom")
)

dev.off()

pdf("4.model/kaps/immue.pan.compare/stat.20200407/pan.TCR.TCR.4.pdf",width = 10,height = 5)
print(
  ggarrange(b2,b3,b4,b5,
            labels = c( "A","B","C","D"),nrow = 1,ncol = 4,
            common.legend = T, legend = "right")
)
dev.off()
########
save(CB.data.pan.kaps, file="4.model/kaps/immue.pan.compare/stat.20200407/immunPan.variables.RData")
