#######self.genesets.score
###########
head(self.genesets.score)
###
self.id.kaps=intersect(rownames(self.genesets.score), rownames(msisubsur))
#####
CB.data.self.kaps=cbind(msisubsur[self.id.kaps,], self.genesets.score[self.id.kaps,])
colnames(CB.data.self.kaps)
dim(CB.data.self.kaps)
######
index.self.kaps=colnames(CB.data.self.kaps)[c(43:231,309:311)]
#my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
#CB.data.self.kaps$kaps.group.agg=factor(CB.data.self.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.self.kaps$kaps.group=factor(CB.data.self.kaps$kaps.group,levels = c("set4","set3","set2","set1"))

for (i in 1:length(index.self.kaps)) {
  pdf(file = paste0("4.model/kaps/immune.self.geneset/stat.20200407/indi.compare/",index.self.kaps[i],".self.genesets.kaps.four.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y =index.self.kaps[i] , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+xlab("kaps.group")+ylab(index.self.kaps[i]) +
          ggtitle(paste0(index.self.kaps[i],"CRC.self"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
########
for (i in 1:length(index.self.kaps)) {
  pdf(file = paste0("4.model/kaps/immune.self.geneset/stat.20200407/total.compare/",index.self.kaps[i],".self.genesets.kaps.four.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y =index.self.kaps[i] , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
          stat_compare_means(label = "p.signif",method = "kruskal.test")+xlab("kaps.group")+ylab(index.self.kaps[i]) +
          ggtitle(paste0(index.self.kaps[i],"CRC.self"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
######################
#####################
#####################
b1=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y ="leukocyte.infiltration" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("leukocyte.infiltration") +
  ggtitle(paste0("leukocyte.infiltration","CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b2=ggviolin(CB.data.pan.kaps[,-(75:78)],x = "kaps.group", y ="Lymphocyte.Infiltration.Signature.Score" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("Lymphocyte.Infiltration") +
  ggtitle(paste0("Lymphocyte.Infiltration"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b3=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y ="IFN.gamma.signature.18.genes.Ayers.etal" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("IFN.gamma response") +
  ggtitle(paste0("IFN.gamma response","CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b4=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y ="hot.tumor.signautre" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("hot.tumor.signautre") +
  ggtitle(paste0("hot.tumor.signautre","CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b5=ggviolin(GEP.stats.kaps,x = "kaps.group", y ="GEP" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
         palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("GEP") +
  ggtitle(paste0("GEP"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

pdf("4.model/kaps/immune.self.geneset/stat.20200407/pan.immue.cell.infilition.panel1.pdf",width = 12.5,height = 5)
print(
  ggarrange(b5,b1,b2,b3,b4,
            labels = c( "A","B","C","D","E"),nrow = 1,ncol = 5,
            common.legend = T, legend = "right")
)
dev.off()
#########################
#########################
b1=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "TcClassII_score", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("TcClassII_score") +
  ggtitle(paste0("TcClassII_score","CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b2=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "TAMsurr_score", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("TAMsurr_score") +
  ggtitle(paste0("TAMsurr_score","CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b3=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "TAMsurr_TcClassII_ratio", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("TAMsurr_TcClassII_ratio") +
  ggtitle(paste0("TAMsurr_TcClassII_ratio","CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b4=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "Tcell_infiltration_1", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("Tcell_infiltration_1") +
  ggtitle(paste0("Tcell_infiltration_1","CRC.pan"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
b5=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "CD4.T.cells.exhuasted.zhang.zeminscRNA", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("CD4.T.cells.exhuasted.zhang.zeminscRNA") +
  ggtitle(paste0("CD4.T.cells.exhuasted.zhang.zeminscRNA"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b6=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "CD8.T.cells.exhuasted.zhang.zeminscRNA", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("CD8.T.cells.exhuasted.zhang.zeminscRNA") +
  ggtitle(paste0("CD8.T.cells.exhuasted.zhang.zeminscRNA"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

pdf("4.model/kaps/immune.self.geneset/stat.20200407/pan.immue.cell.infilition.panel2.pdf",width = 15,height = 5)
print(
  ggarrange(b4,b5,b6,b1,b2,b3,
            labels = c( "A","B","C","D","E","F"),nrow = 1,ncol = 6,
            common.legend = T, legend = "right")
)
dev.off()
####################
####################
b1=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "TGFB.response.wolf", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("TGFB.response.wolf") +
  ggtitle(paste0("TGFB.response.wolf"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
b2=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "up.ECM.signature", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("up.ECM.signature") +
  ggtitle(paste0("up.ECM.signature"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b3=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y = "WESTON_VEGFA_TARGETS_12HR", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("WESTON_VEGFA_TARGETS_12HR") +
  ggtitle(paste0("WESTON_VEGFA_TARGETS_12HR"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b4=ggviolin(CB.data.self.kaps[,-c(238)],x = "kaps.group", y ="EMT_UP", fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("EMT_UP") +
  ggtitle(paste0("EMT_UP"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

b5=ggviolin(CB.data.ipres.kaps,x = "kaps.group", y ="z.mean.IPRES" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
         palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+xlab("TE.cluster")+ylab("z.mean.IPRES") +ggtitle("IPRES")

pdf("4.model/kaps/immune.self.geneset/stat.20200407/pan.immue.evasion.panel3.pdf",width = 12.5,height = 5)
print(
  ggarrange(b1,b2,b3,b4,b5,
            labels = c( "A","B","C","D","E"),nrow = 1,ncol = 5,
            common.legend = T, legend = "right")
)
dev.off()
#####
save(CB.data.self.kaps,CB.data.ipres.kaps,file="4.model/kaps/immune.self.geneset/stat.20200407/pan.immue.evasion.panels.RData")
############################
############################
#################draw heatmap
########

##############################