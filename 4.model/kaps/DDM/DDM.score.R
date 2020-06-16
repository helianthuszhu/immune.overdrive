#####DNA damage signature
####
load("~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/DDM.9.genesets.score.RData")
DDMscore=DDM.9.genesets
head(kaps.td)
#
length(intersect(rownames(kaps.td), rownames(DDMscore)))
###
stat.DDM=cbind(kaps.td$kaps.group, DDMscore[rownames(kaps.td),])
colnames(stat.DDM)[1]="kaps.group"
dim(stat.DDM)
#
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
stat.DDM$kaps.group=factor(stat.DDM$kaps.group,levels = c("set4", "set3","set2","set1"))

for (i in 2:ncol(stat.DDM)) {
  pdf(paste0("4.model/kaps/DDM/indi/",colnames(stat.DDM)[i],".kaps.group.pdf"))
  print(ggviolin(stat.DDM,x = "kaps.group", y = colnames(stat.DDM)[i], fill = "kaps.group",alpha = 1,size = 0.3,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+
          theme_bw()+xlab("kaps.group")+ylab(colnames(stat.DDM)[i])+ggtitle(paste0(colnames(stat.DDM)[i]))
  )
  dev.off()
}
save(stat.DDM, file="4.model/kaps/DDM/DDM.score.RData")
