#################piwi signature
#####
load("~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/ppiwiscore.score.RData")
head(ppiwiscore)
length(intersect(rownames(kaps.td), rownames(ppiwiscore)))
###
stat.piwi=cbind(kaps.td, ppiwiscore[rownames(kaps.td),])
stat.piwi$kaps.group=factor(stat.piwi$kaps.group,levels = c("set4", "set3","set2","set1"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
#######
piwiidx=colnames(ppiwiscore)
for (i in 1:length(piwiidx)) {
  pdf(paste0("4.model/kaps/piwi/",piwiidx[i],".pdf"), width = 5,height = 8)
  print(
    ggviolin(stat.piwi,
             x = "kaps.group", y = piwiidx[i], fill = "kaps.group",alpha = 1,size = 0.3,
             #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
             palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
             add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
      stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
      theme_bw()+xlab("kaps.group")+ylab(paste0(piwiidx[i]))+ggtitle(paste0(piwiidx[i]))
  )
  dev.off()
}
save(stat.piwi, file = "4.model/kaps/piwi/piwi.RData")
