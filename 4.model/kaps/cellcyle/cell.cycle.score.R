######cell cycle
#######
load("4.model/kaps/cellcyle/cell.cycle.score.RData")
dim(cellcycle.score)
head(cellcycle.score)

length(intersect(rownames(cellcycle.score), rownames(kaps.td)))

stat.cellcycle=cbind(kaps.td$kaps.group,cellcycle.score[rownames(kaps.td),])
colnames(stat.cellcycle)[1]="kaps.group"
head(stat.cellcycle)
stat.cellcycle$kaps.group=factor(stat.cellcycle$kaps.group,levels = c("set4", "set3","set2","set1"))

#
for (i in 2:ncol(stat.cellcycle)) {
  pdf(paste0("4.model/kaps/cellcyle/",colnames(stat.cellcycle)[i],".kaps.group.pdf"))
  print(ggviolin(stat.cellcycle,x = "kaps.group", y = colnames(stat.cellcycle)[i], fill = "kaps.group",alpha = 1,size = 0.3,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group.agg")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
          theme_bw()+xlab("kaps.group")+ylab(colnames(stat.cellcycle)[i])+ggtitle(paste0(colnames(stat.cellcycle)[i]))
  )
  dev.off()
}
