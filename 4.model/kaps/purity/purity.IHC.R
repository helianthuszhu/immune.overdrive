purityid=intersect(rownames(purdata), rownames(kaps.td))
stat.purity=cbind(kaps.td[purityid,], purdata[purityid,])

my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
stat.purity$kaps.group=factor(stat.purity$kaps.group,levels = c("set4", "set3","set2","set1"))


pit1=ggviolin(stat.purity,x = "kaps.group", y = "IHC", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means( comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("tumor purity (IHC)")+ggtitle("tumor purity (IHC)")
pit2=ggviolin(stat.purity,x = "kaps.group", y = "IHC", fill = "kaps.group",alpha = 1,size = 0.3,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means( label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("tumor purity (IHC)")+ggtitle("tumor purity (IHC)")

generate.PDF <- function(fig) {
  pdf("4.model/kaps/purity/purity.IHC.pdf",width = 5,height = 8)
  print(pit1)
  print(pit2)
  dev.off()
}
generate.PDF(fig)
save(stat.purity,file="4.model/kaps/purity/purity.IHC.RData")
