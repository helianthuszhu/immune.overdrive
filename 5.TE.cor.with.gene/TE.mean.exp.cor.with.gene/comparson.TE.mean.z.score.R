#######compare the TE mean expression, z.score differences
############
head(mean.stats)
summary(mean.stats$mean.exp)
library(ggridges)
ggplot(mean.stats, aes(x = mean.exp, y = TE.cluster.agg)) +
  geom_density_ridges(aes(fill = TE.cluster.agg)) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
ggdensity(mean.stats, x = "mean.exp",
          add = "mean", rug = TRUE,
          color = "TE.cluster.agg", fill = "TE.cluster.agg",
          palette = c("#0073C2FF", "#FC4E07"))
z.mean.statsM=as.data.frame(scale(mean.stats[,c(3:12)]))
z.mean.statsM=cbind(mean.stats$TE.cluster.agg, z.mean.statsM)
colnames(z.mean.statsM)[1]="TE.cluster.agg"
colnames(z.mean.statsM)[11]="z.of.mean.exp"

z.mean.statsM$z.of.mean.exp.from.indi.TE=rowSums(z.mean.statsM[,2:10])/9
z.mean.statsM$mean.exp=mean.stats$mean.exp
z.mean.statsM$TE.cluster=mean.stats$TE.cluster
head(z.mean.statsM)
my_comparisons <- list( c("cluster_1", "cluster_2"), c("cluster_1", "cluster_3"), c("cluster_2", "cluster_3") )

gs1=ggviolin(z.mean.statsM,x = "TE.cluster", y ="mean.exp" , fill = "TE.cluster",alpha = 1,size = 0.3,
         palette = c("#d73027","#00AFBB", "#E69F00"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("mean.expression.of.9.TEs")

gs2=ggviolin(z.mean.statsM,x = "TE.cluster", y ="z.of.mean.exp" , fill = "TE.cluster",alpha = 1,size = 0.3,
         palette = c("#d73027","#00AFBB", "#E69F00"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("z.of.mean.exp")
gs3=ggviolin(z.mean.statsM,x = "TE.cluster", y ="z.of.mean.exp.from.indi.TE" , fill = "TE.cluster",alpha = 1,size = 0.3,
         palette = c("#d73027","#00AFBB", "#E69F00"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("z.of.mean.exp.from.indi.TE")

generate.PDF <- function(fig) {
  pdf("5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/comparson.TE.mean.z.score.pdf",width = 5,height = 8)
  print(gs1)
  print(gs2)
  print(gs3)
  dev.off()
}
generate.PDF(fig)
save(z.mean.statsM,file = "5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/comparson.TE.mean.z.score.RData")
########
sur.Mexp=cbind(surstats, z.mean.statsM[rownames(surstats),])
sur.Mexp$group=ifelse(sur.Mexp$z.of.mean.exp.from.indi.TE> median(sur.Mexp$z.of.mean.exp.from.indi.TE),"2high","1low")
ss=Surv(sur.Mexp$OS.time, sur.Mexp$OS)

#cox = summary(coxph(ss~factor(surstats$TE.culster.agg,levels = c("TE.high","TE.low"))))
cox = summary(coxph(ss~group,data = sur.Mexp))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
##########





