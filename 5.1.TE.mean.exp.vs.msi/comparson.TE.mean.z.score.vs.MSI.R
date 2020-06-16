head(mean.stats)
msi.com=mean.stats[!is.na(mean.stats$MSI.status),]
#table(msi.com$MSI.status, msi.com$TE.cluster.agg)
msi.com$MSI.status.bin=ifelse(msi.com$MSI.status=="MSI-H","MSI","MSS")
msi.com$msi.TE.group=paste(msi.com$TE.cluster.agg, msi.com$MSI.status,sep = "_")
msi.com$msi.TE.group.bin=paste(msi.com$TE.cluster.agg, msi.com$MSI.status.bin,sep = "_")
msi.com=cbind(msi.com,z.mean.statsM[rownames(msi.com), 11:12])
dim(msi.com)
head(msi.com)

######
my_comparisons <- list( c("TE.high_MSI", "TE.low_MSI"), c("TE.high_MSS", "TE.low_MSS"),
                        c("TE.high_MSI", "TE.high_MSS"), c("TE.low_MSI", "TE.low_MSS"))
msi.com$msi.TE.group.bin=factor(msi.com$msi.TE.group.bin, c("TE.high_MSI", "TE.high_MSS","TE.low_MSI", "TE.low_MSS"))
#
gmsi1=ggviolin(msi.com,x = "msi.TE.group.bin", y ="mean.exp" , fill = "msi.TE.group.bin",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("mean.expression.of.9.TEs")+ggtitle("TE.mean.exp vs MSI.status")
gmsi2=ggviolin(msi.com,x = "msi.TE.group.bin", y ="z.of.mean.exp" , fill = "msi.TE.group.bin",alpha = 1,size = 0.3,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("z.of.mean.exp")+ggtitle("z.of.mean.exp vs MSI.status")
gmsi3=ggviolin(msi.com,x = "msi.TE.group.bin", y ="z.of.mean.exp.from.indi.TE" , fill = "msi.TE.group.bin",alpha = 1,size = 0.3,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("z.of.mean.exp.from.indi.TE")+ggtitle("z.of.mean.exp.from.indi.TE vs MSI.status")

generate.PDF <- function(fig) {
  pdf("5.1.TE.mean.exp.vs.msi/comparson.TE.mean.z.score.vs.MSI.pdf",width = 5,height = 6)
  print(gmsi1)
  print(gmsi2)
  print(gmsi3)
  dev.off()
}
generate.PDF(fig)
######
save(msi.com, file = "5.1.TE.mean.exp.vs.msi/comparson.TE.mean.z.score.vs.MSI.RData")
#############
#############
