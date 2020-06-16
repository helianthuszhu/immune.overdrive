#####correlation of TE mean expression with GEP
load("4.model/immune/GEP/GEP.stats.RData")
load("4.model/Mean.expression.of.9.TE.model.RData")
#######
length(intersect(rownames(GEP.stats),rownames(mean.stats)))
cor.GEP.te.mean=cbind(GEP.stats[,1:2],mean.stats[rownames(GEP.stats),])
head(cor.GEP.te.mean)

gg1=ggscatter(cor.GEP.te.mean, x = "mean.exp", y = "GEP", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "mean.TE.score", ylab = "GEP")+
  ggtitle("mean.TE.exp vs GEP in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#annotate(geom="text", x=quantile(cd$score)[1], y=quantile(cd$exp)[4], label=paste0("pcor =",coef, "\n", "pvalue =",pvalue),color="red")
generate.PDF <- function(fig) {
  pdf("4.model/immune/GEP.vs.mean.TE.exp/cor.mean.TE.with.GEP.pdf",width = 7,height = 7)
  print(gg1)
  dev.off()
}
generate.PDF(fig)
save(cor.GEP.te.mean, file="4.model/immune/GEP.vs.mean.TE.exp/cor.mean.TE.with.GEP.RData")
