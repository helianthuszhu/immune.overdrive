load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/for.kaps.classification.RData")
library(kaps)

fitkaps <- kaps(survival::Surv(OS.time, OS) ~ z.of.mean.exp, data = subset(cor.immune.drive.sur, OS.time >=0), K = 2:4)
save.image(file="~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/for.kaps.classification.output.k2-4.RData")
