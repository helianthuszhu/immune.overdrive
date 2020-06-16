load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/1.tescore.survival.kaps.validation/for.kaps.classification.RData")
library(kaps)
fitkaps.kirc <- kaps(survival::Surv(OS.time, OS) ~ z.of.mean.exp, data = stat.kirc.kaps.vali, K = 4)
save.image(file="~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/1.tescore.survival.kaps.validation/for.kaps.classification.output.k4.kirc.RData")
