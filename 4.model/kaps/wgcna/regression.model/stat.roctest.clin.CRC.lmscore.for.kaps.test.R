load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/wgcna/regression.model/stat.roctest.clin.CRC.lmscore.for.kaps.test.RData")
library(kaps)
stat.roctest.clin.lm.kaps <- kaps(survival::Surv(OS.time, OS) ~ lmscore, data = stat.roctest.clin, K = 4)
save(stat.roctest.clin.lm.kaps, file="~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/wgcna/regression.model/stat.roctest.clin.CRC.lmscore.for.kaps.test.output.RData")

