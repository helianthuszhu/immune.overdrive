########PRJEB23709
clin.PRJEB23709=read.table("6.immnue.treat.data/PRJEB23709/PRJEB23709.clinicaldata.txt",header=T,sep="\t")
clin.PRJEB23709$id=paste(clin.PRJEB23709$drug, clin.PRJEB23709$ptid,clin.PRJEB23709$RNA.Sequencing,sep="_")
head(clin.PRJEB23709)
rownames(clin.PRJEB23709)=clin.PRJEB23709$id
###
table(clin.PRJEB23709$RNA.Sequencing)
########TE expression
teexp.PRJEB23709=readRDS("6.immnue.treat.data/PRJEB23709/RE.output/RE_intergenic_2_counts_normalized.RDS")
teexp.PRJEB23709.signature=as.data.frame(exprs(teexp.PRJEB23709))[rownames(cndi.rep),]
teexp.PRJEB23709.signature=as.data.frame(t(teexp.PRJEB23709.signature))
teexp.PRJEB23709.signature$mean.exp=rowMeans(teexp.PRJEB23709.signature)
teexp.PRJEB23709.signature$z.of.mean.exp=(teexp.PRJEB23709.signature$mean.exp - mean(teexp.PRJEB23709.signature$mean.exp))/sd(teexp.PRJEB23709.signature$mean.exp)
head(teexp.PRJEB23709.signature)
#
length(intersect(rownames(teexp.PRJEB23709.signature), clin.PRJEB23709$id))
idshPRJEB23709=intersect(rownames(clin.PRJEB23709), rownames(teexp.PRJEB23709.signature))
##
stat.PRJEB23709=cbind(teexp.PRJEB23709.signature[idshPRJEB23709,], clin.PRJEB23709[idshPRJEB23709,])
stat.PRJEB23709$response2=gsub(" ","",stat.PRJEB23709$response)
stat.PRJEB23709$response3=ifelse(stat.PRJEB23709$response2=="PD","PD","PRCRSD")
head(stat.PRJEB23709)
table(stat.PRJEB23709$response2)
####
pdf("6.immnue.treat.data/PRJEB23709/TE.score.comparison.pdf",width = 5,height = 6)
ggboxplot(stat.PRJEB23709, x = "response3", y = "z.of.mean.exp",
          color = "response3", palette = c("PD"="#d01c8b","PRCRSD"="#4dac26"),
          add = "jitter")+stat_compare_means()+ggtitle("PRJEB23709")
ggboxplot(stat.PRJEB23709, x = "response3", y = "MER57F",
          color = "response3", palette = c("PD"="#d01c8b","PRCRSD"="#4dac26"),
          add = "jitter")+stat_compare_means()+ggtitle("PRJEB23709")
ggboxplot(stat.PRJEB23709, x = "response3", y = "MER65C",
          color = "response3", palette = c("PD"="#d01c8b","PRCRSD"="#4dac26"),
          add = "jitter")+stat_compare_means()+ggtitle("PRJEB23709")
dev.off()
