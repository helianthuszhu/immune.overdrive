#####pole id
polesampleid=read.csv("4.model/kaps/ploe.compare/POLE_CRC_UCEC.sampleid.csv",header = T)
head(polesampleid)
rownames(polesampleid)=polesampleid$sample
polesampleid$pole.mut=rep("yes",nrow(polesampleid))
#
polesampleid.CRC=subset(polesampleid, cancer_type=="CRC")
dim(polesampleid.CRC)
head(polesampleid.CRC)
#
head(kaps.td)
aa=data.frame(kaps.group=kaps.td$kaps.group)
rownames(aa)=rownames(kaps.td)
aa$sample=substr(rownames(aa),1,12)
head(aa)
aa$pole.mut=ifelse(aa$sample %in% polesampleid.CRC$sample,"pole.mut","pole.non.mut")
table(aa$pole.mut,aa$kaps.group)
chisq.test(aa$pole.mut,aa$kaps.group,correct = T)
write.csv(aa,"4.model/kaps/ploe.compare/kaps.group.pole.mut.csv")
#######
#######
