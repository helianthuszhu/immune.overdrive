load("~/nas/Xiaoqiang/R.pacakages/TE/REdiscoverTE/rollup_annotation/rep.sequence.index.REdis.RData")
head(rep.sequence.index)
head(cndi.rep)

rep.sequence.index.9tes=rep.sequence.index[rep.sequence.index$repName %in% rownames(cndi.rep),]
dim(rep.sequence.index.9tes)
table(rep.sequence.index.9tes$selected_feature)
#write.csv(rep.sequence.index.9tes, "7.TE.sequences/rep.sequence.index.9tes.csv")
for (i  in 1:length(cndi.id)) {
  reps.seq=rep.sequence.index[rep.sequence.index$repName %in% cndi.id[i],]
  write.table(reps.seq[,1], paste0("7.TE.sequences/index/index.of.",cndi.id[i],".txt"),quote = F,sep = "\t",row.names = F,col.names = F)
}
reps.seq=rep.sequence.index[rep.sequence.index$repName %in% cndi.id[1],]
write.table(reps.seq[1:2,1],"7.TE.sequences/test.bed.file.txt",quote = F,sep = "\t",row.names = F,col.names = F)
head(reps.seq)
