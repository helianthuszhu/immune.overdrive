########l1 insertion
#####
l1insertion=read.csv("4.model/kaps/l1.insertion.NG/l1.inseration.Ngenetics.csv",header = T)
head(l1insertion)
l1insertion=l1insertion[,c(1,2,13)]
l1insertion$id=paste(l1insertion$icgc_donor_id, l1insertion$tumor_wgs_icgc_sample_id,sep = "_")
l1insertion=l1insertion[,3:4]
#
l1insertion.agg= l1insertion %>% group_by(id) %>% summarise_all(mean)
l1insertion.agg= as.data.frame(l1insertion.agg)
l1insertion.agg=separate(l1insertion.agg, col = "id",into=c("icgc.donor.id","icgc.sample.id"),sep = "_")
l1insertion.agg=l1insertion.agg[!(duplicated(l1insertion.agg$icgc.donor.id)),]
head(l1insertion.agg)
rownames(l1insertion.agg)=l1insertion.agg$icga.donor.id
colnames(l1insertion.agg)=gsub("[.]","_",colnames(l1insertion.agg))
dim(l1insertion.agg)
write.csv(l1insertion.agg,"4.model/kaps/l1.insertion.NG/l1insertion.tmps.csv")
#################
icgc.sampleid=read.table("4.model/kaps/l1.insertion.NG/specimen.all_projects.tsv",header = T,sep = "\t")
icgc.sampleid=icgc.sampleid[,1:7]
icgc.sampleid=icgc.sampleid[grepl("TCGA",icgc.sampleid$submitted_donor_id),]
icgc.sampleid=icgc.sampleid[!(duplicated(icgc.sampleid$submitted_donor_id)),]
head(icgc.sampleid)
dim(icgc.sampleid)
length(unique(icgc.sampleid$icgc_donor_id))
#
length(intersect(l1insertion.agg$icga.donor.id, icgc.sampleid$icgc_donor_id))
#####
l1sertion.tcga=merge(icgc.sampleid, l1insertion.agg, by="icgc_donor_id",all=FALSE)
l1sertion.tcga$id=substr(l1sertion.tcga$submitted_specimen_id, 1, 15)
l1sertion.tcga$control=substr(l1sertion.tcga$submitted_specimen_id, 14, 15)
l1sertion.tcga.sel=subset(l1sertion.tcga, control=="01")
head(l1sertion.tcga.sel)
#####
table(l1sertion.tcga.sel$control)
##########
##########
head(pan9TE.stat)
aa=pan9TE.stat
aa$id=rownames(aa)
###
l1sertion.tcga.stat=merge(aa, l1sertion.tcga.sel, by="id",all=FALSE)
dim(l1sertion.tcga.stat)
####
pdf("4.model/kaps/l1.insertion.NG/cor.l1insertion.vs.TE.score.pdf")
ggscatter(l1sertion.tcga.stat, x = "z.of.mean.exp", y = "score", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "TE score", ylab = "L1.insertion.score")+
  ggtitle("Cor.TE score vs L1 insertion score")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
