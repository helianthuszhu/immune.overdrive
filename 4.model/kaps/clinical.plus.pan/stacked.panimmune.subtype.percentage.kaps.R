load("4.model/kaps/clinical.plus.pan/paltform.cluster.pancancer.RData")
##########
rownames(cb3)=cb3$id1
head(cb3)
length(intersect(rownames(cb3), rownames(kaps.td)))
head(kaps.td)
tttd=kaps.td
tttd$id1=rownames(tttd)
head(tttd)
#
pc.crc=subset(cb3, Cancer.Type.x=="COAD"|Cancer.Type.x=="READ")
dim(pc.crc)
#####
stat.pan.cluster=merge(tttd, pc.crc,all.x=TRUE, by="id1")
head(stat.pan.cluster)
dim(stat.pan.cluster)
table(stat.pan.cluster$kaps.group)
stat.pan.cluster$kaps.group=factor(stat.pan.cluster$kaps.group,levels = c("set4", "set3","set2","set1"))
#####
table(stat.pan.cluster$cim, stat.pan.cluster$kaps.group)
#########
#####################################
library(reshape2)
library(ggplot2)
library(scales)
x=as.data.frame.matrix(table(stat.pan.cluster$Subtype_Immune_Model_Based, stat.pan.cluster$kaps.group))
x$group=rownames(x)
head(x)
write.csv(x, "4.model/kaps/clinical.plus.pan/immune.subtype.vs.kaps.group.csv")

pval.kaps.immu <- chisq.test(stat.pan.cluster$Subtype_Immune_Model_Based, stat.pan.cluster$kaps.group,correct = T)$p.value

stackdata.kaps.immu=as.data.frame.matrix(table(stat.pan.cluster$Subtype_Immune_Model_Based, stat.pan.cluster$kaps.group))
datm.kaps.immu <- melt(cbind(stackdata.kaps.immu, ind = rownames(stackdata.kaps.immu)), id.vars = c('ind'))

datm.kaps.immu$ind=factor(datm.kaps.immu$ind,levels = rev(unique(datm.kaps.immu$ind)))

gc1=ggplot(datm.kaps.immu,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.immu,4)),size=5)+
  ggtitle("pan.cancer.immune.subtype.vs.TE.kaps.group")
#annotate("text", x=2, y=0.9, label=paste0("p-value: ","\n" ,signif(pval2,4)),size=3)

generate.PDF <- function(fig) {
  pdf(paste0("4.model/kaps/clinical.plus.pan/","stacked.panimmune.subtype.percentage.kaps",".pdf"),width = 5,height = 5)
  print(gc1)
  dev.off()
}
generate.PDF(fig)
####################
###########TCGA cluster
######
clusterTCGA=read.table("4.model/kaps/clinical.plus.pan/cNMF/gdac.broadinstitute.org_COADREAD-TP.Aggregate_Molecular_Subtype_Clusters.Level_4.2016012800.0.0/COADREAD-TP.mergedcluster.txt",header = T,row.names = 1,sep = "\t")
head(clusterTCGA)

######combine immune subtype
tpp=stat.pan.cluster
rownames(tpp)=tpp$id1
head(tpp)
length(intersect(rownames(tpp), rownames(clin.data.CRC.selected)))
######
stat.mole.cluster=cbind(tpp, clusterTCGA[rownames(tpp),],clin.data.CRC.selected[rownames(tpp),])
head(stat.mole.cluster)
table(stat.mole.cluster$kaps.group, stat.mole.cluster$TCGA_subtypes)
save(stat.mole.cluster,file = "4.model/kaps/clinical.plus.pan/stat.mole.cluster.RData")
###################

