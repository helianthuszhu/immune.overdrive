######compare stemness
#########
stemness.rna=read.table("4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/StemnessScores_RNAexp_20170127.2.tsv.gz",sep = "\t",header = T,row.names = 1)
head(stemness.rna)
stemness.rna=as.data.frame(t(stemness.rna))
rownames(stemness.rna)=gsub("[.]","-",rownames(stemness.rna))
####
stemid.rna=intersect(rownames(kaps.td), rownames(stemness.rna))
####
stat.stem.rna=cbind(kaps.td[stemid.rna,], stemness.rna[stemid.rna,])
stat.stem.rna$kaps.group=factor(stat.stem.rna$kaps.group,levels = c("set4","set3","set2","set1"))
head(stat.stem.rna)
#####
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
sg1=ggviolin(stat.stem.rna,x = "kaps.group", y ="RNAss" , fill = "kaps.group",alpha = 1,size = 0.3,
             palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
             add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("RNAss") +ggtitle("CRC.cohort")

sg2=ggviolin(stat.stem.rna,x = "kaps.group", y ="EREG.EXPss" , fill = "kaps.group",alpha = 1,size = 0.3,
         palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("EREG.EXPss") +ggtitle("CRC.cohort")
#####stemness DNA
######
stemness.dna=read.table("4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/StemnessScores_DNAmeth_20170210.tsv.gz",sep = "\t",header = T,row.names = 1)
head(stemness.dna)
stemness.dna=as.data.frame(t(stemness.dna))
rownames(stemness.dna)=gsub("[.]","-",rownames(stemness.dna))
colnames(stemness.dna)[2]="EREG.METHss"
#
####
stemid.dna=intersect(rownames(kaps.td), rownames(stemness.dna))
####
stat.stem.dna=cbind(kaps.td[stemid.dna,], stemness.dna[stemid.dna,])
stat.stem.dna$kaps.group=factor(stat.stem.dna$kaps.group,levels = c("set4","set3","set2","set1"))
head(stat.stem.dna)
#####
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
sg3=ggviolin(stat.stem.dna,x = "kaps.group", y ="DNAss" , fill = "kaps.group",alpha = 1,size = 0.3,
         palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("DNAss") +ggtitle("CRC.cohort")

sg4=ggviolin(stat.stem.dna,x = "kaps.group", y ="EREG.METHss" , fill = "kaps.group",alpha = 1,size = 0.3,
         palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("EREG-METHss") +ggtitle("CRC.cohort")
####
generate.PDF <- function(fig) {
  pdf("4.model/kaps/stemness/stemnetss.pdf",height = 3,width = 5)
  print(sg1)
  print(sg2)
  print(sg3)
  print(sg4)
  dev.off()
}
generate.PDF(fig)
save(stat.stem.dna,stat.stem.rna,file = "4.model/kaps/stemness/stemnetss.RData")
