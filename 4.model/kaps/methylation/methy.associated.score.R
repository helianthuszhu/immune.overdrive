#############methylation associated
##########
load("4.model/kaps/clinical.plus.pan/stat.mole.cluster.RData")
load("4.model/kaps/methylation/stat.mole.cluster.methy.RData")
########
###########cimpNature
x=as.data.frame.matrix(table(stat.mole.cluster$cimpNature, stat.mole.cluster$kaps.group))
x$group=rownames(x)
head(x)
write.csv(x, "4.model/kaps/methylation/cimp.proportion.kaps.group.pdf.nature.csv")

pval.kaps.cimp <- chisq.test(stat.mole.cluster$cimpNature, stat.mole.cluster$kaps.group,correct = T)$p.value

stackdata.kaps.cimp=as.data.frame.matrix(table(stat.mole.cluster$cimpNature, stat.mole.cluster$kaps.group))
datm.kaps.cimp <- melt(cbind(stackdata.kaps.cimp, ind = rownames(stackdata.kaps.cimp)), id.vars = c('ind'))

datm.kaps.cimp$ind=factor(datm.kaps.cimp$ind,levels = rev(unique(datm.kaps.cimp$ind)))

#
pdf("4.model/kaps/methylation/cimp.proportion.kaps.group.pdf",width = 4,height = 6)
ggplot(datm.kaps.cimp,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.cimp,4)),size=5)+
  ggtitle("CIMP.vs.TE.kaps.group")
dev.off()
#
#########cimp1
x=as.data.frame.matrix(table(stat.mole.cluster$cimp1, stat.mole.cluster$kaps.group))
x$group=rownames(x)
head(x)
write.csv(x, "4.model/kaps/methylation/cimp.proportion.kaps.group.pdf.three.csv")

pval.kaps.cimp <- chisq.test(subset(stat.mole.cluster, !(cimp1=="not_available"))$cimp1, subset(stat.mole.cluster, !(cimp1=="not_available"))$kaps.group,correct = T)$p.value

stackdata.kaps.cimp=as.data.frame.matrix(table(stat.mole.cluster$cimp1, stat.mole.cluster$kaps.group))
stackdata.kaps.cimp=stackdata.kaps.cimp[1:3,]
datm.kaps.cimp <- melt(cbind(stackdata.kaps.cimp, ind = rownames(stackdata.kaps.cimp)), id.vars = c('ind'))

datm.kaps.cimp$ind=factor(datm.kaps.cimp$ind,levels = rev(unique(datm.kaps.cimp$ind)))
pdf("4.model/kaps/methylation/cimp.proportion.kaps.group.three.pdf",width = 4,height = 6)
ggplot(datm.kaps.cimp,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#1b9e77","#7570b3","#e7298a","#d95f02","#e7298a","#66a61e","#e6ab02"))+
  ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.cimp,4)),size=5)+
  ggtitle("CIMP.vs.TE.kaps.group")
dev.off()

##########
###########
###########global DNA methylation level
#
globalDNA=read.table("4.model/kaps/methylation/global.DNA.methylation.levels.txt",header = T,row.names = 2,sep = "\t")
globalDNA$id1=rownames(globalDNA)
head(globalDNA)
stat.mole.cluster.methy=merge(stat.mole.cluster, globalDNA, by="id1", all.x=TRUE)
head(stat.mole.cluster.methy)
save(stat.mole.cluster.methy,file="4.model/kaps/methylation/stat.mole.cluster.methy.RData")
#
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
stat.mole.cluster.methy$kaps.group=factor(stat.mole.cluster.methy$kaps.group,levels = c("set4", "set3","set2","set1"))


dnap=ggviolin(stat.mole.cluster.methy,x = "kaps.group", y = "Global.methylation.level", fill = "kaps.group",alpha = 1,size = 0.3,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means( comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Global.methylation.level")+ggtitle("Global.methylation.level")

prol=ggviolin(stat.mole.cluster.methy,x = "kaps.group", y = "Proliferation.score", fill = "kaps.group",alpha = 1,size = 0.3,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means( comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Proliferation.score")+ggtitle("Proliferation.score")

stat.mole.cluster.methy[["cimp1"]][is.na(stat.mole.cluster.methy[["cimp1"]])] ="not_available"
#stat.mole.cluster.methy$Global.methylation.level=as.numeric(stat.mole.cluster.methy$Global.methylation.level)

com.dnaM <- list( c("CIMP.High", "CIMP.Low"),c("CIMP.High", "CIMP.Neg"),c("CIMP.Low", "CIMP.Neg")
)
prol2=ggviolin(subset(stat.mole.cluster.methy, !(cimp1=="not_available")),
               x = "cimp1", y = "Global.methylation.level", fill = "cimp1",alpha = 1,size = 0.3,
               #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
               palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
               add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons =com.dnaM ,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Global.methylation.level")+ggtitle("Global.methylation.level")

#write.csv(stat.mole.cluster.methy,"4.model/kaps/methylation/global.DNA.methylation.level.csv")

generate.PDF <- function(fig) {
  pdf("4.model/kaps/methylation/global.DNA.methylation.level.pdf",width = 5,height = 8)
  print(dnap)
  print(prol)
  print(prol2)
  dev.off()
}
generate.PDF(fig)
############correlation with TE mean
ffc=ggscatter(stat.mole.cluster.methy, x = "z.of.mean.exp", y = "Global.methylation.level", 
              add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
              cor.coef = T, cor.method = "spearman",
              xlab = "z.mean.TE.score", ylab = "Global.methylation.level")+
  ggtitle("mean.TE.exp vs Global.methylation.level in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
generate.PDF <- function(fig) {
  pdf("4.model/kaps/methylation/global.DNA.methylation.level.vs.TE.mean.exp.pdf",width = 6,height = 6)
  print(ffc)
  dev.off()
}
generate.PDF(fig)
############
#######hypermutation####
###########
table(stat.mole.cluster.methy$kaps.group, stat.mole.cluster.methy$hypermutation)

x=as.data.frame.matrix(table(stat.mole.cluster.methy$kaps.group, stat.mole.cluster.methy$hypermutation))
x$group=rownames(x)
head(x)
write.csv(x, "4.model/kaps/methylation/hypermutation.vs.kaps.group.csv")

pval.kaps.hyperM <- chisq.test(subset(stat.mole.cluster.methy, !(hypermutation=="not_available"))$kaps.group, 
                               subset(stat.mole.cluster.methy, !(hypermutation=="not_available"))$hypermutation,correct = T)$p.value

stackdata.kaps.hyperM=as.data.frame.matrix(table( stat.mole.cluster.methy$hypermutation,stat.mole.cluster.methy$kaps.group))
stackdata.kaps.hyperM=stackdata.kaps.hyperM[-2,]
datm.kaps.hyperM <- melt(cbind(stackdata.kaps.hyperM, ind = rownames(stackdata.kaps.hyperM)), id.vars = c('ind'))

datm.kaps.hyperM$ind=factor(datm.kaps.hyperM$ind,levels = rev(unique(datm.kaps.hyperM$ind)))
pdf("4.model/kaps/methylation/hypermutation.vs.kaps.group.pdf",width = 4,height = 6)
ggplot(datm.kaps.hyperM,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.hyperM,4)),size=5)+
  ggtitle("hypermutation.vs.TE.kaps.group")
dev.off()
###########################hypermethylation freuqncy
hyperf=read.table("4.model/kaps/methylation/hyperM.freuencyscore.txt",header = T,row.names = 1,sep = "\t")
head(hyperf)
hyperf$id1=gsub("[.]","-",rownames(hyperf))
stat.mole.cluster.methy.hyperf=merge(stat.mole.cluster.methy, hyperf, by="id1",all.x=TRUE)
#
pdf("4.model/kaps/methylation/hyperM.ROS.pdf",width = 4,height = 6)
ggviolin(stat.mole.cluster.methy.hyperf,
         x = "kaps.group", y = "Hypomethylation.Frequency", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Hypomethylation.Frequency")+ggtitle("Hypomethylation.Frequency")

ggviolin(stat.mole.cluster.methy.hyperf,
         x = "kaps.group", y = "Hypermethylation.Frequency", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Hypermethylation.Frequency")+ggtitle("Hypermethylation.Frequency")

ggviolin(stat.mole.cluster.methy.hyperf,
         x = "kaps.group", y = "DMI_score", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("DMI_score")+ggtitle("DMI_score")

ggviolin(stat.mole.cluster.methy.hyperf,
         x = "kaps.group", y = "ROS.Signature.score", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("ROS.Signature.score")+ggtitle("ROS.Signature.score")
ggviolin(stat.mole.cluster.methy.hyperf,
         x = "kaps.group", y = "H2O2.Signature.score", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("H2O2.Signature.score")+ggtitle("H2O2.Signature.score")
dev.off()
#
#
#
####aneuploidy score
#####
aneuscore=read.delim2("4.model/kaps/methylation/aneuploidy.score.txt",header = T,sep = "\t",na.strings = F,row.names = 1)
aneuscore$id1=rownames(aneuscore)
head(aneuscore) 
#
stat.mole.cluster.methy.aneu=merge(stat.mole.cluster.methy, aneuscore[,-1], by="id1",all.x=TRUE)

#
pdf("4.model/kaps/methylation/aneuploidy.score.pdf",width = 4,height = 6)
ggviolin(stat.mole.cluster.methy.aneu,
         x = "kaps.group", y = "AneuploidyScore", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons ,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("AneuploidyScore")+ggtitle("AneuploidyScore")
dev.off()
#

table(stat.mole.cluster.methy.aneu$kaps.group, stat.mole.cluster.methy.aneu$X4q)
chisq.test(stat.mole.cluster.methy.aneu$kaps.group, stat.mole.cluster.methy.aneu$X,correct = T)
#
datalist=list()
for (i in 102:157) {
  tmppvalue <- chisq.test(stat.mole.cluster.methy.aneu[,39], stat.mole.cluster.methy.aneu[,i],correct = T)$p.value
  tmppvalue=as.data.frame(tmppvalue)
  rownames(tmppvalue)=colnames(stat.mole.cluster.methy.aneu)[i]
  datalist[[i]]=tmppvalue
}
#
armpvalue.2=do.call(rbind, datalist)
armpvalue.2$tmppvalue=round(armpvalue.2$tmppvalue,digits = 4)
#
pdf("4.model/kaps/methylation/aneuploidy.score.cor.TE.mean.pdf")
ggscatter(stat.mole.cluster.methy.aneu, x = "z.of.mean.exp", y = "AneuploidyScore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "AneuploidyScore")+
  ggtitle("mean.TE.exp vs AneuploidyScore in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
#########
#####scna score
#########
scanscore=read.table("4.model/kaps/methylation/cSCNA.score.txt",header = T,sep = "\t",row.names = 1)
scanscore$idd=rownames(scanscore)
stat.mole.cluster.methy.scan=merge(stat.mole.cluster.methy, scanscore[,-c(1:7)], by="idd",all.x=TRUE)
#######
table(stat.mole.cluster.methy.scan$kaps.group,stat.mole.cluster.methy.scan$TP53Mutation)
pdf("4.model/kaps/methylation/cSCNA.score.pdf",width = 4,height = 6)
ggviolin(stat.mole.cluster.methy.scan,
         x = "kaps.group", y = "CellCycle.Signature.Score", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("CellCycle.Signature.Score")+ggtitle("CellCycle.Signature.Score")

ggviolin(stat.mole.cluster.methy.scan,
         x = "kaps.group", y = "SCNA.Level", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("SCNA.Level")+ggtitle("SCNA.Level")

ggviolin(stat.mole.cluster.methy.scan,
         x = "kaps.group", y = "Chrom.SCNA.Level", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Chrom.SCNA.Level")+ggtitle("Chrom.SCNA.Level")

ggviolin(stat.mole.cluster.methy.scan,
         x = "kaps.group", y = "Arm.SCNA.Level", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Arm.SCNA.Level")+ggtitle("Arm.SCNA.Level")

ggviolin(stat.mole.cluster.methy.scan,
         x = "kaps.group", y = "Focal.SCNA.Level", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Focal.SCNA.Level")+ggtitle("Focal.SCNA.Level")

ggviolin(stat.mole.cluster.methy.scan,
         x = "kaps.group", y = "Chrom.Arm.SCNA.Level", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Chrom.Arm.SCNA.Level")+ggtitle("Chrom.Arm.SCNA.Level")

ggviolin(stat.mole.cluster.methy.scan,
         x = "kaps.group", y = "SCNA.Level.normalized.by.size", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("SCNA.Level.normalized.by.size")+ggtitle("SCNA.Level.normalized.by.size")
dev.off()
################
###############correlation
#####
pdf("4.model/kaps/methylation/scnascore.cor.TE.mean.pdf")

ggscatter(stat.mole.cluster.methy.scan, x = "z.of.mean.exp", y = "SCNA.Level", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "SCNA.Level")+
  ggtitle("mean.TE.exp vs SCNA.Level in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(stat.mole.cluster.methy.scan, x = "z.of.mean.exp", y = "Chrom.SCNA.Level", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "Chrom.SCNA.Level")+
  ggtitle("mean.TE.exp vs Chrom.SCNA.Level in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")


ggscatter(stat.mole.cluster.methy.scan, x = "z.of.mean.exp", y = "Arm.SCNA.Level", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "Arm.SCNA.Level")+
  ggtitle("mean.TE.exp vs Arm.SCNA.Level in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(stat.mole.cluster.methy.scan, x = "z.of.mean.exp", y = "Focal.SCNA.Level", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "Focal.SCNA.Level")+
  ggtitle("mean.TE.exp vs Focal.SCNA.Level in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(stat.mole.cluster.methy.scan, x = "z.of.mean.exp", y = "Chrom.Arm.SCNA.Level", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "Chrom.Arm.SCNA.Level")+
  ggtitle("mean.TE.exp vs Chrom.Arm.SCNA.Level in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(stat.mole.cluster.methy.scan, x = "z.of.mean.exp", y = "SCNA.Level.normalized.by.size", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "SCNA.Level.normalized.by.size")+
  ggtitle("mean.TE.exp vs SCNA.Level.normalized.by.size in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

dev.off()
#####################
########PMD score
####
pmdscore=read.table("4.model/kaps/methylation/TCGA_HM450_PMD_HMD_depths.tsv",header = T,sep = '\t')
pmdscore$id1=substr(pmdscore$barcode, 1, 15)
head(pmdscore)
dim(pmdscore)
pmdscore.agg=pmdscore[,4:7 ] %>% group_by(id1) %>% summarise_all(mean)
pmdscore.agg=as.data.frame(pmdscore.agg)
head(pmdscore.agg)
stat.mole.cluster.methy.pmdscore=merge(stat.mole.cluster.methy, pmdscore.agg, by="id1",all.x=TRUE)
head(stat.mole.cluster.methy.pmdscore)
#####
pdf("4.model/kaps/methylation/PMD.HMD.score.pdf",width = 4,height = 6)
ggviolin(stat.mole.cluster.methy.pmdscore,
         x = "kaps.group", y = "commonPMD", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("commonPMD")+ggtitle("commonPMD")

ggviolin(stat.mole.cluster.methy.pmdscore,
         x = "kaps.group", y = "commonHMD", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("commonHMD")+ggtitle("commonHMD")

ggviolin(stat.mole.cluster.methy.pmdscore,
         x = "kaps.group", y = "PMD_minus_HMD", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("PMD_minus_HMD")+ggtitle("PMD_minus_HMD")
dev.off()
#####
pdf("4.model/kaps/methylation/PMD.HMD.score.cor.TE.mean.pdf")
ggscatter(stat.mole.cluster.methy.pmdscore, x = "z.of.mean.exp", y = "commonPMD", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "commonPMD")+
  ggtitle("mean.TE.exp vs commonPMD in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#
ggscatter(stat.mole.cluster.methy.pmdscore, x = "z.of.mean.exp", y = "commonHMD", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "commonHMD")+
  ggtitle("mean.TE.exp vs commonHMD in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#
ggscatter(stat.mole.cluster.methy.pmdscore, x = "z.of.mean.exp", y = "PMD_minus_HMD", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "PMD_minus_HMD")+
  ggtitle("mean.TE.exp vs PMD_minus_HMD in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#
ggscatter(stat.mole.cluster.methy.pmdscore, x = "commonPMD", y = "commonHMD", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "PMD_vs_HMD")+
  ggtitle("PMD_vs_HMD in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
###################
pdf("4.model/kaps/methylation/TET.genes.pdf")
ggscatter(genestats, x = "z.of.mean.exp", y = "TET1", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "TET1")+
  ggtitle("TET1 in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
ggscatter(genestats, x = "z.of.mean.exp", y = "TET2", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "TET2")+
  ggtitle("TET2 in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(genestats, x = "z.of.mean.exp", y = "TET3", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "TET3")+
  ggtitle("TET3 in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
########
