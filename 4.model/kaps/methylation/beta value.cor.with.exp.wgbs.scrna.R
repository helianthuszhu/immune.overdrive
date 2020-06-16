##########get the postion of TE
#######
hg19rep=read.delim2("~/nas/Xiaoqiang/1.retrotransposon/data/RepeatMasker.TE.ref/TE.anno.from.repeatmasker.web/hg19.fa.out",header = T,sep = "\t")
colnames(hg19rep)="id"
#hg19rep$id=gsub(" ","",hg19rep$id)
head(hg19rep)
nchar(hg19rep$id)

hg19rep <- hg19rep %>%
  rowwise() %>%
  mutate_all(funs(str_squish(.))) %>%
  ungroup()

hg19repsep=separate(hg19rep, col = "id",into = c("id1","id2","id3","id4","id5","id6","id7","id8",
                                    "id9","id10","id11","id12","id13","id14","id15"),sep = " ")
hg19repsep=hg19repsep[-1,]
head(hg19repsep)
write.table(hg19repsep[,-c(1:5)],"7.TE.sequences/hg19repsep.txt",quote = F,sep = "\t",row.names = F,col.names = F)
#########
########
#####methylation matrix
#######
mMatrix=read.table("4.model/kaps/methylation/dataxena/HumanMethylation450.gz", header = T,row.names = 1)
mMatrix[1:8,1:4]
dim(mMatrix)
######filter NA
countNAs=as.data.frame(rowSums(is.na(mMatrix)))
colnames(countNAs)="count"
head(countNAs)
countNAs.filter=subset(countNAs, count < 443)
head(countNAs.filter)
#
mMatrix.filter=mMatrix[rownames(countNAs.filter),]
mMatrix.filter[1:4,1:4]
######################
########methylation probe annotation
#####
#meprob=read.csv("4.model/kaps/methylation/dataxena/GPL13534_HumanMethylation450_15017482_v.1.1.clean.csv",header = T,sep = ",")
meprob<-meprob[,c("IlmnID",
          "Infinium_Design_Type",
          "CHR",
          "MAPINFO",
          "UCSC_RefGene_Name",
          "UCSC_RefGene_Group",
          "UCSC_CpG_Islands_Name",
          "Relation_to_UCSC_CpG_Island")]
meprob$CHR=paste0("chr",meprob$CHR)
meprob$start=meprob$MAPINFO-1
meprob$end=meprob$MAPINFO+1
meprob=meprob[,-2]
colnames(meprob)=c("id","chr", "possite","gene", "genepos","islandNane","island","start","end")
head(meprob)
##################
meprob=read.table("4.model/kaps/methylation/dataxena/illuminaMethyl450_hg19_GPL16304_TCGAlegacy",header = F,sep = "\t")
meprob=na.omit(meprob, cols = c("V3"))
colnames(meprob)=c("proName","gene","chr","start","end","strand")
head(meprob)
length(intersect(rownames(mMatrix.filter), meprob$proName))
#
mMatrix.filter=mMatrix.filter[rownames(mMatrix.filter) %in% meprob$proName,]
mMatrix.filter[1:4,1:4]
dim(mMatrix.filter)
##########################
##########################
###########################
#####check the mean beta value of the matrix
#######
global.mean.beta=as.data.frame(colMeans(mMatrix.filter,na.rm = T))
colnames(global.mean.beta)="mena.beta"
global.mean.beta$control=substr(rownames(global.mean.beta),14,15)
head(global.mean.beta)
table(global.mean.beta$control)
global.mean.beta.sel=subset(global.mean.beta, control=="01"|control=="11")

pdata=as.data.frame(rownames(global.mean.beta.sel))
head(pdata)
rownames(pdata)=pdata$`rownames(global.mean.beta)`
colnames(pdata)="id"
pdata$idd=substr(pdata$id,1,12)
pdata$iddd=substr(pdata$id,14,15)
head(pdata)
paired_p=names(table(pdata$idd)[table(pdata$idd)==2])
ppdata=pdata[pdata$idd %in% paired_p,]
ppdata=ppdata[order(ppdata$idd,ppdata$iddd,decreasing = T),]
head(ppdata)
global.mean.beta.parired=global.mean.beta[rownames(global.mean.beta) %in% ppdata$id,]
head(global.mean.beta.parired)

#
pdf("4.model/kaps/methylation/dataxena/mean.beta.compare.pdf",width = 3,height = 5)
ggviolin(global.mean.beta.parired,x = "control", y = "mena.beta", fill = "control",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "control")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("control")+ylab("Global.methylation.level")+ggtitle("mean.beta.value")
dev.off()
###########################################################################
#### compare the mean beta value among kaps group
######
head(global.mean.beta)
rownames(global.mean.beta)=gsub("[.]","-",rownames(global.mean.beta))
#
length(intersect(rownames(global.mean.beta), rownames(kaps.td)))
global.id=intersect(rownames(global.mean.beta), rownames(kaps.td))
######
stat.methy.xena=cbind(kaps.td[global.id,], global.mean.beta[global.id,])
head(stat.methy.xena)
save(stat.methy.xena,file="4.model/kaps/methylation/dataxena/kaps.vs.mean.beta.RData")
#
pdf("4.model/kaps/methylation/dataxena/kaps.vs.mean.beta.pdf")
ggviolin(stat.methy.xena,x = "kaps.group", y = "mena.beta", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("control")+ylab("Global.methylation.level")+ggtitle("mean.beta.value")
##
ggscatter(stat.methy.xena, x = "z.of.mean.exp", y = "mena.beta", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "mena.beta")+
  ggtitle("mean.TE.exp vs mena.beta in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
####################################
###################################
##################################
######calculate individual TE methylation level
########
library(GenomicRanges)
library(IRanges)
library(data.table)
fi<-list.files(full.names=T)
#fi
df1 <- data.frame(chr=hg19repsep$id5,
                  start=hg19repsep$id6,
                  end=hg19repsep$id7,
                  subfamily=hg19repsep$id10,
                  class=hg19repsep$id11
)
df1=df1%>% separate(class,c("class","family"),"/")
df1$family=ifelse(is.na(df1$family), "Simple_repeat", paste(df1$family))
df1$start=as.numeric(paste(df1$start))
df1$end=as.numeric(paste(df1$end))
head(df1)

#
#df2=meprob[grep("cg",meprob$id),c(2,8,9,1,4:7)]
#colnames(df2)=c("chr","start","end","probeN","gene")
df2=meprob[,c(3:5,1:2)]
head(df2)
df2=df2[order(df2$chr,df2$start,df2$end),]
df2=df2[grep("chr",df2$chr),]
#
idokm=intersect(rownames(mMatrix.filter), df2$proName)
mMatrix.filter.ok=mMatrix.filter[idokm,]
df2=df2[df2$proName %in% idokm,]
head(df2)
dim(mMatrix.filter.ok)
mMatrix.filter.ok[1:4,1:4]
#
#start at subfamily
TE.list=cndi.rep
head(TE.list)
length(intersect(unique(df1$subfamily), rownames(TE.list)))
#
datalist=list()
te.list=c(rownames(TE.list),"L1HS","L1PA")

for (k in 1:length(te.list)){
  sub.df1=as.data.frame(subset(df1,df1$subfamily==te.list[k]))
  # next, create IRanges objects
  bed1=with(sub.df1, GRanges(chr, IRanges(start+1, end), subfamily, class, family,strand = NULL))
  head(bed1)
  #length(intersect(rownames(mMatrix.filter), meprob$id))
  
  bed2=with(df2, GRanges(chr, IRanges(start+1, end), proName,gene,strand = NULL))
  head(bed2)
  length(bed2)
  # now find the overlaps
  df3=as.data.frame(subsetByOverlaps(bed2, bed1))
  head(df3)
  table(df3$island)
  #df3$percentage=as.numeric(paste(df3$percentage))
  #df3=na.omit(df3)
  calmmatrix=mMatrix.filter.ok[rownames(mMatrix.filter.ok) %in% df3$proName, ]
  
  calmean=as.data.frame(colMeans(calmmatrix,na.rm = T))
  colnames(calmean)=paste("mean.beta",te.list[k],sep = ".")
  head(calmean)
  datalist[[k]]=calmean
  
}
te.9.mean.beta=do.call(cbind, datalist)
rownames(te.9.mean.beta)=gsub("[.]","-",rownames(te.9.mean.beta))
class(te.9.mean.beta)
head(te.9.mean.beta)
##########
###########
idmeanok=intersect(rownames(kaps.td), rownames(te.9.mean.beta))
#
stat.meanbeta=cbind(kaps.td[idmeanok,], te.9.mean.beta[idmeanok,])
head(stat.meanbeta)
#####
stat.meanbeta.cor=stat.meanbeta[,c(3:11,41:49)]
stat.meanbeta.cor=stat.meanbeta.cor[,-c(2,3,5,11,12,14)]

pdf("4.model/kaps/methylation/cor.exp.vs.methy.9.TEs.pdf",width = 10,height = 10)
ggcorr(stat.meanbeta.cor,nbreaks = 4,geom = "circle",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
ggcorr(stat.meanbeta.cor,nbreaks = 4,geom = "text",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
ggcorr(stat.meanbeta.cor,nbreaks = 4,method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,label = T,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)

dev.off()
save(stat.meanbeta.cor,file = "4.model/kaps/methylation/cor.exp.vs.methy.9.TEs.RData")
#########################
#############wgbs
##########
load("~/nas/Xiaoqiang/opti.data/TCGA.data/pooled.samples.with.WGBS/fq/proced.data/for.TE.9.analysis.expM.RData")
wgbsM=read.table("~/nas/Xiaoqiang/opti.data/TCGA.data/pooled.samples.with.WGBS/output/CB.summarized.subfamily.data/subfamily.methylation.matrix.47.1241.txt",header = T,row.names = 1,sep = "\t")
wgbsM=as.data.frame(t(wgbsM))
rownames(wgbsM)=gsub("[.]","-",rownames(wgbsM))
wgbsgene=as.data.frame(t(exp.33s))
length(intersect(rownames(wgbsgene), rownames(wgbsM)))
wgbsM=wgbsM[rownames(wgbsgene),]
#
wgbs9M=wgbsM[, colnames(wgbsM) %in% rownames(cndi.rep)]
colnames(wgbs9M)=paste("methy",colnames(wgbs9M),sep = ".")
wgbs9g=wgbsgene[, colnames(wgbsgene) %in% rownames(cndi.rep)]
colnames(wgbs9g)=paste("exp",colnames(wgbs9g),sep = ".")
#
wgbs9g[1:4,1:4]
wgbs9M[1:4,1:4]
####
wgbs.stat=cbind(wgbs9g,wgbs9M)
colnames(wgbs.stat)=gsub("-",".",colnames(wgbs.stat))
head(wgbs.stat)
##########
ggscatter(wgbs.stat, x = "exp.MER57F", y = "methy.MER57F", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.mean.TE.score", ylab = "Global.methylation.level")+
  ggtitle("mean.TE.exp vs Global.methylation.level in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
#####
library(GGally)
pdf("4.model/kaps/methylation/tcga.wgbs.33/exp.vs.methy.cor.pdf",width = 10,height = 10)
ggcorr(wgbs.stat,nbreaks = 4,geom = "circle",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
            min_size=4,max_size = 8,
            low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
ggcorr(wgbs.stat,nbreaks = 4,geom = "text",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)

ggcorr(wgbs.stat,nbreaks = 4,method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,label = T,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
dev.off()
save(wgbs.stat, file = "4.model/kaps/methylation/tcga.wgbs.33/exp.vs.methy.cor.RData")
#############################
###########at single cell level
#######
#singleM=read.table("~/nas/Xiaoqiang/opti.data/scRNA/GSE97693.CRC/methy.CRC01.pooled.all.NC.cells/CB.summarized.subfamily.data/cor.matrix.sh.241.cells.csv",header = T,row.names = 1,sep = "\t")
#dim(singleM)
load("~/nas/Xiaoqiang/opti.data/scRNA/GSE97693.CRC/methy.CRC01.pooled.all.NC.cells/CB.summarized.subfamily.data/data.for.subfamily.Me.normalization.RData")
singleM=subfamily.M.non.norm
singleTE=read.table("~/nas/Xiaoqiang/opti.data/scRNA/GSE97693.CRC/salmon.CRC01.pooled.NC/proced.data/outplot/subfamily.TE.expMatrix.logged.can.be.properly.exp.estimated.1194.384.final.txt",
                    header = T,row.names = 1,sep = "\t")
cell.anno.385.texp=read.table("~/nas/Xiaoqiang/opti.data/scRNA/GSE97693.CRC/methy.CRC01.pooled.all.NC.cells/cell.anno.385.texp.txt",header = T,row.names = 1,sep = "\t")
head(cell.anno.385.texp)
head(cell.anno.486.m)
dim(cell.anno.385.texp)
dim(cell.anno.486.m)
singleTE=as.data.frame(t(singleTE))
singleM=as.data.frame(t(singleM))
#
dim(singleTE)
singleTE[1:4,1:4]
dim(singleM)
singleM[1:4,1:4]
#"AluSq4","HERV1_LTRd","LTR21B","MER57F","MER65C","MER92-int","SVA_C","SVA_F","Tigger12A"
#####
singleidsh=merge(cell.anno.385.texp, cell.anno.486.m, by="full.id",all=FALSE)
head(singleidsh)
dim(singleidsh)
######
length(intersect(rownames(singleTE), singleidsh$SRR.id))
length(intersect(rownames(singleM), singleidsh$GSM.id.y))
####
singleTE.sh=singleTE[rownames(singleTE) %in% singleidsh$SRR.id,]
rownames(singleidsh)=singleidsh$SRR.id
singleTE.sh=singleTE.sh[rownames(singleidsh),]
rownames(singleTE.sh)=singleidsh$full.id

singleM.sh=singleM[rownames(singleM) %in% singleidsh$GSM.id.y,]
rownames(singleidsh)=singleidsh$GSM.id.y
singleM.sh=singleM.sh[rownames(singleidsh),]
rownames(singleM.sh)=singleidsh$full.id

singleTE.sh[1:4,1:4]
singleM.sh[1:4,1:4]
####
singleTE.sh.9s=singleTE.sh[, colnames(singleTE.sh) %in% rownames(cndi.rep)]
colnames(singleTE.sh.9s)=paste("exp",colnames(singleTE.sh.9s),sep = ".")
singleM.sh.9s=singleM.sh[, colnames(singleM.sh) %in% rownames(cndi.rep)]
colnames(singleM.sh.9s)=paste("methy",colnames(singleM.sh.9s),sep = ".")
head(singleM.sh.9s)
head(singleTE.sh.9s)
########
stat.single.cor=cbind(singleTE.sh.9s,singleM.sh.9s)

pdf("4.model/kaps/methylation/scRNA/exp.vs.methy.cor.scRNA.pdf",width = 10,height = 10)
ggcorr(stat.single.cor,nbreaks = 4,geom = "circle",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
ggcorr(stat.single.cor,nbreaks = 4,geom = "text",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
ggcorr(stat.single.cor,nbreaks = 4,method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,label = T,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
dev.off()

stat.single.cor.tumor=stat.single.cor
stat.single.cor.tumor$control=substr(rownames(stat.single.cor.tumor),7,9)
stat.single.cor.tumor$pat=substr(rownames(stat.single.cor.tumor),1,5)
stat.single.cor.tumor.pt=stat.single.cor.tumor[grepl("PT",stat.single.cor.tumor$control),]
head(stat.single.cor.tumor.pt)
table(stat.single.cor.tumor.pt$control)
#
pdf("4.model/kaps/methylation/scRNA/exp.vs.methy.cor.scRNA.pt.site.pdf",width = 10,height = 10)
ggcorr(stat.single.cor.tumor.pt[,1:18],nbreaks = 4,geom = "circle",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
ggcorr(stat.single.cor.tumor.pt[,1:18],nbreaks = 4,geom = "text",method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
ggcorr(stat.single.cor.tumor.pt[,1:18],nbreaks = 4,method = c("pairwise", "spearman"),name = "spearman coefs",size=2,
       min_size=4,max_size = 8,label = T,
       low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
dev.off()
#
save(stat.single.cor, stat.single.cor.tumor, stat.single.cor.tumor.pt, file="4.model/kaps/methylation/scRNA/exp.vs.methy.cor.scRNA.pt.site.RData")
########################
########################
##############pan cancer wgbs samples correlation of TE and immune
#########expM
gexp.pancancer=read.table("~/nas/Xiaoqiang/opti.data/TCGA.data/TCGA.pancancer/mRNA.pan/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.txt",header = T,sep = "\t")
gexp.pancancer=gexp.pancancer[!(duplicated(gexp.pancancer$sample)),]
rownames(gexp.pancancer)=gexp.pancancer$sample
gexp.pancancer=gexp.pancancer[,-1]
gexp.pancancer=as.data.frame(t(gexp.pancancer))
gexp.pancancer[1:4,1:4]
rownames(gexp.pancancer)=gsub("[.]","-",rownames(gexp.pancancer))
########global Methylation pan cancer
####
globalM.wgbs=read.table("~/nas/Xiaoqiang/opti.data/TCGA.data/pooled.samples.with.WGBS/output/CB.summarized.data/global.DNA.elevel.47.samples.txt",header = T,row.names = 1,sep = "\t")
head(globalM.wgbs)

rownames(globalM.wgbs)=substr(rownames(globalM.wgbs),1,15)
####
length(intersect(rownames(gexp.pancancer), rownames(globalM.wgbs)))
idwgbsM=intersect(rownames(gexp.pancancer), rownames(globalM.wgbs))
#########
stat.globalM.immune=cbind(globalM.wgbs[idwgbsM,], gexp.pancancer[idwgbsM,])
stat.globalM.immune[1:3,1:15]
dim(stat.globalM.immune)
aa=wgbs9g
rownames(aa)=substr(rownames(aa),1,15)
aa$mean.exp=rowMeans(aa)
length(intersect(rownames(aa), rownames(stat.globalM.immune)))
idr=intersect(rownames(aa), rownames(stat.globalM.immune))
stat.globalM.immune.TEexp=cbind(aa[idr,], stat.globalM.immune[idr,])
stat.globalM.immune.TEexp[1:4,1:20]
######
pdf("4.model/kaps/methylation/tcga.wgbs.33/cor.globalM.level.vs.immune.pdf")
ggscatter(subset(stat.globalM.immune), x = "CD8A", y = "av.beta.not.filtered", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "CD8A", ylab = "Global.methylation.level of WGBS (n=36)")+
  ggtitle("CD8A vs Global.methylation.level")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
######
save(stat.globalM.immune.TEexp,stat.globalM.immune,file = "4.model/kaps/methylation/tcga.wgbs.33/cor.globalM.level.vs.immune.RData" )


