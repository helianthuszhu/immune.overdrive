library(maftools)
#path to TCGA CRC MAF file
laml = read.maf(maf = "4.model/kaps/mutation/snv/data_mutations_extended.txt")
#
getFields(laml)
getClinicalData(laml) 
getSampleSummary(laml)
write.mafSummary(maf = laml, basename = '4.model/kaps/mutation/tcga.maf/maf.CRC.528')
#write.csv(kaps.td,"4.model/kaps/mutation/tcga.maf/kaps.td.csv")
#get the information
cm1=read.table("4.model/kaps/mutation/tcga.maf/maf.CRC.528_sampleSummary.txt",header = T,sep = "\t")
#cm1$id=substr(cm1$Tumor_Sample_Barcode,1,15)
#cm1=cm1[order(cm1$id,cm1$total, decreasing = T),]
#cm1.agg=cm1[!(duplicated(cm1$id)),]
cm1$id=cm1$Tumor_Sample_Barcode
rownames(cm1)=cm1$Tumor_Sample_Barcode
head(cm1)
dim(cm1)
#
cm2=kaps.td
cm2$id=rownames(cm2)
head(cm2)
mdfid=intersect(rownames(cm1), rownames(cm2))
clin.maf=cbind(cm1[mdfid,], cm2[mdfid,])
#
table(clin.maf$kaps.group)
dim(clin.maf)
########subset mafs
kaps.samples=rownames(clin.maf)
kaps.maf = subsetMaf(maf = laml, tsb = kaps.samples, mafObj = TRUE)
getClinicalData(kaps.maf)
dim(getClinicalData(kaps.maf))
kaps.maf
write.mafSummary(maf = kaps.maf, basename = '4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/maf.with.kaps.annoted.492s')
write.table(clin.maf,"4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/clin.maf.with.kaps.group.492s.txt",quote = F,sep = "\t")
#########
kaps.maf.clean=read.maf(maf = "4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/maf.with.kaps.annoted.492s_maftools.maf",
                        clinicalData = clin.maf
                          )
maf.set1=subsetMaf(maf = kaps.maf.clean, tsb = subset(clin.maf, kaps.group=="set1")$Tumor_Sample_Barcode, mafObj = TRUE)
maf.set2=subsetMaf(maf = kaps.maf.clean, tsb = subset(clin.maf, kaps.group=="set2")$Tumor_Sample_Barcode, mafObj = TRUE)
maf.set3=subsetMaf(maf = kaps.maf.clean, tsb = subset(clin.maf, kaps.group=="set3")$Tumor_Sample_Barcode, mafObj = TRUE)
maf.set4=subsetMaf(maf = kaps.maf.clean, tsb = subset(clin.maf, kaps.group=="set4")$Tumor_Sample_Barcode, mafObj = TRUE)
#####
#exclusive/co-occurance event analysis on top 10 mutated genes. 
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/set1.exclusive.co.occurance.gene.pdf")
print(somaticInteractions(maf = maf.set1, top = 25, pvalue = c(0.05, 0.1),fontSize = 0.6)
      )
dev.off()
#

laml.sig.set4= oncodrive(maf = maf.set4,  minMut = 5, pvalMethod = 'zscore')
laml.sig.set3= oncodrive(maf = maf.set3,  minMut = 5, pvalMethod = 'zscore')
laml.sig.set2= oncodrive(maf = maf.set2,  minMut = 5, pvalMethod = 'zscore')
laml.sig.set1= oncodrive(maf = maf.set1,  minMut = 5, pvalMethod = 'zscore')
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/set1.driver.pdf")
print(plotOncodrive(res = laml.sig.set1, fdrCutOff = 0.1, useFraction = TRUE)
      )
dev.off()
#
#Survival analysis based on grouping of DNMT3A mutation status
cgene=c("APC","TTN","TP53","KRAS","FLG","MUC16","SYNE1","RYR2","FAT4","RYR1")
idx.time=c("DSS.time","OS.time","DFI.time","PDF.time")
idx.status=c("DSS","OS","DFI","PDF")
for (i in 1:4) {
  for (j in 1: length(cgene)) {
    pdf(paste0("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/sur.gene.based/",idx.status[i]),"-",cgene[j],".pdf")
    print(mafSurvival(maf =maf.set1 , genes = cgene[j], time = idx.time[i], Status = idx.status[i], isTCGA = TRUE)
    )
    dev.off()
  }
}

#########
#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt.set34 <- mafCompare(m1 = maf.set3, m2 = maf.set4, m1Name = 'set3', m2Name = 'set4', minMut = 5)
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/set3.vs.set4.forest.pdf")
forestPlot(mafCompareRes = pt.vs.rt.set34, pVal = 0.05, color = c("#E69F00","#d73027"), geneFontSize = 0.8)
dev.off()

pt.vs.rt14 <- mafCompare(m1 = maf.set1, m2 = maf.set4, m1Name = 'set1', m2Name = 'set4', minMut = 5)
print(pt.vs.rt14)
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/set1.vs.set4.forest.pdf")
forestPlot(mafCompareRes = pt.vs.rt14, pVal = 0.01, color = c("#E69F00","#d73027"), geneFontSize = 0.3)
dev.off()
#########
#genes = c("PML", "RARA", "RUNX1", "ARID1B", "FLT3")
#coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', genes = genes, removeNonMutated = TRUE)
#
#Clinical enrichment analysis
kaps.ce = clinicalEnrichment(maf = kaps.maf.clean, clinicalFeature = 'kaps.group')
#Results are returned as a list. Significant associations p-value < 0.05
kaps.ce$groupwise_comparision[p_value < 0.05]
######drugable
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/drugable/set4.drug.pdf")
drugInteractions(maf = maf.set4, fontSize = 0.75)
dev.off()
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/drugable/set3.drug.pdf")
drugInteractions(maf = maf.set3, fontSize = 0.75)
dev.off()
##
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/drugable/set2.drug.pdf")
drugInteractions(maf = maf.set2, fontSize = 0.75)
dev.off()
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/drugable/set1.drug.pdf")
drugInteractions(maf = maf.set1, fontSize = 0.75)
dev.off()
###Oncogenic Signaling Pathways
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/onco.paths/set4.oncogenic.pathways.pdf")
OncogenicPathways(maf = maf.set4)
PlotOncogenicPathways(maf = maf.set4, pathways = "RTK-RAS")
dev.off()
####
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/onco.paths/set3.oncogenic.pathways.pdf")
OncogenicPathways(maf = maf.set3)
PlotOncogenicPathways(maf = maf.set3, pathways = "RTK-RAS")
dev.off()
#
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/onco.paths/set2.oncogenic.pathways.pdf")
OncogenicPathways(maf = maf.set2)
PlotOncogenicPathways(maf = maf.set2, pathways = "RTK-RAS")
dev.off()
#
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/onco.paths/set1.oncogenic.pathways.pdf")
OncogenicPathways(maf = maf.set1)
PlotOncogenicPathways(maf = maf.set1, pathways = "RTK-RAS")
dev.off()
####
#Heterogeneity in samples
kaps.het = inferHeterogeneity(maf =kaps.maf.clean,top = nrow(clin.maf))
kaps.het.agg= kaps.het$clusterData[,c(5,8:9)] %>% group_by(Tumor_Sample_Barcode) %>% summarise_all(mean)
kaps.het.agg=as.data.frame(kaps.het.agg)
rownames(kaps.het.agg)=kaps.het.agg$Tumor_Sample_Barcode
head(kaps.het.agg)
rownames(clin.maf)=clin.maf$Tumor_Sample_Barcode
stat.kaps.het=cbind(clin.maf, kaps.het.agg[rownames(clin.maf), 2:3])
#
ggviolin(stat.kaps.het[,-12],x = "kaps.group", y ="MATH" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("MATH")+ggtitle("MATH")
######APOBEC
#Requires BSgenome object
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
#
laml.tnm = trinucleotideMatrix(maf = kaps.maf.clean, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = laml.tnm, maf = kaps.maf.clean, pVal = 0.2)
head(laml.tnm $APOBEC_scores)
apscore=as.data.frame(laml.tnm$APOBEC_scores)
rownames(apscore)=apscore$Tumor_Sample_Barcode
head(apscore)
stat.apscore=cbind(apscore[rownames(clin.maf),]$APOBEC_Enrichment, clin.maf)
colnames(stat.apscore)[1]="APOBEC_Enrichment"
head(stat.apscore)
#
ggviolin(stat.apscore[,-13],x = "kaps.group", y ="APOBEC_Enrichment" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("APOBEC_Enrichment")+ggtitle("APOBEC_Enrichment")
#########
###
library('NMF')
laml.tnm.set4 = trinucleotideMatrix(maf = maf.set4, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
laml.sign.set4 = estimateSignatures(mat = laml.tnm.set4, nTry = 6)
maftools::plotCophenetic(res = laml.sign.set4)
laml.sig = extractSignatures(mat = laml.tnm, n = 2)
#Compate against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
#Compate against updated version3 60 signatures 
#laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")
library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
maftools::plotSignatures(nmfRes = laml.sig, title_size = 0.8)
#####
########mutation signature
####
Poulos.data=read.table("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/mutation.signature/Poulos_et_al_FigS2_data.txt",header = T,sep = "\t")
rownames(Poulos.data)=substr(Poulos.data$TCGA.ID, 1,15)
head(Poulos.data)
pouid=intersect(rownames(kaps.td), rownames(Poulos.data))
Poulos.data.sel=Poulos.data[pouid,-c(1:2)]
head(Poulos.data.sel)
Poulos.data.anno=kaps.td[pouid,]
head(Poulos.data.anno)
####

col_ha.top.poul=columnAnnotation(z.of.mean.exp.score=anno_lines(Poulos.data.anno$z.of.mean.exp),
                             TE.cluster = Poulos.data.anno$TE.cluster,
                             TE.cluster.agg=Poulos.data.anno$TE.cluster.agg,
                             MSI.status.bin=Poulos.data.anno$MSI.status.bin,
                             z.of.mean.exp=Poulos.data.anno$z.of.mean.exp,
                             kaps.group=Poulos.data.anno$kaps.group,
                             col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                      TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                      MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)
puht=Heatmap(t(Poulos.data.sel), name = "30.mutation.signatures", na_col = "black",
              #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
              col = rev(viridis(2)),border = F,
              show_column_names = F,show_row_names = T,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              top_annotation = col_ha.top.poul, 
              #right_annotation = row_ha_right#,
              #show_row_dend = T,show_column_dend = T,
              #row_names_side = "left",
              #left_annotation = row_ha.left
)
generate.PDF <- function(fig) {
  pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/mutation.signature/signature.ht.pdf",height = 7,width =10)
  print(puht)
  dev.off()
}
generate.PDF(fig)
################
polug=cbind(Poulos.data.anno,Poulos.data.sel )
polug$kaps.group=factor(polug$kaps.group,levels = c("set4", "set3","set2","set1"))
polugd1=ggviolin(polug,x = "kaps.group", y ="Signature.1" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Signature.1")+ggtitle("Signature.1")
polugd2=ggviolin(polug,x = "kaps.group", y ="Signature.6" , fill = "kaps.group",alpha = 1,size = 0.3,
                #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
                add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("Signature.6")+ggtitle("Signature.6")
generate.PDF <- function(fig) {
  pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/mutation.signature/signature.1.pdf",height = 4,width =3)
  print(polugd1)
  print(polugd2)
  dev.off()
}
generate.PDF(fig)
###################
###################input cibioportal result
musiggene=read.table("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/res.from.cbioportal/table.snv.among.kaps.group.tsv",header = T,sep = "\t")
musiggene=musiggene[order(musiggene$X.D..set4,decreasing = T),]

musiggene$set4.new=gsub("[(%)]","",musiggene$X.D..set4)
musiggene=separate(musiggene, col = set4.new,into = c("set4.abs","set4.prop"),sep = " ")

musiggene$set3.new=gsub("[(%)]","",musiggene$X.C..set3)
musiggene=separate(musiggene, col = set3.new,into = c("set3.abs","set3.prop"),sep = " ")

musiggene$set2.new=gsub("[(%)]","",musiggene$X.B..set2)
musiggene=separate(musiggene, col = set2.new,into = c("set2.abs","set2.prop"),sep = " ")

musiggene$set1.new=gsub("[(%)]","",musiggene$X.A..set1)
musiggene=separate(musiggene, col = set1.new,into = c("set1.abs","set1.prop"),sep = " ")
######
musiggene.con <- as.data.frame(lapply(as.data.frame(musiggene[,-c(1:9)]), function(x) as.numeric(as.character(x))))
musiggene=cbind(musiggene[,1:9], musiggene.con)
head(musiggene)
####
musiggene.set4.sel=subset(musiggene, p.Value <0.05 & Most.enriched.in =="(D) set4")
musiggene.set4.sel=musiggene.set4.sel[order(musiggene.set4.sel$set4.prop,decreasing = T),]
head(musiggene.set4.sel)
####
musiggene.set3.sel=subset(musiggene, p.Value <0.05& Most.enriched.in =="(C) set3")
musiggene.set3.sel=musiggene.set3.sel[order(musiggene.set3.sel$set3.prop,decreasing = T),]
head(musiggene.set3.sel)
#
musiggene.set2.sel=subset(musiggene, p.Value <0.05& Most.enriched.in =="(B) set2")
musiggene.set2.sel=musiggene.set2.sel[order(musiggene.set2.sel$set2.prop,decreasing = T),]
head(musiggene.set2.sel)
#
musiggene.set1.sel=subset(musiggene, Most.enriched.in =="(A) set1")
musiggene.set1.sel=musiggene.set1.sel[order(musiggene.set1.sel$set1.prop,decreasing = T),]
head(musiggene.set1.sel)
#
#One can use any colors, here in this example color palette from RColorBrewer package is used
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

#Color coding
fabcolors = list(kaps.group=c("set4"="#d73027","set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00"),
                      MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e")
)
#set4
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/set4.sig.top20.pdf")
print(oncoplot(
  maf = maf.set4,genes = musiggene.set4.sel$Gene[1:20],colors = vc_cols,bgCol = "#d9d9d9",
  annotationColor = fabcolors,removeNonMutated = F,fontSize = 0.5,
  #draw_titv = TRUE,
  clinicalFeatures = c('kaps.group', 'MSI.status.bin'),
  #additionalFeature = c("Tumor_Seq_Allele2", "C"),
  sortByAnnotation = TRUE,
  #mutsig = laml.mutsig,
  #exprsTbl = exprs_tbl,
  logColBar = TRUE
)
)
dev.off()
#set3
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/set3.sig.top20.pdf")
print(oncoplot(
  maf = maf.set3,genes = musiggene.set3.sel$Gene[1:21],colors = vc_cols,bgCol = "#d9d9d9",genesToIgnore = "NEXMIF",
  annotationColor = fabcolors,removeNonMutated = F,fontSize = 0.5,
  #draw_titv = TRUE,
  clinicalFeatures = c('kaps.group', 'MSI.status.bin'),
  #additionalFeature = c("Tumor_Seq_Allele2", "C"),
  sortByAnnotation = TRUE,
  #mutsig = laml.mutsig,
  #exprsTbl = exprs_tbl,
  logColBar = TRUE
)
)
dev.off()
####set2
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/set2.sig.top20.pdf")
print(oncoplot(
  maf = maf.set2,genes = musiggene.set2.sel$Gene[1:21],colors = vc_cols,bgCol = "#d9d9d9",genesToIgnore = "ICE1",
  annotationColor = fabcolors,removeNonMutated = F,fontSize = 0.5,
  #draw_titv = TRUE,
  clinicalFeatures = c('kaps.group', 'MSI.status.bin'),
  #additionalFeature = c("Tumor_Seq_Allele2", "C"),
  sortByAnnotation = TRUE,
  #mutsig = laml.mutsig,
  #exprsTbl = exprs_tbl,
  logColBar = TRUE
)
)
dev.off()
##set1
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/set1.sig.top20.pdf")
print(oncoplot(
  maf = maf.set1,genes = musiggene.set1.sel$Gene[1:20],colors = vc_cols,bgCol = "#d9d9d9",
  annotationColor = fabcolors,removeNonMutated = F,fontSize = 0.5,
  #draw_titv = TRUE,
  clinicalFeatures = c('kaps.group', 'MSI.status.bin'),
  #additionalFeature = c("Tumor_Seq_Allele2", "C"),
  sortByAnnotation = TRUE,
  #mutsig = laml.mutsig,
  #exprsTbl = exprs_tbl,
  logColBar = TRUE
)
)
dev.off()
########
#####gradually decreased
musiggene.grad=subset(musiggene, set4.prop <= set3.prop)
musiggene.grad=subset(musiggene.grad, set3.prop <= set2.prop)
musiggene.grad=subset(musiggene.grad, set2.prop <= set1.prop)
musiggene.grad=musiggene.grad[order(musiggene.grad$set1.prop,decreasing = T),]
write.csv(musiggene.grad, "4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/gradually.decrased.mut.genes.csv")
head(musiggene.grad)
#
fabcolors.total = list(kaps.group=c("set4"="#d73027","set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00"),
                 MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e")
)

pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/gradually.decrased.mut.genes.top20.pdf")
print(oncoplot(
  maf = kaps.maf.clean,genes = musiggene.grad$Gene[1:20],colors = vc_cols,bgCol = "#d9d9d9",
  annotationColor = fabcolors,removeNonMutated = F,fontSize = 0.5,
  #draw_titv = TRUE,
  clinicalFeatures = c('kaps.group', 'MSI.status.bin'),
  #additionalFeature = c("Tumor_Seq_Allele2", "C"),
  sortByAnnotation = TRUE,
  #mutsig = laml.mutsig,
  #exprsTbl = exprs_tbl,
  logColBar = TRUE
)
)
dev.off()
#### gradually increase
musiggene.grad2=subset(musiggene, set1.prop <= set2.prop)
musiggene.grad2=subset(musiggene.grad2, set2.prop <= set3.prop)
musiggene.grad2=subset(musiggene.grad2, set3.prop <= set4.prop)
musiggene.grad2=musiggene.grad2[order(musiggene.grad2$set4.prop,decreasing = T),]
write.csv(musiggene.grad2, "4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/gradually.incrased.mut.genes.csv")

head(musiggene.grad2)
#
pdf("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/oncoplot/gradually.increased.mut.genes.top20.pdf")
print(oncoplot(
  maf = kaps.maf.clean,genes = musiggene.grad2$Gene[1:20],colors = vc_cols,bgCol = "#d9d9d9",
  annotationColor = fabcolors,removeNonMutated = F,fontSize = 0.5,
  #draw_titv = TRUE,
  clinicalFeatures = c('kaps.group', 'MSI.status.bin'),
  #additionalFeature = c("Tumor_Seq_Allele2", "C"),
  sortByAnnotation = TRUE,
  #mutsig = laml.mutsig,
  #exprsTbl = exprs_tbl,
  logColBar = TRUE
)
)
dev.off()

#####set4 vs set3
dim(getClinicalData(kaps.maf.clean))
