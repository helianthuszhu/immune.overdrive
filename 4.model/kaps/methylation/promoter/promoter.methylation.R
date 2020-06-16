########compare the methylation level of genes those correlated with TE.mean.exp
#######
mMatrix.filter[1:4,1:4]
head(meprob)
#
meprob.genemodif=meprob %>% separate(gene,c("gene1","gene2"),",")
meprob.genemodif$gene2=ifelse(is.na(meprob.genemodif$gene2),  paste(meprob.genemodif$gene1), paste(meprob.genemodif$gene2))
head(meprob.genemodif)
#eg = bitr(meprob.genemodif$gene1, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
###############
####get the promoter region of gene
library(RnBeads)
library(RnBeads.hg19)
promoter19anno=rnb.annotation2data.frame(rnb.get.annotation("promoters",assembly = "hg19"))
promoter19anno=separate(promoter19anno, col = "symbol",into = c("gene1","gene2"),sep = ";")
head(promoter19anno)
###########
########select out probes in the selected gene promoter region
length(intersect(subset(TE.mean.corM.with.genes, set=="set1")$symbol, promoter19anno$gene1))
table(TE.mean.corM.with.genes$set)
corgenesel=intersect(subset(TE.mean.corM.with.genes, !(set=="set3"))$symbol, promoter19anno$gene1)
#
promoter19anno.sel=promoter19anno[promoter19anno$gene1 %in% corgenesel, ]
dim(promoter19anno.sel)
head(promoter19anno.sel)
#
bed1=with(promoter19anno.sel, GRanges(Chromosome, IRanges(Start+1, End), gene1, entrezID, CpG ,strand = Strand))
head(bed1)
#length(intersect(rownames(mMatrix.filter), meprob$id))
bed2=with(df2, GRanges(chr, IRanges(start+1, end), proName,gene,strand = NULL))
head(bed2)
length(bed2)
# now find the overlaps
df3=as.data.frame(subsetByOverlaps(bed2, bed1))
head(df3)
dim(df3)
#
mMatrix.filter.sel=mMatrix.filter[rownames(mMatrix.filter) %in% df3$proName,]
mMatrix.filter.sel=as.data.frame(t(mMatrix.filter.sel))
rownames(mMatrix.filter.sel)=gsub("[.]","-",rownames(mMatrix.filter.sel))
mMatrix.filter.sel[1:4,1:3]
#
idset=intersect(rownames(kaps.td), rownames(mMatrix.filter.sel))
mMatrix.filter.sel=mMatrix.filter.sel[idset,]
mMatrix.filter.sel=as.data.frame(t(mMatrix.filter.sel))
meanno=kaps.td[idset,]
head(meanno)
save(meanno, mMatrix.filter.sel, promoter19anno,meprob.genemodif,df2,file="4.model/kaps/methylation/promoter/promoter.methylation.RData")
#
col_ha.top.me=columnAnnotation(z.of.mean.exp.score=anno_lines(meanno$z.of.mean.exp),
                             TE.cluster = meanno$TE.cluster,
                             TE.cluster.agg=meanno$TE.cluster.agg,
                             MSI.status.bin=meanno$MSI.status.bin,
                             z.of.mean.exp=meanno$z.of.mean.exp,
                             kaps.group=meanno$kaps.group,
                             col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                      TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                      MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)
#
ind = order(rowVars(as.matrix(mMatrix.filter.sel), na.rm = TRUE), decreasing = TRUE)[1:2000]

mMatrix.filter.sel2=mMatrix.filter.sel[ind,]

pdf("4.model/kaps/methylation/promoter/proM.pdf")
Heatmap(as.matrix(mMatrix.filter.sel), name = "methylation", 
        col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
        #col = rev(viridis(10)),border = F,
        show_column_names = F,show_row_names = F,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha.top.me#, 
        #right_annotation = row_ha.right,
        #show_row_dend = T,show_column_dend = T,
        #row_names_side = "left",
        #left_annotation = row_ha.left
)
dev.off()
###
#BiocManager::install("ChAMP")
#library(ChAMP)
