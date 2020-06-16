########CCLE data
########https://www.synapse.org/#!Synapse:syn18689330
########https://www.synapse.org/#!Synapse:syn18685536/files/
###########
###########
#CCLE gene expression
cclegeneexp=read.table("4.model/kaps/ccleCRC/CCLE_normalized_expression.Syn19689330.txt",header = T,row.names = 1)
cclegeneexp=as.data.frame(t(cclegeneexp))
dim(cclegeneexp)
colnames(cclegeneexp)
#CCLE te expression
ccleteexp.CRC=readRDS("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.plot/consensusClusterPlus/ccle/RE.output/RE_intergenic_2_counts_normalized.RDS")
ccleteexp.CRC=exprs(ccleteexp.CRC)
ccleteexp.CRC=as.data.frame(ccleteexp.CRC)
ccleteexp.CRC=ccleteexp.CRC[rownames(cdrep.sel),]
ccleteexp.CRC=as.data.frame(t(ccleteexp.CRC))
ccleteexp.CRC[1:4,1:4]
dim(ccleteexp.CRC)
#####
ccleCRC.label=read.table("~/nas/Xiaoqiang/opti.data/CCLE/CRC.ccle.id.list",header = F,sep="\t")
colnames(ccleCRC.label)=c("SRRid","cellName")
ccleCRC.label$cellName2=gsub("[-]","",ccleCRC.label$cellName)
ccleCRC.label$cellName2=gsub(" ","",ccleCRC.label$cellName2)
ccleCRC.label$cellName2=gsub("[.]","",ccleCRC.label$cellName2)
#ccleCRC.label$cellName2=gsub("Hs","HS",ccleCRC.label$cellName2)
#
ccleCRC.label$cellName2=toupper(ccleCRC.label$cellName2)
rownames(ccleCRC.label)=ccleCRC.label$SRRid
head(ccleCRC.label)
#write.csv(ccleCRC.label,"4.model/kaps/ccleCRC/cellid.csv")
###ccle CRC msi
ccleCRC.msi=read.table("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.plot/consensusClusterPlus/ccle/cell.crc.cell.id.txt",header=T,sep = "\t")
ccleCRC.msi=ccleCRC.msi[,11:14]
ccleCRC.msi$cell.abb=gsub("[-]","",ccleCRC.msi$cell.abb)
ccleCRC.msi$cell.abb=gsub(" ","",ccleCRC.msi$cell.abb)
ccleCRC.msi$cell.abb=gsub("[.]","",ccleCRC.msi$cell.abb)
ccleCRC.msi$cell.abb=toupper(ccleCRC.msi$cell.abb)
head(ccleCRC.msi)
dim(ccleCRC.msi)
######
######change TEexpM SRR ID into cellName
dim(ccleteexp.CRC)
length(intersect(rownames(ccleteexp.CRC), ccleCRC.label$SRRid))
ccleCRC.label=ccleCRC.label[rownames(ccleteexp.CRC),]
head(ccleCRC.label)
rownames(ccleteexp.CRC)=ccleCRC.label$cellName2
rownames(ccleCRC.label)=ccleCRC.label$cellName2
###### overlapped samples Geneexpression and TEexpression
cellid1=intersect(rownames(ccleteexp.CRC), rownames(cclegeneexp))
###
ccleteexp.CRC.overlapped=ccleteexp.CRC[cellid1,]
cclegeneexp.CRC.overlapped=cclegeneexp[cellid1,]
ccleteexp.CRC.overlapped[1:4,1:4]
cclegeneexp.CRC.overlapped[1:4,1:4]
ccleCRC.label.overlapped=ccleCRC.label[cellid1,]
ccleCRC.label.overlapped$cell.abb=ccleCRC.label.overlapped$cellName2
rownames(ccleCRC.label.overlapped)=ccleCRC.label.overlapped$cell.abb
head(ccleCRC.label.overlapped)
dim(ccleteexp.CRC.overlapped)
ccleCRC.label.overlapped=merge(ccleCRC.label.overlapped, ccleCRC.msi, by="cell.abb",all.x=TRUE)
#########
save(ccleCRC.label.overlapped, cclegeneexp.CRC.overlapped,ccleteexp.CRC.overlapped,file="4.model/kaps/ccleCRC/cellCRC.label.geneexp.teexp.RData")
####
stat.ccleCRC=cbind(ccleteexp.CRC.overlapped[, colnames(ccleteexp.CRC.overlapped) %in% rownames(cndi.rep)],ccleCRC.label.overlapped[rownames(ccleteexp.CRC.overlapped),])
stat.ccleCRC$mean.exp=rowMeans(stat.ccleCRC[,1:9],na.rm = T)
stat.ccleCRC$z.of.mean.exp=(stat.ccleCRC$mean.exp - mean(stat.ccleCRC$mean.exp))/sd(stat.ccleCRC$mean.exp)
stat.ccleCRC=stat.ccleCRC[order(stat.ccleCRC$z.of.mean.exp,decreasing = T),]
#
stat.ccleCRC[["CCLE.MSI.call"]][is.na(stat.ccleCRC[["CCLE.MSI.call"]])] <- "not_available"
head(stat.ccleCRC)
########
summary(stat.ccleCRC$z.of.mean.exp)
#
########heatmap of CCLE
######TE annotation
TE.ann.data=cndi.rep[rownames(cndi.rep),c(2,3)]

TE.an.ha=rowAnnotation(repClass=TE.ann.data$repClass,
                       repFamily=TE.ann.data$repFamily,
                       col=list(repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")),show_annotation_name = FALSE
)
#
col_ha.ccle=columnAnnotation(TE.score=anno_lines(stat.ccleCRC$z.of.mean.exp),
                             z.of.mean.exp=stat.ccleCRC$z.of.mean.exp,
                             msi.status=stat.ccleCRC$CCLE.MSI.call,
                                          col=list(z.of.mean.exp=colorRamp2(c(-6,-2,0,1,2), 
                                                                            c("#ffffcc","#d9f0a3","#addd8e","#78c679","#006837")),
                                                   msi.status=c("inferred-MSI"="#e7298a","inferred-MSS"="#66a61e","not_available"="white")
                                          ),show_annotation_name = TRUE)
#########
aadraw.ccle=as.data.frame(stat.ccleCRC[,1:9])
aadraw.ccle=as.data.frame(scale(aadraw.ccle))
aadraw.ccle[aadraw.ccle< -3] <- -3
aadraw.ccle[aadraw.ccle> 3] <- 3

min_cor = min(as.vector(aadraw.ccle))
max_cor = max(as.vector(aadraw.ccle))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#00F5FF", "white","#FF3E96"))(50))
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
ht.ccle = Heatmap(t(aadraw.ccle),#col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256)
             col=col.pal_cor,
             top_annotation = col_ha.ccle,
             #bottom_annotation = ha,
             show_column_names = F,left_annotation = TE.an.ha,
             cluster_columns = F,cluster_rows = T, 
             column_title = "clinical.comparison among TE clusters")


pdf("4.model/kaps/ccleCRC/ht.9TEs.ccle.pdf",height = 5,width = 10)
draw(ht.ccle, padding = unit(c(22, 80, 20,60), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     #,heatmap_legend_side = "right"
)
dev.off()
####################
####################select the TE inhibitor expression
##########
inhibition.markers=read.table("4.model/kaps/marker.TEs.inhibition/TEs.inhibition.markers.txt",header = T,sep = "\t")
inhibition.markers$symbol=gsub("IL8","CXCL8",inhibition.markers$symbol)
inhibition.markers$symbol=gsub("DAK","TKFC",inhibition.markers$symbol)
rownames(inhibition.markers)=inhibition.markers$symbol
head(inhibition.markers)
dim(inhibition.markers)
#
dim(cclegeneexp.CRC.overlapped)
length(intersect(colnames(cclegeneexp.CRC.overlapped), rownames(inhibition.markers)))  #108s
#write.csv(as.data.frame(t(cclegeneexp.CRC.overlapped[1:2,])),"4.model/kaps/ccleCRC/tmp.geneid.csv")
#write.csv(inhibition.markers, "4.model/kaps/ccleCRC/tmp.marker.csv")
######
ccleinhi.TEmarkersM=cclegeneexp.CRC.overlapped[, colnames(cclegeneexp.CRC.overlapped) %in% rownames(inhibition.markers) ]
ccleinhi.TEmarkersM[1:4,1:4]
######
stat.ccle.inhi.correlation=cbind(stat.ccleCRC$z.of.mean.exp, ccleinhi.TEmarkersM[rownames(stat.ccleCRC),])
colnames(stat.ccle.inhi.correlation)[1]="z.of.mean.exp"
dim(stat.ccle.inhi.correlation)
stat.ccle.inhi.correlation[1:4,1:4]
######calculate correlation 
datalist=list()
for (j in 2:ncol(stat.ccle.inhi.correlation)) {
  res <- cor.test(stat.ccle.inhi.correlation[,1], stat.ccle.inhi.correlation[,j],  method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(stat.ccle.inhi.correlation)[j]
  datalist[[j]] <- res
}
res.cor.ccle.inhi.marker=do.call(rbind, datalist)
res.cor.ccle.inhi.marker=cbind(res.cor.ccle.inhi.marker, inhibition.markers[rownames(res.cor.ccle.inhi.marker),])
head(res.cor.ccle.inhi.marker)
summary(res.cor.ccle.inhi.marker$cor)
write.csv(res.cor.ccle.inhi.marker, "4.model/kaps/ccleCRC/res.cor.ccle.inhi.marker.csv")
######################
######################
save(stat.ccleCRC, inhibition.markers,ccleinhi.TEmarkersM,res.cor.ccle.inhi.marker,file="4.model/kaps/ccleCRC/res.cor.ccle.inhi.marker.RData")
write.csv(stat.ccleCRC, "4.model/kaps/ccleCRC/TE.9s.expM.ccle.CRC.csv")
