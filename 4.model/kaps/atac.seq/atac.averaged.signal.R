#########
atac.seq.probe=read.table("4.model/kaps/atac.seq/TCGA_ATAC_peak.all.probeMap",header = T,sep = "\t")
head(atac.seq.probe)
idd=separate(atac.seq.probe, "id",into = c("cancerName","loci"))
head(idd)
atac.seq.probe=cbind(idd[,1:2], atac.seq.probe)
head(atac.seq.probe)
table(atac.seq.probe$cancerName)
######
######
atac.seq.probe.CRC=subset(atac.seq.probe, cancerName=="COAD")
head(atac.seq.probe.CRC)
rownames(atac.seq.probe.CRC)=atac.seq.probe.CRC$id
######
atac.signal=read.table("4.model/kaps/atac.seq/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz",header = T,row.names = 1,sep = "\t")
atac.signal[1:4,1:4]
####
atac.signal.id=as.data.frame(colnames(atac.signal))
colnames(atac.signal.id)="id"
head(atac.signal.id)
atac.signal.id$idd=substr(atac.signal.id$id,1,15)
atac.signal.id$idd=gsub("[.]","-",atac.signal.id$idd)
rownames(atac.signal.id)=atac.signal.id$idd
#
#length(intersect(rownames(kaps.td), atac.signal.id$idd))
idsss=intersect(rownames(kaps.td), rownames(atac.signal.id))
head(kaps.td)
#####
atac.signal.id.CRC=cbind(atac.signal.id[idsss,], kaps.td[idsss, c(17,38)])
head(atac.signal.id.CRC)
table(atac.signal.id.CRC$kaps.group)
##########
atac.signal.CRC.matrix=atac.signal[rownames(atac.seq.probe.CRC), colnames(atac.signal) %in% atac.signal.id.CRC$id]
atac.signal.CRC.matrix[1:4,1:4]
dim(atac.signal.CRC.matrix)
###################
###################
########overlap the location of TE with signal probes
#
head(df1)
dim(df1)
length(unique(df1$subfamily))

head(cdrep)
#dim(cdrep.sel)
length(intersect(cdrep.sel$repName, unique(df1$subfamily)))
atac.ids=intersect(rownames(cdrep), unique(df1$subfamily))
cdrep.atac=cdrep[atac.ids,]
dim(cdrep.atac)
#write.csv(cdrep, "4.model/kaps/atac.seq/cdrep.csv")
#write.csv(as.data.frame(unique(df1$subfamily)), "4.model/kaps/atac.seq/df1.csv")
#
#te.list=cdrep.atac$repName
te.list=rownames(cndi.rep)
#
df2=atac.seq.probe.CRC
head(df2)

#
datalist=list()
for (k in 1:length(te.list)){
  sub.df1=as.data.frame(subset(df1,df1$subfamily==te.list[k]))
  # next, create IRanges objects
  bed1=with(sub.df1, GRanges(chr, IRanges(start+1, end), subfamily, class, family,strand = NULL))
  head(bed1)
  #length(intersect(rownames(mMatrix.filter), meprob$id))
  
  bed2=with(df2, GRanges(chrom, IRanges(chromStart+1, chromEnd), id,strand = strand))
  head(bed2)
  length(bed2)
  # now find the overlaps
  df3=as.data.frame(subsetByOverlaps(bed2, bed1))
  head(df3)
  
  
  
  calmmatrix=atac.signal.CRC.matrix[df3$id, ]
  calmmatrix[1:4,1:7]
  #
  #calmean=as.data.frame(colMeans(calmmatrix))
  #colnames(calmean)=te.list[k]
  calmean=cbind(TEname=rep(te.list[k], times=nrow(calmmatrix)),
                calmmatrix)
  

  datalist[[k]]=calmean
  
}
#atac.agg.signal=do.call(cbind, datalist)
#rownames(atac.agg.signal)=substr(rownames(atac.agg.signal),1,15)
#rownames(atac.agg.signal)=gsub("[.]","-",rownames(atac.agg.signal))
#
atac.agg.signal=do.call(rbind, datalist)
colnames(atac.agg.signal)=substr(colnames(atac.agg.signal),1,15)
colnames(atac.agg.signal)=gsub("[.]","-",colnames(atac.agg.signal))
atac.agg.signal[1:4,1:6]
dim(atac.agg.signal)
####
head(atac.signal.id.CRC)
dim(atac.signal.id.CRC)
length(intersect(atac.signal.id.CRC$idd, colnames(atac.agg.signal)))
####
#atac.stat=cbind(atac.signal.id.CRC, atac.agg.signal[rownames(atac.signal.id.CRC),])
#atac.stat[1:4,1:8]
#write.csv(atac.stat, "4.model/kaps/atac.seq/ht.atac.averaged.signal.csv")
####
####
#heatmap
#
#atac.drawmatrix=as.data.frame(t(atac.stat[,-c(1:4)]))
atac.drawmatrix=atac.agg.signal[,-1]
atac.drawmatrix[1:6,1:4]
#atac.drawmatrix=na.omit(atac.drawmatrix)
##
head(cdrep.atac)
#####
#row_ha_left.atac= rowAnnotation(
#  TE.class=cdrep.atac[rownames(atac.drawmatrix),]$repClass.new2,
#  col=list(TE.class=c(DNA="#1f78b4",LINE="#d95f02",LTR="#7570b3",Retroposon="#e7298a",Satellite="#66a61e",SINE="#e6ab02","other.repeats"="gray")
#  ),
#  show_annotation_name = FALSE)

#col_ha_top.atac= columnAnnotation(
#  TE.cluster=atac.stat$kaps.group,
#  col=list(TE.cluster=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )
#  ),
#  show_annotation_name = FALSE)
##
row_ha_left.atac= rowAnnotation(
  TE.class=atac.agg.signal$TEname,
  col=list(TE.class=c("AluSq4"="#1f78b4","HERV1_LTRd"="#d95f02","LTR21B"="#7570b3","MER57F"="#e7298a",
                      "MER65C"="#66a61e","MER92-int"="#e6ab02","SVA_C"="#e41a1c","SVA_F"="#984ea3","Tigger12A"="#ff7f00")
  ),
  show_annotation_name = FALSE)
col_ha_top.atac= columnAnnotation(
  TE.cluster=atac.signal.id.CRC$kaps.group,
  col=list(TE.cluster=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )
  ),
  show_annotation_name = FALSE)

hmht.atac=Heatmap(atac.drawmatrix, name = "signal", 
                 col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                 #col = rev(viridis(10)),
                 width = unit(5, "cm"),
                 height = unit(12, "cm"),
                 border = F,
                 #col=col.pal_cor,
                 show_column_names = F,show_row_names = F,
                 cluster_columns = F,
                 row_names_gp = gpar(fontsize = 5),
                 column_names_gp = gpar(fontsize = 5),
                 top_annotation = col_ha_top.atac,
                 #right_annotation = row_ha.right,
                 show_row_dend = F,show_column_dend = F,
                 #row_names_side = "left",
                 left_annotation = row_ha_left.atac,
                 column_title=paste0("average signal of promoter/enhancer inside TE location"),
                 column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/atac.seq/ht.atac.averaged.signal.9.tes.pdf",width = 8,height = 6)
draw(hmht.atac, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
##############
#######uprsteam and downstream location
##############
head(df1)
df1new=df1
df1new$start.leftup1kb=df1new$start-1000
df1new$start.leftdown1kb=df1new$start+1000
df1new$start.rightup1kb=df1new$end-1000
df1new$start.rightdown1kb=df1new$end+1000
head(df1new)
##############
#####left
datalist=list()
for (k in 1:length(te.list)){
  #left
  sub.df1left=as.data.frame(subset(df1new,df1new$subfamily==te.list[k]))   #change
  # next, create IRanges objects
  bed1.left=with(sub.df1left, GRanges(chr, IRanges(start.leftup1kb+1, start), subfamily, class, family,strand = NULL))
  head(bed1.left)
  #right
  sub.df1right=as.data.frame(subset(df1new,df1new$subfamily==te.list[k]))  #change
  # next, create IRanges objects
  bed1.right=with(sub.df1right, GRanges(chr, IRanges(end+1, start.rightdown1kb), subfamily, class, family,strand = NULL))
  head(bed1.right)
  ######
  bed2=with(df2, GRanges(chrom, IRanges(chromStart+1, chromEnd), id,strand = strand))
  head(bed2)
  length(bed2)
  ######
  ######
  # now find the overlaps
  #left
  df3.left=as.data.frame(subsetByOverlaps(bed2, bed1.left))
  head(df3.left)
  #right
  df3.right=as.data.frame(subsetByOverlaps(bed2, bed1.right))
  head(df3.right)
  #########
  df3=rbind(df3.left, df3.right)
  #
  calmmatrix=atac.signal.CRC.matrix[df3$id, ]
  calmmatrix[1:4,1:7]
  #
  calmean=as.data.frame(colMeans(calmmatrix))
  colnames(calmean)=te.list[k]
  
  datalist[[k]]=calmean
  
}
atac.agg.signal.stream=do.call(cbind, datalist)
rownames(atac.agg.signal.stream)=substr(rownames(atac.agg.signal.stream),1,15)
rownames(atac.agg.signal.stream)=gsub("[.]","-",rownames(atac.agg.signal.stream))
atac.agg.signal.stream[1:4,1:6]
#
####
head(atac.signal.id.CRC)
dim(atac.signal.id.CRC)
length(intersect(atac.signal.id.CRC$idd, rownames(atac.agg.signal.stream)))
####
atac.stat.stream=cbind(atac.signal.id.CRC, atac.agg.signal.stream[rownames(atac.signal.id.CRC),])
atac.stat.stream[1:4,1:8]
write.csv(atac.stat.stream, "4.model/kaps/atac.seq/stream/ht.atac.averaged.signal.stream.csv")
####
####
#heatmap
#
atac.drawmatrix.stream=as.data.frame(t(atac.stat.stream[,-c(1:4)]))
atac.drawmatrix.stream[1:6,1:4]
atac.drawmatrix.stream=na.omit(atac.drawmatrix.stream)
##
head(cdrep.atac)
#####
row_ha_left.atac= rowAnnotation(
  TE.class=cdrep.atac[rownames(atac.drawmatrix.stream),]$repClass.new2,
  col=list(TE.class=c(DNA="#1f78b4",LINE="#d95f02",LTR="#7570b3",Retroposon="#e7298a",Satellite="#66a61e",SINE="#e6ab02","other.repeats"="gray")
  ),
  show_annotation_name = FALSE)

col_ha_top.atac= columnAnnotation(
  TE.cluster=atac.stat.stream$kaps.group,
  col=list(TE.cluster=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )
  ),
  show_annotation_name = FALSE)
##
hmht.atac.stream=Heatmap(atac.drawmatrix.stream, name = "signal", 
                  col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                  #col = rev(viridis(10)),
                  width = unit(5, "cm"),
                  height = unit(12, "cm"),
                  border = F,
                  #col=col.pal_cor,
                  show_column_names = F,show_row_names = F,
                  cluster_columns = F,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  top_annotation = col_ha_top.atac,
                  #right_annotation = row_ha.right,
                  show_row_dend = F,show_column_dend = F,
                  #row_names_side = "left",
                  left_annotation = row_ha_left.atac,
                  column_title=paste0("average signal of promoter/enhancer up/down dtream TE location"),
                  column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/atac.seq/stream/ht.atac.averaged.signal.stream.pdf",width = 8,height = 6)
draw(hmht.atac.stream, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#####
