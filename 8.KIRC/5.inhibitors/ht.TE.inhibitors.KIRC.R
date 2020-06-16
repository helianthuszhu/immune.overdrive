##########TE inhibition markers in KIRC
###get expmatrx and kaps.group kirc
genestats.KIRC=as.data.frame(t(kirc.expfireh))
genestats.KIRC[1:4,1:4]
###get kaps group kirc
head(clin.heat.kirc.var)
table(clin.heat.kirc.var$kaps.group.kirc)
length(intersect(rownames(clin.heat.kirc.var), rownames(genestats.KIRC)))
###
genestats.KIRC=cbind(kaps.group=clin.heat.kirc.var$kaps.group.kirc, genestats.KIRC[rownames(clin.heat.kirc.var),])
genestats.KIRC[1:4,1:4]
##########
inhibition.markers=read.table("4.model/kaps/marker.TEs.inhibition/TEs.inhibition.markers.txt",header = T,sep = "\t")
rownames(inhibition.markers)=inhibition.markers$symbol
dim(inhibition.markers)
#
####gene exp level compare
inhibition.markers.exp.KIRC=genestats.KIRC[, colnames(genestats.KIRC) %in% rownames(inhibition.markers)]
inhibition.markers.exp.KIRC[1:4,1:10]
dim(inhibition.markers.exp.KIRC)
#count the na 
inhinacount.KIRC=as.data.frame(colSums(is.na(inhibition.markers.exp.KIRC)))
colnames(inhinacount.KIRC)="count"
inhinacount.sel.KIRC=subset(inhinacount.KIRC, count< 400)
inhibition.markers.sel.KIRC=inhibition.markers[rownames(inhinacount.sel.KIRC),]
head(inhibition.markers.sel.KIRC)
dim(inhibition.markers.sel.KIRC)
#
inhibition.markers.exp.KIRC=genestats.KIRC[, colnames(genestats.KIRC) %in% rownames(inhinacount.sel.KIRC)]
dim(inhibition.markers.exp.KIRC)

inhibition.markers.exp.KIRC=as.data.frame(scale(inhibition.markers.exp.KIRC))
inhibition.markers.exp.KIRC[1:4,1:10]

inhibition.markers.exp.KIRC=cbind(genestats.KIRC$kaps.group,inhibition.markers.exp.KIRC[rownames(genestats.KIRC),])

colnames(inhibition.markers.exp.KIRC)[1]="kaps.group"
inhibition.markers.exp.KIRC[1:4,1:10]

class(inhibition.markers.exp.KIRC$kaps.group)
inhibition.markers.exp.KIRC$kaps.group=as.factor(inhibition.markers.exp.KIRC$kaps.group)   ##should be factor

dim(inhibition.markers.exp.KIRC)
table(inhibition.markers.exp.KIRC$kaps.group)
#########
######pvalue
datalist <- list()
for(i in names(inhibition.markers.exp.KIRC[,2:ncol(inhibition.markers.exp.KIRC)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = inhibition.markers.exp.KIRC,na.action=na.omit)
}
pvalue.ht.inhites.KIRC=do.call(rbind, datalist)
pvalue.ht.inhites.KIRC=as.data.frame(pvalue.ht.inhites.KIRC)
########################
######matrix
ht.stat.inhi.KIRC= inhibition.markers.exp.KIRC %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.inhi.KIRC=as.data.frame(ht.stat.inhi.KIRC)
rownames(ht.stat.inhi.KIRC)=ht.stat.inhi.KIRC$kaps.group
ht.stat.inhi.KIRC=ht.stat.inhi.KIRC[,-1]
ht.stat.inhi.KIRC=as.data.frame(t(ht.stat.inhi.KIRC))
ht.stat.inhi.KIRC=ht.stat.inhi.KIRC[,c(4:1)]
head(ht.stat.inhi.KIRC)
#
save(pvalue.ht.inhites.KIRC, ht.stat.inhi.KIRC,genestats.KIRC,inhibition.markers,file="8.KIRC/5.inhibitors/ht.TE.inhibitors.KIRC.RData")
#######
###########show individually
############
###########APOBECs
APOBEC.matrix.KIRC=ht.stat.inhi.KIRC[grep("APOBEC",rownames(ht.stat.inhi.KIRC)),]
pvalue.ht.inhites.apobec.KIRC=pvalue.ht.inhites.KIRC[rownames(APOBEC.matrix.KIRC),]
#colore set
min_cor = min(as.vector(APOBEC.matrix.KIRC))
max_cor = max(as.vector(APOBEC.matrix.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.inhi.KIRC = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.apobec.KIRC$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-3,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi.KIRC = columnAnnotation(
  kaps.group=colnames(APOBEC.matrix.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht.KIRC=Heatmap(APOBEC.matrix.KIRC, name = "median.of.z.score", 
               #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),
               width = unit(2, "cm"),
               #height = unit(12, "cm"),
               border = F,
               col=col.pal_cor,
               show_column_names = T,show_row_names = T,
               cluster_columns = F,cluster_rows = T,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               top_annotation = col_ha_top.inhi.KIRC,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi.KIRC,
               column_title="APOBECs.KIRC",
               column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/ht.TE.inhibitors.apobec.KIRC.pdf",width = 5,height = 6)
draw(inhiht.KIRC, padding = unit(c(50, 5, 50,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
##########
#########
#########chromatin
head(inhibition.markers.sel.KIRC)
chromatin.matrix.KIRC=ht.stat.inhi.KIRC[rownames(subset(inhibition.markers.sel.KIRC, TEs.suppression.mechanism=="chromatin.remodeling.proteins.repressing.TE")),]
pvalue.ht.inhites.chromatin.KIRC=pvalue.ht.inhites.KIRC[rownames(chromatin.matrix.KIRC),]
#colore set
min_cor = min(as.vector(chromatin.matrix.KIRC))
max_cor = max(as.vector(chromatin.matrix.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.inhi.KIRC = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin.KIRC$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-3,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi.KIRC = columnAnnotation(
  kaps.group=colnames(chromatin.matrix.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht.KIRC=Heatmap(chromatin.matrix.KIRC, name = "median.of.z.score", 
               #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),
               width = unit(2, "cm"),
               height = unit(1, "cm"),
               border = F,
               col=col.pal_cor,
               show_column_names = T,show_row_names = T,
               cluster_columns = F,cluster_rows = T,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               top_annotation = col_ha_top.inhi.KIRC,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi.KIRC,
               column_title="chromatin.remodeling.proteins.KIRC",
               column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/ht.TE.inhibitors.chromatin.remodeling.proteins.KIRC.pdf",width = 5,height = 6)
draw(inhiht.KIRC, padding = unit(c(50, 5, 50,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
############################
##########OLi
##########
head(inhibition.markers.sel.KIRC)
table(inhibition.markers.sel.KIRC$TEs.suppression.mechanism)
chromatin.matrix.KIRC=ht.stat.inhi.KIRC[rownames(subset(inhibition.markers.sel.KIRC, TEs.suppression.mechanism=="Oligoadenylate.synthetases")),]
pvalue.ht.inhites.chromatin.KIRC=pvalue.ht.inhites.KIRC[rownames(chromatin.matrix.KIRC),]
#colore set
min_cor = min(as.vector(chromatin.matrix.KIRC))
max_cor = max(as.vector(chromatin.matrix.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.inhi.KIRC = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin.KIRC$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.1,0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi.KIRC = columnAnnotation(
  kaps.group=colnames(chromatin.matrix.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht.KIRC=Heatmap(chromatin.matrix.KIRC, name = "median.of.z.score", 
               #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),
               width = unit(2, "cm"),
               #height = unit(12, "cm"),
               border = F,
               col=col.pal_cor,
               show_column_names = T,show_row_names = T,
               cluster_columns = F,cluster_rows = T,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               top_annotation = col_ha_top.inhi.KIRC,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi.KIRC,
               column_title="Oligoadenylate.synthetases.KIRC",
               column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/ht.TE.inhibitors.Oligoadenylate.synthetases.KIRC.pdf",width = 5,height = 6)
draw(inhiht.KIRC, padding = unit(c(55, 5, 55,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#############
############
############TLRs
head(inhibition.markers.sel.KIRC)
table(inhibition.markers.sel.KIRC$TEs.suppression.mechanism)
chromatin.matrix.KIRC=ht.stat.inhi.KIRC[rownames(subset(inhibition.markers.sel.KIRC, TEs.suppression.mechanism=="RNA.sensor.TLRs")),]
pvalue.ht.inhites.chromatin.KIRC=pvalue.ht.inhites.KIRC[rownames(chromatin.matrix.KIRC),]
#colore set
min_cor = min(as.vector(chromatin.matrix.KIRC))
max_cor = max(as.vector(chromatin.matrix.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.inhi.KIRC = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin.KIRC$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.1,0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi.KIRC = columnAnnotation(
  kaps.group=colnames(chromatin.matrix.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht.KIRC=Heatmap(chromatin.matrix.KIRC, name = "median.of.z.score", 
               #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),
               width = unit(2, "cm"),
               #height = unit(12, "cm"),
               border = F,
               col=col.pal_cor,
               show_column_names = T,show_row_names = T,
               cluster_columns = F,cluster_rows = T,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               top_annotation = col_ha_top.inhi.KIRC,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi.KIRC,
               column_title="TLRs.KIRC",
               column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/ht.TE.inhibitors.TLRs.KIRC.pdf",width = 5,height = 6)
draw(inhiht.KIRC, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
################
################RIG-like
################
head(inhibition.markers.sel.KIRC)
table(inhibition.markers.sel.KIRC$TEs.suppression.mechanism)
chromatin.matrix.KIRC=ht.stat.inhi.KIRC[rownames(subset(inhibition.markers.sel.KIRC, TEs.suppression.mechanism=="RNA.sensor.RIG.like.receptor")),]
pvalue.ht.inhites.chromatin.KIRC=pvalue.ht.inhites.KIRC[rownames(chromatin.matrix.KIRC),]
class(pvalue.ht.inhites.chromatin.KIRC)
####only sig
pvalue.ht.inhites.chromatin.KIRC=subset(pvalue.ht.inhites.chromatin.KIRC, p.value<0.05)
chromatin.matrix.KIRC=chromatin.matrix.KIRC[rownames(pvalue.ht.inhites.chromatin.KIRC),]
#colore set
min_cor = min(as.vector(chromatin.matrix.KIRC))
max_cor = max(as.vector(chromatin.matrix.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.inhi.KIRC = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin.KIRC$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#ffffcc","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi.KIRC = columnAnnotation(
  kaps.group=colnames(chromatin.matrix.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht.KIRC=Heatmap(chromatin.matrix.KIRC, name = "median.of.z.score", 
               #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),
               width = unit(2, "cm"),
               #height = unit(12, "cm"),
               border = F,
               col=col.pal_cor,
               show_column_names = T,show_row_names = T,
               cluster_columns = F,cluster_rows = T,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               top_annotation = col_ha_top.inhi.KIRC,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi.KIRC,
               column_title="RIG-like.KIRC",
               column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/ht.TE.inhibitors.RIG-like.KIRC.only.sig.pdf",width = 5,height = 6)
draw(inhiht.KIRC, padding = unit(c(5, 5, 5,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
##################
##################
##################
############other inhibitors
head(inhibition.markers.sel.KIRC)
table(inhibition.markers.sel.KIRC$TEs.suppression.mechanism)
chromatin.matrix1.KIRC=ht.stat.inhi.KIRC[rownames(subset(inhibition.markers.sel.KIRC, TEs.suppression.mechanism=="RNA.sensor"|TEs.suppression.mechanism=="SVA.repression")),]
chromatin.matrix2.KIRC=ht.stat.inhi.KIRC[rownames(subset(inhibition.markers.sel.KIRC, TEs.suppression.mechanism=="L1 inhibitor")),][-c(2:9),]
chromatin.matrix.KIRC=rbind(chromatin.matrix1.KIRC,chromatin.matrix2.KIRC)
pvalue.ht.inhites.chromatin.KIRC=pvalue.ht.inhites.KIRC[rownames(chromatin.matrix.KIRC),]
#colore set
min_cor = min(as.vector(chromatin.matrix.KIRC))
max_cor = max(as.vector(chromatin.matrix.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.inhi.KIRC = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin.KIRC$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.1,0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi.KIRC = columnAnnotation(
  kaps.group=colnames(chromatin.matrix.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht.KIRC=Heatmap(chromatin.matrix.KIRC, name = "median.of.z.score", 
               #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),
               width = unit(2, "cm"),
               #height = unit(12, "cm"),
               border = F,
               col=col.pal_cor,
               show_column_names = T,show_row_names = T,
               cluster_columns = F,cluster_rows = T,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               top_annotation = col_ha_top.inhi.KIRC,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi.KIRC,
               column_title="other RNA sensors.KIRC",
               column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/ht.TE.inhibitors.other.RNA.sensors.KIRC.pdf",width = 5,height = 6)
draw(inhiht.KIRC, padding = unit(c(35, 5, 40,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###############

