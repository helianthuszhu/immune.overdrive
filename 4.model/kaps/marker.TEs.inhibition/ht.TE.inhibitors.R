##########TE inhibition markers
##########
inhibition.markers=read.table("4.model/kaps/marker.TEs.inhibition/TEs.inhibition.markers.txt",header = T,sep = "\t")
rownames(inhibition.markers)=inhibition.markers$symbol
head(inhibition.markers)
####gene exp level compare
inhibition.markers.exp=genestats.orderd[, colnames(genestats.orderd) %in% rownames(inhibition.markers)]
inhibition.markers.exp[1:4,1:10]
#count the na 
inhinacount=as.data.frame(colSums(is.na(inhibition.markers.exp)))
colnames(inhinacount)="count"
inhinacount.sel=subset(inhinacount, count< 500)
inhibition.markers.sel=inhibition.markers[rownames(inhinacount.sel),]
head(inhibition.markers.sel)
dim(inhibition.markers.sel)
#
inhibition.markers.exp=genestats.orderd[, colnames(genestats.orderd) %in% rownames(inhinacount.sel)]

inhibition.markers.exp=as.data.frame(scale(inhibition.markers.exp))
inhibition.markers.exp=cbind(genestats.orderd$kaps.group,inhibition.markers.exp)
colnames(inhibition.markers.exp)[1]="kaps.group"
inhibition.markers.exp[1:4,1:10]
dim(inhibition.markers.exp)
class(inhibition.markers.exp$kaps.group)
#########
######pvalue
datalist <- list()
for(i in names(inhibition.markers.exp[,2:ncol(inhibition.markers.exp)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = inhibition.markers.exp,na.action=na.omit)
}
pvalue.ht.inhites=do.call(rbind, datalist)
pvalue.ht.inhites=as.data.frame(pvalue.ht.inhites)
########################
######matrix
ht.stat.inhi= inhibition.markers.exp %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.inhi=as.data.frame(ht.stat.inhi)
rownames(ht.stat.inhi)=ht.stat.inhi$kaps.group
ht.stat.inhi=ht.stat.inhi[,-1]
ht.stat.inhi=as.data.frame(t(ht.stat.inhi))
ht.stat.inhi=ht.stat.inhi[,c(4,3,2,1)]
head(ht.stat.inhi)
#
save(pvalue.ht.inhites, ht.stat.inhi,genestats.orderd,inhibition.markers,file="4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.RData")
#######
#colore set
min_cor = min(as.vector(ht.stat.inhi))
max_cor = max(as.vector(ht.stat.inhi))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.inhi = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites$p.value),
                                 category=inhibition.markers.sel$TEs.suppression.mechanism,
                                col=list(
                                  kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                            c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837")),
                                  category=c("chromatin.remodeling.proteins.repressing.TE"="#a6cee3",
                                             "L1 inhibitor"="#1f78b4",
                                             "Oligoadenylate.synthetases"="#b2df8a",
                                             "RNA.sensor"="#33a02c",
                                             "RNA.sensor.RIG.like.receptor"="#fb9a99",
                                             "RNA.sensor.TLRs"="#e31a1c",
                                             "SVA.repression"="#fdbf6f"
                                  )
                                ),show_annotation_name = FALSE)
col_ha_top.inhi = columnAnnotation(
  kaps.group=colnames(ht.stat.inhi),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

inhiht=Heatmap(ht.stat.inhi, name = "median.of.z.score", 
              #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
              #col = rev(viridis(10)),
              width = unit(2, "cm"),
              #height = unit(12, "cm"),
              border = F,
              col=col.pal_cor,
              show_column_names = T,show_row_names = T,
              cluster_columns = F,cluster_rows = F,
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5),
              top_annotation = col_ha_top.inhi,
              #right_annotation = row_ha.right,
              show_row_dend = F,show_column_dend = F,
              #row_names_side = "left",
              left_annotation = row_ha.left.inhi,
              column_title="TE inhibitiors",
              column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.pdf",width = 6,height = 7)
draw(inhiht, padding = unit(c(10, 5, 10,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###########show individually
############
##########
##########
###########APOBECs
APOBEC.matrix=ht.stat.inhi[grep("APOBEC",rownames(ht.stat.inhi)),]
pvalue.ht.inhites.apobec=pvalue.ht.inhites[rownames(APOBEC.matrix),]
#colore set
min_cor = min(as.vector(APOBEC.matrix))
max_cor = max(as.vector(APOBEC.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.inhi = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.apobec$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-3,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                  
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi = columnAnnotation(
  kaps.group=colnames(APOBEC.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht=Heatmap(APOBEC.matrix, name = "median.of.z.score", 
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
               top_annotation = col_ha_top.inhi,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi,
               column_title="APOBECs",
               column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.apobec.pdf",width = 5,height = 6)
draw(inhiht, padding = unit(c(50, 5, 50,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
##########
#########
#########chromatin
head(inhibition.markers.sel)
chromatin.matrix=ht.stat.inhi[rownames(subset(inhibition.markers.sel, TEs.suppression.mechanism=="chromatin.remodeling.proteins.repressing.TE")),]
pvalue.ht.inhites.chromatin=pvalue.ht.inhites[rownames(chromatin.matrix),]
#colore set
min_cor = min(as.vector(ht.stat.inhi))
max_cor = max(as.vector(ht.stat.inhi))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.inhi = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.1,0.05,0.01), 
                                                             c("#4d4d4d","#ffffcc","#d9f0a3"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi = columnAnnotation(
  kaps.group=colnames(chromatin.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht=Heatmap(chromatin.matrix, name = "median.of.z.score", 
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
               top_annotation = col_ha_top.inhi,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi,
               column_title="chromatin.remodeling.proteins",
               column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.chromatin.remodeling.proteins.new.pdf",width = 5,height = 6)
draw(inhiht, padding = unit(c(60, 5, 65,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
############################
##########OLi
##########
head(inhibition.markers.sel)
table(inhibition.markers.sel$TEs.suppression.mechanism)
chromatin.matrix=ht.stat.inhi[rownames(subset(inhibition.markers.sel, TEs.suppression.mechanism=="Oligoadenylate.synthetases")),]
pvalue.ht.inhites.chromatin=pvalue.ht.inhites[rownames(chromatin.matrix),]
#colore set
min_cor = min(as.vector(chromatin.matrix))
max_cor = max(as.vector(chromatin.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.inhi = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.1,0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi = columnAnnotation(
  kaps.group=colnames(chromatin.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht=Heatmap(chromatin.matrix, name = "median.of.z.score", 
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
               top_annotation = col_ha_top.inhi,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi,
               column_title="Oligoadenylate.synthetases",
               column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.Oligoadenylate.synthetases.pdf",width = 5,height = 6)
draw(inhiht, padding = unit(c(55, 5, 55,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#############
############
############TLRs
head(inhibition.markers.sel)
table(inhibition.markers.sel$TEs.suppression.mechanism)
chromatin.matrix=ht.stat.inhi[rownames(subset(inhibition.markers.sel, TEs.suppression.mechanism=="RNA.sensor.TLRs")),]
pvalue.ht.inhites.chromatin=pvalue.ht.inhites[rownames(chromatin.matrix),]
#colore set
min_cor = min(as.vector(chromatin.matrix))
max_cor = max(as.vector(chromatin.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.inhi = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.1,0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi = columnAnnotation(
  kaps.group=colnames(chromatin.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht=Heatmap(chromatin.matrix, name = "median.of.z.score", 
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
               top_annotation = col_ha_top.inhi,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi,
               column_title="TLRs",
               column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.TLRs.pdf",width = 5,height = 6)
draw(inhiht, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
################
################RIG-like
################
head(inhibition.markers.sel)
table(inhibition.markers.sel$TEs.suppression.mechanism)
chromatin.matrix=ht.stat.inhi[rownames(subset(inhibition.markers.sel, TEs.suppression.mechanism=="RNA.sensor.RIG.like.receptor")),]
pvalue.ht.inhites.chromatin=pvalue.ht.inhites[rownames(chromatin.matrix),]
class(pvalue.ht.inhites.chromatin)
pvalue.ht.inhites.chromatin=subset(pvalue.ht.inhites.chromatin, p.value<0.05)
chromatin.matrix=chromatin.matrix[rownames(pvalue.ht.inhites.chromatin),]
#colore set
min_cor = min(as.vector(chromatin.matrix))
max_cor = max(as.vector(chromatin.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.inhi = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi = columnAnnotation(
  kaps.group=colnames(chromatin.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht=Heatmap(chromatin.matrix, name = "median.of.z.score", 
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
               top_annotation = col_ha_top.inhi,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi,
               column_title="RIG-like",
               column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.RIG-like.pdf",width = 5,height = 6)
draw(inhiht, padding = unit(c(5, 5, 5,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
##################
##################
##################
############other inhibitors
head(inhibition.markers.sel)
table(inhibition.markers.sel$TEs.suppression.mechanism)
chromatin.matrix1=ht.stat.inhi[rownames(subset(inhibition.markers.sel, TEs.suppression.mechanism=="RNA.sensor"|TEs.suppression.mechanism=="SVA.repression")),]
chromatin.matrix2=ht.stat.inhi[rownames(subset(inhibition.markers.sel, TEs.suppression.mechanism=="L1 inhibitor")),][-c(2:9),]
chromatin.matrix=rbind(chromatin.matrix1,chromatin.matrix2)
pvalue.ht.inhites.chromatin=pvalue.ht.inhites[rownames(chromatin.matrix),]
#colore set
min_cor = min(as.vector(chromatin.matrix))
max_cor = max(as.vector(chromatin.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.inhi = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.inhites.chromatin$p.value),
                                 
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.1,0.05,0.001,10e-5,10e-7,10e-15), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   
                                 ),show_annotation_name = FALSE)
col_ha_top.inhi = columnAnnotation(
  kaps.group=colnames(chromatin.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
inhiht=Heatmap(chromatin.matrix, name = "median.of.z.score", 
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
               top_annotation = col_ha_top.inhi,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.inhi,
               column_title="other RNA sensors",
               column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/ht.TE.inhibitors.other.RNA.sensors.pdf",width = 5,height = 6)
draw(inhiht, padding = unit(c(35, 5, 40,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
