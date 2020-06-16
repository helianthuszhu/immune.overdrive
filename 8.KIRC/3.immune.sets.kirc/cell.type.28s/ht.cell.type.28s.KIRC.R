#####################28 cell types
######
#######
load("8.KIRC/3.immune.sets.kirc/genesets.cell.type.score.kircfh.RData")
dim(genesets.cell.type.score.kircfh)
newmanscore.kirc=genesets.cell.type.score.kircfh
#
newmanscore.kirc=as.data.frame(scale(newmanscore.kirc))
head(newmanscore)
ids.kirc=intersect(rownames(stat.kirc.kaps.vali), rownames(newmanscore.kirc))

ht.newmanscore.kirc=cbind(stat.kirc.kaps.vali$kaps.group.kirc, newmanscore.kirc[rownames(stat.kirc.kaps.vali),])
colnames(ht.newmanscore.kirc)[1]="kaps.group.kirc"
head(ht.newmanscore.kirc)
dim(ht.newmanscore.kirc)
####

######pvalue
datalist <- list()
for(i in names(ht.newmanscore.kirc[,2:ncol(ht.newmanscore.kirc)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group.kirc")), data = ht.newmanscore.kirc)
}
pvalue.ht.newman.kirc=do.call(rbind, datalist)
pvalue.ht.newman.kirc=as.data.frame(pvalue.ht.newman.kirc)
######matrix
ht.stat.newman.kirc= ht.newmanscore.kirc %>% group_by(kaps.group.kirc) %>% summarise_all(median,na.rm = TRUE)
ht.stat.newman.kirc=as.data.frame(ht.stat.newman.kirc)
rownames(ht.stat.newman.kirc)=ht.stat.newman.kirc$kaps.group.kirc
ht.stat.newman.kirc=ht.stat.newman.kirc[,-1]
ht.stat.newman.kirc=as.data.frame(t(ht.stat.newman.kirc))
#ht.stat.newman.kirc=ht.stat.newman.kirc[,c(4,3,2,1)]
head(ht.stat.newman.kirc)
save(pvalue.ht.newman.kirc,ht.stat.newman.kirc,file="8.KIRC/3.immune.sets.kirc/cell.type.28s/ht.cell.type.28s.KIRC..RData" )
#######
#colore set
min_cor = min(as.vector(ht.stat.newman.kirc))
max_cor = max(as.vector(ht.stat.newman.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#00F5FF", "white","#dd1c77"))(50))

col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.newman.kirc = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.newman.kirc$p.value),
                                   col=list(
                                     kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                               c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   ),show_annotation_name = FALSE)
col_ha_top.newman.kirc = columnAnnotation(
  kaps.group=colnames(ht.stat.newman.kirc),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

newmanht.kirc=Heatmap(ht.stat.newman.kirc, name = "median.of.z.score", 
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
                 top_annotation = col_ha_top.newman.kirc,
                 #right_annotation = row_ha.right,
                 show_row_dend = F,show_column_dend = F,
                 #row_names_side = "left",
                 left_annotation = row_ha.left.newman.kirc,
                 column_title="28 cell type fraction.KIRC",
                 column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/3.immune.sets.kirc/cell.type.28s/ht.cell.type.28s.KIRC.pdf",width = 5,height = 6)
draw(newmanht.kirc, padding = unit(c(5, 5, 5,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
