#####################28 cell types
######
head(rds.stats.gse170422)
#######
load("vali.GSE107422/immunesets/genesets.cell.type.score.GSE017422.RData")
dim(genesets.cell.type.score.GSE107422)
newmanscore.GSE107422=genesets.cell.type.score.GSE107422
#
newmanscore.GSE107422=as.data.frame(scale(newmanscore.GSE107422))
head(newmanscore.GSE107422)
ids.gse107422=intersect(rownames(rds.stats.gse170422), rownames(newmanscore.GSE107422))

ht.newmanscore.GSE107422=cbind(rds.stats.gse170422$group, newmanscore.GSE107422[rownames(rds.stats.gse170422),])
colnames(ht.newmanscore.GSE107422)[1]="kaps.group.gse107422"
head(ht.newmanscore.GSE107422)
dim(ht.newmanscore.GSE107422)
####

######pvalue
datalist <- list()
for(i in names(ht.newmanscore.GSE107422[,2:ncol(ht.newmanscore.GSE107422)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group.gse107422")), data = ht.newmanscore.GSE107422)
}
pvalue.ht.newman.gse107422=do.call(rbind, datalist)
pvalue.ht.newman.gse107422=as.data.frame(pvalue.ht.newman.gse107422)
######matrix
ht.stat.newman.gse107422= ht.newmanscore.GSE107422 %>% group_by(kaps.group.gse107422) %>% summarise_all(median,na.rm = TRUE)
ht.stat.newman.gse107422=as.data.frame(ht.stat.newman.gse107422)
rownames(ht.stat.newman.gse107422)=ht.stat.newman.gse107422$kaps.group.gse107422
ht.stat.newman.gse107422=ht.stat.newman.gse107422[,-1]
ht.stat.newman.gse107422=as.data.frame(t(ht.stat.newman.gse107422))
#ht.stat.newman.kirc=ht.stat.newman.kirc[,c(4,3,2,1)]
head(ht.stat.newman.gse107422)
save(pvalue.ht.newman.gse107422,ht.stat.newman.gse107422,file="vali.GSE107422/immunesets/celltype/ht.cell.type.28s.GSE107422.RData" )
#######
#colore set
min_cor = min(as.vector(ht.stat.newman.gse107422))
max_cor = max(as.vector(ht.stat.newman.gse107422))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#00F5FF", "white","#dd1c77"))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#5e4fa2", "white","#9e0142"))(50))
#
row_ha.left.newman.gse107422 = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.newman.gse107422$p.value),
                                        col=list(
                                          kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                                    c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                        ),show_annotation_name = FALSE)
col_ha_top.newman.gse107422 = columnAnnotation(
  kaps.group=colnames(ht.stat.newman.gse107422),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

newmanht.gse107422=Heatmap(ht.stat.newman.gse107422, name = "median.of.z.score", 
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
                      top_annotation = col_ha_top.newman.gse107422,
                      #right_annotation = row_ha.right,
                      show_row_dend = F,show_column_dend = F,
                      #row_names_side = "left",
                      left_annotation = row_ha.left.newman.gse107422,
                      column_title="28 cell type fraction.GSE107422",
                      column_title_gp = gpar(fontsize = 8)
)
pdf("vali.GSE107422/immunesets/celltype/ht.cell.type.28s.GSE107422.pdf",width = 5,height = 6)
draw(newmanht.gse107422, padding = unit(c(5, 5, 5,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
