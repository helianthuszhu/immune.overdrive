#####################
#####################28 cell types
######
load("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/cell.type.28s.gsva.RData")

newmanscore=genesets17s.score
#
newmanscore=as.data.frame(scale(newmanscore))
head(newmanscore)
ids=intersect(rownames(kaps.td), rownames(newmanscore))

ht.newmanscore=cbind(kaps.td$kaps.group, newmanscore[rownames(kaps.td),])
colnames(ht.newmanscore)[1]="kaps.group"
head(ht.newmanscore)
####

######pvalue
datalist <- list()
for(i in names(ht.newmanscore[,2:ncol(ht.newmanscore)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = ht.newmanscore)
}
pvalue.ht.newman=do.call(rbind, datalist)
pvalue.ht.newman=as.data.frame(pvalue.ht.newman)
######matrix
ht.stat.newman= ht.newmanscore %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.newman=as.data.frame(ht.stat.newman)
rownames(ht.stat.newman)=ht.stat.newman$kaps.group
ht.stat.newman=ht.stat.newman[,-1]
ht.stat.newman=as.data.frame(t(ht.stat.newman))
ht.stat.newman=ht.stat.newman[,c(4,3,2,1)]
head(ht.stat.newman)
save(pvalue.ht.newman,ht.stat.newman,file="4.model/kaps/cell.type.ssgsea/stat.2020408/ht.cell.type.28s.RData" )
#######
#colore set
min_cor = min(as.vector(ht.stat.newman))
max_cor = max(as.vector(ht.stat.newman))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.newman = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.newman$p.value),
                                   col=list(
                                     kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                               c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                   ),show_annotation_name = FALSE)
col_ha_top.newman = columnAnnotation(
  kaps.group=colnames(ht.stat.newman),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

newmanht=Heatmap(ht.stat.newman, name = "median.of.z.score", 
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
                 top_annotation = col_ha_top.newman,
                 #right_annotation = row_ha.right,
                 show_row_dend = F,show_column_dend = F,
                 #row_names_side = "left",
                 left_annotation = row_ha.left.newman,
                 column_title="cell type fraction from MCPcounter",
                 column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/cell.type.ssgsea/stat.2020408/ht.cell.type.28s.pdf",width = 5,height = 6)
draw(newmanht, padding = unit(c(5, 5, 5,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()