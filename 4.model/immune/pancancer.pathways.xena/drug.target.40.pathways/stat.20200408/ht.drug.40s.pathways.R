stat.pathyp.drug[1:4,1:4]
######
g40s.drug=stat.pathyp.drug[,-c(1:2)]
g40s.drug=as.data.frame(scale(g40s.drug))
g40s.drug=cbind(stat.pathyp.drug$kaps.group,g40s.drug)
colnames(g40s.drug)[1]="kaps.group"
colnames(g40s.drug)=gsub("[/]",".",colnames(g40s.drug))
colnames(g40s.drug)=gsub("-",".",colnames(g40s.drug))
colnames(g40s.drug)=gsub("[(]",".",colnames(g40s.drug))
colnames(g40s.drug)=gsub("[)]",".",colnames(g40s.drug))
######pvalue
datalist <- list()
for(i in names(g40s.drug[,2:ncol(g40s.drug)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = g40s.drug)
}
pvalue.ht.40sdrug=do.call(rbind, datalist)
pvalue.ht.40sdrug=as.data.frame(pvalue.ht.40sdrug)
pvalue.ht.40sdrug.sig=subset(pvalue.ht.40sdrug, p.value<0.05)
######matrix
ht.stat.40sdrug= g40s.drug %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.40sdrug=as.data.frame(ht.stat.40sdrug)
rownames(ht.stat.40sdrug)=ht.stat.40sdrug$kaps.group
ht.stat.40sdrug=ht.stat.40sdrug[,-1]
ht.stat.40sdrug=as.data.frame(t(ht.stat.40sdrug))
ht.stat.40sdrug=ht.stat.40sdrug[,c(4,3,2,1)]
#
ht.stat.40sdrug.sig=ht.stat.40sdrug[rownames(pvalue.ht.40sdrug.sig),]
ht.stat.40sdrug.sig=(ht.stat.40sdrug.sig - rowMeans(ht.stat.40sdrug.sig))/apply(ht.stat.40sdrug.sig,1,sd)
#
head(ht.stat.40sdrug.sig)

save(pvalue.ht.40sdrug, ht.stat.40sdrug,ht.stat.40sdrug.sig,file="4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/stat.20200408/ht.drug.40s.pathways.RData")
#######
#colore set
min_cor = min(as.vector(ht.stat.40sdrug.sig))
max_cor = max(as.vector(ht.stat.40sdrug.sig))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.40sdrug = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.40sdrug$p.value),
                                col=list(
                                  kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                            c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                ),show_annotation_name = FALSE)
col_ha_top.40sdrug = columnAnnotation(
  kaps.group=colnames(ht.stat.40sdrug),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))
#
row_ha.left.40sdrug = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.40sdrug.sig$p.value),
                                    col=list(
                                      kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                                c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                    ),show_annotation_name = FALSE)
col_ha_top.40sdrug = columnAnnotation(
  kaps.group=colnames(ht.stat.40sdrug.sig),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))


####

drught=Heatmap(ht.stat.40sdrug.sig, name = "median.of.z.score", 
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
              top_annotation = col_ha_top.40sdrug,
              #right_annotation = row_ha.right,
              show_row_dend = F,show_column_dend = F,
              #row_names_side = "left",
              left_annotation = row_ha.left.40sdrug,
              column_title=paste0("gene program and","\n","canonical targetable pathways"),
              column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/stat.20200408/ht.drug.40s.pathways.scale.only.sig.pdf",width = 5,height = 6)
draw(drught, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
