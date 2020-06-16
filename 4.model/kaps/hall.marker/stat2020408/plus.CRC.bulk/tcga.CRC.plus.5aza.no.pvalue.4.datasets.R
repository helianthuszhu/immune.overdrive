min_cor = min(as.vector(ht.hallM))
max_cor = max(as.vector(ht.hallM))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))

row_ha.left.hm = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.hm$p.value),
                               col=list(
                                 kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                           c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                               ),show_annotation_name = FALSE)
#
ht.hallM.or=ht.hallM[,c(4:1)]
col_ha_top.hm = columnAnnotation(
  kaps.group=colnames(ht.hallM.or),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))


ht0=Heatmap(ht.hallM.or, name = "TCGA.CRC", 
            #col = viridis(10),
            col=col.pal_cor,
            #width = unit(2, "cm"),height = unit(12, "cm"),
            border = F,
            show_column_names = F,show_row_names = T,cluster_columns = F,
            row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
            top_annotation = col_ha_top.hm,
            #right_annotation = row_ha.right,
            show_row_dend = F,show_column_dend = F#,
            #row_names_side = "left",
            #left_annotation = row_ha.left.hm
            #column_title=paste0("hall.marker.genesets"," GSE41586 CRC"),
            #column_title_gp = gpar(fontsize = 8)
)


pcb.2=ht0+ht1+ht2+ht3+ht4+ht5
#
pdf("4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/tcga.CRC.plus.5aza.no.pvalue.pdf",width = 12,height = 9)
draw(pcb.2, padding = unit(c(20, 15, 20,15), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()

pcb.3=ht0+ht1+ht2+ht3
pcb.4=ht1+ht2+ht3+ht0
#
pdf("4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/tcga.CRC.plus.5aza.no.pvalue.4.datasets.2.colorchanged.pdf",width = 12,height = 9)
draw(pcb.4, padding = unit(c(20, 15, 20,15), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()

#########################
#########################
#########################
CB.cell.line.ssgsea.pvalue.gse80317.tr <- apply(CB.cell.line.ssgsea.pvalue.gse80317,2,as.character)
CB.cell.line.ssgsea.pvalue.gse22250.tr<- apply(CB.cell.line.ssgsea.pvalue.gse22250,2,as.character)
CB.cell.line.ssgsea.pvalue.gse51811.tr <- apply(CB.cell.line.ssgsea.pvalue.gse51811,2,as.character)

fwrite(CB.cell.line.ssgsea.pvalue.gse5816,"4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/ssgsea.pvalue.gse5816.csv")
write.csv(CB.cell.line.ssgsea.pvalue.gse80317.tr,"4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/ssgsea.pvalue.gse80137.csv")
fwrite(CB.cell.line.ssgsea.pvalue.gse22250.tr,"4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/ssgsea.pvalue.gse22250.csv")
fwrite(CB.cell.line.ssgsea.pvalue.gse51811.tr,"4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/ssgsea.pvalue.gse51811.csv")
fwrite(CB.cell.line.ssgsea.pvalue.gse41586,"4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/ssgsea.pvalue.gse41586.csv")
fwrite(pvalue.hm,"4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/ssgsea.pvalue.CRC.bulk.csv")


