######draw heatmap of MCPcounter 10 cell types
######
load("8.KIRC/3.immune.sets.kirc/mcpCounter/socre.Mcpcounter.kirc.RData")

mcpscore.kirc=as.data.frame(t(score.mcp.kirc))
mcpscore.kirc=as.data.frame(scale(mcpscore.kirc))
head(mcpscore.kirc)
#ids=intersect(rownames(kaps.td), rownames(mcpscore))

ht.mcpscore.kirc=cbind(stat.kirc.kaps.vali$kaps.group.kirc, mcpscore.kirc[rownames(stat.kirc.kaps.vali),])
colnames(ht.mcpscore.kirc)[1]="kaps.group.kirc"
head(ht.mcpscore.kirc)
######
######pvalue
datalist <- list()
for(i in names(ht.mcpscore.kirc[,2:ncol(ht.mcpscore.kirc)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group.kirc")), data = ht.mcpscore.kirc)
}
pvalue.ht.mcp.kirc=do.call(rbind, datalist)
pvalue.ht.mcp.kirc=as.data.frame(pvalue.ht.mcp.kirc)
######matrix
ht.stat.mcp.kirc= ht.mcpscore.kirc %>% group_by(kaps.group.kirc) %>% summarise_all(median,na.rm = TRUE)
ht.stat.mcp.kirc=as.data.frame(ht.stat.mcp.kirc)
rownames(ht.stat.mcp.kirc)=ht.stat.mcp.kirc$kaps.group.kirc
ht.stat.mcp.kirc=ht.stat.mcp.kirc[,-1]
ht.stat.mcp.kirc=as.data.frame(t(ht.stat.mcp.kirc))
#ht.stat.mcp.kirc=ht.stat.mcp.kirc[,c(4,3,2,1)]
head(ht.stat.mcp.kirc)

save(pvalue.ht.mcp.kirc, ht.stat.mcp.kirc,file="8.KIRC/3.immune.sets.kirc/mcpCounter/ht.cell.type.MCPcounter.kirc.RData")
#######
#colore set
min_cor = min(as.vector(ht.stat.mcp.kirc))
max_cor = max(as.vector(ht.stat.mcp.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.mcp.kirc = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.mcp.kirc$p.value),
                                col=list(
                                  kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                            c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                ),show_annotation_name = FALSE)
col_ha_top.mcp.kirc = columnAnnotation(
  kaps.group=colnames(ht.stat.mcp.kirc),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

mcpht.kirc=Heatmap(ht.stat.mcp.kirc, name = "mean.of.z.score", 
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
              top_annotation = col_ha_top.mcp.kirc,
              #right_annotation = row_ha.right,
              show_row_dend = F,show_column_dend = F,
              #row_names_side = "left",
              left_annotation = row_ha.left.mcp.kirc,
              column_title="cell type fraction from MCPcounter.KIRC",
              column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/3.immune.sets.kirc/mcpCounter/ht.cell.type.MCPcounter.kirc.pdf",width = 5,height = 6)
draw(mcpht.kirc, padding = unit(c(30, 5, 30,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()


