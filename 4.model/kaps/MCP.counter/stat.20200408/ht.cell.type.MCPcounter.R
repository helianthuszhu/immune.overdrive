######draw heatmap of MCPcounter 10 cell types
######
load("4.model/kaps/MCP.counter/socre.Mcpcounter.CRC.RData")
mcpscore=as.data.frame(t(score.mcp.CRC))
mcpscore=as.data.frame(scale(mcpscore))
head(mcpscore)
ids=intersect(rownames(kaps.td), rownames(mcpscore))

ht.mcpscore=cbind(kaps.td$kaps.group, mcpscore[rownames(kaps.td),])
colnames(ht.mcpscore)[1]="kaps.group"
head(ht.mcpscore)
######
######pvalue
datalist <- list()
for(i in names(ht.mcpscore[,2:ncol(ht.mcpscore)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = ht.mcpscore)
}
pvalue.ht.mcp=do.call(rbind, datalist)
pvalue.ht.mcp=as.data.frame(pvalue.ht.mcp)
######matrix
ht.stat.mcp= ht.mcpscore %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.mcp=as.data.frame(ht.stat.mcp)
rownames(ht.stat.mcp)=ht.stat.mcp$kaps.group
ht.stat.mcp=ht.stat.mcp[,-1]
ht.stat.mcp=as.data.frame(t(ht.stat.mcp))
ht.stat.mcp=ht.stat.mcp[,c(4,3,2,1)]
head(ht.stat.mcp)

save(pvalue.ht.mcp, ht.stat.mcp,file="4.model/kaps/MCP.counter/stat.20200408/ht.cell.type.MCPcounter.RData")
#######
#colore set
min_cor = min(as.vector(ht.stat.mcp))
max_cor = max(as.vector(ht.stat.mcp))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.mcp = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.mcp$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
col_ha_top.mcp = columnAnnotation(
  kaps.group=colnames(ht.stat.mcp),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

mcpht=Heatmap(ht.stat.mcp, name = "mean.of.z.score", 
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
             top_annotation = col_ha_top.mcp,
             #right_annotation = row_ha.right,
             show_row_dend = F,show_column_dend = F,
             #row_names_side = "left",
             left_annotation = row_ha.left.mcp,
             column_title="cell type fraction from MCPcounter",
             column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/MCP.counter/stat.20200408/ht.cell.type.MCPcounter.pdf",width = 5,height = 6)
draw(mcpht, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()


