#######heatmap of IPRES
########
ht.ipres=CB.data.ipres.kaps[order(CB.data.ipres.kaps$z.of.mean.exp,decreasing = T),c(45:68)]


col_ha.ipres=columnAnnotation(z.mean.IPRES=CB.data.ipres.kaps[order(CB.data.ipres.kaps$z.of.mean.exp,decreasing = T),]$z.mean.IPRES,
                             kaps.group=CB.data.ipres.kaps[order(CB.data.ipres.kaps$z.of.mean.exp,decreasing = T),]$kaps.group,
                             col=list(
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)
#
ann.ipres=CB.data.ipres.kaps[order(CB.data.ipres.kaps$z.of.mean.exp,decreasing = T),c(38,43)]

ada=as.data.frame(t(ht.ipres))
z_ada=(ada - rowMeans(ada))/apply(ada,1,sd)
z_ada[z_ada< -3] <- -3
z_ada[z_ada>3] <- 3

min_cor = min(as.vector(z_ada))
max_cor = max(as.vector(z_ada))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=100)
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "Spectral")))(100))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#7fcdbb", "white","#ce1256"))(100))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#4393c3", "white","#ce1256"))(100))


pdf("4.model/kaps/ipres/htmap.ipres.pdf",width = 8,height = 3)
Heatmap(z_ada, name = "ipres", 
              #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
              #col = rev(viridis(10)),border = F,
              col=col.pal_cor,
              show_column_names = F,show_row_names = T,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 8),
              top_annotation = col_ha.ipres#, 
              #right_annotation = row_ha.right,
              #show_row_dend = T,show_column_dend = T,
              #row_names_side = "left",
              #left_annotation = row_ha.left
)
dev.off()

save(CB.data.ipres.kaps,ht.ipres, ann.ipres,file="4.model/kaps/ipres/htmap.ipres.RData")






#
pdf("4.model/kaps/ipres/htmap.ipres.pdf",width = 8,height = 2)
print(
  pheatmap(ada,show_rownames = T,border_color = NA,main = "ipres",
           show_colnames = F,fontsize = 4,#na_col = "black",
           #cluster_rows = hc,
           cluster_cols =F,
           annotation_col = ann.ipres,
           #annotation_row = cndi.rep[rownames(ppd),c(2,3)],
           #scale="row",
           color=col.pal_cor,
           #color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8), 
           #color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
           #annotation_legend = T,
           annotation_colors  = list( 
             kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")#,
             #repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
             #repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
           )
  )
)
dev.off()
