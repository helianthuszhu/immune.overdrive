############show the sig correlation matrx of 14 TEs
########
head(res.cor.29s)
sig.14.corM=res.cor.29s[res.cor.29s$repName %in% rownames(te.counts),]
head(sig.14.corM)
sigid=as.character(rownames(te.counts))
datalist=list()
for (i in 1:length(sigid)) {
  aa=as.data.frame(sig.14.corM[which(sig.14.corM$repName==sigid[i]),])
  aa2=as.data.frame(aa$cor)
  rownames(aa2)=aa$geneset
  colnames(aa2)=sigid[i]
  datalist[[i]]=aa2
}
sig.14.corM=do.call(cbind, datalist)
sig.14.corM=sig.14.corM[c(1:3,5:13,4,14:29),]
head(sig.14.corM)
dim(sig.14.corM)
#
min_cor = min(as.vector(sig.14.corM))
max_cor = max(as.vector(sig.14.corM))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=100)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
#####
head(te.counts)
col_ha_top.sig.14TEs=columnAnnotation(
  repClass=te.counts$repClass,
  col=list(repClass=c("DNA"="#1f78b4","LINE"="#d95f02","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02")
  ),
  show_annotation_name = FALSE,gp = gpar(col = "black"))
#
#pdf("4.model/kaps/piwi/piRNA/ht.cor.pdf",width = 6,height = 10)
ht.14tecounts=Heatmap(sig.14.corM,border = F,
        #col = col.pal_cor,
        cluster_columns = F,cluster_rows = F,
        #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
        #col=rev(viridis(10)),
        width = unit(6, "cm"),
        height = unit(7, "cm"),
        col=col.pal_cor,
        #rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        heatmap_legend_param = list(title = "Cor"),
        column_title = paste0("correlation of TEs with immune sets","\n",
                              "significantly correlated with","\n",
                              "at least one immune set (n=14, cor >=0.4)"),
        column_title_gp = gpar(fontsize = 10),
        column_title_side = "top",
        show_row_dend = F,show_column_dend = F,
        top_annotation = col_ha_top.sig.14TEs
)
pdf("3.immu.correlated/stat.20200405/cor.14s.TEs.vs.immune.sets.new2.pdf",width = 10)
draw(ht.14tecounts, padding = unit(c(20, 20, 20,20), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "left"
)
dev.off()
######
save(sig.14.corM,res.cor.29s,te.counts,file="3.immu.correlated/stat.20200405/cor.14s.TEs.vs.immune.sets.RData")
write.csv(sig.14.corM, "3.immu.correlated/stat.20200405/sig.14.corMatrix.csv")
