######piRNA
######
stat.piRNAm.sel[1:4,1:4]
normd.pirna=as.data.frame(scale(stat.piRNAm.sel[,-c(1:2)]))
normd.pirna=cbind(stat.piRNAm.sel$kaps.group, normd.pirna)
colnames(normd.pirna)[1]="kaps.group"
normd.pirna[1:4,1:4]
#pvalue
#
datalist <- list()
for(i in names(normd.pirna[,2:ncol(normd.pirna)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = normd.pirna,na.action=na.omit)
}
pvalue.ht.piRNA=do.call(rbind, datalist)
pvalue.ht.piRNA=as.data.frame(pvalue.ht.piRNA)
pvalue.ht.piRNA.sig=subset(pvalue.ht.piRNA, p.value<0.05)
dim(pvalue.ht.piRNA.sig)
#matrix
#
ht.stat.piRNA= normd.pirna %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.piRNA =as.data.frame(ht.stat.piRNA)
rownames(ht.stat.piRNA)=ht.stat.piRNA$kaps.group
ht.stat.piRNA=ht.stat.piRNA[,-1]
ht.stat.piRNA=as.data.frame(t(ht.stat.piRNA))
#ht.stat.piRNA=ht.stat.piRNA[,c(4,3,2,1)]
ht.stat.piRNA.sel=ht.stat.piRNA[rownames(pvalue.ht.piRNA.sig),]
ht.stat.piRNA.sel.z=(ht.stat.piRNA.sel - rowMeans(ht.stat.piRNA.sel))/apply(ht.stat.piRNA.sel,1,sd)
ht.stat.piRNA.sel.z=na.omit(ht.stat.piRNA.sel.z)
#
ht.stat.piRNA.sel=ht.stat.piRNA.sel[rownames(ht.stat.piRNA.sel.z),]
head(ht.stat.piRNA.sel)
#
pvalue.ht.piRNA.sig.sel=pvalue.ht.piRNA.sig[rownames(ht.stat.piRNA.sel.z),]

#################
##########
#colore set
min_cor = min(as.vector(ht.stat.piRNA.sel))
max_cor = max(as.vector(ht.stat.piRNA.sel))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.piRNA = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.piRNA.sig.sel$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,0.01,10e-3,10e-4,10e-5,10e-6), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
col_ha_top.piRNA = columnAnnotation(
  kaps.group=colnames(ht.stat.piRNA.sel),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

pirnaht=Heatmap(ht.stat.piRNA.sel, name = "median.of.z.score", 
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
               top_annotation = col_ha_top.piRNA,
               #right_annotation = row_ha.right,
               show_row_dend = F,show_column_dend = F,
               #row_names_side = "left",
               left_annotation = row_ha.left.piRNA,
               column_title=paste("piRNA differentially expressed ", "\n","19 out of 153"),
               column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/piwi/piRNA/stat.20200409/ht.piRNA.19s.2.pdf",width = 5,height = 6)
draw(pirnaht, padding = unit(c(40, 5, 40,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
####
save(stat.piRNAm.sel, normd.pirna, pvalue.ht.piRNA,file="4.model/kaps/piwi/piRNA/stat.20200409/ht.piRNA.19s.RData")
