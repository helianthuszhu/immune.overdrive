#################
#######
isgs.IFN=read.table("4.model/kaps/marker.TEs.inhibition/additional/TEs.inhibition.markers.ISG.txt",header = T,sep = "\t")
isgs.IFN=isgs.IFN[!(duplicated(isgs.IFN$symbol)),]
rownames(isgs.IFN)=isgs.IFN$symbol

head(isgs.IFN)
dim(isgs.IFN)
table(isgs.IFN$TEs.suppression.mechanism)
#######
genestats.KIRC[1:4,1:4]

isgs.IFN.expMatrix.KIRC=genestats.KIRC[, colnames(genestats.KIRC) %in% rownames(isgs.IFN)]
isgs.IFN.expMatrix.KIRC=as.data.frame(scale(isgs.IFN.expMatrix.KIRC))
isgs.IFN.expMatrix.KIRC=cbind(kaps.group=genestats.KIRC$kaps.group, isgs.IFN.expMatrix.KIRC[rownames(genestats.KIRC),])
isgs.IFN.expMatrix.KIRC$kaps.group=as.factor(isgs.IFN.expMatrix.KIRC$kaps.group)
isgs.IFN.expMatrix.KIRC[1:3,1:4]
######
######pvalue
datalist <- list()
for(i in names(isgs.IFN.expMatrix.KIRC[,2:ncol(isgs.IFN.expMatrix.KIRC)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = isgs.IFN.expMatrix.KIRC,na.action=na.omit)
}
pvalue.ht.isgs.KIRC=do.call(rbind, datalist)
pvalue.ht.isgs.KIRC=as.data.frame(pvalue.ht.isgs.KIRC)
pvalue.ht.isgs.sig.KIRC=subset(pvalue.ht.isgs.KIRC, p.value<0.05)
#####
######matrix
ht.stat.isgs.KIRC= isgs.IFN.expMatrix.KIRC %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.isgs.KIRC=as.data.frame(ht.stat.isgs.KIRC)
rownames(ht.stat.isgs.KIRC)=ht.stat.isgs.KIRC$kaps.group
ht.stat.isgs.KIRC=ht.stat.isgs.KIRC[,-1]
ht.stat.isgs.KIRC=as.data.frame(t(ht.stat.isgs.KIRC))
ht.stat.isgs.KIRC=ht.stat.isgs.KIRC[,c(4:1)]
ht.stat.isgs.KIRC[1:4,1:4]
#
ht.stat.isgs.KIRC=ht.stat.isgs.KIRC[rownames(pvalue.ht.isgs.sig.KIRC),]
#####draw
#colore set
min_cor = min(as.vector(ht.stat.isgs.KIRC))
max_cor = max(as.vector(ht.stat.isgs.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))

#
row_ha.left.isgs.KIRC = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.isgs.sig.KIRC$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
col_ha_top.isgs.KIRC = columnAnnotation(
  kaps.group=colnames(ht.stat.isgs.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

ht.isgs.KIRC=Heatmap(ht.stat.isgs.KIRC, name = "median.of.z.score", 
                #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                #col = rev(viridis(10)),
                width = unit(2, "cm"),
                height = unit(12, "cm"),
                border = F,
                col=col.pal_cor,
                show_column_names = T,show_row_names = T,
                cluster_columns = F,cluster_rows = F,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 5),
                top_annotation = col_ha_top.isgs.KIRC,
                #right_annotation = row_ha.right,
                show_row_dend = F,show_column_dend = F,
                #row_names_side = "left",
                left_annotation = row_ha.left.isgs.KIRC,
                column_title="interferon-stimulated genes (ISGs).KIRC",
                column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/additional/ht.interferon.stimulated.genes.significant.KIRC.pdf",width = 5,height = 7)
draw(ht.isgs.KIRC, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###############
save(ht.stat.isgs.KIRC,pvalue.ht.isgs.KIRC, file="8.KIRC/5.inhibitors/additional/ht.interferon.stimulated.genes.all.KIRC.RData")
