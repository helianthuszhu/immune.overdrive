#################
#######
isgs.IFN=read.table("4.model/kaps/marker.TEs.inhibition/additional/TEs.inhibition.markers.ISG.txt",header = T,sep = "\t")
isgs.IFN=isgs.IFN[!(duplicated(isgs.IFN$symbol)),]
rownames(isgs.IFN)=isgs.IFN$symbol

head(isgs.IFN)
dim(isgs.IFN)
table(isgs.IFN$TEs.suppression.mechanism)
#######
isgs.IFN.expMatrix=genestats.orderd[, colnames(genestats.orderd) %in% rownames(isgs.IFN)]
isgs.IFN.expMatrix=as.data.frame(scale(isgs.IFN.expMatrix))
isgs.IFN.expMatrix=cbind(kaps.group=genestats.orderd$kaps.group, isgs.IFN.expMatrix[rownames(genestats.orderd),])
isgs.IFN.expMatrix[1:3,1:4]
######
######pvalue
datalist <- list()
for(i in names(isgs.IFN.expMatrix[,2:ncol(isgs.IFN.expMatrix)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = isgs.IFN.expMatrix,na.action=na.omit)
}
pvalue.ht.isgs=do.call(rbind, datalist)
pvalue.ht.isgs=as.data.frame(pvalue.ht.isgs)
pvalue.ht.isgs.sig=subset(pvalue.ht.isgs, p.value<0.05)
#####
######matrix
ht.stat.isgs= isgs.IFN.expMatrix %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.isgs=as.data.frame(ht.stat.isgs)
rownames(ht.stat.isgs)=ht.stat.isgs$kaps.group
ht.stat.isgs=ht.stat.isgs[,-1]
ht.stat.isgs=as.data.frame(t(ht.stat.isgs))
ht.stat.isgs=ht.stat.isgs[,c(4:1)]
#
ht.stat.isgs.sig=ht.stat.isgs[rownames(pvalue.ht.isgs.sig),]
#####draw
#colore set
min_cor = min(as.vector(ht.stat.isgs.sig))
max_cor = max(as.vector(ht.stat.isgs.sig))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.isgs = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.isgs.sig$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
col_ha_top.isgs = columnAnnotation(
  kaps.group=colnames(ht.stat.isgs.sig),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

ht.isgs=Heatmap(ht.stat.isgs.sig, name = "median.of.z.score", 
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
                top_annotation = col_ha_top.isgs,
                #right_annotation = row_ha.right,
                show_row_dend = F,show_column_dend = F,
                #row_names_side = "left",
                left_annotation = row_ha.left.isgs,
                column_title="interferon-stimulated genes (ISGs)",
                column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/additional/ht.interferon.stimulated.genes.only.significant.pdf",width = 5,height = 7)
draw(ht.isgs, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###############
save(ht.stat.isgs,pvalue.ht.isgs, file="4.model/kaps/marker.TEs.inhibition/additional/ht.interferon.stimulated.genes.RData")
