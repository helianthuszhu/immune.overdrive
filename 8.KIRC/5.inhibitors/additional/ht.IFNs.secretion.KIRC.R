##########################################KIRC
#####additional TE inhibitors
#####
###get kaps group kirc
head(clin.heat.kirc.var)
table(clin.heat.kirc.var$kaps.group.kirc)
###get IFNs secretion ssgsea score
load("8.KIRC/5.inhibitors/additional/genesets.IFNs.kircfh.RData")
head(genesets.IFNs..kircfh)
###
CB.data.self.kaps.IFNs=cbind(kaps.td[self.id.kaps.IFNs,c(17,38)], IFNtypes.score[self.id.kaps.IFNs,])

CB.data.self.kaps.IFNs.KIRC=cbind(kaps.group=clin.heat.kirc.var$kaps.group.kirc, genesets.IFNs..kircfh[rownames(clin.heat.kirc.var),])
####draw IFN secrection difference among four groups
####
ht.IFNs.KIRC=data.frame(kaps.group=CB.data.self.kaps.IFNs.KIRC$kaps.group,
                   IFN_ALPHA_SECRETION=CB.data.self.kaps.IFNs.KIRC$GO_INTERFERON_ALPHA_SECRETION,
                   IFN_BETA_SECRETION=CB.data.self.kaps.IFNs.KIRC$GO_INTERFERON_BETA_SECRETION,
                   IFN_GAMMA_SECRETION=CB.data.self.kaps.IFNs.KIRC$GO_INTERFERON_GAMMA_SECRETION,
                   IFN_ALPHA_PRODUCTION=CB.data.self.kaps.IFNs.KIRC$GO_INTERFERON_ALPHA_PRODUCTION,
                   IFN_BETA_PRODUCTION=CB.data.self.kaps.IFNs.KIRC$GO_INTERFERON_BETA_PRODUCTION,
                   IFN_GAMMA_PRODUCTION=CB.data.self.kaps.IFNs.KIRC$GO_INTERFERON_GAMMA_PRODUCTION
)
rownames(ht.IFNs.KIRC)=rownames(CB.data.self.kaps.IFNs.KIRC)
head(ht.IFNs.KIRC)
############normalize score
ht.IFNs.KIRC=cbind(kaps.group=ht.IFNs.KIRC$kaps.group, as.data.frame(scale(ht.IFNs.KIRC[,-1])))
ht.IFNs.KIRC$kaps.group=as.factor(ht.IFNs.KIRC$kaps.group)
ht.IFNs.KIRC[1:4,1:4]
############
######pvalue
datalist <- list()
for(i in names(ht.IFNs.KIRC[,2:ncol(ht.IFNs.KIRC)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = ht.IFNs.KIRC)
}
pvalue.ht.IFNs.KIRC=do.call(rbind, datalist)
pvalue.ht.IFNs.KIRC=as.data.frame(pvalue.ht.IFNs.KIRC)
#
######matrix
ht.stat.IFNs.KIRC= ht.IFNs.KIRC %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.IFNs.KIRC=as.data.frame(ht.stat.IFNs.KIRC)
rownames(ht.stat.IFNs.KIRC)=ht.stat.IFNs.KIRC$kaps.group
ht.stat.IFNs.KIRC=ht.stat.IFNs.KIRC[,-1]
ht.stat.IFNs.KIRC=as.data.frame(t(ht.stat.IFNs.KIRC))
ht.stat.IFNs.KIRC=ht.stat.IFNs.KIRC[,c(4:1)]
######draw heatmap
######

#colore set
min_cor = min(as.vector(ht.stat.IFNs.KIRC))
max_cor = max(as.vector(ht.stat.IFNs.KIRC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))

#
row_ha.left.IFNs.KIRC= rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.IFNs.KIRC$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
col_ha_top.IFNs.KIRC = columnAnnotation(
  kaps.group=colnames(ht.stat.IFNs.KIRC),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

ht.IFNs.KIRC=Heatmap(ht.stat.IFNs.KIRC, name = "median.of.z.score", 
                #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                #col = rev(viridis(10)),
                width = unit(2, "cm"),
                height = unit(3, "cm"),
                border = F,
                col=col.pal_cor,
                show_column_names = T,show_row_names = T,
                cluster_columns = F,cluster_rows = F,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 5),
                top_annotation = col_ha_top.IFNs.KIRC,
                #right_annotation = row_ha.right,
                show_row_dend = F,show_column_dend = F,
                #row_names_side = "left",
                left_annotation = row_ha.left.IFNs.KIRC,
                column_title="IFNs.secretion.KIRC",
                column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/5.inhibitors/additional/ht.IFNs.secretion.KIRC.pdf",width = 5,height = 6)
draw(ht.IFNs.KIRC, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
##################################
#########################################
save(ht.stat.IFNs.KIRC,pvalue.ht.IFNs.KIRC, file="8.KIRC/5.inhibitors/additional/ht.IFNs.secretion.KIRC.RData")
