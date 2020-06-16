#IFN three types response signatures
####
load("4.model/kaps/marker.TEs.inhibition/additional/IFNtype.response.score.RData")
########
head(kaps.td)
head(IFNtypes.score)
self.id.kaps.IFNs=intersect(rownames(IFNtypes.score), rownames(kaps.td))
#####
CB.data.self.kaps.IFNs=cbind(kaps.td[self.id.kaps.IFNs,c(17,38)], IFNtypes.score[self.id.kaps.IFNs,])
colnames(CB.data.self.kaps.IFNs)
dim(CB.data.self.kaps.IFNs)
######
index.self.kaps.IFNs=colnames(CB.data.self.kaps.IFNs)[-c(1:2)]
#my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),c("set3", "set4"))
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
#CB.data.self.kaps$kaps.group.agg=factor(CB.data.self.kaps$kaps.group.agg,levels = c("set4","set3","set_1&2"))
CB.data.self.kaps.IFNs$kaps.group=factor(CB.data.self.kaps.IFNs$kaps.group,levels = c("set4","set3","set2","set1"))

for (i in 1:length(index.self.kaps.IFNs)) {
  pdf(file = paste0("4.model/kaps/marker.TEs.inhibition/additional/indi.com/",index.self.kaps.IFNs[i],".IFNs.kaps.four.pdf"),height  = 5,width = 4)
  print(ggviolin(CB.data.self.kaps.IFNs,x = "kaps.group", y =index.self.kaps.IFNs[i] , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+xlab("kaps.group")+ylab(index.self.kaps.IFNs[i]) +
          ggtitle(paste0(index.self.kaps.IFNs[i],"CRC"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
############
############
####draw IFN secrection difference among four groups
####
ht.IFNs=data.frame(kaps.group=CB.data.self.kaps.IFNs$kaps.group,
                   IFN_ALPHA_SECRETION=CB.data.self.kaps.IFNs$GO_INTERFERON_ALPHA_SECRETION,
                   IFN_BETA_SECRETION=CB.data.self.kaps.IFNs$GO_INTERFERON_BETA_SECRETION,
                   IFN_GAMMA_SECRETION=CB.data.self.kaps.IFNs$GO_INTERFERON_GAMMA_SECRETION,
                   IFN_ALPHA_PRODUCTION=CB.data.self.kaps.IFNs$GO_INTERFERON_ALPHA_PRODUCTION,
                   IFN_BETA_PRODUCTION=CB.data.self.kaps.IFNs$GO_INTERFERON_BETA_PRODUCTION,
                   IFN_GAMMA_PRODUCTION=CB.data.self.kaps.IFNs$GO_INTERFERON_GAMMA_PRODUCTION
                   )
rownames(ht.IFNs)=rownames(CB.data.self.kaps.IFNs)
head(ht.IFNs)
############normalize score
ht.IFNs=cbind(kaps.group=ht.IFNs$kaps.group, as.data.frame(scale(ht.IFNs[,-1])))
############
######pvalue
datalist <- list()
for(i in names(ht.IFNs[,2:ncol(ht.IFNs)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = ht.IFNs)
}
pvalue.ht.IFNs=do.call(rbind, datalist)
pvalue.ht.IFNs=as.data.frame(pvalue.ht.IFNs)
#
######matrix
ht.stat.IFNs= ht.IFNs %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.IFNs=as.data.frame(ht.stat.IFNs)
rownames(ht.stat.IFNs)=ht.stat.IFNs$kaps.group
ht.stat.IFNs=ht.stat.IFNs[,-1]
ht.stat.IFNs=as.data.frame(t(ht.stat.IFNs))
######draw heatmap
######

#colore set
min_cor = min(as.vector(ht.stat.IFNs))
max_cor = max(as.vector(ht.stat.IFNs))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.IFNs = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.ht.IFNs$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
col_ha_top.IFNs = columnAnnotation(
  kaps.group=colnames(ht.stat.IFNs),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

ht.IFNs=Heatmap(ht.stat.IFNs, name = "median.of.z.score", 
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
             top_annotation = col_ha_top.IFNs,
             #right_annotation = row_ha.right,
             show_row_dend = F,show_column_dend = F,
             #row_names_side = "left",
             left_annotation = row_ha.left.IFNs,
             column_title="IFNs.secretion",
             column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/marker.TEs.inhibition/additional/ht/ht.IFNs.secretion.pdf",width = 5,height = 6)
draw(ht.IFNs, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
################
########
save(ht.stat.IFNs,pvalue.ht.IFNs, file="4.model/kaps/marker.TEs.inhibition/additional/ht/ht.IFNs.secretion.RData")
