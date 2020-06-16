######hall marker genesets
#######
load("4.model/kaps/hall.marker/genesets.hallmarker.score.RData")
length(intersect(rownames(genesets.hall.marker.score), rownames(kaps.td)))

stat.hallm=cbind(kaps.td$kaps.group,genesets.hall.marker.score[rownames(kaps.td),])
colnames(stat.hallm)[1]="kaps.group"
head(stat.hallm)
stat.hallm$kaps.group=factor(stat.hallm$kaps.group,levels = c("set4", "set3","set2","set1"))

#
for (i in 2:ncol(stat.hallm)) {
  pdf(paste0("4.model/kaps/hall.marker/",colnames(stat.hallm)[i],".kaps.group.pdf"),width = 4,height = 6)
  print(ggviolin(stat.hallm,x = "kaps.group", y = colnames(stat.hallm)[i], fill = "kaps.group",alpha = 1,size = 0.3,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group.agg")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
          theme_bw()+xlab("kaps.group")+ylab(colnames(stat.hallm)[i])+ggtitle(paste0(colnames(stat.hallm)[i]))
  )
  dev.off()
}
####################
####################
###################
for (i in 2:ncol(stat.hallm)) {
  pdf(paste0("4.model/kaps/hall.marker/stat2020408/indi.compare/",colnames(stat.hallm)[i],".kaps.group.pdf"),height  = 5,width = 4)
  print(ggviolin(stat.hallm,x = "kaps.group", y = colnames(stat.hallm)[i], fill = "kaps.group",alpha = 1,size = 0.3,width = 1,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
          xlab("kaps.group")+ylab(colnames(stat.hallm)[i])+ggtitle(paste0(colnames(stat.hallm)[i]))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
######
for (i in 2:ncol(stat.hallm)) {
  pdf(paste0("4.model/kaps/hall.marker/stat2020408/total.compare/",colnames(stat.hallm)[i],".kaps.group.pdf"),height  = 5,width = 4)
  print(ggviolin(stat.hallm,x = "kaps.group", y = colnames(stat.hallm)[i], fill = "kaps.group",alpha = 1,size = 0.3,width = 1,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(label = "p.signif",method = "kruskal.test")+
          xlab("kaps.group")+ylab(colnames(stat.hallm)[i])+ggtitle(paste0(colnames(stat.hallm)[i]))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  )
  dev.off()
}
##############
##################draw heatmap
#####
ht.stat.hallm=cbind(stat.hallm$kaps.group, as.data.frame(scale(stat.hallm[,-1])))
colnames(ht.stat.hallm)[1]="kaps.group"
head(ht.stat.hallm)

ht.hallM=ht.stat.hallm %>% group_by(kaps.group) %>% summarise_all(median)
ht.hallM=as.data.frame(ht.hallM)
rownames(ht.hallM)=ht.hallM$kaps.group
ht.hallM=ht.hallM[,-1]
ht.hallM=as.data.frame(t(ht.hallM))
head(ht.hallM)

#######calculate p value
#######

datalist <- list()
for(i in names(stat.hallm[,2:ncol(stat.hallm)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = stat.hallm)
}
pvalue.hm=do.call(rbind, datalist)
pvalue.hm=as.data.frame(pvalue.hm)
pvalue.hm.sig=subset(pvalue.hm, p.value< 0.05)

#
ht.hallM.sig=ht.hallM[rownames(pvalue.hm.sig),]
min_cor = min(as.vector(ht.hallM.sig))
max_cor = max(as.vector(ht.hallM.sig))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))



###
row_ha.left.hm = rowAnnotation(kruskal.pvalue=as.numeric(pvalue.hm.sig$p.value),
                               col=list(
                                 kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                               ),show_annotation_name = FALSE)


col_ha_top.hm = columnAnnotation(
                                   kaps.group=colnames(ht.hallM.sig),
                                   col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
                                   show_annotation_name = FALSE,gp = gpar(col = "black"))

####
hmht=Heatmap(ht.hallM.sig, name = "hall.marker", 
             #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             #col = rev(viridis(10)),
             width = unit(2, "cm"),
             height = unit(12, "cm"),
             border = F,
             col=col.pal_cor,
             show_column_names = T,show_row_names = T,
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             top_annotation = col_ha_top.hm,
             #right_annotation = row_ha.right,
             show_row_dend = F,show_column_dend = F,
             #row_names_side = "left",
             left_annotation = row_ha.left.hm,
             column_title=paste0("hall.marker.genesets","\n","(significant 32 out of 50)"),
             column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/hall.marker/stat2020408/ht.hallM.pdf",width = 5,height = 6)
draw(hmht, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()

save(ht.hallM,pvalue.hm,pvalue.hm.sig,ht.hallM.sig,file="4.model/kaps/hall.marker/stat2020408/ht.hallM.RData")
