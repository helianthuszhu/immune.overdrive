#####DNA damage signature
####
load("~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/hypoxia.genesets.score.RData")
hypoxiascore=hypoxia.genesets
head(kaps.td)
#
length(intersect(rownames(kaps.td), rownames(hypoxiascore)))
###
stat.hypoxia=cbind(kaps.td$kaps.group, hypoxiascore[rownames(kaps.td),])
colnames(stat.hypoxia)[1]="kaps.group"
dim(stat.hypoxia)
#
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
stat.hypoxia$kaps.group=factor(stat.hypoxia$kaps.group,levels = c("set4", "set3","set2","set1"))

for (i in 2:ncol(stat.hypoxia)) {
  pdf(paste0("4.model/kaps/hypoxia/",colnames(stat.hypoxia)[i],".kaps.group.pdf"))
  print(ggviolin(stat.hypoxia,x = "kaps.group", y = colnames(stat.hypoxia)[i], fill = "kaps.group",alpha = 1,size = 0.3,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+
          theme_bw()+xlab("kaps.group")+ylab(colnames(stat.hypoxia)[i])+ggtitle(paste0(colnames(stat.hypoxia)[i]))
  )
  dev.off()
}
save(stat.hypoxia,file="4.model/kaps/hypoxia/hypoxia.score.RData")
#############
#############

sel_gmt.hypoxia=read.table("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/hypoxia.signature.txt",header = T,na.strings = c(""),sep = "\t")
dim(sel_gmt.hypoxia)
sets.hypoxia=as.list(sel_gmt.hypoxia)
sets.hypoxia=lapply(sets.hypoxia, function(x) x[!is.na(x)])
sets.hypoxia[1]$Ye.hypoxia.signature


panelM.hypoxia=genestats.orderd[, colnames(genestats.orderd) %in% sets.hypoxia[4]$Hu.signature]
#panelM=genestats.orderd[, colnames(genestats.orderd) %in% subset(markerpanel, category=="CD.8.T.exhuasted")$symbol]
colnames(panelM.hypoxia)
dim(panelM.hypoxia)
#
col_ha.top=genestats.orderd[rownames(panelM),c(1,13,14,17,27)]
class(col_ha.top)
head(col_ha.top)
#
pk1=pheatmap(t(panelM.hypoxia),show_rownames = T,border_color = NA,main = "Th-1.signatures",
             show_colnames = F,fontsize = 8,#na_col = "black",
             #cluster_rows = hc,
             cluster_cols =F,
             annotation_col = col_ha.top,
             #annotation_row = cndi.rep[rownames(ppd),c(2,3)],
             scale="row",
             #color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8), 
             color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             #annotation_legend = T,
             annotation_colors  = list( TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                        TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                        #MSI.status=c("MSI"="#a65628", "MSI"="#ffd92f","MSS"="#8da0cb"),
                                        MSI.status.bin=c("2MSI"="#e7298a","1MSS"="#66a61e"),
                                        kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")#,
                                        #repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                        #repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
             )
)



col_ha.top1=columnAnnotation(z.of.mean.exp.score=anno_lines(col_ha.top$z.of.mean.exp),
                             TE.cluster = col_ha.top$TE.cluster,
                             TE.cluster.agg=col_ha.top$TE.cluster.agg,
                             MSI.status.bin=col_ha.top$MSI.status.bin,
                             z.of.mean.exp=col_ha.top$z.of.mean.exp,
                             kaps.group=col_ha.top$kaps.group,
                             col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                      TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                      MSI.status.bin=c("2MSI"="#e7298a","1MSS"="#66a61e"),
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)

comp=Heatmap(t(scale(panelM.hypoxia)), name = "Th-1.signatures", 
             col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             #col = rev(viridis(10)),border = F,
             show_column_names = F,show_row_names = T,
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 8),
             top_annotation = col_ha.top1#, 
             #right_annotation = row_ha.right,
             #show_row_dend = T,show_column_dend = T,
             #row_names_side = "left",
             #left_annotation = row_ha.left
)