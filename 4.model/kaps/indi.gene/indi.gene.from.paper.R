ggviolin(genestats,
         x = "kaps.group", y = "AZI2", fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")
#+ylab("Global.methylation.level")+ggtitle("Global.methylation.level")
##########

apobec.gene=c("APOBEC1","APOBEC3A","APOBEC3B","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G","APOBEC3H",
                            "MOV10","SAMHD1","AICDA","RNASEL","TREX1","ZC3HAV1")

apobec.gene=read.table("4.model/kaps/indi.gene/genemarker.txt",header = F,sep = "\t")
apobec.gene=read.csv("4.model/kaps/indi.gene/TLRs.csv",header = F,sep = "\t")
apobec.gene=apobec.gene$V1

panelM.hypoxia=genestats.orderd[, colnames(genestats.orderd) %in% apobec.gene]
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

Heatmap(t(scale(panelM.hypoxia)), name = "Th-1.signatures", 
             col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             #col = rev(viridis(10)),border = F,
             show_column_names = F,show_row_names = T,
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             top_annotation = col_ha.top1#, 
             #right_annotation = row_ha.right,
             #show_row_dend = T,show_column_dend = T,
             #row_names_side = "left",
             #left_annotation = row_ha.left
)
##################################
############show the logFC
########
head(kaps.four.degs)
APOBEC.logFC=kaps.four.degs[rownames(kaps.four.degs) %in% apobec.gene,]
APOBEC.logFC=APOBEC.logFC[,grepl("logFC",colnames(APOBEC.logFC))]
APOBEC.logFC
#########

min_cor = min(as.vector(APOBEC.logFC))
max_cor = max(as.vector(APOBEC.logFC))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=100)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
#pdf("4.model/kaps/piwi/piRNA/ht.cor.pdf",width = 6,height = 10)
Heatmap(APOBEC.logFC,
        #col = col.pal_cor,
        cluster_columns = F,
        #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
        #col=rev(viridis(10)),
        col=col.pal_cor,
        #rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        heatmap_legend_param = list(title = "logFC"),
        column_title = paste0("Histone methyltransferase","\n",
                              "logFC compared with the remaining three for each set"),
        column_title_gp = gpar(fontsize = 10),
        column_title_side = "top",
        show_row_dend = F,show_column_dend = F
)


