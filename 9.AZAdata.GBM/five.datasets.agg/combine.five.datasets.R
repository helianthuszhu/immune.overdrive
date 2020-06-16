#################
#################combine four data
#
#GSE80317 GMB 3 cell lines
CB.cell.line.expMatrix.gse80137=res.deg.aza
CB.cell.line.ssgsea.gse80317=stat.hallmarkerscore.aza
CB.cell.line.ssgsea.pvalue.gse80317=pvalue.hm.aza
CB.cell.line.degs.gse.gse80137=res.deg.aza[,1:9]
#GSE51811 CRC HCT116
CB.cell.line.expMatrix.gse51811=expdata.gse51811.neqcnormed
CB.cell.line.ssgsea.gse51811=stat.hallmarkerscore.aza.gse51811.necqnormed
CB.cell.line.ssgsea.pvalue.gse51811=pvalue.hm.aza.gse51811.necqnormed
CB.cell.line.degs.gse51811=res.deg.aza.gse51811
#GSE5816 lungs
CB.cell.line.expMatrix.gse5816=combat_edata.azaexp.GEOdata
CB.cell.line.ssgsea.gse5816=stat.hallmarkerscore.aza.GEOdata
CB.cell.line.ssgsea.pvalue.gse5816=pvalue.hm.aza.GEOdata
CB.cell.line.degs.gse5816=res.deg.aza.GEOdata.GSE5816
#GSE22250 BRCA 7 cell lines
load("9.AZAdata.GBM/GSE22250/combat_edata.azaexp.GEOdata.GSE22250.RData")
dim(combat_edata.azaexp.GEOdata)
#######
CB.cell.line.expMatrix.gse22250=combat_edata.azaexp.GEOdata
CB.cell.line.ssgsea.gse22250=stat.hallmarkerscore.aza.GEOdata
CB.cell.line.ssgsea.pvalue.gse22250=pvalue.hm.aza.GEOdata
CB.cell.line.degs.gse22250=read.csv("9.AZAdata.GBM/GSE22250/degs.genes.AZA.gse22250.BRCA.csv",header = T,row.names = 1)
#GSE41586  CRC HT29
#
CB.cell.line.expMatrix.gse41586=GEO.data_exp.agg
CB.cell.line.ssgsea.gse41586=stat.hallmarkerscore.aza.gse41586.arraydata
CB.cell.line.ssgsea.pvalue.gse41586=pvalue.hm.aza.gse41586.arraydata
CB.cell.line.degs.gse41586=res.deg.aza.gse41586.array
##########
##########
save(CB.cell.line.degs.gse22250,CB.cell.line.degs.gse41586,CB.cell.line.degs.gse51811,CB.cell.line.degs.gse5816,CB.cell.line.degs.gse.gse80137,
     CB.cell.line.expMatrix.gse22250,CB.cell.line.expMatrix.gse41586,CB.cell.line.expMatrix.gse51811,CB.cell.line.expMatrix.gse5816,CB.cell.line.expMatrix.gse80137,
     CB.cell.line.ssgsea.gse22250,CB.cell.line.ssgsea.gse41586,CB.cell.line.ssgsea.gse51811,CB.cell.line.ssgsea.gse5816,CB.cell.line.ssgsea.gse80317,
     CB.cell.line.ssgsea.pvalue.gse22250,CB.cell.line.ssgsea.pvalue.gse41586,CB.cell.line.ssgsea.pvalue.gse51811,CB.cell.line.ssgsea.pvalue.gse5816,CB.cell.line.ssgsea.pvalue.gse80317,
     file = "9.AZAdata.GBM/five.datasets.agg/five.datasets.agg.expMatrix.stat.cell.line.ssgsea.score.pvalue.degs.RData")
##########
##########
##########
##########draw heatmap
##########
head(CB.cell.line.ssgsea.gse5816)

ascore.gse5816=CB.cell.line.ssgsea.gse5816
ascore.gse5816=cbind(time.point=rep("D6",nrow(ascore.gse5816)), GEOid=rep("GSE5816",nrow(ascore.gse5816)), ascore.gse5816)
ascore.gse5816[1:4,1:5]
#
ascore.gse80137=CB.cell.line.ssgsea.gse80317[,-2]
ascore.gse80137=cbind(time.point=rep("D3",nrow(ascore.gse80137)), cell.type=rep("GBM",nrow(ascore.gse80137)),GEOid=rep("GSE80137",nrow(ascore.gse80137)),  ascore.gse80137 )
ascore.gse80137=ascore.gse80137[order(ascore.gse80137$treatment,decreasing = T),]
ascore.gse80137[1:4,1:6]
#
ascore.gse22250=CB.cell.line.ssgsea.gse22250
ascore.gse22250=cbind(time.point=rep("D4",nrow(ascore.gse22250)), cell.type=rep("BRCA",nrow(ascore.gse22250)), GEOid=rep("GSE22250",nrow(ascore.gse22250)),  ascore.gse22250 )
ascore.gse22250=ascore.gse22250[order(ascore.gse22250$treatment,decreasing = T),]
ascore.gse22250[1:4,1:6]
#
ascore.gse51811=CB.cell.line.ssgsea.gse51811[,-c(3,4)]
ascore.gse51811=cbind(cell.line=rep("HCT.116",nrow(ascore.gse51811)), cell.type=rep("CRC",nrow(ascore.gse51811)), GEOid=rep("GSE51811",nrow(ascore.gse51811)),  ascore.gse51811 )
ascore.gse51811=ascore.gse51811[c(1,6,2,7,3,8,4,9,5,10),]
ascore.gse51811[1:4,1:6]
#
ascore.gse41586=CB.cell.line.ssgsea.gse41586[,-2]
ascore.gse41586=cbind(time.point=rep("D5",nrow(ascore.gse41586)), cell.type=rep("CRC",nrow(ascore.gse41586)), GEOid=rep("GSE41586",nrow(ascore.gse41586)),  ascore.gse41586 )
###########combine ssgsea score 
####
ascore.five.CB=rbind(ascore.gse5816, ascore.gse80137,ascore.gse22250,ascore.gse51811,ascore.gse41586)
ascore.five.CB[1:4,1:6]
table(ascore.five.CB$GEOid)
####
mat1=subset(ascore.five.CB, GEOid=="GSE5816")[,-c(1:5)]
mat1=as.data.frame(scale(mat1))
mat1[mat1< -2] <- -2
mat1[mat1>2] <- 2

mat2=subset(ascore.five.CB, GEOid=="GSE80137")[,-c(1:5)]
mat2=as.data.frame(scale(mat2))
mat2[mat2< -2] <- -2
mat2[mat2>2] <- 2

mat3=subset(ascore.five.CB, GEOid=="GSE22250")[,-c(1:5)]
mat3=as.data.frame(scale(mat3))
mat3[mat3< -2] <- -2
mat3[mat3>2] <- 2

mat4=subset(ascore.five.CB, GEOid=="GSE51811")[,-c(1:5)]
mat4=as.data.frame(scale(mat4))
mat4[mat4< -2] <- -2
mat4[mat4>2] <- 2

mat5=subset(ascore.five.CB, GEOid=="GSE41586")[,-c(1:5)]
mat5=as.data.frame(scale(mat5))
mat5[mat5< -2] <- -2
mat5[mat5>2] <- 2

####annotate column in each dataset
#
baseline #ffffe5
D3 #d9f0a3
D4 #addd8e
D5 #78c679
D6 #41ab5d
D14#238443
D24 #006837
D42 #004529

col_ha_top.mat1 = columnAnnotation(
  GEOID=ascore.gse5816$GEOid,
  cell.line=ascore.gse5816$cell.line,
  cell.type=ascore.gse5816$cell.type,
  time.point=ascore.gse5816$time.point,
  treatment=ascore.gse5816$treatment,
  col=list(treatment=c("3High"="#cc4c02","2Low"="#E7B800", "1Control"="#FC4E07"),
           time.point=c("D6"= "#41ab5d"),
           cell.type=c("Breat.cancer"="#8dd3c7","Bronchial.epithelial"="#ffffb3","Colon.cancer"="#a65628","Lung.cancer"="#fb8072"),
           cell.line=c("A549"="#8dd3c7","H1299"="#ffffb3","H157"="#bebada","H1819"="#80b1d3","H1993"="#fdb462","H2347"="#b3de69",
                       "H460"="#fccde5","H526"="#d9d9d9","HBEC2"="#fb8072", "HBEC2.Rep2"="#e31a1c", "HBEC3"="#bc80bd","HBEC3.Rep2"="#6a3d9a", 
                       "HBEC4"="#f1b6da", "HBEC4.Rep2"="#c51b7d", "HCT116"="#ccebc5","MCF7"="#ffed6f"),
           GEOID=c("GSE5816"="#e41a1c")
  ),
  show_annotation_name = F,gp = gpar(col = "white"))
#
col_ha_top.mat2 = columnAnnotation(
  GEOID=ascore.gse80137$GEOid,
  cell.line=ascore.gse80137$cell.line,
  cell.type=ascore.gse80137$cell.type,
  time.point=ascore.gse80137$time.point,
  treatment=ascore.gse80137$treatment,
  col=list(treatment=c("aza"="#E7B800", "untreatment"="#FC4E07"),
           time.point=c("D3"= "#d9f0a3"),
           cell.line=c("LNT.229"="#00AFBB","T98G"="#e41a1c","U.87"="#377eb8"),
           cell.type=c("GBM"="#6a51a3"),
           GEOID=c("GSE80137"="#d95f02")
  ),
  show_annotation_name = F,gp = gpar(col = "white"))
#
col_ha_top.mat3= columnAnnotation(
  GEOID=ascore.gse22250$GEOid,
  cell.line=ascore.gse22250$cell.line,
  cell.type=ascore.gse22250$cell.type,
  time.point=ascore.gse22250$time.point,
  treatment=ascore.gse22250$treatment,
  col=list(treatment=c("5-AZA"="#E7B800", "WT"="#FC4E07"),
           cell.line=c("MCF.7"="#1b9e77","T47D"="#d95f02","SKBR3"="#7570b3","BT20"="#e7298a","MDA.MB.231"="#66a61e","MDA.MB.361"="#e6ab02","ZR.75.1"="#a6761d"),
           time.point=c("D4"= "#addd8e"),
           cell.type=c("BRCA"="#67a9cf"),
           GEOID=c("GSE22250"="#7570b3")
  ),
  show_annotation_name = F,gp = gpar(col = "white"))
#
col_ha_top.mat4 = columnAnnotation(
  GEOID=ascore.gse51811$GEOid,
  cell.line=ascore.gse51811$cell.line,
  cell.type=ascore.gse51811$cell.type,
  time.point=ascore.gse51811$time.point,
  treatment=ascore.gse51811$treatment,
  col=list(treatment=c("aza"="#E7B800", "untreatment"="#FC4E07"),
           time.point=c("baseline"="#c51b7d","D5"="#78c679","D14"="#238443","D24"="#006837","D42"="#004529"),
           cell.type=c("CRC"="#a65628"),
           cell.line=c("HCT.116"="#ccebc5"),
           GEOID=c("GSE51811"="#e7298a")
  ),
  show_annotation_name = F,gp = gpar(col = "white"))
#
col_ha_top.mat5 = columnAnnotation(
  GEOID=ascore.gse41586$GEOid,
  cell.line=ascore.gse41586$cell.line,
  cell.type=ascore.gse41586$cell.type,
  time.point=ascore.gse41586$time.point,
  treatment=ascore.gse41586$treatment,
  col=list(treatment=c("High"="#662506","Low"="#E7B800", "Control"="#FC4E07"),
           cell.line=c("HT29"="#8dd3c7"),
           cell.type=c("CRC"="#a65628"),
           time.point=c("D5"= "#78c679"),
           GEOID=c("GSE41586"="#66a61e")
  ),
  show_annotation_name = TRUE,gp = gpar(col = "white"))
######p value annotation
row_ha.left.mat1 = rowAnnotation(kruskal.pvalue=as.numeric(CB.cell.line.ssgsea.pvalue.gse5816$p.value),
                               col=list(
                                 kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                           c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                               ),show_annotation_name = FALSE)
row_ha.left.mat2 = rowAnnotation(kruskal.pvalue=as.numeric(CB.cell.line.ssgsea.pvalue.gse80317$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
row_ha.left.mat3 = rowAnnotation(kruskal.pvalue=as.numeric(CB.cell.line.ssgsea.pvalue.gse22250$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
row_ha.left.mat4 = rowAnnotation(kruskal.pvalue=as.numeric(CB.cell.line.ssgsea.pvalue.gse51811$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
row_ha.left.mat5 = rowAnnotation(kruskal.pvalue=as.numeric(CB.cell.line.ssgsea.pvalue.gse41586$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
#
#color setting
#mat1
min_cor = min(as.vector(mat1))
max_cor = max(as.vector(mat1))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor.mat1 = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#mat2
min_cor = min(as.vector(mat2))
max_cor = max(as.vector(mat2))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor.mat2 = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#mat3
min_cor = min(as.vector(mat3))
max_cor = max(as.vector(mat3))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor.mat3 = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))


ht1=Heatmap(t(mat1), name = "GSE5816", 
            #col = viridis(10),
            col=col.pal_cor.mat1,
            #width = unit(4, "cm"),height = unit(12, "cm"),
            border = F,
        show_column_names = F,show_row_names = F,cluster_columns = F,
        row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
        top_annotation = col_ha_top.mat1,
        #right_annotation = row_ha.right,
        show_row_dend = F,show_column_dend = F
        #row_names_side = "left",
        #left_annotation = row_ha.left.mat1
        #column_title=paste0("hall.marker.genesets"," GSE5816 lung cancer"),
        #column_title_gp = gpar(fontsize = 8)
)
ht2=Heatmap(t(mat2), name = "GSE80137", 
            #col = viridis(10),
            col=col.pal_cor.mat2,
            #width = unit(2, "cm"),height = unit(12, "cm"),
            border = F,
            show_column_names = F,show_row_names = F,cluster_columns = F,
            row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
            top_annotation = col_ha_top.mat2,
            #right_annotation = row_ha.right,
            show_row_dend = F,show_column_dend = F#,
            #row_names_side = "left",
            #left_annotation = row_ha.left.mat2
            #column_title=paste0("hall.marker.genesets"," GSE80137 GBM"),
            #column_title_gp = gpar(fontsize = 8)
)
ht3=Heatmap(t(mat3), name = "GSE22250", 
            #col = viridis(10),
            col=col.pal_cor.mat3,
            #width = unit(2, "cm"),height = unit(12, "cm"),
            border = F,
            show_column_names = F,show_row_names = T,cluster_columns = F,
            row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
            top_annotation = col_ha_top.mat3,
            #right_annotation = row_ha.right,
            show_row_dend = F,show_column_dend = F#,
            #row_names_side = "left",
            #left_annotation = row_ha.left.mat3
            #column_title=paste0("hall.marker.genesets"," GSE22250 BRCA"),
            #column_title_gp = gpar(fontsize = 8)
)
ht4=Heatmap(t(mat4), name = "GSE51811", col = viridis(10),
            #width = unit(2, "cm"),height = unit(12, "cm"),
            border = F,
            show_column_names = F,show_row_names = F,cluster_columns = F,
            row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
            top_annotation = col_ha_top.mat4,
            #right_annotation = row_ha.right,
            show_row_dend = F,show_column_dend = F#,
            #row_names_side = "left",
            #left_annotation = row_ha.left.mat4
            #column_title=paste0("hall.marker.genesets"," GSE51811 CRC"),
            #column_title_gp = gpar(fontsize = 8)
)
ht5=Heatmap(t(mat5), name = "GSE41586", col = viridis(10),
            #width = unit(2, "cm"),height = unit(12, "cm"),
            border = F,
            show_column_names = F,show_row_names = T,cluster_columns = F,
            row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
            top_annotation = col_ha_top.mat5,
            #right_annotation = row_ha.right,
            show_row_dend = F,show_column_dend = F#,
            #row_names_side = "left",
            #left_annotation = row_ha.left.mat5
            #column_title=paste0("hall.marker.genesets"," GSE41586 CRC"),
            #column_title_gp = gpar(fontsize = 8)
)

#row_ha.right = rowAnnotation(foo2 = anno_mark(at = c(1:4), labels = rownames(type4.mat)[1:4]))
pcb=ht1+ht2+ht3+ht4+ht5
#
pdf("9.AZAdata.GBM/five.datasets.agg/ht.cb.five.pdf",width = 12,height = 9)
draw(pcb, padding = unit(c(20, 15, 20,15), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#########################################################
#########################################################
######check the DEGs among five datasets
######
head(CB.cell.line.degs.gse.gse80137)
head(CB.cell.line.degs.gse5816)
head(CB.cell.line.degs.gse22250)
head(CB.cell.line.degs.gse51811)
head(CB.cell.line.degs.gse41586)
#########################
length(intersect(CB.cell.line.degs.gse.gse80137$Gene.name, rownames(CB.cell.line.degs.gse22250)))
length(intersect(rownames(CB.cell.line.degs.gse51811), rownames(CB.cell.line.degs.gse22250)))
length(unique(CB.cell.line.degs.gse.gse80137$Gene.name))
##########regular way
degs.sel.gse80137=subset(CB.cell.line.degs.gse.gse80137, adj.P.Val <0.05 &  abs(logFC) >=1 )
degs.sel.gse5816=subset(CB.cell.line.degs.gse5816, adj.P.Val <0.05 &  abs(logFC) >=1 )
degs.sel.gse22250=subset(CB.cell.line.degs.gse22250, adj.P.Val <0.05 &  abs(logFC) >=1 )
degs.sel.gse51811=subset(CB.cell.line.degs.gse51811, adj.P.Val <0.05 &  abs(logFC) >=1 )
degs.sel.gse41586=subset(CB.cell.line.degs.gse41586, adj.P.Val <0.05 &  abs(logFC) >=1 )
#########
write.csv()
dim(degs.sel.gse80137)
######################
######################
