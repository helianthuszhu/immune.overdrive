dim(CB.data.self.kaps)
dim(GEP.stats.kaps)
dim(CB.data.ipres.kaps)
###########
#########self genesets
#########
#########
colnames(CB.data.self.kaps)
ht.infiltration=data.frame(leukocyte.infiltration=CB.data.self.kaps$leukocyte.infiltration,
                           IFN.gamma.signature.18.genes.Ayers.etal=CB.data.self.kaps$IFN.gamma.signature.18.genes.Ayers.etal,
                           hot.tumor.signautre=CB.data.self.kaps$hot.tumor.signautre,
                           MHC.II=CB.data.self.kaps$MHC.II,
                           MHC.I=CB.data.self.kaps$MHC.I,
                           Tcell_infiltration_1=CB.data.self.kaps$Tcell_infiltration_1,
                           CD4.T.cells.exhuasted=CB.data.self.kaps$CD4.T.cells.exhuasted.zhang.zeminscRNA,
                           CD8.T.cells.exhuasted=CB.data.self.kaps$CD8.T.cells.exhuasted.zhang.zeminscRNA,
                           TcClassII_score=CB.data.self.kaps$TcClassII_score,
                           TAMsurr_score=CB.data.self.kaps$TAMsurr_score,
                           TAMsurr_TcClassII_ratio=CB.data.self.kaps$TAMsurr_TcClassII_ratio,
                      #
                      TGFB.response.wolf=CB.data.self.kaps$TGFB.response.wolf,
                      up.ECM.signature=CB.data.self.kaps$up.ECM.signature,
                      WESTON_VEGFA_TARGETS_12HR=CB.data.self.kaps$WESTON_VEGFA_TARGETS_12HR,
                      EMT_UP=CB.data.self.kaps$EMT_UP
                      )
rownames(ht.infiltration)=rownames(CB.data.self.kaps)
ht.infiltration=cbind(GEP.stats.kaps[rownames(ht.infiltration),]$GEP, 
                      CB.data.ipres.kaps[rownames(ht.infiltration),]$z.mean.IPRES,
                      stat.th1sig[rownames(ht.infiltration),]$mean.Th1.sig,
                      ht.infiltration)

colnames(ht.infiltration)[1:3]=c("GEP","IPRES","Th1.signature")

ht.infiltration=as.data.frame(scale(ht.infiltration))
ht.infiltration=cbind(CB.data.self.kaps$kaps.group,ht.infiltration)
colnames(ht.infiltration)[1]="kaps.group"
head(ht.infiltration)
######pvalue
datalist <- list()
for(i in names(ht.infiltration[,2:ncol(ht.infiltration)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = ht.infiltration)
}
pvalue.ht.inf=do.call(rbind, datalist)
pvalue.ht.inf=as.data.frame(pvalue.ht.inf)
######matrix
ht.stat.inf= ht.infiltration %>% group_by(kaps.group) %>% summarise_all(median)
ht.stat.inf=as.data.frame(ht.stat.inf)
rownames(ht.stat.inf)=ht.stat.inf$kaps.group
ht.stat.inf=ht.stat.inf[,-1]
ht.stat.inf=as.data.frame(t(ht.stat.inf))
#############
######in pan
########
##########
dim(CB.data.pan.kaps)
ht.BCR=data.frame(TCR.Shannon=CB.data.pan.kaps$TCR.Shannon,
                  TCR.Richness=CB.data.pan.kaps$TCR.Richness,
                  TCR.Evenness=CB.data.pan.kaps$TCR.Evenness,
                  BCR.Shannon=CB.data.pan.kaps$BCR.Shannon,
                  BCR.Richness=CB.data.pan.kaps$BCR.Richness,
                  BCR.Evenness=CB.data.pan.kaps$BCR.Evenness,
                  Wound.Healing=CB.data.pan.kaps$Wound.Healing,
                  Lymphocyte.Infiltration.Signature.Score=CB.data.pan.kaps$Lymphocyte.Infiltration.Signature.Score,
                  Aneuploidy.Score=CB.data.pan.kaps$Aneuploidy.Score,
                  Homologous.Recombination.Defects=CB.data.pan.kaps$Homologous.Recombination.Defects,
                  Number.of.Segments=CB.data.pan.kaps$Number.of.Segments,
                  Fraction.Altered=CB.data.pan.kaps$Fraction.Altered,
                  SNV.Neoantigens=CB.data.pan.kaps$SNV.Neoantigens,
                  Indel.Neoantigens=CB.data.pan.kaps$Indel.Neoantigens,
                  Intratumor.Heterogeneity=CB.data.pan.kaps$Intratumor.Heterogeneity
                  )
rownames(ht.BCR)=rownames(CB.data.pan.kaps)
#######
ht.BCR=as.data.frame(scale(ht.BCR))
ht.BCR=cbind(CB.data.pan.kaps$kaps.group,ht.BCR)
colnames(ht.BCR)[1]="kaps.group"
head(ht.BCR)
######pvalue
datalist <- list()
for(i in names(ht.BCR[,2:ncol(ht.BCR)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group")), data = ht.BCR)
}
pvalue.ht.BCR=do.call(rbind, datalist)
pvalue.ht.BCR=as.data.frame(pvalue.ht.BCR)
######matrix
ht.stat.BCR= ht.BCR %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
ht.stat.BCR=as.data.frame(ht.stat.BCR)
rownames(ht.stat.BCR)=ht.stat.BCR$kaps.group
ht.stat.BCR=ht.stat.BCR[,-1]
ht.stat.BCR=as.data.frame(t(ht.stat.BCR))
head(ht.stat.BCR)
###########
###########combine matrix and p valie
####
ht.cb.matrix=rbind(ht.stat.inf, ht.stat.BCR)
pvalue.cb=rbind(pvalue.ht.inf,pvalue.ht.BCR)
save(ht.cb.matrix,pvalue.cb,file="4.model/kaps/immune.self.geneset/stat.20200407/ht.show/ht.immune.infiltration.RData")
#############draw now
###innmune infilitation
###
rownames(ht.cb.matrix)
infi.matrix=ht.cb.matrix[c(1,26,4:6,3,9:14,7:8),]
#infi.matrix=(infi.matrix - rowMeans(infi.matrix))/apply(infi.matrix,1,sd)
infi.pvalue=pvalue.cb[rownames(infi.matrix),]
rownames(infi.matrix)
#colore set
min_cor = min(as.vector(infi.matrix))
max_cor = max(as.vector(infi.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.infi = rowAnnotation(kruskal.pvalue=as.numeric(infi.pvalue$p.value),
                               col=list(
                                 kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                 c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                               ),show_annotation_name = FALSE)
col_ha_top.infi = columnAnnotation(
  kaps.group=colnames(infi.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

inht=Heatmap(infi.matrix, name = "mean.of.z.score", 
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
        top_annotation = col_ha_top.infi,
        #right_annotation = row_ha.right,
        show_row_dend = F,show_column_dend = F,
        #row_names_side = "left",
        left_annotation = row_ha.left.infi,
        column_title="Immune.infiltration.signature",
        column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/immune.self.geneset/stat.20200407/ht.show/ht.immune.infiltration.pdf",width = 5,height = 6)
draw(inht, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###########
###########BCR
###########
rownames(ht.cb.matrix)
BCR.matrix=ht.cb.matrix[c(18,21,19,22,20,23),]
BCR.pvalue=pvalue.cb[rownames(BCR.matrix),]
rownames(BCR.matrix)
#colore set
min_cor = min(as.vector(BCR.matrix))
max_cor = max(as.vector(BCR.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.BCR = rowAnnotation(kruskal.pvalue=as.numeric(BCR.pvalue$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#4d4d4d","#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354"))
                                 ),show_annotation_name = FALSE)
col_ha_top.BCR = columnAnnotation(
  kaps.group=colnames(BCR.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

BCRht=Heatmap(BCR.matrix, name = "mean.of.z.score", 
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
             top_annotation = col_ha_top.BCR,
             #right_annotation = row_ha.right,
             show_row_dend = F,show_column_dend = F,
             #row_names_side = "left",
             left_annotation = row_ha.left.BCR,
             column_title="T/B cell receptor score",
             column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/immune.self.geneset/stat.20200407/ht.show/ht.TBcR.pdf",width = 5,height = 6)
draw(BCRht, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###########
###########immune evasion
##########
rownames(ht.cb.matrix)
eva.matrix=ht.cb.matrix[c(14:17,2),]
#eva.matrix=(eva.matrix - rowMeans(eva.matrix))/apply(eva.matrix,1,sd)
eva.pvalue=pvalue.cb[rownames(eva.matrix),]
rownames(eva.matrix)
#colore set
min_cor = min(as.vector(eva.matrix))
max_cor = max(as.vector(eva.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.eva = rowAnnotation(kruskal.pvalue=as.numeric(eva.pvalue$p.value),
                                col=list(
                                  kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                            c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                ),show_annotation_name = FALSE)
col_ha_top.eva = columnAnnotation(
  kaps.group=colnames(eva.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

evaht=Heatmap(eva.matrix, name = "mean.of.z.score", 
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
              top_annotation = col_ha_top.eva,
              #right_annotation = row_ha.right,
              show_row_dend = F,show_column_dend = F,
              #row_names_side = "left",
              left_annotation = row_ha.left.eva,
              column_title="Immune evasion genesets",
              column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/immune.self.geneset/stat.20200407/ht.show/ht.immune.evasion2.pdf",width = 5,height = 6)
draw(evaht, padding = unit(c(50, 5, 50,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
####################
####################genetic changes
######
rownames(ht.cb.matrix)
genetic.matrix=ht.cb.matrix[c(27:33),]
genetic.matrix=(genetic.matrix - rowMeans(genetic.matrix))/apply(genetic.matrix,1,sd)
genetic.pvalue=pvalue.cb[rownames(genetic.matrix),]
rownames(genetic.matrix)
#colore set
min_cor = min(as.vector(genetic.matrix))
max_cor = max(as.vector(genetic.matrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.genetic = rowAnnotation(kruskal.pvalue=as.numeric(genetic.pvalue$p.value),
                                col=list(
                                  kruskal.pvalue=colorRamp2(c(0.05,0.01,0.001,0.0001,10e-5,10e-10,10e-15), 
                                                            c("#4d4d4d","#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                ),show_annotation_name = FALSE)

col_ha_top.genetic = columnAnnotation(
  kaps.group=colnames(genetic.matrix),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

geneticht=Heatmap(genetic.matrix, name = "mean.of.z.score", 
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
              top_annotation = col_ha_top.genetic,
              #right_annotation = row_ha.right,
              show_row_dend = F,show_column_dend = F,
              #row_names_side = "left",
              left_annotation = row_ha.left.genetic,
              column_title="Genetic changes",
              column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/immune.self.geneset/stat.20200407/ht.show/ht.genetic.changes.pdf",width = 5,height = 6)
draw(geneticht, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
#######################
######################
###############
##########
dim(CB.data.pan.kaps)
ht.genetic.change.loged=data.frame(
                  kaps.group=CB.data.pan.kaps$kaps.group,
                  Aneuploidy.Score=log2(CB.data.pan.kaps$Aneuploidy.Score+1),
                  Homologous.Recombination.Defects=log2(CB.data.pan.kaps$Homologous.Recombination.Defects+1),
                  Number.of.Segments=log2(CB.data.pan.kaps$Number.of.Segments),
                  Fraction.Altered=CB.data.pan.kaps$Fraction.Altered,
                  SNV.Neoantigens=log2(CB.data.pan.kaps$SNV.Neoantigens),
                  Indel.Neoantigens=log2(CB.data.pan.kaps$Indel.Neoantigens+1),
                  Intratumor.Heterogeneity=CB.data.pan.kaps$Intratumor.Heterogeneity+1
)
rownames(ht.genetic.change.loged)=rownames(CB.data.pan.kaps)
#ht.genetic.change.loged=as.data.frame(t(ht.genetic.change.loged))
head(ht.genetic.change.loged)
#ht.genetic.change.loged.scaled=(ht.genetic.change.loged - rowMeans(ht.genetic.change.loged,na.rm = T))/apply(ht.genetic.change.loged,1,sd,na.rm=TRUE)
#ht.genetic.change.loged.scaled= scale(ht.genetic.change.loged, center = FALSE, scale = apply(ht.genetic.change.loged, 2, sd, na.rm = TRUE))
#stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
#ht.genetic.change.loged.scaled = lapply(ht.genetic.change.loged , stdize, na.rm = T)
#ht.genetic.change.loged.scaled=scale(ht.genetic.change.loged,center = F)
#ht.genetic.change.loged.scaled=as.data.frame(t(ht.genetic.change.loged.scaled))
datalist=list()
for (i in 1:7) {
  aa=ht.genetic.change.loged[,c(1,i+1)]
  aa=na.omit(aa)
  aa[,3]=(aa[,2] - mean(aa[,2],na.rm=T))/sd(aa[,2],na.rm = T)
  colnames(aa)[3]=paste(colnames(aa)[2], "new",sep="_")
  bb=aa  %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
  bb=as.data.frame(bb)
  rownames(bb)=bb$kaps.group
  bb=bb[,-1]
  datalist[[i]]=bb
}
ht.genetic.change.loged.scaled.median=do.call(cbind, datalist)
ht.genetic.change.loged.scaled.median=as.data.frame(t(ht.genetic.change.loged.scaled.median))
ht.genetic.change.loged.scaled.median=ht.genetic.change.loged.scaled.median[grepl(pattern = "new",rownames(ht.genetic.change.loged.scaled.median)),]
rownames(ht.genetic.change.loged.scaled.median)=substr(rownames(ht.genetic.change.loged.scaled.median),1,nchar(rownames(ht.genetic.change.loged.scaled.median))-4)
head(ht.genetic.change.loged.scaled.median)
#######
ht.genetic.change.loged.pvalue=pvalue.cb[rownames(pvalue.cb) %in% colnames(ht.genetic.change.loged),]
ht.genetic.change.loged.pvalue
######
#ht.genetic.change.loged.scaled=scale(ht.genetic.change.loged)
#head(ht.genetic.change.loged.scaled)
#colore set
min_cor = min(as.vector(ht.genetic.change.loged.scaled.median))
max_cor = max(as.vector(ht.genetic.change.loged.scaled.median))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
#
row_ha.left.genetic = rowAnnotation(kruskal.pvalue=as.numeric(ht.genetic.change.loged.pvalue$p.value),
                                    col=list(
                                      kruskal.pvalue=colorRamp2(c(0.05,0.01,0.001,0.0001,10e-5,10e-10,10e-15), 
                                                                c("#4d4d4d","#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                    ),show_annotation_name = FALSE)

col_ha_top.genetic = columnAnnotation(
  kaps.group=colnames(ht.genetic.change.loged.scaled.median),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
htgeneic=Heatmap(ht.genetic.change.loged.scaled.median, name = "mean.of.z.score", 
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
                  top_annotation = col_ha_top.genetic,
                  #right_annotation = row_ha.right,
                  show_row_dend = F,show_column_dend = F,
                  #row_names_side = "left",
                  left_annotation = row_ha.left.genetic,
                  column_title="Genetic changes",
                  column_title_gp = gpar(fontsize = 8)
)
pdf("4.model/kaps/immune.self.geneset/stat.20200407/ht.show/ht.genetic.changes..ew.pdf",width = 5,height = 6)
draw(htgeneic, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
