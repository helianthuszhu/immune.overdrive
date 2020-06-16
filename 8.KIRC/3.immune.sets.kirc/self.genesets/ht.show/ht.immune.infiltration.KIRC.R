####self score
dim(CB.data.self.kaps.kirc)
colnames(CB.data.self.kaps.kirc)
#load self genesets score from wolf
load("8.KIRC/3.immune.sets.kirc/self.genesets/pan.108s/CB.data.self.kaps.kirc.108s.RData")   ###containing GEP score
dim(CB.data.self.kaps.kirc.108s)
colnames(CB.data.self.kaps.kirc.108s)
length(intersect(rownames(CB.data.self.kaps.kirc),rownames(CB.data.self.kaps.kirc.108s)))
##combine two score
CB.data.score.all.self.wolf.kaps.kirc=cbind(CB.data.self.kaps.kirc,CB.data.self.kaps.kirc.108s[rownames(CB.data.self.kaps.kirc),-c(1:59)] )
CB.data.score.all.self.wolf.kaps.kirc=CB.data.score.all.self.wolf.kaps.kirc[,!(duplicated(colnames(CB.data.score.all.self.wolf.kaps.kirc)))]

#load the EMT score
load("8.KIRC/3.immune.sets.kirc/strom.fib.emt.ecm.sigs.kircfh.RData")
dim(strom.fib.emt.ecm.sigs.score.kircfh)
head(strom.fib.emt.ecm.sigs.score.kircfh)
CB.data.score.all.self.wolf.kaps.kirc=cbind(CB.data.score.all.self.wolf.kaps.kirc, strom.fib.emt.ecm.sigs.score.kircfh[rownames(CB.data.score.all.self.wolf.kaps.kirc),])
#load ipres score
load("8.KIRC/3.immune.sets.kirc/ipres/IPRES.score.KIRC.RData")
dim(ipres.score.kirc)
head(ipres.score.kirc)
CB.data.score.all.self.wolf.kaps.kirc=cbind(CB.data.score.all.self.wolf.kaps.kirc, ipres.score.kirc[rownames(CB.data.score.all.self.wolf.kaps.kirc),c(1:2)])
#Th1 signature
head(stat.th1sig.kirc)
dim(stat.th1sig.kirc)
#
CB.data.score.all.self.wolf.kaps.kirc=cbind(CB.data.score.all.self.wolf.kaps.kirc, stat.th1sig.kirc[rownames(CB.data.score.all.self.wolf.kaps.kirc),c(1,61)])
#change the GEP score back
###
CB.data.score.all.self.wolf.kaps.kirc$GEP=stat.kirc.kaps.vali.GEP[rownames(CB.data.score.all.self.wolf.kaps.kirc),]$GEP
CB.data.score.all.self.wolf.kaps.kirc[1:4,1:4]
#kaps.group.kirc
############
###########
#########self genesets
#########
#########
colnames(CB.data.score.all.self.wolf.kaps.kirc)
ht.infiltration.kirc=data.frame(
                                GEP=CB.data.score.all.self.wolf.kaps.kirc$GEP,
                           leukocyte.infiltration=CB.data.score.all.self.wolf.kaps.kirc$leukocyte.infiltration,
                           IFN.gamma.signature.18.genes.Ayers.etal=CB.data.score.all.self.wolf.kaps.kirc$IFN.gamma.signature.18.genes.Ayers.etal,
                           hot.tumor.signautre=CB.data.score.all.self.wolf.kaps.kirc$hot.tumor.signautre,
                           Th1.signature=CB.data.score.all.self.wolf.kaps.kirc$mean.Th1.sig,
                           MHC.II=CB.data.score.all.self.wolf.kaps.kirc$MHC.II,
                           MHC.I=CB.data.score.all.self.wolf.kaps.kirc$MHC.I,
                           Tcell_infiltration_1=CB.data.score.all.self.wolf.kaps.kirc$Tcell_infiltration_1,
                           CD4.T.cells.exhuasted=CB.data.score.all.self.wolf.kaps.kirc$CD4.T.cells.exhuasted.zhang.zeminscRNA,
                           CD8.T.cells.exhuasted=CB.data.score.all.self.wolf.kaps.kirc$CD8.T.cells.exhuasted.zhang.zeminscRNA,
                           TcClassII_score=CB.data.score.all.self.wolf.kaps.kirc$TcClassII_score,
                           TAMsurr_score=CB.data.score.all.self.wolf.kaps.kirc$TAMsurr_score,
                           TAMsurr_TcClassII_ratio=CB.data.score.all.self.wolf.kaps.kirc$TAMsurr_TcClassII_ratio,
                           #
                           TGFB.response.wolf=CB.data.score.all.self.wolf.kaps.kirc$TGFB.response.wolf,
                           up.ECM.signature=CB.data.score.all.self.wolf.kaps.kirc$up.ECM.signature,
                           WESTON_VEGFA_TARGETS_12HR=CB.data.score.all.self.wolf.kaps.kirc$WESTON_VEGFA_TARGETS_12HR,
                           EMT_UP=CB.data.score.all.self.wolf.kaps.kirc$EMT_UP,
                           IPRES=CB.data.score.all.self.wolf.kaps.kirc$z.mean.IPRES
)
rownames(ht.infiltration.kirc)=rownames(CB.data.score.all.self.wolf.kaps.kirc)
#scale
ht.infiltration.kirc=as.data.frame(scale(ht.infiltration.kirc))
ht.infiltration.kirc=cbind(CB.data.score.all.self.wolf.kaps.kirc$kaps.group.kirc,ht.infiltration.kirc )
colnames(ht.infiltration.kirc)[1]="kaps.group.kirc"

head(ht.infiltration.kirc)
######pvalue
datalist <- list()
for(i in names(ht.infiltration.kirc[,2:ncol(ht.infiltration.kirc)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group.kirc")), data = ht.infiltration.kirc)
}
pvalue.ht.inf.kirc=do.call(rbind, datalist)
pvalue.ht.inf.kirc=as.data.frame(pvalue.ht.inf.kirc)
######matrix
ht.stat.inf.kirc= ht.infiltration.kirc %>% group_by(kaps.group.kirc) %>% summarise_all(median,na.rm = TRUE)
ht.stat.inf.kirc=as.data.frame(ht.stat.inf.kirc)
rownames(ht.stat.inf.kirc)=ht.stat.inf.kirc$kaps.group.kirc
ht.stat.inf.kirc=ht.stat.inf.kirc[,-1]
ht.stat.inf.kirc=as.data.frame(t(ht.stat.inf.kirc))
#############
######in pan
########
##########
dim(CB.data.pan.kaps.kirc)
ht.BCR.kirc=data.frame(
                  TCR.Shannon=CB.data.pan.kaps.kirc$TCR.Shannon,
                  TCR.Richness=CB.data.pan.kaps.kirc$TCR.Richness,
                  TCR.Evenness=CB.data.pan.kaps.kirc$TCR.Evenness,
                  BCR.Shannon=CB.data.pan.kaps.kirc$BCR.Shannon,
                  BCR.Richness=CB.data.pan.kaps.kirc$BCR.Richness,
                  BCR.Evenness=CB.data.pan.kaps.kirc$BCR.Evenness,
                  Wound.Healing=CB.data.pan.kaps.kirc$Wound.Healing,
                  Lymphocyte.Infiltration.Signature.Score=CB.data.pan.kaps.kirc$Lymphocyte.Infiltration.Signature.Score,
                  Aneuploidy.Score=CB.data.pan.kaps.kirc$Aneuploidy.Score,
                  Homologous.Recombination.Defects=CB.data.pan.kaps.kirc$Homologous.Recombination.Defects,
                  Number.of.Segments=CB.data.pan.kaps.kirc$Number.of.Segments,
                  Fraction.Altered=CB.data.pan.kaps.kirc$Fraction.Altered,
                  SNV.Neoantigens=CB.data.pan.kaps.kirc$SNV.Neoantigens,
                  Indel.Neoantigens=CB.data.pan.kaps.kirc$Indel.Neoantigens,
                  Intratumor.Heterogeneity=CB.data.pan.kaps.kirc$Intratumor.Heterogeneity
)
rownames(ht.BCR.kirc)=rownames(CB.data.pan.kaps.kirc)
#######
ht.BCR.kirc=as.data.frame(scale(ht.BCR.kirc))
ht.BCR.kirc=cbind(CB.data.pan.kaps.kirc$kaps.group.kirc,ht.BCR.kirc )

colnames(ht.BCR.kirc)[1]="kaps.group.kirc"

head(ht.BCR.kirc)
######pvalue
datalist <- list()
for(i in names(ht.BCR.kirc[,2:ncol(ht.BCR.kirc)])){  
  datalist[[i]] <- kruskal.test(formula(paste(i, "~ kaps.group.kirc")), data = ht.BCR.kirc)
}
pvalue.ht.BCR.kirc=do.call(rbind, datalist)
pvalue.ht.BCR.kirc=as.data.frame(pvalue.ht.BCR.kirc)
######matrix
ht.stat.BCR.kirc= ht.BCR.kirc %>% group_by(kaps.group.kirc) %>% summarise_all(median,na.rm = TRUE)
ht.stat.BCR.kirc=as.data.frame(ht.stat.BCR.kirc)
rownames(ht.stat.BCR.kirc)=ht.stat.BCR.kirc$kaps.group.kirc
ht.stat.BCR.kirc=ht.stat.BCR.kirc[,-1]
ht.stat.BCR.kirc=as.data.frame(t(ht.stat.BCR.kirc))
head(ht.stat.BCR.kirc)
###########
###########combine matrix and p valie
####
ht.cb.matrix.kirc=rbind(ht.stat.inf.kirc, ht.stat.BCR.kirc)
pvalue.cb.kirc=rbind(pvalue.ht.inf.kirc,pvalue.ht.BCR.kirc)
save(ht.cb.matrix.kirc,pvalue.cb.kirc,file="8.KIRC/3.immune.sets.kirc/self.genesets/ht.show/ht.immune.infiltration.KIRC.RData")
#############draw now
###innmune infilitation
###
rownames(ht.cb.matrix.kirc)
infi.matrix.kirc=ht.cb.matrix.kirc[c(1,26,2:5,8:13,6,7),]
#infi.matrix=(infi.matrix - rowMeans(infi.matrix))/apply(infi.matrix,1,sd)
infi.pvalue.kirc=pvalue.cb.kirc[rownames(infi.matrix.kirc),]
rownames(infi.matrix.kirc)
#colore set
min_cor = min(as.vector(infi.matrix.kirc))
max_cor = max(as.vector(infi.matrix.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.infi.kirc = rowAnnotation(kruskal.pvalue=as.numeric(infi.pvalue.kirc$p.value),
                                 col=list(
                                   kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                             c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                 ),show_annotation_name = FALSE)
col_ha_top.infi.kirc = columnAnnotation(
  kaps.group=colnames(infi.matrix.kirc),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

inht.kirc=Heatmap(infi.matrix.kirc, name = "median.of.z.score", 
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
             top_annotation = col_ha_top.infi.kirc,
             #right_annotation = row_ha.right,
             show_row_dend = F,show_column_dend = F,
             #row_names_side = "left",
             left_annotation = row_ha.left.infi.kirc,
             column_title="Immune.infiltration.signature",
             column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/3.immune.sets.kirc/self.genesets/ht.show/ht.immune.infiltration.KIRC.GEP.changed.pdf",width = 5,height = 6)
draw(inht.kirc, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###########
###########BCR
###########
rownames(ht.cb.matrix.kirc)
BCR.matrix.kirc=ht.cb.matrix.kirc[c(19,22,20,23,21,24),]
BCR.pvalue.kirc=pvalue.cb.kirc[rownames(BCR.matrix.kirc),]
rownames(BCR.matrix.kirc)
#colore set
min_cor = min(as.vector(BCR.matrix.kirc))
max_cor = max(as.vector(BCR.matrix.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.BCR.kirc = rowAnnotation(kruskal.pvalue=as.numeric(BCR.pvalue.kirc$p.value),
                                col=list(
                                  kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                            c("#4d4d4d","#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354"))
                                ),show_annotation_name = FALSE)
col_ha_top.BCR.kirc = columnAnnotation(
  kaps.group=colnames(BCR.matrix.kirc),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

BCRht.kirc=Heatmap(BCR.matrix.kirc, name = "median.of.z.score", 
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
              top_annotation = col_ha_top.BCR.kirc,
              #right_annotation = row_ha.right,
              show_row_dend = F,show_column_dend = F,
              #row_names_side = "left",
              left_annotation = row_ha.left.BCR.kirc,
              column_title="T/B cell receptor score",
              column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/3.immune.sets.kirc/self.genesets/ht.show/ht.TBcR.KIRC.pdf",width = 5,height = 6)
draw(BCRht.kirc, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
###########
###########immune evasion
##########
rownames(ht.cb.matrix.kirc)
eva.matrix.kirc=ht.cb.matrix.kirc[c(14:18),]
#eva.matrix.kirc=(eva.matrix.kirc - rowMeans(eva.matrix.kirc))/apply(eva.matrix.kirc,1,sd)
eva.pvalue.kirc=pvalue.cb.kirc[rownames(eva.matrix.kirc),]
rownames(eva.matrix.kirc)
#colore set
min_cor = min(as.vector(eva.matrix.kirc))
max_cor = max(as.vector(eva.matrix.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.eva.kirc = rowAnnotation(kruskal.pvalue=as.numeric(eva.pvalue.kirc$p.value),
                                col=list(
                                  kruskal.pvalue=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), 
                                                            c("#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                ),show_annotation_name = FALSE)
col_ha_top.eva.kirc = columnAnnotation(
  kaps.group=colnames(eva.matrix.kirc),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

evaht.kirc=Heatmap(eva.matrix.kirc, name = "median.of.z.score", 
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
              top_annotation = col_ha_top.eva.kirc,
              #right_annotation = row_ha.right,
              show_row_dend = F,show_column_dend = F,
              #row_names_side = "left",
              left_annotation = row_ha.left.eva.kirc,
              column_title="Immune evasion genesets",
              column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/3.immune.sets.kirc/self.genesets/ht.show/ht.immune.evasion.KIRC.scaled.pdf",width = 5,height = 6)
draw(evaht.kirc, padding = unit(c(50, 5, 50,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
####################
####################genetic changes
######
rownames(ht.cb.matrix.kirc)
genetic.matrix.kirc=ht.cb.matrix.kirc[c(27:33),]
genetic.matrix.kirc=(genetic.matrix.kirc - rowMeans(genetic.matrix.kirc))/apply(genetic.matrix.kirc,1,sd)
genetic.pvalue.kirc=pvalue.cb.kirc[rownames(genetic.matrix.kirc),]
rownames(genetic.matrix.kirc)
#colore set
min_cor = min(as.vector(genetic.matrix.kirc))
max_cor = max(as.vector(genetic.matrix.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.genetic.kirc = rowAnnotation(kruskal.pvalue=as.numeric(genetic.pvalue.kirc$p.value),
                                    col=list(
                                      kruskal.pvalue=colorRamp2(c(0.05,0.01,0.001,0.0001,10e-5,10e-10,10e-15), 
                                                                c("#4d4d4d","#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                    ),show_annotation_name = FALSE)

col_ha_top.genetic.kirc = columnAnnotation(
  kaps.group=colnames(genetic.matrix.kirc),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####

geneticht.kirc=Heatmap(genetic.matrix.kirc, name = "median.of.z.score", 
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
                  top_annotation = col_ha_top.genetic.kirc,
                  #right_annotation = row_ha.right,
                  show_row_dend = F,show_column_dend = F,
                  #row_names_side = "left",
                  left_annotation = row_ha.left.genetic.kirc,
                  column_title="Genetic changes",
                  column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/3.immune.sets.kirc/self.genesets/ht.show/ht.genetic.changes.KIRC.pdf",width = 5,height = 6)
draw(geneticht.kirc, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
####################
###################
##################genetic changes
##################
##########
dim(CB.data.pan.kaps.kirc)
ht.genetic.change.loged.kirc=data.frame(
  kaps.group=CB.data.pan.kaps.kirc$kaps.group,
  Aneuploidy.Score=log2(CB.data.pan.kaps.kirc$Aneuploidy.Score+1),
  Homologous.Recombination.Defects=log2(CB.data.pan.kaps.kirc$Homologous.Recombination.Defects+1),
  Number.of.Segments=log2(CB.data.pan.kaps.kirc$Number.of.Segments),
  Fraction.Altered=CB.data.pan.kaps.kirc$Fraction.Altered,
  SNV.Neoantigens=log2(CB.data.pan.kaps.kirc$SNV.Neoantigens),
  Indel.Neoantigens=log2(CB.data.pan.kaps.kirc$Indel.Neoantigens+1),
  Intratumor.Heterogeneity=CB.data.pan.kaps.kirc$Intratumor.Heterogeneity+1
)
rownames(ht.genetic.change.loged.kirc)=rownames(CB.data.pan.kaps.kirc)
#ht.genetic.change.loged=as.data.frame(t(ht.genetic.change.loged))
head(ht.genetic.change.loged.kirc)
#ht.genetic.change.loged.scaled=(ht.genetic.change.loged - rowMeans(ht.genetic.change.loged,na.rm = T))/apply(ht.genetic.change.loged,1,sd,na.rm=TRUE)
#ht.genetic.change.loged.scaled= scale(ht.genetic.change.loged, center = FALSE, scale = apply(ht.genetic.change.loged, 2, sd, na.rm = TRUE))
#stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
#ht.genetic.change.loged.scaled = lapply(ht.genetic.change.loged , stdize, na.rm = T)
#ht.genetic.change.loged.scaled=scale(ht.genetic.change.loged,center = F)
#ht.genetic.change.loged.scaled=as.data.frame(t(ht.genetic.change.loged.scaled))
datalist=list()
for (i in 1:7) {
  aa=ht.genetic.change.loged.kirc[,c(1,i+1)]
  aa=na.omit(aa)
  aa[,3]=(aa[,2] - mean(aa[,2],na.rm=T))/sd(aa[,2],na.rm = T)
  colnames(aa)[3]=paste(colnames(aa)[2], "new",sep="_")
  bb=aa  %>% group_by(kaps.group) %>% summarise_all(median,na.rm = TRUE)
  bb=as.data.frame(bb)
  rownames(bb)=bb$kaps.group
  bb=bb[,-1]
  datalist[[i]]=bb
}
ht.genetic.change.loged.scaled.median.kirc=do.call(cbind, datalist)
ht.genetic.change.loged.scaled.median.kirc=as.data.frame(t(ht.genetic.change.loged.scaled.median.kirc))
ht.genetic.change.loged.scaled.median.kirc=ht.genetic.change.loged.scaled.median.kirc[grepl(pattern = "new",rownames(ht.genetic.change.loged.scaled.median.kirc)),]
rownames(ht.genetic.change.loged.scaled.median.kirc)=substr(rownames(ht.genetic.change.loged.scaled.median.kirc),1,nchar(rownames(ht.genetic.change.loged.scaled.median.kirc))-4)
head(ht.genetic.change.loged.scaled.median.kirc)
ht.genetic.change.loged.scaled.median.kirc=ht.genetic.change.loged.scaled.median.kirc[,4:1]
#######
ht.genetic.change.loged.pvalue.kirc=pvalue.cb.kirc[rownames(pvalue.cb.kirc) %in% colnames(ht.genetic.change.loged.kirc),]
ht.genetic.change.loged.pvalue.kirc
######
#ht.genetic.change.loged.scaled=scale(ht.genetic.change.loged)
#head(ht.genetic.change.loged.scaled)
#colore set
min_cor = min(as.vector(ht.genetic.change.loged.scaled.median.kirc))
max_cor = max(as.vector(ht.genetic.change.loged.scaled.median.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))
#
row_ha.left.genetic.kirc = rowAnnotation(kruskal.pvalue=as.numeric(ht.genetic.change.loged.pvalue.kirc$p.value),
                                    col=list(
                                      kruskal.pvalue=colorRamp2(c(0.05,0.01,0.001,0.0001,10e-5,10e-10,10e-15), 
                                                                c("#4d4d4d","#ffffcc","#d9f0a3","#addd8e","#78c679","#31a354","#006837"))
                                    ),show_annotation_name = FALSE)

col_ha_top.genetic.kirc = columnAnnotation(
  kaps.group=colnames(ht.genetic.change.loged.scaled.median.kirc),
  col=list(kaps.group=c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" )),
  show_annotation_name = FALSE,gp = gpar(col = "black"))

####
htgeneic.kirc=Heatmap(ht.genetic.change.loged.scaled.median.kirc, name = "mean.of.z.score", 
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
                 top_annotation = col_ha_top.genetic.kirc,
                 #right_annotation = row_ha.right,
                 show_row_dend = F,show_column_dend = F,
                 #row_names_side = "left",
                 left_annotation = row_ha.left.genetic.kirc,
                 column_title="Genetic changes",
                 column_title_gp = gpar(fontsize = 8)
)
pdf("8.KIRC/3.immune.sets.kirc/self.genesets/ht.show/ht.genetic.changes.KIRC.new.pdf",width = 5,height = 6)
draw(htgeneic.kirc, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()
