########### compare pan cancer pathways activities
##############
pathypancancer=read.table("4.model/immune/pancancer.pathways.xena/PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level_Z.txt.gz",header = T,row.names = 1,sep = "\t")
pathypancancer=as.data.frame(t(pathypancancer))
rownames(pathypancancer)=gsub("[.]","-",rownames(pathypancancer))
pathypancancer[1:4,1:4]
##############
pathyid=intersect(rownames(pathypancancer), rownames(kaps.td))
###
stat.pathyp=cbind(kaps.td[pathyid,]$kaps.group,pathypancancer[pathyid,])
stat.pathyp=cbind(z.mean.statsM[rownames(stat.pathyp),]$z.of.mean.exp,stat.pathyp)
colnames(stat.pathyp)[1:2]=c("z.of.mean.exp","kaps.group")
stat.pathyp[1:4,1:4]
######1 calculate the cor value
########
#####
datalist=list()
for (i in 3:ncol(stat.pathyp)) {
  aa=cbind(stat.pathyp[,i], stat.pathyp[,1])
  colnames(aa)=c(colnames(stat.pathyp)[i],colnames(stat.pathyp)[1])
  aa=na.omit(aa)
  res <- cor.test(aa[,2], aa[,1],method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(stat.pathyp)[i]
  rownames(res)=paste(colnames(stat.pathyp)[1], rownames(res),sep = "&")
  datalist[[i]] <- res
}
cor.res.pathyp=do.call(rbind, datalist)
head(cor.res.pathyp)
summary(cor.res.pathyp$cor)
######2 calculate Anova p value
datalist=list()
for (i in 3:ncol(stat.pathyp)) {
  stat.pathyp$kaps.group=as.factor(stat.pathyp$kaps.group)
  res.aov <- aov(stat.pathyp[,i] ~ kaps.group, data = stat.pathyp)
  # Summary of the analysis
  pv=data.frame(summary(res.aov)[[1]][["Pr(>F)"]][1])
  colnames(pv)='anova.p.value'
  rownames(pv)=colnames(stat.pathyp)[i]
  datalist[[i]] <- pv
}
anova.p.pathyp = do.call(rbind, datalist)
head(anova.p.pathyp)
anova.p.pathyp$pathyName=rownames(anova.p.pathyp)
######combine two stats
resCB.pathyp=cbind(cor.res.pathyp, anova.p.pathyp)
dim(subset(resCB.pathyp, abs(cor)>= 0.4))
head(resCB.pathyp)
#####
#resCB.pathyp.sig=subset(resCB.pathyp, pvalue < 0.05 & abs(cor)>= 0.3 &  anova.p.value < 0.05)
resCB.pathyp.sig=subset(resCB.pathyp, anova.p.value < 0.05)
rownames(resCB.pathyp.sig)=resCB.pathyp.sig$pathyName
head(resCB.pathyp.sig)
dim(resCB.pathyp.sig)
#####
ann.col.pathyp=kaps.td[rownames(stat.pathyp),]
ann.col.pathyp=ann.col.pathyp[order(ann.col.pathyp$z.of.mean.exp,decreasing = T),]
dim(ann.col.pathyp)
head(ann.col.pathyp)
#
stat.pathypMatrix=as.data.frame(t(stat.pathyp[rownames(ann.col.pathyp),-c(1:2)]))[rownames(resCB.pathyp.sig),]
dim(stat.pathypMatrix)
stat.pathypMatrix[1:4,1:4]
#
col_ha.top1.pathyp=columnAnnotation(z.of.mean.exp.score=anno_lines(ann.col.pathyp$z.of.mean.exp),
                             TE.cluster = ann.col.pathyp$TE.cluster,
                             TE.cluster.agg=ann.col.pathyp$TE.cluster.agg,
                             MSI.status.bin=ann.col.pathyp$MSI.status.bin,
                             z.of.mean.exp=ann.col.pathyp$z.of.mean.exp,
                             kaps.group=ann.col.pathyp$kaps.group,
                             col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                      TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                      MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)
#
row_ha_left.pathyp = rowAnnotation(Spearman.Rho=resCB.pathyp.sig$cor,
                                   Anova.p.value=resCB.pathyp.sig$anova.p.value,
                                   col=list(Spearman.Rho=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
                                            Anova.p.value=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), c("#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#91003f"))
                                   ),
  annotation_name_side = "bottom")
#######sig cor >=0.4 to show
########
####choose some specifi tes to show
pathshow=subset(resCB.pathyp.sig, abs(cor)>= 0.4)
dim(pathshow)
head(pathshow)
showingp=rownames(pathshow)
idxpath=match(showingp,rownames(stat.pathypMatrix))
#
row_ha.right = rowAnnotation(foo2 = anno_mark(at = idxpath, labels = showingp,
                                              link_width = unit(6, "mm")),show_annotation_name = FALSE)

pathpp=Heatmap(stat.pathypMatrix, 
               name = paste0("significant among kaps group (anova < 0.05)","\n",
                             "847 sig and 36 ( abs.cor >=0.4)"), 
               na_col = "gray",
               col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),border = F,
               show_column_names = F,show_row_names = F,
               cluster_columns = F,
               cluster_rows = T,
               row_names_gp = gpar(fontsize = 3),
               column_names_gp = gpar(fontsize = 8),
               top_annotation = col_ha.top1.pathyp, 
               right_annotation = row_ha.right,
               #show_row_dend = T,show_column_dend = T,
               #row_names_side = "left",
               left_annotation = row_ha_left.pathyp
)
sdssp=t(scale(t(stat.pathypMatrix)))
sdssp[sdssp< -3] <- -3
sdssp[sdssp>3] <- 3
dim(sdssp)
dim(resCB.pathyp)
pathpp2=Heatmap(sdssp, 
                name = paste0("significant among kaps group (anova < 0.05)","\n",
                              "847 sig and 36 ( abs.cor >=0.4)"), 
               na_col = "gray",
               col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = viridis(10),border = F,
               show_column_names = F,show_row_names = F,
               cluster_columns = F,
               cluster_rows = T,
               row_names_gp = gpar(fontsize = 3),
               column_names_gp = gpar(fontsize = 8),
               top_annotation = col_ha.top1.pathyp, 
               right_annotation = row_ha.right,
               #show_row_dend = T,show_column_dend = T,
               #row_names_side = "left",
               left_annotation = row_ha_left.pathyp
)

generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/heatmp.256.sig.pathys.annoted.with.pvalue.2.only.anova.0.05.showing.cor.0.4.pdf",height = 9,width = 20)
  print(pathpp)
  print(pathpp2)
  dev.off()
}
generate.PDF(fig)

save(resCB.pathyp,resCB.pathyp.sig,ann.col.pathyp, stat.pathypMatrix,file = "4.model/immune/pancancer.pathways.xena/heatmp.256.sig.pathys.annoted.with.pvalue.RData")
write.csv(resCB.pathyp.sig, "4.model/immune/pancancer.pathways.xena/heatmp.256.sig.pathys.csv")

pathpp=Heatmap(stat.pathypMatrix, 
               name = paste0("positively correlated with","\n",
                              "mean.TE.exp (Rho >= 0.3)", "\n",
                               " and significant among kaps group"), 
        na_col = "gray",
        col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
        #col = rev(viridis(10)),border = F,
        show_column_names = F,show_row_names = T,
        cluster_columns = F,
        cluster_rows = T,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha.top1.pathyp
        #right_annotation = row_ha_right#,
        #show_row_dend = T,show_column_dend = T,
        #row_names_side = "left",
        #left_annotation = row_ha_left.possig
)
generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/heatmp.256.sig.pathys.pdf",height = 5,width = 10)
  print(pathpp)
  dev.off()
}
generate.PDF(fig)
#
generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/heatmp.256.sig.pathys.long.pdf",height = 50,width = 10)
  print(pathpp)
  dev.off()
}
generate.PDF(fig)
######


