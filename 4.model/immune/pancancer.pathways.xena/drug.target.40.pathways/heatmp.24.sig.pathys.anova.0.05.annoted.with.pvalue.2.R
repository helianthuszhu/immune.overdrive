########### compare pan cancer pathways activities
##############
pathypancancer.drug=read.table("4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/Pancan12_GenePrograms_drugTargetCanon_in_Pancan33.tsv.gz",header = T,row.names = 1,sep = "\t")
pathypancancer.drug=as.data.frame(t(pathypancancer.drug))
rownames(pathypancancer.drug)=gsub("[.]","-",rownames(pathypancancer.drug))
pathypancancer.drug[1:4,1:4]
##############
pathyid.drug=intersect(rownames(pathypancancer.drug), rownames(kaps.td))
###
stat.pathyp.drug=cbind(kaps.td[pathyid.drug,]$kaps.group,pathypancancer.drug[pathyid.drug,])
stat.pathyp.drug=cbind(z.mean.statsM[rownames(stat.pathyp.drug),]$z.of.mean.exp,stat.pathyp.drug)
colnames(stat.pathyp.drug)[1:2]=c("z.of.mean.exp","kaps.group")
stat.pathyp.drug[1:4,1:4]
######1 calculate the cor value
########
#####
datalist=list()
for (i in 3:ncol(stat.pathyp.drug)) {
  aa=cbind(stat.pathyp.drug[,i], stat.pathyp.drug[,1])
  colnames(aa)=c(colnames(stat.pathyp.drug)[i],colnames(stat.pathyp.drug)[1])
  aa=na.omit(aa)
  res <- cor.test(aa[,2], aa[,1],method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(stat.pathyp.drug)[i]
  rownames(res)=paste(colnames(stat.pathyp.drug)[1], rownames(res),sep = "&")
  datalist[[i]] <- res
}
cor.res.pathyp.drug=do.call(rbind, datalist)
head(cor.res.pathyp.drug)
summary(cor.res.pathyp$cor)
######2 calculate Anova p value
datalist=list()
for (i in 3:ncol(stat.pathyp.drug)) {
  stat.pathyp.drug$kaps.group=as.factor(stat.pathyp.drug$kaps.group)
  res.aov <- aov(stat.pathyp.drug[,i] ~ kaps.group, data = stat.pathyp.drug)
  # Summary of the analysis
  pv=data.frame(summary(res.aov)[[1]][["Pr(>F)"]][1])
  colnames(pv)='anova.p.value'
  rownames(pv)=colnames(stat.pathyp.drug)[i]
  datalist[[i]] <- pv
}
anova.p.pathyp.drug = do.call(rbind, datalist)
head(anova.p.pathyp.drug)
anova.p.pathyp.drug$pathyName=rownames(anova.p.pathyp.drug)
######combine two stats
resCB.pathyp.drug=cbind(cor.res.pathyp.drug, anova.p.pathyp.drug)
head(resCB.pathyp.drug)
dim(subset(resCB.pathyp.drug, anova.p.value < 0.05))
#####consider cor value
resCB.pathyp.sig.drug=subset(resCB.pathyp.drug, anova.p.value < 0.05)
rownames(resCB.pathyp.sig.drug)=resCB.pathyp.sig.drug$pathyName
head(resCB.pathyp.sig.drug)
dim(resCB.pathyp.sig.drug)
#####
ann.col.pathyp.drug=kaps.td[rownames(stat.pathyp.drug),]
ann.col.pathyp.drug=ann.col.pathyp.drug[order(ann.col.pathyp.drug$z.of.mean.exp,decreasing = T),]
dim(ann.col.pathyp.drug)
head(ann.col.pathyp.drug)
#
stat.pathypMatrix.drug=as.data.frame(t(stat.pathyp.drug[rownames(ann.col.pathyp.drug),-c(1:2)]))[rownames(resCB.pathyp.sig.drug),]


dim(stat.pathypMatrix.drug)
stat.pathypMatrix.drug[1:4,1:4]
#
col_ha.top1.pathyp.drug=columnAnnotation(z.of.mean.exp.score=anno_lines(ann.col.pathyp.drug$z.of.mean.exp),
                                    TE.cluster = ann.col.pathyp.drug$TE.cluster,
                                    TE.cluster.agg=ann.col.pathyp.drug$TE.cluster.agg,
                                    MSI.status.bin=ann.col.pathyp.drug$MSI.status.bin,
                                    z.of.mean.exp=ann.col.pathyp.drug$z.of.mean.exp,
                                    kaps.group=ann.col.pathyp.drug$kaps.group,
                                    col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                             TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                             MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                             kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                                    ),show_annotation_name = TRUE)
#
row_ha_left.pathyp.drug = rowAnnotation(Spearman.Rho=resCB.pathyp.sig.drug$cor,
                                   Anova.p.value=resCB.pathyp.sig.drug$anova.p.value,
                                   col=list(Spearman.Rho=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
                                            Anova.p.value=colorRamp2(c(0.05,10e-5,10e-5,10e-10,10e-20,10e-30), c("#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#91003f"))
                                   ),
                                   annotation_name_side = "bottom")
sdss=t(scale(t(stat.pathypMatrix.drug)))
sdss[sdss< -3] <- -3
sdss[sdss>3] <- 3
pathpp1=Heatmap(sdss, 
               name = paste0("significant among kaps group","\n","(anova < 0.05)","\n",
                             "24 out of 40 pathways"), 
               na_col = "gray",
               col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               #col = rev(viridis(10)),border = F,
               show_column_names = F,show_row_names = T,
               cluster_columns = F,
               cluster_rows = T,
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               top_annotation = col_ha.top1.pathyp.drug, 
               #right_annotation = row_ha_right#,
               #show_row_dend = T,show_column_dend = T,
               #row_names_side = "left",
               left_annotation = row_ha_left.pathyp.drug
)
#
#
generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/heatmp.24.sig.pathys.anova.0.05.annoted.with.pvalue.2.pdf",height = 5,width = 15)
  print(pathpp1)
  dev.off()
}
generate.PDF(fig)

save(resCB.pathyp.drug,resCB.pathyp.sig.drug,ann.col.pathyp.drug, stat.pathypMatrix.drug,
     file = "4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/heatmp.24.sig.pathys.anova.0.05.annoted.with.pvalue.2.RData")
write.csv(resCB.pathyp.sig.drug, "4.model/immune/pancancer.pathways.xena/drug.target.40.pathways/heatmp.24.sig.pathys.anova.0.05.annoted.with.pvalue.2.csv")
######


