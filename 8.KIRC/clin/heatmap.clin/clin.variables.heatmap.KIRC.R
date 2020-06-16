######
######
########heatmap of TE expression and clinical variables in KIRC
########
head(stat.clin.kirc)
tmp.kirc.clin=stat.clin.kirc
rownames(tmp.kirc.clin)=tmp.kirc.clin$idss
head(tmp.kirc.clin)
####TMB
head(totoal.M)
colnames(totoal.M)[5]="idss"
####combine
stat.clin.kirc.TMB=merge(tmp.kirc.clin, totoal.M, by="idss",all.x=TRUE)
rownames(stat.clin.kirc.TMB)=stat.clin.kirc.TMB$idss
head(stat.clin.kirc.TMB)
###deal with the TE expression columns AluSq4 is missed, take from stat.kirc.kaps.vali.GEP
###
length(intersect(rownames(stat.clin.kirc.TMB), rownames(stat.kirc.kaps.vali.GEP)))
stat.clin.kirc.TMB.full=cbind(stat.kirc.kaps.vali.GEP[rownames(stat.clin.kirc.TMB), c(1:10)],stat.clin.kirc.TMB[,-c(1:11)])
stat.clin.kirc.TMB.full$TMB=stat.clin.kirc.TMB.full$total_mut/50
###mRNA.microRNA cluster change
head(clin.KIRC.sel.full)
head(stat.clin.kirc.TMB.full)
length(intersect(rownames(clin.KIRC.sel.full), rownames(stat.clin.kirc.TMB.full)))
stat.clin.kirc.TMB.full$mRNA_cluster=clin.KIRC.sel.full[rownames(stat.clin.kirc.TMB.full),]$mRNA_cluster
stat.clin.kirc.TMB.full$microRNA_cluster=clin.KIRC.sel.full[rownames(stat.clin.kirc.TMB.full),]$microRNA_cluster
#
stat.clin.kirc.TMB.full$mRNA_cluster=paste("cluster_",stat.clin.kirc.TMB.full$mRNA_cluster)
stat.clin.kirc.TMB.full$microRNA_cluster=paste("cluster_",stat.clin.kirc.TMB.full$microRNA_cluster)
stat.clin.kirc.TMB.full$mRNA_cluster=gsub("cluster_ NA","not_available",stat.clin.kirc.TMB.full$mRNA_cluster)
stat.clin.kirc.TMB.full$microRNA_cluster=gsub("cluster_ NA","not_available",stat.clin.kirc.TMB.full$microRNA_cluster)
table(stat.clin.kirc.TMB.full$microRNA_cluster)
##immune subtype
panimmunesubtype=read.table("8.KIRC/clin/Subtype_Immune_Model_Based.txt.gz",header = T,row.names = 1,sep = "\t")
panimmunesubtype$id=rownames(panimmunesubtype)
head(panimmunesubtype)
panmethysubtype=read.delim2("8.KIRC/clin/TCGASubtype.20170308.tsv.gz",header = T)
colnames(panmethysubtype)[1]="id"
head(panmethysubtype)
#combine
tcgapan.subtype=merge(panimmunesubtype, panmethysubtype, by="id",all.x=TRUE)
colnames(tcgapan.subtype)=paste0(colnames(tcgapan.subtype), ".pancancer")
tcgapan.subtype[is.na(tcgapan.subtype)] <- "not_available"
colnames(tcgapan.subtype)[1]="id"
head(tcgapan.subtype)
#save(tcgapan.subtype,file="8.KIRC/clin/tcgapancancer.subtype.RData")
###################
#######
stat.clin.kirc.TMB.full=merge(stat.clin.kirc.TMB.full, tcgapan.subtype, by="id",all.x=TRUE)
rownames(stat.clin.kirc.TMB.full)=stat.clin.kirc.TMB.full$id
head(stat.clin.kirc.TMB.full)
###
#deal wtih NA value
#
clin.heat.kirc.var <- sapply(stat.clin.kirc.TMB.full[,-c(3:11)], as.character)
clin.heat.kirc.var[is.na(clin.heat.kirc.var)] <- "not_available"
clin.heat.kirc.var=as.data.frame(clin.heat.kirc.var)
rownames(clin.heat.kirc.var)=clin.heat.kirc.var$id
clin.heat.kirc.var$z.of.mean.exp=as.numeric(paste(clin.heat.kirc.var$z.of.mean.exp))
clin.heat.kirc.var$TMB=as.numeric(paste(clin.heat.kirc.var$TMB))
clin.heat.kirc.var$TMB=log2(clin.heat.kirc.var$TMB)
clin.heat.kirc.var=clin.heat.kirc.var[order(clin.heat.kirc.var$z.of.mean.exp,decreasing = T),]
save(clin.heat.kirc.TE,clin.heat.kirc.var, stat.clin.kirc.TMB.full,file="8.KIRC/clin/heatmap.clin/clin.variables.heatmap.KIRC.RData")
###

####
head(clin.heat.kirc.var)
#chisq p value calculation
#########calculate the pvalue
age.pvalue.kirc=chisq.test(clin.heat.kirc.var$age.bi, clin.heat.kirc.var$kaps.group.kirc,correct = T)$p.value
gender.pvalue.kirc=chisq.test(clin.heat.kirc.var$gender.x, clin.heat.kirc.var$kaps.group.kirc,correct = T)$p.value

location.pvalue.kirc=chisq.test(table(clin.heat.kirc.var$laterality, clin.heat.kirc.var$kaps.group.kirc),correct = T)$p.value

AJCC.pvalue.kirc=chisq.test(table(clin.heat.kirc.var$AJCC.stage, clin.heat.kirc.var$kaps.group.kirc)[-1,],correct = T)$p.value

hypermutation.pvalue.kirc=chisq.test(table(clin.heat.kirc.var$Hyper.methylated.subtype.cluster.7, clin.heat.kirc.var$kaps.group.kirc)[-2,],correct = T)$p.value

mRNA.subtype.pvalue.kirc=chisq.test(table(clin.heat.kirc.var$mRNA_cluster, clin.heat.kirc.var$kaps.group.kirc)[-5,],correct = T)$p.value
microRNA.subtype.pvalue.kirc=chisq.test(table(clin.heat.kirc.var$microRNA_cluster, clin.heat.kirc.var$kaps.group.kirc)[-5,],correct = T)$p.value

Immune.subtype.pvalue.kirc=chisq.test(table(clin.heat.kirc.var$Subtype_Immune_Model_Based.pancancer, clin.heat.kirc.var$kaps.group.kirc)[-5,],correct = T)$p.value

TMB.pvalue.kirc=as.numeric(paste(compare_means(TMB ~ kaps.group.kirc, data = clin.heat.kirc.var,method = "kruskal.test")[2]))


#color set
######set color
kaps_col.kirc = c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
age_col.kirc=c("1low"="#c2e699","2high"="#238443")
gender_col.kirc=c("MALE"="#f03b20","FEMALE"="#feb24c")
tumor.location_col.kirc=c("Right"="#8856a7","Left"="#9ebcda","not_available"="white")

AJCC.stage_col.kirc=c("StageI"="#fcae91","StageII"="#fb6a4a","StageIII"="#de2d26","StageIV"="#a50f15","not_available"="white")

hypermutation_col.kirc=c("yes"="#1c9099","no"="#a6bddb","not_available"="white")


mRNA.subtype_col.kirc=c("cluster_ 1"="#386cb0","cluster_ 2"="#ff7f00","cluster_ 3"="#a6d854","cluster_ 4"="#f0027f","not_available"="white")

microRNA.subtype_col.kirc=c("cluster_ 1"="#e41a1c","cluster_ 2"="#377eb8","cluster_ 3"="#4daf4a","cluster_ 4"="#984ea3","not_available"="white")

Immune.subtype_col.kirc=c("IFN-gamma Dominant (Immune C2)"="#e6ab02","Immunologically Quiet (Immune C5)"="#66a61e","Inflammatory (Immune C3)"="#e7298a",
                     "Lymphocyte Depleted (Immune C4)"="#7570b3","TGF-beta Dominant (Immune C6)"="#d95f02","Wound Healing (Immune C1)"="#1b9e77","not_available"="white")
summary(clin.heat.kirc.var$TMB)
summary(clin.heat.kirc.var$z.of.mean.exp)

z.of.mean.exp.col.kirc=colorRamp2(c(-4,-2,0,2,4), 
                         c("#ffffcc","#d9f0a3","#addd8e","#78c679","#006837"))
TMB.col.kirc=colorRamp2(c(-3,-2,0,2,4), 
                                  c("#ffffcc","#d9f0a3","#addd8e","#78c679","#006837"))

#TMB.col.kirc=colorRamp2(c(-3,-2,0,2,4),    c("#ffffd4","#fed98e","#fe9929","#d95f0e","#993404"))







########
#
ha.kirc = HeatmapAnnotation(
  kaps.group=clin.heat.kirc.var$kaps.group.kirc,
  z.TE.score=clin.heat.kirc.var$z.of.mean.exp,
  age=clin.heat.kirc.var$age.bi,
  gender=clin.heat.kirc.var$gender.x,
  tumor.location=clin.heat.kirc.var$laterality,
  AJCC.stage=clin.heat.kirc.var$AJCC.stage,
  TMB=clin.heat.kirc.var$TMB,
  hypermutation=clin.heat.kirc.var$Hyper.methylated.subtype.cluster.7,
  mRNA.subtype=clin.heat.kirc.var$mRNA_cluster,
  microRNA.subtype=clin.heat.kirc.var$microRNA_cluster,
  Immune.subtype=clin.heat.kirc.var$Subtype_Immune_Model_Based.pancancer,
  #
  col = list(
    kaps.group=kaps_col.kirc,
    z.TE.score=z.of.mean.exp.col.kirc,
    age=age_col.kirc,
    gender=gender_col.kirc,
    tumor.location=tumor.location_col.kirc,
    TMB=TMB.col.kirc,
    AJCC.stage=AJCC.stage_col.kirc,
    hypermutation=hypermutation_col.kirc,
    mRNA.subtype=  mRNA.subtype_col.kirc,
    microRNA.subtype=microRNA.subtype_col.kirc,
    Immune.subtype=Immune.subtype_col.kirc
  ),
  na_col = "white", border = TRUE,
  show_legend = c(TRUE, TRUE, TRUE, TRUE,TRUE,TRUE, TRUE, TRUE, TRUE,TRUE,TRUE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    kaps.group=list(title="TE.cluster"),
    z.TE.score=list(title="z.TE.score"),
    age=list(title="age"),
    gender=list(title="gender"),
    tumor.location=list(title="tumor.location"),
    AJCC.stage=list(title="AJCC.stage"),
    TMB=list(title="TMB"),
    hypermutation=list(title="hypermutation"),
    mRNA.subtype=  list(title="mRNA.subtype"),
    microRNA.subtype=list(title="microRNA.subtype"),
    Immune.subtype=list(title="Immune.subtype")
  )
)
######TE annotation
TE.ann.data=cndi.rep[rownames(cndi.rep),c(2,3)]

TE.an.ha=rowAnnotation(repClass=TE.ann.data$repClass,
                       repFamily=TE.ann.data$repFamily,
                       col=list(repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")),show_annotation_name = FALSE
)
#
#####TE expression
clin.heat.kirc.TE <- stat.clin.kirc.TMB.full[rownames(clin.heat.kirc.var),3:11]
head(clin.heat.kirc.TE)
#aadraw.kirc=as.data.frame(t(as.data.frame(scale(clin.heat.kirc.TE))))
aadraw.kirc=as.data.frame(t(scale(clin.heat.kirc.TE)))

#aadraw.kirc=(aadraw.kirc - rowMeans(aadraw.kirc))/apply(aadraw.kirc,1,sd)
#aadraw.kirc[aadraw.kirc< -3] <- -3
#aadraw.kirc[aadraw.kirc> 3] <- 3

min_cor = min(as.vector(aadraw.kirc))
max_cor = max(as.vector(aadraw.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#00F5FF", "white","#FF3E96"))(50))
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
ht.kirc.clin = Heatmap((aadraw.kirc),#col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256)
             col=col.pal_cor,
             bottom_annotation = ha.kirc,
             show_column_names = F,left_annotation = TE.an.ha,
             cluster_columns = F,cluster_rows = T, column_title = "clinical.comparison among TE clusters.KIRC")


pdf("8.KIRC/clin/heatmap.clin/clin.variables.heatmap.KIRC.pdf",height = 10,width = 15)
draw(ht.kirc.clin, padding = unit(c(25, 80, 25,60), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     #,heatmap_legend_side = "right"
)


annotation_titles2.kirc = c(
  kaps.group=paste("TE.cluster"),
  z.TE.score=paste("z.TE.score"),
  age=paste0("age"),
  gender=paste0("gender"),
  tumor.location=paste0("tumor.location"),
  AJCC.stage=paste0("AJCC.stage"),
  TMB=paste0("log2(TMB)"),
  hypermutation=paste0("hypermutation"),
  mRNA.subtype=  paste0("mRNA.subtype"),
  microRNA.subtype=paste0("microRNA.subtype"),
  Immune.subtype=paste0("Immune.subtype")
)
annotation_titles3.kirc = c(
  kaps.group=paste("TE.cluster"),
  z.TE.score=paste("z.TE.score"),
  age=paste0("p=", round(age.pvalue.kirc,digits = 4)),
  gender=paste0("p=", round(gender.pvalue.kirc,digits = 4)),
  tumor.location=paste0("p=", round(location.pvalue.kirc,digits = 4)),
  AJCC.stage=paste0("p=", round(AJCC.pvalue.kirc,digits = 4)),
  TMB=paste0("p=", round(TMB.pvalue.kirc,digits = 4)),
  hypermutation=paste0("p=", round(hypermutation.pvalue.kirc,digits = 4)),
  mRNA.subtype=paste0("p<0.0001"),
  microRNA.subtype=paste0("p=", round(microRNA.subtype.pvalue.kirc,digits = 4)),
  Immune.subtype=paste0("p<0.0001")
)
for(an in names(annotation_titles3.kirc)) {
  decorate_annotation(an, {
    grid.text(annotation_titles3.kirc[an], unit(165, "mm"), just = "left")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
for(an in names(annotation_titles2.kirc)) {
  decorate_annotation(an, {
    grid.text(annotation_titles2.kirc[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

dev.off()



