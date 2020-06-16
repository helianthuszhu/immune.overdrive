head(stat.clin)
head(totoal.M)
stat.clin$id=rownames(stat.clin)
totoal.M$id=rownames(totoal.M)
head(cb3)
dim(cb3)
cb3$id=rownames(cb3)
########
clin.heat=merge(stat.clin, totoal.M, by="id",all.x=TRUE)
clin.heat=merge(clin.heat, cb3, by="id",all.x=TRUE)
clin.heat$TMB=clin.heat$total_mut/50
clin.heat$TMB=log2(clin.heat$TMB)
clin.heat[["TMB"]][is.na(clin.heat[["TMB"]])] <- "not_available"
clin.heat$TMB=as.numeric(paste(clin.heat$TMB))
class(clin.heat$TMB)
dim(clin.heat)
head(clin.heat)
table(clin.heat$cimpNature)
#deal wtih NA value
clin.heat <- sapply(clin.heat, as.character)
clin.heat[is.na(clin.heat)] <- "not_available"
clin.heat=as.data.frame(clin.heat)
clin.heat$z.of.mean.exp=as.numeric(paste(clin.heat$z.of.mean.exp))
clin.heat$TMB=as.numeric(paste(clin.heat$TMB))
clin.heat=clin.heat[order(clin.heat$z.of.mean.exp,decreasing = T),]
rownames(clin.heat)=clin.heat$id
save(clin.heat, file="4.model/kaps/fur.cluster/clin.variables.heatmap.RData")
########
#########calculate the pvalue
age.pvalue=chisq.test(clin.heat$age.bi, clin.heat$kaps.group,correct = T)$p.value
gender.pvalue=chisq.test(clin.heat$gender, clin.heat$kaps.group,correct = T)$p.value
location.pvalue=chisq.test(table(clin.heat$location.bi, clin.heat$kaps.group)[-2,],correct = T)$p.value
AJCC.pvalue=chisq.test(table(clin.heat$AJCC.stage, clin.heat$kaps.group)[-1,],correct = T)$p.value
MSI.pvalue=chisq.test(table(clin.heat$MSI.status.bin, clin.heat$kaps.group),correct = T)$p.value
MLH1.pvalue=chisq.test(table(clin.heat$MLH1_silencing, clin.heat$kaps.group)[-3,],correct = T)$p.value
hypermutation.pvalue=chisq.test(table(clin.heat$hypermutation, clin.heat$kaps.group)[-2,],correct = T)$p.value
CIMP1.pvalue=chisq.test(table(clin.heat$cimp1, clin.heat$kaps.group)[-4,],correct = T)$p.value
CIMP2.pvalue=chisq.test(table(clin.heat$cimpNature, clin.heat$kaps.group)[-5,],correct = T)$p.value
TCGA.subtype.pvalue=chisq.test(table(clin.heat$TCGA_subtypes, clin.heat$kaps.group)[-4,],correct = T)$p.value
CRCassigner.subtype.pvalue=chisq.test(table(clin.heat$mole_subtype, clin.heat$kaps.group)[-3,],correct = T)$p.value
CMS.subtype.pvalue=chisq.test(table(clin.heat$CMS.cluster, clin.heat$kaps.group)[-5,],correct = T)$p.value
Immune.subtype.pvalue=chisq.test(table(clin.heat$Subtype_Immune_Model_Based, clin.heat$kaps.group)[-4,],correct = T)$p.value
TMB.pvalue=as.numeric(paste(compare_means(TMB ~ kaps.group, data = clin.heat,method = "kruskal.test")[2]))

######set color
kaps_col = c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
age_col=c("1low"="#c2e699","2high"="#238443")
gender_col=c("MALE"="#f03b20","FEMALE"="#feb24c")
tumor.location_col=c("proximal"="#8856a7","distal"="#9ebcda","not_available"="white")
AJCC.stage_col=c("stageI"="#fcae91","stageII"="#fb6a4a","stageIII"="#de2d26","stageIV"="#a50f15","not_available"="white")
MSI.status_col=c("MSI"="#e7298a","MSS"="#66a61e","not_available"="white")
MLH1._col=c("1"="#2b8cbe","0"="#a6bddb","not_available"="white")
hypermutation_col=c("TRUE"="#1c9099","FALSE"="#a6bddb","not_available"="white")
CIMP1_col=c("CIMP.High"="#7a0177","CIMP.Low"="#c51b8a","CIMP.Neg"="#f768a1","not_available"="white")
CIMP2_col=c("CIMP.H"="#7a0177","CIMP.L"="#c51b8a","Cluster3"="#f768a1","Cluster4"="#fbb4b9","not_available"="white")

TCGA.subtype_col=c("CIN"="#377EB8","Invasive"="#1a9850","MSI/CIMP"="#E41A1C","not_available"="white")

CRCassigner.subtype_col=c("Enterocyte/Goblet-like"="#386cb0","Inflammatory"="#ff7f00","Stem-like"="#a6d854","TA"="#f0027f","not_available"="white")

CMS.subtype_col=c("CMS1"="#e41a1c","CMS2"="#377eb8","CMS3"="#4daf4a","CMS4"="#984ea3","not_available"="white")

Immune.subtype_col=c("IFN-gamma Dominant (Immune C2)"="#e6ab02","Immunologically Quiet (Immune C5)"="#66a61e","Inflammatory (Immune C3)"="#e7298a",
                     "Lymphocyte Depleted (Immune C4)"="#7570b3","TGF-beta Dominant (Immune C6)"="#d95f02","Wound Healing (Immune C1)"="#1b9e77","not_available"="white")

summary(clin.heat$TMB)
TMB_col=colorRamp2(c(-4,-2,0,2,4,6,8), 
                        c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02"))


#
ha = HeatmapAnnotation(
  kaps.group=clin.heat$kaps.group,
  z.TE.score=clin.heat$z.of.mean.exp,
  age=clin.heat$age.bi,
  gender=clin.heat$gender,
  tumor.location=clin.heat$location,
  AJCC.stage=clin.heat$AJCC.stage,
  MSI.status=clin.heat$MSI.status.bin,
  MLH1.silencing=clin.heat$MLH1_silencing,
  TMB=clin.heat$TMB,
  hypermutation=clin.heat$hypermutation,
  CIMP1=clin.heat$cimp1,
  CIMP2=clin.heat$cimpNature,
  TCGA.subtype=clin.heat$TCGA_subtypes,
  CRCassigner.subtype=clin.heat$mole_subtype,
  CMS.subtype=clin.heat$CMS.cluster,
  Immune.subtype=clin.heat$Subtype_Immune_Model_Based,
#
  col = list(
    kaps.group=kaps_col,
    age=age_col,
    gender=gender_col,
    tumor.location=tumor.location_col,
    AJCC.stage=AJCC.stage_col,
    MSI.status=MSI.status_col,
    MLH1.silencing=MLH1._col,
    hypermutation=hypermutation_col,
    TMB=TMB_col,
    CIMP1=CIMP1_col,
    CIMP2=CIMP2_col,
    TCGA.subtype=TCGA.subtype_col,
    CRCassigner.subtype=CRCassigner.subtype_col,
    CMS.subtype=CMS.subtype_col,
    Immune.subtype=Immune.subtype_col
    ),
  na_col = "white", border = TRUE,
  show_legend = c(TRUE, TRUE, TRUE, TRUE,TRUE,TRUE, TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE,TRUE, TRUE,TRUE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    kaps.group=list(title="TE.cluster"),
    z.TE.score=list(title="z.TE.score"),
    age=list(title="age"),
    gender=list(title="gender"),
    tumor.location=list(title="tumor.location"),
    AJCC.stage=list(title="AJCC.stage"),
    MSI.status=list(title="MSI.status"),
    MLH1.silencing=list(title="MLH1.silencing"),
    TMB=list(title="TMB"),
    hypermutation=list(title="hypermutation"),
    CIMP1=list(title="CIMP1"),
    CIMP2=list(title="CIMP2"),
    TCGA.subtype=list(title="TCGA.subtype"),
    CRCassigner.subtype=list(title="mRNA.subtype"),
    CMS.subtype=list(title="CMS.subtype"),
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

#zero_row_mat = matrix(nrow = 0, ncol = 590)
#ht = Heatmap(zero_row_mat, top_annotation = ha,cluster_columns = F,cluster_rows = F, column_title = "clinical.comparison among TE clusters")

aadraw=as.data.frame(t(kaps.td[,3:11]))
aadraw=(aadraw - rowMeans(aadraw))/apply(aadraw,1,sd)
aadraw[aadraw< -2] <- -3
aadraw[aadraw> 2] <- 3

min_cor = min(as.vector(aadraw))
max_cor = max(as.vector(aadraw))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#00F5FF", "white","#FF3E96"))(50))
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
ht = Heatmap((aadraw),#col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256)
             col=col.pal_cor,
             bottom_annotation = ha,
             show_column_names = F,left_annotation = TE.an.ha,
             cluster_columns = F,cluster_rows = T, column_title = "clinical.comparison among TE clusters")


pdf("4.model/kaps/fur.cluster/clin.variables.heatmap.2.new.pdf",height = 10,width = 15)
draw(ht, padding = unit(c(22, 80, 20,60), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     #,heatmap_legend_side = "right"
)


annotation_titles = c(
                      kaps.group=paste("TE.cluster"),
                      z.TE.score=paste("z.TE.score"),
                      age=paste0("age"," ","p=", round(age.pvalue,digits = 4)),
                      gender=paste0("gender"," ","p=", round(gender.pvalue,digits = 4)),
                      tumor.location=paste0("tumor.location"," ","p=", round(location.pvalue,digits = 4)),
                      AJCC.stage=paste0("AJCC.stage"," ","p=", round(AJCC.pvalue,digits = 4)),
                      MSI.status=paste0("MSI.status"," ","p=", round(MSI.pvalue,digits = 4)),
                      MLH1.silencing=paste0("MLH1.silencing"," ","p=", round(MLH1.pvalue,digits = 4)),
                      TMB=paste0("log2(TMB)"," ","p=", round(TMB.pvalue,digits = 4)),
                      hypermutation=paste0("hypermutation"," ","p=", round(hypermutation.pvalue,digits = 4)),
                      CIMP1=paste0("CIMP1"," ","p=", round(CIMP1.pvalue,digits = 4)),
                      CIMP2=paste0("CIMP2"," ","p=", round(CIMP2.pvalue,digits = 4)),
                      TCGA.subtype=paste0("TCGA.subtype"," ","p=", round(TCGA.subtype.pvalue,digits = 4)),
                      CRCassigner.subtype=paste0("molecular.subtype"," ","p=", round(CRCassigner.subtype.pvalue,digits = 4)),
                      CMS.subtype=paste0("CMS.subtype"," ","p=", round(CMS.subtype.pvalue,digits = 4)),
                      Immune.subtype=paste0("Immune.subtype"," ","p=", round(Immune.subtype.pvalue,digits = 4))
                      )

annotation_titles2 = c(
  kaps.group=paste("TE.cluster"),
  z.TE.score=paste("z.TE.score"),
  age=paste0("age"),
  gender=paste0("gender"),
  tumor.location=paste0("tumor.location"),
  AJCC.stage=paste0("AJCC.stage"),
  MSI.status=paste0("MSI.status"),
  MLH1.silencing=paste0("MLH1.silencing"),
  TMB=paste0("log2(TMB)"),
  hypermutation=paste0("hypermutation"),
  CIMP1=paste0("CIMP1"),
  CIMP2=paste0("CIMP2"),
  TCGA.subtype=paste0("TCGA.subtype"),
  CRCassigner.subtype=paste0("molecular.subtype"),
  CMS.subtype=paste0("CMS.subtype"),
  Immune.subtype=paste0("Immune.subtype")
)
annotation_titles3 = c(
  kaps.group=paste("TE.cluster"),
  z.TE.score=paste("z.TE.score"),
  age=paste0("p=", round(age.pvalue,digits = 4)),
  gender=paste0("p=", round(gender.pvalue,digits = 4)),
  tumor.location=paste0("p=", round(location.pvalue,digits = 4)),
  AJCC.stage=paste0("p=", round(AJCC.pvalue,digits = 4)),
  MSI.status=paste0("p=", round(MSI.pvalue,digits = 4)),
  MLH1.silencing=paste0("p=", round(MLH1.pvalue,digits = 4)),
  TMB=paste0("p=", round(TMB.pvalue,digits = 4)),
  hypermutation=paste0("p=", round(hypermutation.pvalue,digits = 4)),
  CIMP1=paste0("p=", round(CIMP1.pvalue,digits = 4)),
  CIMP2=paste0("p=", round(CIMP2.pvalue,digits = 4)),
  TCGA.subtype=paste0("p=", round(TCGA.subtype.pvalue,digits = 4)),
  CRCassigner.subtype=paste0("p=", round(CRCassigner.subtype.pvalue,digits = 4)),
  CMS.subtype=paste0("p<0.0001"),
  Immune.subtype=paste0("p<0.0001")
)
for(an in names(annotation_titles3)) {
  decorate_annotation(an, {
    grid.text(annotation_titles3[an], unit(165, "mm"), just = "left")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
for(an in names(annotation_titles2)) {
  decorate_annotation(an, {
    grid.text(annotation_titles2[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

dev.off()
#####################CMS distribution
########################
x=as.data.frame.matrix(table(clin.heat$CMS.cluster, clin.heat$kaps.group))
x$group=rownames(x)
head(x)
write.csv(x, "4.model/kaps/fur.cluster/stacked.CMS.subtype.vs.kaps.group.csv")

pval.kaps.cms <- chisq.test(table(clin.heat$CMS.cluster, clin.heat$kaps.group)[-5,],correct = T)$p.value

stackdata.kaps.cms=as.data.frame.matrix(table(clin.heat$CMS.cluster, clin.heat$kaps.group)[-5,])
datm.kaps.cms <- melt(cbind(stackdata.kaps.cms, ind = rownames(stackdata.kaps.cms)), id.vars = c('ind'))
datm.kaps.cms$variable=factor(datm.kaps.cms$variable,levels = c("set4","set3","set2","set1"))
datm.kaps.cms$ind=factor(datm.kaps.cms$ind,levels = rev(unique(datm.kaps.cms$ind)))
pdf("4.model/kaps/fur.cluster/stacked.CMS.subtype.vs.kaps.group.pdf",width = 4,height = 6)
ggplot(datm.kaps.cms,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())+#coord_flip() +
  theme(legend.position="right")+
  scale_fill_manual(values= c("#984ea3","#4daf4a","#377eb8","#e41a1c"))+
  ylab("Percentage(%)")+xlab("kaps.group")+
  annotate("text", x=3, y=0.9, label=paste0("p-value: ","\n" ,signif(pval.kaps.cms,4)),size=5)+
  ggtitle("CMS subtype.vs.TE.kaps.group")
dev.off()

####################
head(clin.heat)
##############
##############
##############
#########subgroup analysis on stage
#####
clin.heat$OS.time=as.numeric(paste(clin.heat$OS.time))
clin.heat$DSS.time=as.numeric(paste(clin.heat$DSS.time))
clin.heat$DFI.time=as.numeric(paste(clin.heat$DFI.time))
clin.heat$PFI.time=as.numeric(paste(clin.heat$PFI.time))
#
clin.heat$OS=as.numeric(paste(clin.heat$OS))
clin.heat$DSS=as.numeric(paste(clin.heat$DSS))
clin.heat$DFI=as.numeric(paste(clin.heat$DFI))
clin.heat$PFI=as.numeric(paste(clin.heat$PFI))



subgroup.immune=subset(clin.heat,Subtype_Immune_Model_Based=="Wound Healing (Immune C1)")
subgroup.immune=subset(clin.heat,Subtype_Immune_Model_Based=="IFN-gamma Dominant (Immune C2)")

subgroup.immune=subset(clin.heat,CMS.cluster=="CMS4")
#
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = subgroup.immune )
kp1=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
#####
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group, data = subgroup.immune)
kp2=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
###
fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group, data = subgroup.immune)
kp3=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               ylab="disease-free interval probability",
               xlab="Time(days)"
)
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data =subgroup.immune )
kp4=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
##
generate.PDF <- function(fig) {
  pdf("4.model/kaps/fur.cluster/sub.immune/surC.4.cluster.CMS4.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
###############further compare set4 vs set3
#######
subgroup.immune=subset(clin.heat,Subtype_Immune_Model_Based=="Wound Healing (Immune C1)")
subgroup.immune=subset(clin.heat,Subtype_Immune_Model_Based=="IFN-gamma Dominant (Immune C2)")
subgroup.immune=subset(subgroup.immune, kaps.group=="set4"|kaps.group=="set3")
#
ss=Surv(subgroup.immune$OS.time, subgroup.immune$OS)
cox = summary(coxph(ss~factor(subgroup.immune$kaps.group,levels = c("set3","set4"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = subgroup.immune )
kp1=ggsurvplot(fit2, data =subgroup.immune,
               title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set3","set4"),  # Change legend labels
               palette = c("#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
#####
ss=Surv(subgroup.immune$DSS.time, subgroup.immune$DSS)
cox = summary(coxph(ss~factor(subgroup.immune$kaps.group,levels = c("set3","set4"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group, data = subgroup.immune)
kp2=ggsurvplot(fit2, data =subgroup.immune,
               title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set3","set4"),  # Change legend labels
               palette = c("#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
###
ss=Surv(subgroup.immune$DFI.time, subgroup.immune$DFI)
cox = summary(coxph(ss~factor(subgroup.immune$kaps.group,levels = c("set3","set4"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group, data = subgroup.immune)
kp3=ggsurvplot(fit2, data =subgroup.immune,
               title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               ylab="disease-free interval probability",
               xlab="Time(days)"
)
####
ss=Surv(subgroup.immune$PFI.time, subgroup.immune$PFI)
cox = summary(coxph(ss~factor(subgroup.immune$kaps.group,levels = c("set3","set4"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data =subgroup.immune )
kp4=ggsurvplot(fit2, data =subgroup.immune,
               title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set3","set4"),  # Change legend labels
               palette = c("#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
##
generate.PDF <- function(fig) {
  pdf("4.model/kaps/fur.cluster/sub.immune/surC.4.cluster.C2.set3.vs.set4.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
###################
###################
dim(expM)
eg=bitr(rownames(expM), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg=eg[!(duplicated(eg$SYMBOL)),]
rownames(eg)=eg$SYMBOL
head(eg)
expM.ENTREZID=expM[rownames(eg),]
rownames(expM.ENTREZID)=eg$ENTREZID
expM.ENTREZID=expM.ENTREZID[, colnames(expM.ENTREZID) %in% rownames(kaps.td)]

expM.ENTREZID[1:4,1:4]
#
library(CMSclassifier)
Rfcms <- CMSclassifier::classifyCMS(expM.ENTREZID,method="RF")[[3]]
SScms <- CMSclassifier::classifyCMS(expM.ENTREZID,method="SSP")[[3]]
Rfcms=cbind(Rfcms, clin.heat[rownames(Rfcms),])
SScms=cbind(SScms, clin.heat[rownames(SScms),])
save(Rfcms,SScms, file = "4.model/kaps/fur.cluster/CMS.estimation.output.RData")
#
table(Rfcms$RF.predictedCMS,Rfcms$CMS.cluster)
#table(SScms$SSP.predictedCMS,SScms$CMS.cluster)
########
########
######
subgroup.immune=subset(Rfcms,RF.predictedCMS=="CMS4")
#
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = subgroup.immune )
kp1=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
#####
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group, data = subgroup.immune)
kp2=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
###
fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group, data = subgroup.immune)
kp3=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               ylab="disease-free interval probability",
               xlab="Time(days)"
)
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data =subgroup.immune )
kp4=ggsurvplot(fit2, data =subgroup.immune,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
##
generate.PDF <- function(fig) {
  pdf("4.model/kaps/fur.cluster/sub.immune/CMS.estimated/surC.4.cluster.CMS4.estimation.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
#####################total mutation
####
pdf("4.model/kaps/fur.cluster/TMB.violin.pdf",width = 3.5,height = 5)

ggviolin(clin.heat,x = "kaps.group", y ="TMB" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
            palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
            add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+xlab("TE.cluster")+ylab("TMB") +ggtitle("TMB")

ggviolin(clin.heat,x = "kaps.group", y ="TMB" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
         palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif",method = "kruskal.test")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+xlab("TE.cluster")+ylab("TMB") +ggtitle("TMB")
dev.off()
################TE score correlation with TMB
head(clin.heat)
pdf("4.model/kaps/fur.cluster/correlation.TMB.vs.TE.score.pdf",width = 5,height = 6)
ggscatter(clin.heat, x = "z.of.mean.exp", y = "TMB", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "mean.TE.score", ylab = "TMB")+
  ggtitle("mean.TE.exp vs TMB in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
################
pdf("4.model/kaps/fur.cluster/TE.score.MSIvsMSS.pdf",width = 3,height = 5)
ggboxplot(clin.heat, x = "MSI.status.bin", y = "z.of.mean.exp",
          color = "MSI.status.bin", palette = c("MSI"="#d01c8b","MSS"="#4dac26"),
          add = "jitter")+stat_compare_means()
ggviolin(clin.heat,x = "MSI.status.bin", y ="z.of.mean.exp" , fill = "MSI.status.bin",alpha = 1,size = 0.01,width = 1,
         palette = c("MSI"="#d01c8b","MSS"="#4dac26"),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "MSI.status.bin")+
  stat_compare_means(label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+xlab("MSI.status")+ylab("TE.score") +ggtitle("TE.score.vs.MSI.status")
dev.off()
###############
###############
###############subgroup analysis based on TMB
###############
head(clin.heat)
TMB.subgroupdata=clin.heat
TMB.subgroupdata$total_mut=as.numeric(paste(TMB.subgroupdata$total_mut))
TMB.subgroupdata=subset(TMB.subgroupdata,  total_mut> 0)
TMB.subgroupdata$tmb.group=ifelse(TMB.subgroupdata$TMB> median(TMB.subgroupdata$TMB), "TMB.high","TMB.low")
TMB.subgroupdata$OS.time=as.numeric(paste(TMB.subgroupdata$OS.time))
TMB.subgroupdata$DSS.time=as.numeric(paste(TMB.subgroupdata$DSS.time))
TMB.subgroupdata$DFI.time=as.numeric(paste(TMB.subgroupdata$DFI.time))
TMB.subgroupdata$PFI.time=as.numeric(paste(TMB.subgroupdata$PFI.time))

TMB.subgroupdata$OS=as.numeric(paste(TMB.subgroupdata$OS))
TMB.subgroupdata$DSS=as.numeric(paste(TMB.subgroupdata$DSS))
TMB.subgroupdata$DFI=as.numeric(paste(TMB.subgroupdata$DFI))
TMB.subgroupdata$PFI=as.numeric(paste(TMB.subgroupdata$PFI))
#############
###########
##########
TMB.subgroupdata.sub=subset(TMB.subgroupdata, tmb.group=="TMB.low" )

fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = TMB.subgroupdata.sub)
kp1=ggsurvplot(fit2, data =TMB.subgroupdata.sub,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
#####
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group, data = TMB.subgroupdata.sub)
kp2=ggsurvplot(fit2, data =TMB.subgroupdata.sub,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
###
fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group, data = TMB.subgroupdata.sub)
kp3=ggsurvplot(fit2, data =TMB.subgroupdata.sub,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               #ylab="progression-free interval probability",
               ylab="disease-free interval probability",
               xlab="Time(days)"
)
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data =TMB.subgroupdata.sub )
kp4=ggsurvplot(fit2, data =TMB.subgroupdata.sub,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
               palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),  # Use JCO journal color palette,
               #palette =c("#d73027","#E69F00","#00AFBB"),
               #palette = c("#ca0020","#0571b0","#4daf4a"),
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,             # Hide tables y axis text
               #ylab="Disease specific survival probability",
               #ylab="Overall survival probability",
               ylab="progression-free interval probability",
               #ylab="disease-free interval probability",
               xlab="Time(days)"
)
##
generate.PDF <- function(fig) {
  pdf("4.model/kaps/fur.cluster/sub.TMB/subgroup.TMB.low.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
#################
write.csv(clin.heat, "4.model/kaps/fur.cluster/clin.data.full.CRC.590s.csv")
