#############compare the clinal variable among kaps group in KIRC 
############
head(clin.KIRC.sel.full)
head(stat.kirc.kaps.vali)
stat.kirc.kaps.vali$id1=stat.kirc.kaps.vali$bcr_patient_barcode
###########
dim(stat.kirc.kaps.vali)
length(unique(stat.kirc.kaps.vali$id1))
length(unique(clin.KIRC.sel.full$id1))
#
length(intersect(clin.KIRC.sel.full$id1, stat.kirc.kaps.vali$id1))
stat.clin.kirc=merge(stat.kirc.kaps.vali, clin.KIRC.sel.full, by="id1",all.x=TRUE)
dim(stat.clin.kirc)
head(stat.clin.kirc)
#######################
table(stat.clin.kirc$kaps.group.kirc, stat.clin.kirc$COCA.sybtype)
chisq.test(table(stat.clin.kirc$kaps.group.kirc, stat.clin.kirc$microRNA_cluster)[,1:4],correct = T)$p.value


pdf("8.KIRC/clin/total.DNA.methylation.index.pdf",height  = 5,width = 4)
ggviolin(stat.clin.kirc,x = "kaps.group.kirc", y ="total.DNA.methylation.index" , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
         palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("kaps.group")+ylab("total.DNA.methylation.index") +
  ggtitle(paste0("total.DNA.methylation.index",".KIRC.cutoff"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
###########multiCox analysis
###########
###########
stat.clin.kirc$age.bi=ifelse(stat.clin.kirc$age >= 65,"2high","1low")
stat.clin.kirc$AJCC.stage=gsub(" ","", stat.clin.kirc$pathologic_stage)
stat.clin.kirc$AJCC.stage=ifelse(stat.clin.kirc$AJCC.stage=="StageI"|stat.clin.kirc$AJCC.stage=="StageII"|
                                   stat.clin.kirc$AJCC.stage=="StageIII"|stat.clin.kirc$AJCC.stage=="StageIV",stat.clin.kirc$AJCC.stage,NA)
stat.clin.kirc$AJCC.stage.bi=gsub("StageIV","III&IV",stat.clin.kirc$AJCC.stage)
stat.clin.kirc$AJCC.stage.bi=gsub("StageIII","III&IV",stat.clin.kirc$AJCC.stage.bi)
stat.clin.kirc$AJCC.stage.bi=gsub("StageII","I&II",stat.clin.kirc$AJCC.stage.bi)
stat.clin.kirc$AJCC.stage.bi=gsub("StageI","I&II",stat.clin.kirc$AJCC.stage.bi)
stat.clin.kirc$grade=ifelse(stat.clin.kirc$neoplasm_histologic_grade=="G1"|stat.clin.kirc$neoplasm_histologic_grade=="G2"|
                              stat.clin.kirc$neoplasm_histologic_grade=="G3"|stat.clin.kirc$neoplasm_histologic_grade=="G4",stat.clin.kirc$neoplasm_histologic_grade, NA)
stat.clin.kirc$mRNA_cluster=paste("mRNA.cluster",stat.clin.kirc$mRNA_cluster,sep = ".")
stat.clin.kirc$microRNA_cluster=paste("microRNA.cluster",stat.clin.kirc$microRNA_clusterr,sep = ".")
#
stat.clin.kirc$grade.bi=gsub("G1","G1&2",stat.clin.kirc$grade)
stat.clin.kirc$grade.bi=gsub("G2","G1&2",stat.clin.kirc$grade.bi)
stat.clin.kirc$grade.bi=gsub("G3","G3&4",stat.clin.kirc$grade.bi)
stat.clin.kirc$grade.bi=gsub("G4","G3&4",stat.clin.kirc$grade.bi)
#
table(stat.clin.kirc$grade)
###########
stat.clin.kirc$kaps.group.kirc <- factor(stat.clin.kirc$kaps.group.kirc, levels = c("set3","set1","set2","set4"))
stat.clin.kirc$laterality <- factor(stat.clin.kirc$laterality, levels = c("Right","Left"))
#my.surv=Surv(kaps.td$OS.time,kaps.td$OS)
my.surv=Surv(stat.clin.kirc$DFI.time,stat.clin.kirc$DFI)
multicox=coxph(my.surv ~ kaps.group.kirc+gender.x+ age.bi+laterality+AJCC.stage.bi, data =  stat.clin.kirc)
#####
cc=summary(multicox)
cc
colnames(stat.clin.kirc)
HR <- round(cc$coefficients[,2],2)
P_Value <- round(cc$coefficients[,5],4)
LCI <- round(cc$conf.int[,3],2)
UCI <- round(cc$conf.int[,4],2)
CI95 <- paste(LCI,'-',UCI)
dfm <- cbind(Hazard_Ratio = HR,CI95 = CI95,P_Value = P_Value)
dfm
#rownames(dfm)=c("risk.group","age","gender","MSI.status","AJCC.stage")
colnames(dfm)=paste(colnames(dfm),"DFI",sep = ".")
MulCox=as.data.frame(dfm)
MulCox.DFI=MulCox
write.csv(MulCox.DFI,"8.KIRC/1.tescore.survival.kaps.validation/multiCox/MulCox.DFI.res.kirc.csv")
#####
save(stat.clin.kirc,file="8.KIRC/1.tescore.survival.kaps.validation/multiCox/stat.clin.kirc.RData")
########
########forestplot
######
library(forestplot)
##########################OS
cochrane_from_rmeta <- 
  structure(list(
    HR = c(NA, NA,1.00,1.08,	2.12,	1.96,	0.98,	1.59,	1.36,	3.56), 
    lower = c(NA, NA,NA,0.65,	1.02,	1.14,	0.71,	1.17,	1.01,	2.57),
    upper = c(NA, NA,NA,1.81,	4.40,	3.34,	1.35,	2.17,	1.85,	4.95)),
    .Names = c("HR", "lower", "upper"), 
    row.names = c(NA, -10L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "risk.factors","set3 (ref)", "set1","set2","set4","gender (male vs female)","age (>=65 vs <65)","location (left vs right)","AJCC.stage (III&IV vs I&II)"),
  c("", "HR","1.00","1.08",	"2.12",	"1.96",	"0.98",	"1.59",	"1.36",	"3.56"),
  c("", "95%CI",NA, "(0.65 - 1.81)",	"(1.02 - 4.40)",	"(1.14 - 3.34)",	"(0.71 - 1.35)",	"(1.17 - 2.17)",	"(1.01 - 1.85)",	"(2.57 - 4.95)"),
  c("", "P value",NA,"0.7713",	"0.0442",	"0.0141",	"0.9091",	"0.0033",	"0.0453",	"<0.0001"))

pdf("8.KIRC/1.tescore.survival.kaps.validation/multiCox/forestplot.OS.kirc.pdf",width = 7,height = 4)
forestplot(tabletext, title = paste0("MultiCox for OS.KIRC"),
           cochrane_from_rmeta,
           clip=c(0.1,Inf), lwd.zero=2,ci.vertices.height=0.05,
           xlog=TRUE, boxsize = .15,lwd.ci=2,lty.ci=1,ci.vertices=TRUE,
           col=fpColors(box="#0571b0",line="#ca0020", summary="royalblue"))
dev.off()
###########
###########DSS
##
cochrane_from_rmeta <- 
  structure(list(
    HR = c(NA, NA,1.00,1.25,	2.68,	2.36,	1.18,	1.08,	1.57,	8.71), 
    lower = c(NA, NA,NA,0.63,	0.99,	1.18,	0.77,	0.72,	1.06,	5.19),
    upper = c(NA, NA,NA,2.47,	7.25,	4.74,	1.82,	1.63,	2.32,	14.61)),
    .Names = c("HR", "lower", "upper"), 
    row.names = c(NA, -10L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "risk.factors","set3 (ref)", "set1","set2","set4","gender (male vs female)","age (>=65 vs <65)","location (left vs right)","AJCC.stage (III&IV vs I&II)"),
  c("", "HR","1.00","1.25",	"2.68",	"2.36",	"1.18",	"1.08",	"1.57",	"8.71"),
  c("", "95%CI",NA, "(0.63 - 2.47)",	"(0.99 - 7.25)",	"(1.18 - 4.74)","(0.77 - 1.82)",	"(0.72 - 1.63)",	"(1.06 - 2.32)",	"(5.19 - 14.61)"),
  c("", "P value",NA,"0.5307","0.0529",	"0.0158",	"0.441",	"0.6971",	"0.0236",	"<0.0001"))

pdf("8.KIRC/1.tescore.survival.kaps.validation/multiCox/forestplot.DSS.kirc.pdf",width = 7,height = 4)
forestplot(tabletext, title = paste0("MultiCox for DSS.KIRC"),
           cochrane_from_rmeta,
           clip=c(0.1,Inf), lwd.zero=2,ci.vertices.height=0.05,
           xlog=TRUE, boxsize = .15,lwd.ci=2,lty.ci=1,ci.vertices=TRUE,
           col=fpColors(box="#0571b0",line="#ca0020", summary="royalblue"))
dev.off()
############
###########PFI
cochrane_from_rmeta <- 
  structure(list(
    HR = c(NA, NA,1.00,1.63,	3.20,	2.17,	1.52,	1.10,	1.68,	7.15), 
    lower = c(NA, NA,NA,0.90,	1.37,	1.17,	1.06,	0.79,	1.21,	4.90),
    upper = c(NA, NA,NA,2.94,	7.45,	4.04,	2.19,	1.54,	2.32,	10.43)),
    .Names = c("HR", "lower", "upper"), 
    row.names = c(NA, -10L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "risk.factors","set3 (ref)", "set1","set2","set4","gender (male vs female)","age (>=65 vs <65)","location (left vs right)","AJCC.stage (III&IV vs I&II)"),
  c("", "HR","1.00","1.63",	"3.20",	"2.17",	"1.52",	"1.10",	"1.68",	"7.15"),
  c("", "95%CI",NA, "(0.9 - 2.94)",	"(1.37 - 7.45)",	"(1.17 - 4.04)","(1.06 - 2.19)",	"(0.79 - 1.54)",	"(1.21 - 2.32)",	"(4.9 - 10.43)"),
  c("", "P value",NA,"0.1037",	"0.007",	"0.0144",	"0.0233",	"0.5733",	"0.0018",	"<0.0001"))

pdf("8.KIRC/1.tescore.survival.kaps.validation/multiCox/forestplot.PFI.kirc.pdf",width = 7,height = 4)
forestplot(tabletext, title = paste0("MultiCox for PFI.KIRC"),
           cochrane_from_rmeta,
           clip=c(0.1,Inf), lwd.zero=2,ci.vertices.height=0.05,
           xlog=TRUE, boxsize = .15,lwd.ci=2,lty.ci=1,ci.vertices=TRUE,
           col=fpColors(box="#0571b0",line="#ca0020", summary="royalblue"))
dev.off()
###########

