#########compare four kaps group
#########multiCox
length(intersect(rownames(kaps.td), rownames(clin.data.CRC.selected)))
stat.clin=cbind(kaps.td, clin.data.CRC.selected[rownames(kaps.td),])
stat.clin$age.bi=ifelse(stat.clin$age >= 65,"2high","1low")
stat.clin$AJCC.stage.bi=gsub("stageIV","III&IV",stat.clin$AJCC.stage)
stat.clin$AJCC.stage.bi=gsub("stageIII","III&IV",stat.clin$AJCC.stage.bi)
stat.clin$AJCC.stage.bi=gsub("stageII","I&II",stat.clin$AJCC.stage.bi)
stat.clin$AJCC.stage.bi=gsub("stageI","I&II",stat.clin$AJCC.stage.bi)
#stat.clin=subset(stat.clin, AJCC.stage=="stageII"|AJCC.stage=="stageIII")
stat.clin$kaps.group.agg=gsub("set1","set.other",stat.clin$kaps.group)
stat.clin$kaps.group.agg=gsub("set2","set.other",stat.clin$kaps.group.agg)
stat.clin$location.bi=gsub("not_available",NA,stat.clin$location)
head(stat.clin)
table(stat.clin$location.bi)
#
#stat.clin$kaps.group = relevel(stat.clin$kaps.group, ref = "set3")
stat.clin$kaps.group <- factor(stat.clin$kaps.group, levels = c("set3","set1","set2","set4"))
#my.surv=Surv(kaps.td$OS.time,kaps.td$OS)
my.surv=Surv(stat.clin$DFI.time,stat.clin$DFI)
multicox=coxph(my.surv ~ kaps.group+gender+ age.bi+MSI.status.bin+lymphatic_invasion+location.bi+AJCC.stage.bi, data =  stat.clin)
#####
cc=summary(multicox)
cc
HR <- round(cc$coefficients[,2],2)
P_Value <- round(cc$coefficients[,5],4)
LCI <- round(cc$conf.int[,3],2)
UCI <- round(cc$conf.int[,4],2)
CI95 <- paste(LCI,'-',UCI)
dfm <- cbind(Hazard_Ratio = HR,CI95 = CI95,P_Value = P_Value)
dfm
#rownames(dfm)=c("risk.group","age","gender","MSI.status","AJCC.stage")
colnames(dfm)=paste(colnames(dfm),"PFI",sep = ".")
MulCox=as.data.frame(dfm)
MulCox.PFI=MulCox
write.csv(MulCox.PFI,"4.model/kaps/fur.cluster/multiCox/MulCox.PFI.res.csv")
#########################
########################
#####forestplot
#######
library(forestplot)
##########################OS
cochrane_from_rmeta <- 
  structure(list(
    HR = c(NA, NA,1.00,2.51,	4.19,	3.98,	1.10	,2.45	,1.39	,1.69,	1.74,	2.63), 
    lower = c(NA, NA,NA,0.77,	1.25,	1.09,	0.73,	1.50,	0.75,	1.06,	1.13,	1.63),
    upper = c(NA, NA,NA,8.19,	14.03,	14.57,	1.67,	3.99,	2.55,	2.70,	2.70,	4.23)),
    .Names = c("HR", "lower", "upper"), 
    row.names = c(NA, -12L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "risk.factors","set3 (ref)", "set1","set2","set4","gender (male vs female)","age (>=65 vs <65)","MSI.status (MSS vs MSI)","lymphatic_invasion (Yes vs No)","Location (proximal vs distal)","AJCC.stage (III&IV vs I&II)"),
  c("", "HR","1.00","2.51",	"4.19",	"3.98",	"1.10"	,"2.45"	,"1.39"	,"1.69",	"1.74",	"2.63"),
  c("", "95%CI",NA, "(0.77 - 8.19)",	"(1.25 - 14.03)",	"(1.09 - 14.57)",	"(0.73 - 1.67)",	"(1.50 - 3.99)",	"(0.75 - 2.55)",	"(1.06 - 2.70)",	"(1.13 - 2.70)",	"(1.63 - 4.23)"),
  c("", "P value",NA,"0.126",	"0.02",	"0.037",	"0.6381",	"<0.0001",	"0.2926",	"0.0284",	"0.0125",	"<0.0001"))

pdf("4.model/kaps/fur.cluster/multiCox/forestplot.OS.pdf",width = 7,height = 4)
forestplot(tabletext, title = paste0("MultiCox for OS"),
           cochrane_from_rmeta,
           clip=c(0.1,Inf), lwd.zero=2,ci.vertices.height=0.05,
           xlog=TRUE, boxsize = .15,lwd.ci=2,lty.ci=1,ci.vertices=TRUE,
           col=fpColors(box="#0571b0",line="#ca0020", summary="royalblue"))
dev.off()
######################DSS
cochrane_from_rmeta <- 
  structure(list(
    HR = c(NA, NA,1.00,5.47,	7.96,	9.52,	1.06,	1.83,	1.37,	2.17,	1.69,	4.09), 
    lower = c(NA, NA,NA,0.74,	1.04,	1.18,	0.64,	1.05,	0.60,	1.20,	0.99,	2.11),
    upper = c(NA, NA,NA,40.47,	60.75,	76.54,	1.76,	3.19,	3.10,	3.92,	2.87,	7.93)),
    .Names = c("HR", "lower", "upper"), 
    row.names = c(NA, -12L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "risk.factors","set3 (ref)", "set1","set2","set4","gender (male vs female)","age (>=65 vs <65)","MSI.status (MSS vs MSI)",
    "lymphatic_invasion (Yes vs No)","Location (proximal vs distal)","AJCC.stage (III&IV vs I&II)"),
  c("", "HR","1.00","5.47",	"7.96",	"9.52",	"1.06",	"1.83",	"1.37",	"2.17",	"1.69",	"4.09"),
  c("", "95%CI",NA, "(0.74 - 40.47)",	"(1.04 - 60.75)",	"(1.18 - 76.54)",	"(0.64 - 1.76)","(1.05 - 3.19)",	"(0.60 - 3.10)",	"(1.20 - 3.92)",	"(0.99 - 2.87)",	"(2.11 - 7.93)"),
  c("", "P value",NA,"0.0961"	,"0.0454",	"0.0342",	"0.8154","0.032","0.4559",	"0.0103",	"0.0526",	"<0.0001"))

pdf("4.model/kaps/fur.cluster/multiCox/forestplot.DSS.pdf",width = 7,height = 4)
forestplot(tabletext, title = paste0("MultiCox for DSS"),
           cochrane_from_rmeta,
           clip=c(0.1,Inf), lwd.zero=2,ci.vertices.height=0.05,
           xlog=TRUE, boxsize = .15,lwd.ci=2,lty.ci=1,ci.vertices=TRUE,
           col=fpColors(box="#0571b0",line="#ca0020", summary="royalblue"))
dev.off()
####################PFI
cochrane_from_rmeta <- 
  structure(list(
    HR = c(NA, NA,1.00,2.03,	1.91,	2.79,	1.29,	1.15,	1.36,	1.61,	1.27,	2.44), 
    lower = c(NA, NA,NA,0.87,	0.77,	1.07,	0.90,	0.79,	0.76,	1.08,	0.88,	1.62),
    upper = c(NA, NA,NA,4.71,	4.74,	7.31,	1.84,	1.67,	2.44,	2.40,	1.85,	3.68)),
    .Names = c("HR", "lower", "upper"), 
    row.names = c(NA, -12L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "risk.factors","set3 (ref)", "set1","set2","set4","gender (male vs female)","age (>=65 vs <65)",
    "MSI.status (MSS vs MSI)","lymphatic_invasion (Yes vs No)","Location (proximal vs distal)","AJCC.stage (III&IV vs I&II)"),
  c("", "HR","1.00","2.03",	"1.91",	"2.79",	"1.29",	"1.15",	"1.36",	"1.61",	"1.27",	"2.44"),
  c("", "95%CI",NA, "(0.87 - 4.71)",	"(0.77 - 4.74)",	"(1.07 - 7.31)",	"(0.90 - 1.84)",	"(0.79 - 1.67)",	"(0.76 - 2.44)",	"(1.08 - 2.40)",	"(0.88 - 1.85)",	"(1.62 - 3.68)"),
  c("", "P value",NA,"0.0994",	"0.1639",	"0.0365",	"0.1669",	"0.4592",	"0.2945",	"0.0193",	"0.206",	"<0.0001"))

pdf("4.model/kaps/fur.cluster/multiCox/forestplot.PFI.pdf",width = 7,height = 4)
forestplot(tabletext, title = paste0("MultiCox for PFI"),
           cochrane_from_rmeta,
           clip=c(0.1,Inf), lwd.zero=2,ci.vertices.height=0.05,
           xlog=TRUE, boxsize = .15,lwd.ci=2,lty.ci=1,ci.vertices=TRUE,
           col=fpColors(box="#0571b0",line="#ca0020", summary="royalblue"))
dev.off()
#################DFI
cochrane_from_rmeta <- 
  structure(list(
    HR = c(NA, NA,1.00,1.31,	0.70,	0.95,	2.67,	0.77,	2.50,	1.30,	0.85,	1.12), 
    lower = c(NA, NA,NA,0.17,	0.06,	0.08,	1.04,	0.31,	0.54,	0.51,	0.34,	0.44),
    upper = c(NA, NA,NA,10.18,	7.79,	11.28,	6.87,	1.88,	11.51,	3.31,	2.12,	2.89)),
    .Names = c("HR", "lower", "upper"), 
    row.names = c(NA, -12L), 
    class = "data.frame")

tabletext<-cbind(
  c("", "risk.factors","set3 (ref)", "set1","set2","set4","gender (male vs female)","age (>=65 vs <65)",
    "MSI.status (MSS vs MSI)","lymphatic_invasion (Yes vs No)","Location (proximal vs distal)","AJCC.stage (III&IV vs I&II)"),
  c("", "HR","1.00","1.31",	"0.70",	"0.95",	"2.67",	"0.77",	"2.50",	"1.30",	"0.85",	"1.12"),
  c("", "95%CI",NA, "(0.17 - 10.18)",	"(0.06 - 7.79)",	"(0.08 - 11.28)",	"(1.04 - 6.87)",	"(0.31 - 1.88)",	"(0.54 - 11.51)",	"(0.51 - 3.31)",	"(0.34 - 2.12)",	"(0.44 - 2.89)"),
  c("", "P value",NA,"0.7958",	"0.769",	"0.9699",	"0.0418",	"0.5608",	"0.2404",	"0.5898",	"0.7278",	"0.8073"))

pdf("4.model/kaps/fur.cluster/multiCox/forestplot.DFI.pdf",width = 7,height = 4)
forestplot(tabletext, title = paste0("MultiCox for DFI"),
           cochrane_from_rmeta,
           clip=c(0.1,Inf), lwd.zero=2,ci.vertices.height=0.05,
           xlog=TRUE, boxsize = .15,lwd.ci=2,lty.ci=1,ci.vertices=TRUE,
           col=fpColors(box="#0571b0",line="#ca0020", summary="royalblue"))
dev.off()
###################
#########subgroup analysis on stage
#####
subgroup.ajcc=subset(stat.clin, AJCC.stage=="stageII"| AJCC.stage=="stageIII")

fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = subgroup.ajcc )
kp1=ggsurvplot(fit2, data =subgroup.ajcc,
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
fit2<- survfit(Surv(DSS.time, DSS) ~ kaps.group, data = subgroup.ajcc )
kp2=ggsurvplot(fit2, data =subgroup.ajcc,
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
fit2<- survfit(Surv(DFI.time, DFI) ~ kaps.group, data = subgroup.ajcc )
kp3=ggsurvplot(fit2, data =subgroup.ajcc,
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
fit2<- survfit(Surv(PFI.time, PFI) ~ kaps.group, data =subgroup.ajcc )
kp4=ggsurvplot(fit2, data =subgroup.ajcc,
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
  pdf("4.model/kaps/fur.cluster/sub.stage.4.group/surC.4.cluster.stage2&3.pdf",width = 6,height = 8)
  print(kp1)
  print(kp2)
  print(kp3)
  print(kp4)
  dev.off()
}
generate.PDF(fig)
#######################
###########
save(stat.clin,file="4.model/kaps/fur.cluster/Cox.subgroup.analysis.stat.clin.RData")




