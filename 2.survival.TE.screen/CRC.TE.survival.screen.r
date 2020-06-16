
###get the survival data
CRC.survivaldata=read.table("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/clinical.data/COADREAD_survival.from.pan-cancer.Atlas.txt",header = T,row.names = 1,sep = "\t")

###get the TE expression
load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.plot/tmp.RData")
load('~/nas/Xiaoqiang/opti.data//1.CGC.TE.CRC.REdiscoverTE/stats.20200311/1.TE.used/TE.label.used.RData')

###inclue TE used
head(cdrep)
table(cdrep$repClass.new2)
cdrep.sel=subset(cdrep, !(repClass.new2=="other.repeats"))
dim(cdrep.sel)
###TE matrix selection
te.tumorexp.sel=te.tumorexp[, colnames(te.tumorexp) %in% rownames(cdrep.sel)]
te.tumorexp.sel[1:4,1:4]

#screen for OS
#OS
surOS.tumor=subset(CRC.survivaldata,OS.time > 0 & OS>= 0)
idsur=intersect(rownames(te.tumorexp.sel),rownames(surOS.tumor))
length(idsur)
cbOSdata=cbind(surOS.tumor[idsur,2:3],te.tumorexp.sel[idsur,])

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$OS,  #CHANGE
                             os=cbOSdata$OS.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.os=cox_results[,cox_results[4,]<0.05]
dim(Results.os)

#DSS
surDSS.tumor=subset(CRC.survivaldata,DSS.time > 0 & DSS >=0)
idsur=intersect(rownames(te.tumorexp.sel),rownames(surDSS.tumor))
length(idsur)
cbOSdata=cbind(surDSS.tumor[idsur,4:5],te.tumorexp.sel[idsur,])     ####570 SAMPLES

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$DSS,  #CHANGE
                             os=cbOSdata$DSS.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.dss=cox_results[,cox_results[4,]<0.05]
dim(Results.dss)

#DFI
surDFI.tumor=subset(CRC.survivaldata,DFI.time > 0 & DFI >=0)
idsur=intersect(rownames(te.tumorexp.sel),rownames(surDFI.tumor))
length(idsur)
cbOSdata=cbind(surDFI.tumor[idsur,6:7],te.tumorexp.sel[idsur,])     #####234 SAMPLES

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$DFI,  #CHANGE
                             os=cbOSdata$DFI.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.dfi=cox_results[,cox_results[4,]<0.05]
dim(Results.dfi)

#PFI PFI.time
surPFI.tumor=subset(CRC.survivaldata,PFI.time > 0 & PFI >= 0)
idsur=intersect(rownames(te.tumorexp.sel),rownames(surPFI.tumor))
length(idsur)
cbOSdata=cbind(surPFI.tumor[idsur,8:9],te.tumorexp.sel[idsur,])     #####590 SAMPPLES
##

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$PFI,  #CHANGE
                             os=cbOSdata$PFI.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.pfi=cox_results[,cox_results[4,]<0.05]
dim(Results.pfi)

save(Results.os,Results.dfi,Results.pfi,Results.dss,
     file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/survival.screen.subfamily.level.1072s.RData")

###screen at family level
#######
###############
tefamily.tumorexp.sel=tefamily.tumorexp.filterd[, colnames(tefamily.tumorexp.filterd) %in% cdrep.sel$repFamily ]
dim(tefamily.tumorexp.sel)

#screen for OS
#OS
surOS.tumor=subset(CRC.survivaldata,OS.time > 0 & OS>= 0)
idsur=intersect(rownames(tefamily.tumorexp.sel),rownames(surOS.tumor))
length(idsur)
cbOSdata=cbind(surOS.tumor[idsur,2:3],tefamily.tumorexp.sel[idsur,])

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$OS,  #CHANGE
                             os=cbOSdata$OS.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.os.family=cox_results[,cox_results[4,]<0.05]
dim(Results.os.family)

#DSS
surDSS.tumor=subset(CRC.survivaldata,DSS.time > 0 & DSS >=0)
idsur=intersect(rownames(tefamily.tumorexp.sel),rownames(surDSS.tumor))
length(idsur)
cbOSdata=cbind(surDSS.tumor[idsur,4:5],tefamily.tumorexp.sel[idsur,])     ####570 SAMPLES

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$DSS,  #CHANGE
                             os=cbOSdata$DSS.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.dss.family=cox_results[,cox_results[4,]<0.05]
dim(Results.dss.family)

#DFI
surDFI.tumor=subset(CRC.survivaldata,DFI.time > 0 & DFI >=0)
idsur=intersect(rownames(tefamily.tumorexp.sel),rownames(surDFI.tumor))
length(idsur)
cbOSdata=cbind(surDFI.tumor[idsur,6:7],tefamily.tumorexp.sel[idsur,])     #####234 SAMPLES

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$DFI,  #CHANGE
                             os=cbOSdata$DFI.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.dfi.family=cox_results[,cox_results[4,]<0.05]
dim(Results.dfi.family)

#PFI PFI.time
surPFI.tumor=subset(CRC.survivaldata,PFI.time > 0 & PFI >= 0)
idsur=intersect(rownames(tefamily.tumorexp.sel),rownames(surPFI.tumor))
length(idsur)
cbOSdata=cbind(surPFI.tumor[idsur,8:9],tefamily.tumorexp.sel[idsur,])     #####590 SAMPPLES
##

##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$PFI,  #CHANGE
                             os=cbOSdata$PFI.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.pfi.family=cox_results[,cox_results[4,]<0.05]
dim(Results.pfi.family)

save(Results.os.family,Results.dfi.family,Results.pfi.family,Results.dss.family,
     file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/survival.screen.family.level.47s.RData")

###screen at class level
#######
###############
teclass.tumorexp.sel=teclass.tumorexp.filterd[, colnames(teclass.tumorexp.filterd) %in% cdrep.sel$repClass ]

#screen for OS
#OS
surOS.tumor=subset(CRC.survivaldata,OS.time > 0 & OS>= 0)
idsur=intersect(rownames(teclass.tumorexp.sel),rownames(surOS.tumor))
length(idsur)
cbOSdata=cbind(surOS.tumor[idsur,2:3],teclass.tumorexp.sel[idsur,])
#
##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$OS,  #CHANGE
                             os=cbOSdata$OS.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.os.class=cox_results
dim(Results.os.class)


#DSS
surDSS.tumor=subset(CRC.survivaldata,DSS.time > 0 & DSS >=0)
idsur=intersect(rownames(teclass.tumorexp.sel),rownames(surDSS.tumor))
length(idsur)
cbOSdata=cbind(surDSS.tumor[idsur,4:5],teclass.tumorexp.sel[idsur,])     ####570 SAMPLES
#
##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$DSS,  #CHANGE
                             os=cbOSdata$DSS.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.dss.class=cox_results
dim(Results.dss.class)

#DFI
surDFI.tumor=subset(CRC.survivaldata,DFI.time > 0 & DFI >=0)
idsur=intersect(rownames(teclass.tumorexp.sel),rownames(surDFI.tumor))
length(idsur)
cbOSdata=cbind(surDFI.tumor[idsur,6:7],teclass.tumorexp.sel[idsur,])     #####234 SAMPLES
#
##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$DFI,  #CHANGE
                             os=cbOSdata$DFI.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.dfi.class=cox_results
dim(Results.dfi.class)

#PFI PFI.time
surPFI.tumor=subset(CRC.survivaldata,PFI.time > 0 & PFI >= 0)
idsur=intersect(rownames(),rownames(surPFI.tumor))
length(idsur)
cbOSdata=cbind(surPFI.tumor[idsur,8:9],teclass.tumorexp.sel[idsur,])     #####590 SAMPPLES
##
##########  useful
cox_results <- apply(cbOSdata[,-c(1,2)] , 2, function(values1){
  group=ifelse(values1>median(values1),'Migh','low')
  #group=values1
  survival_dat <- data.frame(group=group,
                             os.status=cbOSdata$PFI,  #CHANGE
                             os=cbOSdata$PFI.time,   #CHANGE
                             #age=clin$age,
                             #gender=gender,
                             #stage=stage,
                             #smoking=smoking,
                             stringsAsFactors = F)
  library(survival)
  my.surv <- Surv(survival_dat$os,survival_dat$os.status)
  #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
  m=coxph(my.surv ~ group, data =  survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['groupMigh',])
  
})
Results.pfi.class=cox_results
dim(Results.pfi.class)

save(Results.os.class,Results.dfi.class,Results.pfi.class,Results.dss.class,
     file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/survival.screen.class.level.6s.RData")

######subfamily level
#########plot the data
###############
#
sursc=as.data.frame(t(Results.dss))
surstat=cbind(sursc, repinfo.sel[rownames(sursc),])
surstat$endpoint=rep("dss", times=nrow(surstat))
surstat.dss=surstat
#
sursc=as.data.frame(t(Results.os))
surstat=cbind(sursc, repinfo.sel[rownames(sursc),])
surstat$endpoint=rep("os", times=nrow(surstat))
surstat.os=surstat
#
sursc=as.data.frame(t(Results.dfi))
surstat=cbind(sursc, repinfo.sel[rownames(sursc),])
surstat$endpoint=rep("dfi", times=nrow(surstat))
surstat.dfi=surstat
#
sursc=as.data.frame(t(Results.pfi))
surstat=cbind(sursc, repinfo.sel[rownames(sursc),])
surstat$endpoint=rep("pfi", times=nrow(surstat))
surstat.pfi=surstat
###
cb.data=rbind(surstat.dfi, surstat.dss, surstat.os, surstat.pfi)
table(cb.data$repClass,cb.data$endpoint)

# stacked bar chart; Barplot of the count
library(easyGgplot2)
#cbPalette= c("#984ea3","#377eb8","#e41a1c","#4daf4a","#ff7f00","grey")
cbPalette= c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")
pc=ggplot2.barplot(data=cb.data, xName="endpoint",
                   #brewerPalette="Blues",
                   groupName="repClass")+
  scale_fill_manual(values= cbPalette)+theme_classic()+
  ylab("No.TEs (log rank p value < 0.05)")+xlab("survival endpoint")
generate.PDF <- function(fig) {
  pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/stacked.plot.endpoint.pdf",width = 4,height = 5)
  print(pc)
  dev.off()
}
generate.PDF(fig)
ggplot2.barplot(data=cb.data, xName="endpoint",
                   #brewerPalette="Blues",
                   groupName="repClass")+
  scale_fill_manual(values= cbPalette)+theme_classic()+
  ylab("No.TEs (log rank p value < 0.05)")+xlab("survival endpoint")

library(ggpubr)
ggplot(cb.data, aes_string(x="repClass", y="HR", group="repClass")) + 
    geom_boxplot(aes(fill = repClass), alpha = 0.5,outlier.colour = "black",outlier.size = 1,notch = F)+
    stat_compare_means(label = "p.signif",method = "kruskal.test")+
    geom_jitter(aes(color = repClass),width = 0.15,size=0.3)+
    facet_grid(. ~ endpoint)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=cbPalette)+ # Boxplot fill color
    scale_color_manual(values = cbPalette)+ylab("Hazard ratio")

ph=ggplot(cb.data, aes_string(x="repClass", y="HR", group="repClass")) + 
    geom_boxplot(aes(fill = repClass), alpha = 0.5,outlier.colour = "black",outlier.size = 1,notch = F)+
    stat_compare_means(label = "p.signif",method = "kruskal.test")+
    geom_jitter(aes(color = repClass),width = 0.15,size=0.3)+
    facet_grid(. ~ endpoint)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=cbPalette)+ # Boxplot fill color
    scale_color_manual(values = cbPalette)+ylab("Hazard ratio")+ggtitle("survival at subfamily level")
generate.PDF <- function(fig) {
  pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/HR.endpoint.subfamily.pdf",width = 6,height = 5)
  print(ph)
  dev.off()
}
generate.PDF(fig)

######family level
#########plot the data
###############
#
repf=cdrep.sel[!(duplicated(cdrep.sel$repFamily)),]
rownames(repf)=repf$repFamily
####
sursc=as.data.frame(t(Results.dss.family))
surstat=cbind(sursc, repf[rownames(sursc),])
surstat$endpoint=rep("dss", times=nrow(surstat))
surstat.dss.family=surstat
#
sursc=as.data.frame(t(Results.os.family))
surstat=cbind(sursc, repf[rownames(sursc),])
surstat$endpoint=rep("os", times=nrow(surstat))
surstat.os.family=surstat
#
sursc=as.data.frame(t(Results.dfi.family))
surstat=cbind(sursc, repf[rownames(sursc),])
surstat$endpoint=rep("dfi", times=nrow(surstat))
surstat.dfi.family=surstat
#
sursc=as.data.frame(t(Results.pfi.family))
surstat=cbind(sursc, repf[rownames(sursc),])
surstat$endpoint=rep("pfi", times=nrow(surstat))
surstat.pfi.family=surstat
###
cb.data.family=rbind(surstat.dfi.family, surstat.dss.family, surstat.os.family, surstat.pfi.family)
#cb.data.family$repClass=gsub(pattern = "[?]","_",cb.data.family$repClass)
#cb.data.family=subset(cb.data.family, !(repClass=="SINE_"))
head(cb.data.family)

library(ggpubr)
cbPalette= c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","gray")


ggplot(cb.data.family, aes_string(x="repClass", y="HR", group="repClass")) + 
    geom_boxplot(aes(fill = repClass), alpha = 0.5,outlier.colour = "black",outlier.size = 1,notch = F)+
    stat_compare_means(label = "p.signif",method = "kruskal.test")+
    geom_jitter(aes(color = repClass),width = 0.15,size=0.3)+
    facet_grid(. ~ endpoint)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=cbPalette)+ # Boxplot fill color
    scale_color_manual(values = cbPalette)+ylab("Hazard ratio")

ph=ggplot(cb.data.family, aes_string(x="repClass", y="HR", group="repClass")) + 
    geom_boxplot(aes(fill = repClass), alpha = 0.5,outlier.colour = "black",outlier.size = 1,notch = F)+
    stat_compare_means(label = "p.signif",method = "kruskal.test")+
    geom_jitter(aes(color = repClass),width = 0.15,size=0.3)+
    facet_grid(. ~ endpoint)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=cbPalette)+ # Boxplot fill color
    scale_color_manual(values = cbPalette)+ylab("Hazard ratio")
generate.PDF <- function(fig) {
  pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/HR.endpoint.family.pdf",width = 6,height = 5)
  print(ph)
  dev.off()
}
generate.PDF(fig)


 ggplot(cb.data.family, aes( repFamily,HR, colour = repClass)) + 
       geom_pointrange(aes(ymin = HRCILL, ymax = HRCIUL))+
       facet_grid(. ~ endpoint)+
       scale_fill_manual(values=cbPalette)+ # Boxplot fill color
       scale_color_manual(values = cbPalette)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+rotate()+ylab("Hazard ratio")
#
 phf=ggplot(cb.data.family, aes( repFamily,HR, colour = repClass)) + 
       geom_pointrange(aes(ymin = HRCILL, ymax = HRCIUL),fatten = 2)+
       facet_grid(. ~ endpoint)+
       scale_fill_manual(values=cbPalette)+ # Boxplot fill color
       scale_color_manual(values = cbPalette)+
       #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
       rotate()+ylab("Hazard ratio")+ggtitle("survival at repFamily level (n = 47)")
generate.PDF <- function(fig) {
  pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/HR.endpoint.family.indi.pdf",
      width = 8,height = 4)
  print(phf)
  dev.off()
}
generate.PDF(fig)

#####class level
####
surstat=as.data.frame(t(Results.dss.class))
surstat$endpoint=rep("dss", times=nrow(surstat))
surstat$repClass=rownames(surstat)
surstat.dss.class=surstat
#
surstat=as.data.frame(t(Results.os.class))
surstat$endpoint=rep("os", times=nrow(surstat))
surstat$repClass=rownames(surstat)
surstat.os.class=surstat
#
surstat=as.data.frame(t(Results.dfi.class))
surstat$endpoint=rep("dfi", times=nrow(surstat))
surstat$repClass=rownames(surstat)
surstat.dfi.class=surstat
#
surstat=as.data.frame(t(Results.pfi.class))
surstat$endpoint=rep("pfi", times=nrow(surstat))
surstat$repClass=rownames(surstat)
surstat.pfi.class=surstat
###
cb.data.class=rbind(surstat.dfi.class, surstat.dss.class, surstat.os.class, surstat.pfi.class)
head(cb.data.class)
table(cb.data.class$repClass)

 ggplot(cb.data.class, aes( repClass,HR, colour = repClass)) + 
       geom_pointrange(aes(ymin = HRCILL, ymax = HRCIUL))+
       facet_grid(. ~ endpoint)+
       scale_fill_manual(values=cbPalette)+ # Boxplot fill color
       scale_color_manual(values = cbPalette)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+rotate()+ylab("Hazard ratio")
#
 phf=ggplot(cb.data.class, aes( repClass,HR, colour = repClass)) + 
       geom_pointrange(aes(ymin = HRCILL, ymax = HRCIUL),fatten = 2)+
       facet_grid(. ~ endpoint)+
       scale_fill_manual(values=cbPalette)+ # Boxplot fill color
       scale_color_manual(values = cbPalette)+
       #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
       rotate()+ylab("Hazard ratio")+ggtitle("survival at class level (n = 6)")
generate.PDF <- function(fig) {
  pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/HR.endpoint.class.pdf",
      width = 8,height = 4)
  print(phf)
  dev.off()
}
generate.PDF(fig)

####check the overlapped te among these four endpoints
####
library(VennDiagram)
OS=subset(cb.data, endpoint=="os")$repName
DSS=subset(cb.data, endpoint=="dss")$repName
DFI=subset(cb.data, endpoint=="dfi")$repName
PFI=subset(cb.data, endpoint=="pfi")$repName
#
venn.diagram(x= list(OS = OS, DSS = DSS,DFI = DFI, PFI=PFI), compression ="lzw",main="sig TEs overlapped among four endpoints",main.cex = 0.45,
             filename = "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/venn.subfamily.tif", height = 450,
             width = 450,resolution =300,  col ="transparent", 
             fill =c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"),alpha = 0.5, 
             cex = 0.45,fontfamily = "serif", fontface = "bold",
             #cat.col =c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"), 
             cat.cex = 0.45,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif", 
             rotation.degree = 0)

####check what these overlapped TE belong to 
#######
teso=Reduce(intersect, list(OS,DSS,DFI,PFI))
sig.te.overlapped=cdrep[rownames(cdrep) %in% teso,]
write.csv(sig.te.overlapped,
          "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/sig.te.four.endpoints.overlapped.csv")

save(cb.data, cb.data.class, cb.data.family, 
    file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen/sur.res.combined.endpoints.RData")

setwd("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.survival.TE.screen")


