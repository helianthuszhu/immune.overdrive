cbOSdata[1:4,1:4]

te.tumorexp.sel=te.tumorexp[, colnames(te.tumorexp) %in% rownames(cdrep.sel)]
te.tumorexp.sel[1:4,1:4]

head(CRC.survivaldata)
#####
da1exp=te.tumorexp.sel
da1sur=CRC.survivaldata
######shuffle expMatrix data
######
datalist.os=list()
datalist.dss=list()
datalist.dfi=list()
datalist.pfi=list()
#
for (i in 1:100) {
  set.seed(i)
  rownames(da1exp)=sample(rownames(da1exp))
  ####
  #screen for OS
  #OS
  surOS.tumor=subset(CRC.survivaldata,OS.time > 0 & OS>= 0)
  idsur=intersect(rownames(da1exp),rownames(surOS.tumor))
  length(idsur)
  cbOSdata=cbind(surOS.tumor[idsur,2:3],da1exp[idsur,])
  
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
  Results.os.shuffle=cox_results[,cox_results[4,]<0.05]
  Results.os.shuffle=as.data.frame(t(Results.os.shuffle))
  Results.os.shuffle$shuffle.id=rep(i, times=nrow(Results.os.shuffle))
  dim(Results.os.shuffle)
  datalist.os[[i]]=Results.os.shuffle
  
  #DSS
  surDSS.tumor=subset(CRC.survivaldata,DSS.time > 0 & DSS >=0)
  idsur=intersect(rownames(da1exp),rownames(surDSS.tumor))
  length(idsur)
  cbOSdata=cbind(surDSS.tumor[idsur,4:5],da1exp[idsur,])     ####570 SAMPLES
  
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
  Results.dss.shuffle=cox_results[,cox_results[4,]<0.05]
  dim(Results.dss.shuffle)
  Results.dss.shuffle=as.data.frame(t(Results.dss.shuffle))
  Results.dss.shuffle$shuffle.id=rep(i, times=nrow(Results.dss.shuffle))
  dim(Results.dss.shuffle)
  datalist.dss[[i]]=Results.dss.shuffle
  
  #DFI
  surDFI.tumor=subset(CRC.survivaldata,DFI.time > 0 & DFI >=0)
  idsur=intersect(rownames(da1exp),rownames(surDFI.tumor))
  length(idsur)
  cbOSdata=cbind(surDFI.tumor[idsur,6:7],da1exp[idsur,])     #####234 SAMPLES
  
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
  Results.dfi.shuffle=cox_results[,cox_results[4,]<0.05]
  Results.dfi.shuffle=as.data.frame(t(Results.dfi.shuffle))
  Results.dfi.shuffle$shuffle.id=rep(i, times=nrow(Results.dfi.shuffle))
  dim(Results.dfi.shuffle)
  datalist.dfi[[i]]=Results.dfi.shuffle
  
  #PFI PFI.time
  surPFI.tumor=subset(CRC.survivaldata,PFI.time > 0 & PFI >= 0)
  idsur=intersect(rownames(da1exp),rownames(surPFI.tumor))
  length(idsur)
  cbOSdata=cbind(surPFI.tumor[idsur,8:9],da1exp[idsur,])     #####590 SAMPPLES
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
  Results.pfi.shuffle=cox_results[,cox_results[4,]<0.05]
  Results.pfi.shuffle=as.data.frame(t(Results.pfi.shuffle))
  Results.pfi.shuffle$shuffle.id=rep(i, times=nrow(Results.pfi.shuffle))
  dim(Results.pfi.shuffle)
  datalist.pfi[[i]]=Results.pfi.shuffle

}
#############
#############combine the result
#
Results.os.shuffle.agg=do.call(rbind,datalist.os)
Results.os.shuffle.agg$repName=rownames(Results.os.shuffle.agg)

Results.dss.shuffle.agg=do.call(rbind,datalist.dss)
Results.dss.shuffle.agg$repName=rownames(Results.dss.shuffle.agg)

Results.dfi.shuffle.agg=do.call(rbind,datalist.dfi)
Results.dfi.shuffle.agg$repName=rownames(Results.dfi.shuffle.agg)

Results.pfi.shuffle.agg=do.call(rbind,datalist.pfi)
Results.pfi.shuffle.agg$repName=rownames(Results.pfi.shuffle.agg)
#
###find the overlapped with ture significant TEs
head(cb.data)
#os
Results.os.shuffle.agg.overlapped=Results.os.shuffle.agg[rownames(Results.os.shuffle.agg) %in% rownames(subset(cb.data, endpoint=="os")),]
#dss
Results.dss.shuffle.agg.overlapped=Results.dss.shuffle.agg[rownames(Results.dss.shuffle.agg) %in% rownames(subset(cb.data, endpoint=="dss")),]
#dfi
Results.dfi.shuffle.agg.overlapped=Results.dfi.shuffle.agg[rownames(Results.dfi.shuffle.agg) %in% rownames(subset(cb.data, endpoint=="dfi")),]
#pfi
Results.pfi.shuffle.agg.overlapped=Results.pfi.shuffle.agg[rownames(Results.pfi.shuffle.agg) %in% rownames(subset(cb.data, endpoint=="pfi")),]
#######

#
shuffle.stat.os=as.data.frame(table(Results.pfi.shuffle.agg.overlapped$shuffle.id))
colnames(shuffle.stat.os)=c("shuffleid","count")
shuffle.stat.os$ratio=shuffle.stat.os$count/length(rownames(subset(cb.data, endpoint=="os")))
shuffle.stat.os$endpoint=rep("os",times=nrow(shuffle.stat.os))
shuffle.stat.os
#
shuffle.stat.dss=as.data.frame(table(Results.dss.shuffle.agg.overlapped$shuffle.id))
colnames(shuffle.stat.dss)=c("shuffleid","count")
shuffle.stat.dss$ratio=shuffle.stat.dss$count/length(rownames(subset(cb.data, endpoint=="dss")))
shuffle.stat.dss$endpoint=rep("dss",times=nrow(shuffle.stat.dss))
shuffle.stat.dss
#
shuffle.stat.dfi=as.data.frame(table(Results.dfi.shuffle.agg.overlapped$shuffle.id))
colnames(shuffle.stat.dfi)=c("shuffleid","count")
shuffle.stat.dfi$ratio=shuffle.stat.dfi$count/length(rownames(subset(cb.data, endpoint=="dfi")))
shuffle.stat.dfi$endpoint=rep("dfi",times=nrow(shuffle.stat.dfi))
shuffle.stat.dfi
#
shuffle.stat.pfi=as.data.frame(table(Results.pfi.shuffle.agg.overlapped$shuffle.id))
colnames(shuffle.stat.pfi)=c("shuffleid","count")
shuffle.stat.pfi$ratio=shuffle.stat.pfi$count/length(rownames(subset(cb.data, endpoint=="pfi")))
shuffle.stat.pfi$endpoint=rep("pfi",times=nrow(shuffle.stat.pfi))
shuffle.stat.pfi
####
shuffle.stat.CB=rbind(shuffle.stat.os,shuffle.stat.dss,shuffle.stat.dfi,shuffle.stat.pfi)
shuffle.stat.CB$ratio=shuffle.stat.CB$ratio*100
head(shuffle.stat.CB)


#############
#BiocManager::install("tidyverse")
library(ggplot2)
library(ggridges)
library(tidyverse)
library(dplyr)
pdf("data/te/stat20200428/accuracy.of.survival.shuffle.pdf",width = 6,height = 3)
ggplot(shuffle.stat.CB, aes(x = ratio, y = endpoint)) +
  geom_density_ridges(aes(fill = endpoint),quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values = c("os"="#66C2A5","dss"="#FC8D62","dfi"="#8DA0CB","pfi"="#E78AC3"))+
  geom_text(data=shuffle.stat.CB %>% group_by(endpoint) %>% 
              summarise(ratio=median(ratio)),
            aes(label=sprintf("%1.1f", ratio)),
            position=position_nudge(y=-0.1), colour="black", size=3.5)+
  theme_classic()+scale_x_continuous(breaks=c(0,1,2,4,8,12))+
  xlab("accuracy ratio")
dev.off()
#
save(shuffle.stat.CB, file = "data/te/stat20200428/accuracy.of.survival.shuffle.RData")
####################
######################
#####################
##################shuffle for immune screen
#################
dim(cor.score)
dim(cor.exp)
cor.score[1:4,1:4]
cor.exp[1:4,1:4]
length(intersect(rownames(cor.score), rownames(cor.exp)))
####shuffle the cor.score sample id
cor.score.shuffled=cor.score
datalist.immune.shuffle=list()
for (k in 1:100){
  set.seed(k)
  rownames(cor.score.shuffled)=sample(rownames(cor.score.shuffled))
  cor.exp.shuffled=cor.exp[rownames(cor.score.shuffled),]
  #
  datalist=list()
  oklist=list()
  for (i in 1:ncol(cor.exp.shuffled)) {
    for (j in 1:ncol(cor.score.shuffled)) {
      res <- cor.test(cor.exp.shuffled[,i], cor.score.shuffled[,j],  method = "spearman")
      pvalue=res$p.value
      cor=res$estimate
      res=data.frame(pvalue, cor)
      rownames(res)=colnames(cor.score.shuffled)[j]
      rownames(res)=paste(colnames(cor.exp.shuffled)[i], rownames(res),sep = "&")
      datalist[[j]] <- res
    }
    oklist[[i]]=do.call(rbind, datalist)
    shuffleresimmune=do.call(rbind, oklist)
    shuffleresimmune$shuffle.id=rep(k, times=nrow(shuffleresimmune))
  }
  datalist.immune.shuffle[[k]]=shuffleresimmune
}

Results.immune.shuffle=do.call(rbind, datalist.immune.shuffle)
save(Results.immune.shuffle, file = "data/te/stat20200428/accuracy.of.immune.shuffle.RData")




head(Results.immune.shuffle)
table(Results.immune.shuffle$shuffle.id)
summary(Results.immune.shuffle$cor)
#
