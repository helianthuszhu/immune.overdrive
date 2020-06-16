#####
####pan cancer TE expression matrix
#######
library(Biobase)
library(dplyr)
library(tidyr)
aa=readRDS("~/nas/Xiaoqiang/R.pacakages/TE/REdiscoverTEdata/inst/Fig4_data/eset_TCGA_TE_intergenic_logCPM.RDS")
rexp=as.data.frame(t(exprs(aa)))
rexp[1:4,1:4]
rclin=pData(aa)
rclin$control=substr(rownames(rclin),14,15)
rclin$control=ifelse(rclin$control=="01","tumor","normal")
rclin$id=substr(rownames(rclin), 1,15)
rclin$id2=substr(rownames(rclin), 16,16)
head(rclin)
table(rclin$id2)
######
bb=cbind(rclin$id, rexp[rownames(rclin),])
colnames(bb)[1]="id"
bb[1:4,1:4]
rexp.agg= bb %>% group_by(id) %>% summarise_all(mean)
rexp.agg=as.data.frame(rexp.agg)
rownames(rexp.agg)=rexp.agg$id
rexp.agg=rexp.agg[,-1]
rexp.agg[1:4,1:4]
dim(rexp.agg)
#######
rclin.agg=rclin[!(duplicated(rclin$id)),]
rownames(rclin.agg)=rclin.agg$id
rclin.agg=rclin.agg[rownames(rexp.agg),]
rclin.agg.tumor=subset(rclin.agg, control=="tumor")
head(rclin.agg.tumor)
###
table(rclin.agg$id2)
#######load TE used
load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/1.TE.used/TE.label.used.RData")
head(cdrep)
te.used=subset(cdrep, !(repClass.new2=="other.repeats"))
dim(te.used)
###############
rexp.agg.T=rexp.agg[rownames(rclin.agg.tumor),colnames(rexp.agg) %in% rownames(te.used)]
rexp.agg.T[1:4,1:4]
dim(rexp.agg.T)
##############
library("FactoMineR")
library("factoextra")
iris.pca.clean <- PCA(as.matrix(rexp.agg.T), graph = FALSE)
ff=fviz_pca_ind(iris.pca.clean,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = rclin.agg.tumor$indication, # color by groups
             #palette = c("#377eb8","#984ea3"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "batch",title="PCA of pan.cancer (N=6274)")
#table(rclin.agg$indication, rclin.agg$control)
generate.PDF <- function(fig) {
  pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.PCA/PCA.total.pdf",width = 12,height = 6)
  print(ff)
  dev.off()
}
generate.PDF(fig)
save(rexp.agg.T, rclin.agg.tumor,file = "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.PCA/pca.data.RData")
###########
#######
#########survival screen
surpan=read.csv("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/clinical.data/pan/pan.sur.clinical.full.info.csv",header = T)
rownames(surpan)=surpan$bcr_patient_barcode
head(surpan)
dim(surpan)
######
cc=rclin.agg.tumor
rownames(cc)=cc$patient_ID
head(cc)
ids=intersect(rownames(cc), rownames(surpan))
###
stat.surclin=cbind(surpan[ids,], cc[ids,])
rownames(stat.surclin)=stat.surclin$id
head(stat.surclin)
###
table(stat.surclin$indication,stat.surclin$OS)
stat.rexp=rexp.agg.T[rownames(stat.surclin),]
####
library(survival)
library(survminer)

#da2=list()
surdata=data.frame(time=stat.surclin$OS.time,
                   status=stat.surclin$OS,
                   type=stat.surclin$indication)
rownames(surdata)=rownames(stat.surclin)
head(surdata)
idx=unique(surdata$type)[-17]
#
datalist=list()
for (j in 1:length(idx)) {
    a1=subset(surdata, type==idx[j],)
    a1$time=as.numeric(paste(a1$time))
    a1$status=as.numeric(paste(a1$status))
    a2=stat.rexp[rownames(a1),]
    #
    cox_results <- apply(a2 , 2, function(values1){
      #group=ifelse(values1>median(values1),'Migh','low')
      group=values1
      survival_dat <- data.frame(exp=group,
                                 time=a1$time,
                                 status=a1$status,
                                 stringsAsFactors = F)
      library(survival)
      my.surv <- Surv(survival_dat$time,survival_dat$status)
      #m=coxph(my.surv ~ group+age+gender+stage+smoking, data =  survival_dat)
      m=coxph(my.surv ~ exp, data =  survival_dat)
      beta <- coef(m)
      se <- sqrt(diag(vcov(m)))
      HR <- exp(beta)
      HRse <- HR * se
      
      #summary(m)
      tmp <- data.frame(coef = beta, 
                        se = se,
                        z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                        HR = HR, HRse = HRse,
                        HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                        HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                        HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)
      )
      return(tmp)
    })
    coef.out=do.call(rbind,cox_results)
    rownames(coef.out)=paste(rownames(coef.out), idx[j],sep = "&")
    datalist[[j]]=coef.out
  }
coef.out.pan=do.call(rbind,datalist)
coef.out.pan$id=rownames(coef.out.pan)
stat.coef=separate(coef.out.pan, col = id, into = c('TE', 'type'),sep = "&")
stat.coef.sig=subset(stat.coef, p < 0.05 & HR >1)
length(unique(stat.coef.sig$TE))
dim(stat.coef.sig)
head(stat.coef.sig)
summary(stat.coef.sig$p)
sss=as.data.frame.matrix(table(stat.coef.sig$TE, stat.coef.sig$type))
sss$count=rowSums(sss)
sss=sss[order(sss$count,decreasing = T),]
sssa=data.frame(counts=colSums(sss))
head(sss)
#

