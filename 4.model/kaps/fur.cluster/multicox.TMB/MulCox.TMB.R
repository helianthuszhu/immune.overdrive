#multiCox for TMB
head(clin.heat)
head(stat.clin)
dim(clin.heat)
dim(stat.clin)
#
length(intersect(rownames(clin.heat), rownames(stat.clin)))
#########
#########
multicox.data.TMB=cbind(TMB=clin.heat$TMB,totalM=clin.heat$total_mut, stat.clin[rownames(clin.heat),])
multicox.data.TMB$totalM=as.numeric(paste(multicox.data.TMB$totalM))
multicox.data.TMB$TMB.log=log2(multicox.data.TMB$totalM)
table(multicox.data.TMB$kaps.group)
#
multicox.data.TMB.sub=subset(multicox.data.TMB, kaps.group=="set4"|kaps.group=="set3")
multicox.data.TMB.sub=subset(multicox.data.TMB.sub, TMB.log>0)
#######
summary(multicox.data.TMB.sub$TMB.log)
table(multicox.data.TMB.sub$kaps.group)
####
chisq.test(multicox.data.TMB.sub$kaps.group, multicox.data.TMB.sub$MSI.status.bin,correct = T)
#multicox.data.TMB=subset(multicox.data.TMB,  totalM> 0)
#
#multicox.data.TMB$TMB.group=ifelse(multicox.data.TMB$TMB> median(multicox.data.TMB$TMB), "TMB.high","TMB.low")
#########
#multicox
my.surv=Surv(multicox.data.TMB.sub$PFI.time,multicox.data.TMB.sub$PFI)
multicox=coxph(my.surv ~ kaps.group+gender+ age.bi+lymphatic_invasion+location.bi+AJCC.stage.bi, data =  multicox.data.TMB.sub)
#####
cc=summary(multicox)
cc
HR <- round(cc$coefficients[,2],2)
P_Value <- round(cc$coefficients[,5],4)
LCI <- round(cc$conf.int[,3],2)
UCI <- round(cc$conf.int[,4],2)
CI95 <- paste(LCI,'-',UCI)
dfmTMB <- cbind(Hazard_Ratio = HR,CI95 = CI95,P_Value = P_Value)
dfmTMB
#rownames(dfm)=c("risk.group","age","gender","MSI.status","AJCC.stage")
colnames(dfmTMB)=paste(colnames(dfmTMB),"PFI",sep = ".")
MulCox.TMB.PFI=as.data.frame(dfmTMB)
write.csv(MulCox.TMB.PFI,"4.model/kaps/fur.cluster/multicox.TMB/MulCox.TMB.PFI.res.csv")
