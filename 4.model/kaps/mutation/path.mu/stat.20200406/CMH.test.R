######
head(caldata.pmu)
aa=caldata.pmu
dim(aa)
aa=data.frame(msi.status=kaps.td[rownames(aa),]$MSI.status.bin, aa)
class(aa$Cell.Cycle)
aa$Cell.Cycle=ifelse(aa$Cell.Cycle=="0" |aa$Cell.Cycle=="1", aa$Cell.Cycle,NA)
head(aa)
#
pathidx=colnames(aa)[-c(1,2)]
datalist=list()
for (i in 1:10) {

#########
cmhtestdata=aa[,c(1,2,i+2)] %>% group_by_all() %>% summarise(COUNT = n()) 
cmhtestdata=na.omit(cmhtestdata)
cmhtestdata
#cmhtestdata$msi.status = factor(cmhtestdata$msi.status,levels=unique(cmhtestdata$msi.status))
#cmhtestdata$kaps.group = factor(cmhtestdata$kaps.group,levels=unique(cmhtestdata$kaps.group))
#cmhtestdata$HIPPO = factor(cmhtestdata$HIPPO,levels=unique(cmhtestdata$HIPPO))
#
#cmhtestf = xtabs(paste0("COUNT ~ msi.status + kaps.group +", pathidx[i]),data=cmhtestdata)

cmhtestf = xtabs(paste0("COUNT ~", pathidx[i] , "+ kaps.group + msi.status") ,data=cmhtestdata)

#cmhtestf = xtabs(paste0("COUNT ~", pathidx[i] , " + msi.status + kaps.group") ,data=cmhtestdata)

########
###  Note that the grouping variable is last in the xtabs function
ftable(cmhtestf)                     # Display a flattened table
cmhtestf.res1=mantelhaen.test(cmhtestf)
library(vcd)
cmhtestf.res2=woolf_test(cmhtestf)
###########
ress=data.frame(cmh.pvalue=cmhtestf.res1$p.value,
                woolf.pvalue=cmhtestf.res2$p.value)
rownames(ress)=paste(pathidx[i])
datalist[[i]]=ress
}
#
cmhtest.result=do.call(rbind, datalist)
###########
cmhtest.result=data.frame(cmhtest.result, chisq.p.value=dfm.pmu$pvalue$pvalue)
#cmhtest.result=as.data.frame(cmhtest.result)
#colnames(cmhtest.result)[3]="chisq.p.value"
cmhtest.result

############
####p value for msi
dfm.msi=NULL
for (k in 3:ncol(aa)){
  p.value=chisq.test(aa[,k],aa$msi.status,correct = T)
  val=as.data.frame(p.value$p.value)
  rownames(val)=colnames(aa)[k]
  dfm.msi=rbind(dfm.msi,val)
}
head(dfm.msi)
colnames(dfm.msi)[1]="chisq.p.value.msi"
########
cmhtest.result=data.frame(cmhtest.result, chisq.p.value.msi=dfm.msi$chisq.p.value.msi)
cmhtest.result=round(cmhtest.result,digits = 3)
write.csv(cmhtest.result,"4.model/kaps/mutation/path.mu/stat.20200406/cmh.test.combined.csv")
