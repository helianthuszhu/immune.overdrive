##########sva survival in pan cancer#######
###########
table(stat.surclin$type)
dim(stat.rexp)
dim(stat.surclin)
###########
svaexp.pan=stat.rexp[, grep("SVA",colnames(stat.rexp))]
####
length(intersect(rownames(stat.rexp), rownames(stat.surclin)))
#
stat.sva.pan=cbind(svaexp.pan[rownames(stat.surclin),], stat.surclin)
stat.sva.pan$OS=as.numeric(paste(stat.sva.pan$OS))
stat.sva.pan$DSS=as.numeric(paste(stat.sva.pan$DSS))
stat.sva.pan$DFI=as.numeric(paste(stat.sva.pan$DFI))
stat.sva.pan$PFI=as.numeric(paste(stat.sva.pan$PFI))
stat.sva.pan$OS.time=as.numeric(paste(stat.sva.pan$OS.time))
stat.sva.pan$DSS.time=as.numeric(paste(stat.sva.pan$DSS.time))
stat.sva.pan$DFI.time=as.numeric(paste(stat.sva.pan$DFI.time))
stat.sva.pan$PFI.time=as.numeric(paste(stat.sva.pan$PFI.time))
save(stat.sva.pan, file = "2.survival.TE.screen/sva.survival/sva.pan.vali/stat.sva.pan.RData")
head(stat.sva.pan)
svaid=colnames(stat.sva.pan)[1:6]
idtype=unique(stat.sva.pan$type)
datalist=list()
datalist2=list()
for (j in 1:length(idtype)) {
  
  for (i in 1:6) {
    aa=subset(stat.sva.pan, type==idtype[j])
    
    ss=Surv(aa$PFI.time, aa$PFI)
    aa$risk.group.sva=ifelse(aa[, colnames(aa) %in% svaid[i]]> median(aa[, colnames(aa) %in% svaid[i]]),"2high","1low")
    cox = summary(coxph(ss~aa$risk.group.sva))
    pvalue=cox$sctest['pvalue']
    hr = round(cox$conf.int[1],2)
    hr_left = round(cox$conf.int[3],2)
    hr_right = round(cox$conf.int[4],2)
    conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep=""); 
    #txt = paste("HR = ", hr, conf_int, "\nlogrank P = ", signif(pvalue, 2), sep="")
    #text(grconvertX(0.98, "npc"), grconvertY(.97, "npc"),labels=txt, adj=c(1, 1))
    res=data.frame(pvalue, hr, hr_left, hr_right)
    colnames(res)=c("pvalue", "hr", "hr_left", "hr_right")
    rownames(res)=paste(svaid[i],"OS",sep="&")
    datalist[[i]]=res
  }
  res2=do.call(rbind, datalist)
  rownames(res2)=paste(rownames(res2),idtype[j],sep = "&")
  datalist2[[j]]=res2
}

sva.sur.os.pan.output=do.call(rbind, datalist2)
subset(sva.sur.os.pan.output, pvalue<0.05)
