##########compare piRNA expression among kaps group
#####
piRNAm=read.table("4.model/kaps/piwi/piwiRNA.expMatrix.txt",header = T,sep = "\t")
##

piRNAm=cbind(substr(piRNAm$TCGA.sample.ID,1,15), piRNAm)
colnames(piRNAm)[1]="id"
piRNAm=piRNAm[,-c(2:5)]
piRNAm[1:4,1:7]
########
piRNAm=piRNAm %>% group_by(id) %>% summarise_all(mean)
piRNAm=as.data.frame(piRNAm)
rownames(piRNAm)=piRNAm$id
piRNAm=piRNAm[,-1]
#
length(intersect(rownames(kaps.td), rownames(piRNAm)))
idpirna=intersect(rownames(kaps.td), rownames(piRNAm))
######
stat.piRNAm=cbind(kaps.td[idpirna,]$kaps.group, kaps.td[idpirna,]$z.of.mean.exp, piRNAm[idpirna,])
colnames(stat.piRNAm)[1:2]=c("kaps.group","z.of.mean.exp")
stat.piRNAm[1:3,1:4]
###filter on piRNA
#####
pirnacount=as.data.frame(colSums(stat.piRNAm[,-(1:2)] != 0))
colnames(pirnacount)="count"
head(pirnacount)
pirnacount.sel=subset(pirnacount, pirnacount[,1] >40)
dim(pirnacount.sel)
head(pirnacount.sel)
######
stat.piRNAm.sel=cbind(stat.piRNAm[,1:2], stat.piRNAm[, colnames(stat.piRNAm) %in% rownames(pirnacount.sel)])
dim(stat.piRNAm.sel)
######correlation of each piRNA with z.of mean exp
datalist=list()
for (i in 3:ncol(stat.piRNAm.sel)) {
  aa=cbind(stat.piRNAm.sel[,i], stat.piRNAm.sel[,2])
  colnames(aa)=c(colnames(stat.piRNAm.sel)[i],colnames(stat.piRNAm.sel)[2])
  #aa=na.omit(aa)
  aa=subset(aa, aa[,1] > 0)
  res <- cor.test(aa[,2], aa[,1],method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(stat.piRNAm.sel)[i]
  rownames(res)=paste(colnames(stat.piRNAm.sel)[2], rownames(res),sep = "&")
  datalist[[i]] <- res
}
TE.mean.corM.with.piRNA=do.call(rbind, datalist)
TE.mean.corM.with.piRNA=subset(TE.mean.corM.with.piRNA, pvalue< 0.05)
TE.mean.corM.with.piRNA=TE.mean.corM.with.piRNA[order(TE.mean.corM.with.piRNA$cor),]

TE.mean.corM.with.piRNA$piName=substr(rownames(TE.mean.corM.with.piRNA),15,nchar(rownames(TE.mean.corM.with.piRNA)))

head(TE.mean.corM.with.piRNA)
summary(TE.mean.corM.with.piRNA$cor)
save(TE.mean.corM.with.piRNA,file="4.model/kaps/piwi/piRNA/TE.mean.corM.with.piRNA.RData")
#######
stat.piRNAm.sel$kaps.group=factor(stat.piRNAm.sel$kaps.group,levels = c("set4", "set3","set2","set1"))
pirnaidx=TE.mean.corM.with.piRNA$piName[1:10]
for (i in 1:length(pirnaidx)) {
  
  aa=cbind(stat.piRNAm.sel[,pirnaidx[i]], stat.piRNAm.sel[,1:2])
  colnames(aa)=c(pirnaidx[i],colnames(stat.piRNAm.sel)[1:2])
  #aa=na.omit(aa)
  aa=subset(aa, aa[,1] > 0)
 
  pdf(paste0("4.model/kaps/piwi/piRNA/top.10.neg.corrlated/",pirnaidx[i],".pdf"), width = 7,height = 7)
  print(
    ggviolin(aa,
             x = "kaps.group", y = pirnaidx[i], fill = "kaps.group",alpha = 1,size = 0.3,
             #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
             palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
             add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
      stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
      theme_bw()+xlab("kaps.group")+ylab(paste0(pirnaidx[i]))+ggtitle(paste0(pirnaidx[i]))
  )
  print(ggscatter(aa, x = "z.of.mean.exp", y = pirnaidx[i], 
                  add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
                  cor.coef = T, cor.method = "spearman",
                  xlab = "z.mean.TE.score", ylab = paste0(pirnaidx[i]))+
                    ggtitle(paste0("z.of.mean.exp ", "vs ",pirnaidx[i]))+
                    geom_point(fill="#dd1c77",color="#dd1c77")
  )
  dev.off()
}
#######################
##############correlate with each 9 TEs
#################
stat.piRNAm.indiTE=cbind(kaps.td[rownames(stat.piRNAm.sel), c(3:11)],stat.piRNAm.sel)
stat.piRNAm.indiTE[1:4,1:15]


datalist1=list()
datalist2=list()
for (i in 1:9) {
  aa=cbind(stat.piRNAm.indiTE[,i], stat.piRNAm.indiTE[,-c(1:11)])
  colnames(aa)[1]=colnames(stat.piRNAm.indiTE)[i]
  aa[1:3,1:10]
  for (j in 2:ncol(aa)) {
    bb=aa[,c(1, j)]
    ccc=subset(bb, bb[,2] > 0)
    res <- cor.test(ccc[,2], ccc[,1],method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    #rownames(res)=colnames(stat.piRNAm.indiTE)[i]
    rownames(res)=paste(colnames(aa)[j], colnames(stat.piRNAm.indiTE)[i],sep = "&")
    datalist1[[j]] <- res
  }
  datalist2[[i]]=do.call(rbind, datalist1)
}
pirna.indicor.res=do.call(rbind, datalist2)
pirna.indicor.res$id=rownames(pirna.indicor.res)
pirna.indicor.res=separate(pirna.indicor.res, "id",into = c("pirnaName","teName"),sep = "&")
head(pirna.indicor.res)
summary(pirna.indicor.res$cor)
###########
idxx=unique(pirna.indicor.res$teName)
for (i in 1:length(idxx)) {
  pdf(paste0("4.model/kaps/piwi/piRNA/cor.indi.TE.with.piRNA/",idxx[i],".cor.distribution.pdf"))
  print(ggplot(subset(pirna.indicor.res, teName==idxx[i]), aes(x=cor))+
          geom_density(color="#92c5de", fill="#92c5de")+
          #geom_vline(xintercept = 0.3,color = "#b2182b", size=0.3)+
          geom_vline(xintercept = -0.2, color = "#b2182b", size=0.3)+
          #annotate(geom="text", x=0.2, y=1, label=paste0("No.piRNAs = ",nrow(subset(pirna.indicor.res, teName==idxx[i] & cor > 0.3) )),color="red")+
          annotate(geom="text", x=-0.2, y=1, label=paste0("No.piRNAs = ",nrow(subset(pirna.indicor.res, teName==idxx[i] & cor < -0.2))),color="red")
        )
  dev.off()
}
save(pirna.indicor.res, file = "4.model/kaps/piwi/piRNA/pirna.indicor.res.RData")
################
################draw hearmap to show
sig.indiTE.cor=subset(pirna.indicor.res,  cor <= -0.2)
#
datalist=list()
tarpirna=unique(sig.indiTE.cor$pirnaName)
for (i in 1:length(tarpirna)) {
  aa=as.data.frame(subset(pirna.indicor.res, pirnaName==tarpirna[i]))
  bb=as.data.frame(aa$cor)
  rownames(bb)=aa$teName
  colnames(bb)=tarpirna[i]
  datalist[[i]]=bb
}
sigTE.corMatrix=do.call(cbind, datalist)
dim(sigTE.corMatrix)
###
datalist=list()
for (i in 1:length(tarpirna)) {
  aa=as.data.frame(subset(pirna.indicor.res, pirnaName==tarpirna[i]))
  bb=as.data.frame(aa$pvalue)
  rownames(bb)=aa$teName
  colnames(bb)=tarpirna[i]
  datalist[[i]]=bb
}
sigTE.pvalueMatrix=do.call(cbind, datalist)
dim(sigTE.pvalueMatrix)
#

min_cor = min(as.vector(sigTE.corMatrix))
max_cor = max(as.vector(sigTE.corMatrix))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=100)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
pdf("4.model/kaps/piwi/piRNA/ht.cor.pdf",width = 6,height = 10)
Heatmap(t(sigTE.corMatrix),col = col.pal_cor,
        #rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        heatmap_legend_param = list(title = "Cor"),
        column_title = "Cor: expression TE vs piRNA (coef <= -0.2 with at least one TE)",
        column_title_gp = gpar(fontsize = 10),
        column_title_side = "top",
        show_row_dend = F,show_column_dend = F
        )
dev.off()

save(sigTE.corMatrix,sigTE.pvalueMatrix,file = "4.model/kaps/piwi/piRNA/ht.cor.RData" )

#####################check the PIWI RNA level 
#########
dim(genestats)
piwimRNA=genestats[, grep("PIWI",colnames(genestats))]
head(piwimRNA)
######
stat.piwimRNA=cbind(piwimRNA, kaps.td[rownames(piwimRNA),])
head(stat.piwimRNA)
##########

colnames(stat.piwimRNA)
datalist1=list()
datalist2=list()
for (i in 1:4) {
  aa=cbind(stat.piwimRNA[,i], stat.piwimRNA[,c(7:15,21)])
  colnames(aa)[1]=colnames(stat.piwimRNA)[i]
  
  for (j in 2:ncol(aa)) {
    bb=aa[,c(1, j)]
    ccc=na.omit(bb)
    res <- cor.test(ccc[,2], ccc[,1],method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    #rownames(res)=colnames(stat.piRNAm.indiTE)[i]
    rownames(res)=paste(colnames(aa)[j], colnames(stat.piwimRNA)[i],sep = "&")
    datalist1[[j]] <- res
  }
  datalist2[[i]]=do.call(rbind, datalist1)
}
piwimRNA.indicor.res=do.call(rbind, datalist2)
piwimRNA.indicor.res$id=rownames(piwimRNA.indicor.res)
piwimRNA.indicor.res=separate(piwimRNA.indicor.res, "id",into = c("piwimNRA","teName"),sep = "&")
head(piwimRNA.indicor.res)

save(piwimRNA.indicor.res, file="4.model/kaps/piwi/piRNA/piwimRNA.indicor.res.RData")













#
########################
######2 calculate Anova p value
datalist=list()
for (i in 3:ncol(stat.piRNAm)) {
  stat.piRNAm$kaps.group=as.factor(stat.piRNAm$kaps.group)
  res.aov <- aov(stat.piRNAm[,i] ~ kaps.group, data = stat.piRNAm)
  # Summary of the analysis
  pv=data.frame(summary(res.aov)[[1]][["Pr(>F)"]][1])
  colnames(pv)='anova.p.value'
  rownames(pv)=colnames(stat.piRNAm)[i]
  datalist[[i]] <- pv
}
anova.p.piRNA = do.call(rbind, datalist)
head(anova.p.piRNA)
anova.p.piRNA$piName=rownames(anova.p.piRNA)
anova.p.piRNA=anova.p.piRNA[order(anova.p.piRNA$anova.p.value),]
#summary(anova.p.piRNA$anova.p.value)
########
pirnaidx=anova.p.piRNA$piName[1:10]
for (i in 1:length(pirnaidx)) {
  pdf(paste0("4.model/kaps/piwi/piRNA/top.anova.piRNA/",pirnaidx[i],".pdf"), width = 5,height = 8)
  print(
    ggviolin(stat.piRNAm,
             x = "kaps.group", y = pirnaidx[i], fill = "kaps.group",alpha = 1,size = 0.3,
             #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
             palette = c(c("#d73027","#E69F00","#756bb1","#00AFBB")),
             add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
      stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
      theme_bw()+xlab("kaps.group")+ylab(paste0(pirnaidx[i]))+ggtitle(paste0(pirnaidx[i]))
  )
  dev.off()
}
##########################################







