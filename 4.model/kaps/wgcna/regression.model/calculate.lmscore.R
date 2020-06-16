#########
#########ROC test for each gene with kaps group
head(group_g)
dim(group_g)
#############
######prepare the taining data
head(group_g.draw)
dim(group_g.draw)
dim(datExpr)
glmnet.data=datExpr[, colnames(datExpr) %in% rownames(group_g.draw) ]
glmnet.data=as.data.frame(scale(glmnet.data))
glmnet.data[1:4,1:4]
######prepare GSE39582 test data
load("4.model/kaps/glmnet/GSE39582_eset.rda")
gse39582exp=as.data.frame(exprs(GSE39582_eset))
gse39582exp[1:3,1:4]
#####
gseid=as.data.frame(rownames(gse39582exp))
colnames(gseid)="id"
gseid=separate(gseid, col = "id",into = c("symbol","other"),sep = "[/]")
gseid$id.new=gseid$symbol
gseid$id.new=gsub("BMS1","LOC96610",gseid$id.new)
gseid$id.new=gsub("THEMIS2","C1orf38",gseid$id.new)
gseid$id.new=gsub("IGFLR1","TMEM149",gseid$id.new)
gseid$id.new=gsub("NXPE3","FAM55C",gseid$id.new)
head(gseid)
###
gse39582exp=cbind(gseid, gse39582exp)
gse39582exp=gse39582exp[,-c(1:2)]
gse39582exp.agg= gse39582exp %>% group_by(id.new) %>% summarise_all(mean)
gse39582exp.agg=as.data.frame(gse39582exp.agg)
rownames(gse39582exp.agg)=gse39582exp.agg$id.new
gse39582exp.agg=gse39582exp.agg[,-1]
gse39582exp.agg=as.data.frame(t(gse39582exp.agg))
dim(gse39582exp.agg)
###
modulegeneshid=intersect(colnames(glmnet.data), colnames(gse39582exp.agg))
######
glmnet.data.module=cbind(datTraits[rownames(glmnet.data),],glmnet.data[,modulegeneshid])
glmnet.data.module=glmnet.data.module[,-c(1)]
glmnet.data.module[1:4,1:4]
dim(glmnet.data.module)
#
gse39582exp.agg.module=gse39582exp.agg[,modulegeneshid]
gse39582exp.agg.module=as.data.frame(scale(gse39582exp.agg.module))
gse39582exp.agg.module[1:4,1:4]
dim(gse39582exp.agg.module)
########prepare kIRC expmatirx
#expM
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/8.KIRC/expMatrix.KIRC.RData")
kirc.expfireh[1:4,1:4]
load("8.KIRC/clin/heatmap.clin/clin.variables.heatmap.KIRC.RData")
head(clin.heat.kirc.var)
#
length(intersect(rownames(kirc.expfireh),modulegeneshid ))
glmnet.data.module.kirc=as.data.frame(t(kirc.expfireh[modulegeneshid,]))[rownames(clin.heat.kirc.var),]
glmnet.data.module.kirc=as.data.frame(scale(glmnet.data.module.kirc))
dim(glmnet.data.module.kirc)
#######calculate ROC value for each gene
########
library(pROC)
#roc.list.module <- roc(kaps.group ~ IGJ +ADAM6+STAT2, auc=TRUE, ci=TRUE,data = glmnet.data.module[,-1])
datalist=list()
#idxmodule=colnames(glmnet.data.module)[-c(1:2)]
for (i in 3:ncol(glmnet.data.module)) {
  aa=multiclass.roc(glmnet.data.module[,2] ~ glmnet.data.module[,i] )
  aucres=as.data.frame(aa$auc)
  colnames(aucres)="Multi.class.auc"
  rownames(aucres)=colnames(glmnet.data.module)[i]
  datalist[[i]]=aucres
}
roc.list.module=do.call(rbind, datalist)
#roc.list.module=multiclass.roc(glmnet.data.module[,2] ~ glmnet.data.module[,4])
#roc.list.module=multiclass.roc( kaps.group ~ IGJ,data=glmnet.data.module)
#as.data.frame(roc.list.module$auc)
head(roc.list.module)
dim(roc.list.module)
summary(roc.list.module$Multi.class.auc)
roc.list.module.sel=subset(roc.list.module,Multi.class.auc >= 0.7 )
##linear model
linearMod=lm( z.of.mean.exp~ JAK3+LAT+NLRC3+C5orf56+BTN3A1+STAT2, data=glmnet.data.module)
linearMod=as.data.frame(linearMod$coefficients)
###
logiMod=glm(kaps.group~ JAK3+LAT+NLRC3+C5orf56+BTN3A1+STAT2,family=binomial(link="logit"),data=glmnet.data.module)
logiMod=as.data.frame(logiMod$coefficients)
#
stat.roctest=glmnet.data.module[, colnames(glmnet.data.module) %in% rownames(roc.list.module.sel) ]
stat.roctest=cbind(stat.roctest, glmnet.data.module[, 1:2])
stat.roctest$lmscore=stat.roctest[,1]*linearMod[2,1]+stat.roctest[,2]*linearMod[3,1]+stat.roctest[,3]*linearMod[4,1]+
                     stat.roctest[,4]*linearMod[5,1]+stat.roctest[,5]*linearMod[6,1]+stat.roctest[,6]*linearMod[7,1]+linearMod[1,1]

stat.roctest$logiscore=stat.roctest[,1]*logiMod[2,1]+stat.roctest[,2]*logiMod[3,1]+stat.roctest[,3]*logiMod[4,1]+
                     stat.roctest[,4]*logiMod[5,1]+stat.roctest[,5]*logiMod[6,1]+stat.roctest[,6]*logiMod[7,1]+logiMod[1,1]



multiclass.roc(stat.roctest$kaps.group ~ stat.roctest$logiscore )
multiclass.roc(stat.roctest$kaps.group ~ stat.roctest$lmscore )
head(stat.roctest)
#####
ggscatter(stat.roctest, x = "z.of.mean.exp", y = "lmscore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.TE.score", ylab = "lmscore")+
  ggtitle("mean.TE.exp vs lmscore in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")


ggscatter(stat.roctest, x = "z.of.mean.exp", y = "logiscore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.TE.score", ylab = "logiscore")+
  ggtitle("mean.TE.exp vs logiscore in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(stat.roctest, x = "lmscore", y = "logiscore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "lmscore", ylab = "logiscore")+
  ggtitle("lmscore vs logiscore in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
####################
####################validation on kirc data

##
stat.roctest.kirc=as.data.frame(t(kirc.expfireh[rownames(roc.list.module.sel),]))
stat.roctest.kirc=as.data.frame(scale(stat.roctest.kirc))
length(intersect(rownames(stat.roctest.kirc), rownames(clin.heat.kirc.var)))
stat.roctest.kirc=cbind(stat.roctest.kirc[rownames(clin.heat.kirc.var),],clin.heat.kirc.var[,c(46,49)])
head(stat.roctest.kirc)
stat.roctest.kirc$lmscore=stat.roctest.kirc[,1]*linearMod[2,1]+stat.roctest.kirc[,2]*linearMod[3,1]+stat.roctest.kirc[,3]*linearMod[4,1]+
  stat.roctest.kirc[,4]*linearMod[5,1]+stat.roctest.kirc[,5]*linearMod[6,1]+stat.roctest.kirc[,6]*linearMod[7,1]+linearMod[1,1]

stat.roctest.kirc$logiscore=stat.roctest.kirc[,1]*logiMod[2,1]+stat.roctest.kirc[,2]*logiMod[3,1]+stat.roctest.kirc[,3]*logiMod[4,1]+
  stat.roctest.kirc[,4]*logiMod[5,1]+stat.roctest.kirc[,5]*logiMod[6,1]+stat.roctest.kirc[,6]*logiMod[7,1]+logiMod[1,1]
#
ggscatter(stat.roctest.kirc, x = "z.of.mean.exp", y = "lmscore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.TE.score", ylab = "lmscore")+
  ggtitle("mean.TE.exp vs lmscore in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")


ggscatter(stat.roctest.kirc, x = "z.of.mean.exp", y = "logiscore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.TE.score", ylab = "logiscore")+
  ggtitle("mean.TE.exp vs logiscore in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(stat.roctest.kirc, x = "lmscore", y = "logiscore", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "lmscore", ylab = "logiscore")+
  ggtitle("lmscore vs logiscore in KIRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
###########################################
####use kapa on lmscore
head(stat.roctest)
length(intersect(rownames(stat.roctest), rownames(kaps.td)))
#
stat.roctest.clin=cbind(stat.roctest[,-c(8,9)], kaps.td[rownames(stat.roctest),])
save(stat.roctest.clin,file="4.model/kaps/wgcna/regression.model/stat.roctest.clin.CRC.lmscore.for.kaps.test.RData")
#####################
