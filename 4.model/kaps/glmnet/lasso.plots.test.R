#############
dim(expM)
dim(datExpr0)
#############
dim(traindata)
x <- as.matrix(traindata[,-1])
#y <- traindata$kaps.group
#x=as.matrix(datExpr0[rownames(kaps.td),])
#x=expM[subset(TE.mean.corM.with.genes, !(set=="set3"))$symbol,colnames(expM) %in% rownames(kaps.td)]
#xx=expM[subset(TE.mean.corM.with.genes, !(set=="set3"))$symbol,colnames(expM) %in% rownames(kaps.td)]
#x=as.data.frame(t(x))
y <- kaps.td[rownames(x),]$kaps.group
y=ifelse(y=="set4",1,2)
dim(x)
#library(dplyr)
#x=x %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))

#
#x=scale(x)
#dim(x)

fit = glmnet(x, y, family = "multinomial", type.multinomial = "grouped")
plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")
cvfit=cv.glmnet(x, y, family="multinomial", type.multinomial = "grouped", parallel = TRUE)
plot(cvfit)
coef(cvfit, cvfit$lambda.min)
coef(cvfit, cvfit$lambda.1se)

length(coef(cvfit, cvfit$lambda.1se)$set3@i)

##############
#############
cv.lasso=cvfit
pdf(file = "4.model/kaps/glmnet/lasso.plots.pdf",width = 14)
# Plot Cross-Validation LASSO model
par(mfrow=c(1,2))
plot(cv.lasso,las=1, main="Lasso fit CV")
abline(v=log(cv.lasso$lambda.min), col="#980043", lwd=2, lty=3)
#dev.off()
# Plot lambda fit
#pdf(file = "4.model/kaps/glmnet/lasso.plots2.pdf",width = 10)
plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")
abline(v=log(cv.lasso$lambda.min), col="#980043", lwd=2,lty=3)
text(log(cv.lasso$lambda.min), 1.4, paste("lambda.min=",round(cv.lasso$lambda.min,4),"\n",
                                          length(coef(cv.lasso, cv.lasso$lambda.min)$set4@i)-1, " TEs" ,sep=""), 
     col="#980043", cex=0.75, pos=4)
dev.off()
################
########final model
x.test <- as.matrix(traindata[,-1])
#x.test=as.matrix(datExpr0[rownames(kaps.td),])
#x.test=x
prob.train=predict(cvfit, newx = x.test, s = "lambda.min", type = "class")
prob.train=as.data.frame(prob.train)
colnames(prob.train)="kaps.group.predict"
prob.train.roc=cbind(prob.train, glmnet.data.module[,1:2])
#prob.train.roc=cbind(prob.train, kaps.td[colnames(xx),c(17,38)])
#prob.train.roc$kaps.group.predict.numeric=as.numeric(paste(substr(prob.train.roc$kaps.group.predict,4,4)))
prob.train.roc$kaps.group.bi=ifelse(prob.train.roc$kaps.group=="set4","set4","set.other")
prob.train.roc$kaps.group.bi.numeric=ifelse(prob.train.roc$kaps.group=="set4",1,2)
prob.train.roc$kaps.group.predict=as.numeric(paste(prob.train.roc$kaps.group.predict))
head(prob.train.roc)
#####
#multiclass.roc(prob.train.roc$kaps.group~prob.train.roc$kaps.group.predict.numeric )
roc(kaps.group.bi.numeric ~ kaps.group.predict, auc=TRUE, ci=TRUE,data = prob.train.roc)
#
table(prob.train.roc$kaps.group.predict,prob.train.roc$kaps.group.bi)

###########################
###############KIRC prediction
dim(glmnet.data.module.kirc)
#x.test.kirc <- as.matrix(glmnet.data.module.kirc)
x.test.kirc=as.data.frame(t(kirc.expfireh[colnames(x),]))[rownames(clin.heat.kirc.var),]


x.test.kirc=x.test.kirc %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
rownames(x.test.kirc)=rownames(clin.heat.kirc.var)
x.test.kirc=scale(x.test.kirc)
class(x.test.kirc)
x.test.kirc[1:4,1:4]

prob.kirc=predict(cvfit, newx = x.test.kirc, s = "lambda.min", type = "class")
prob.kirc=as.data.frame(prob.kirc)
colnames(prob.kirc)="kaps.group.kirc.predict"
#rownames(prob.kirc)=rownames(glmnet.data.module.kirc)
rownames(prob.kirc)=rownames(x.test.kirc)
length(intersect(rownames(prob.kirc), rownames(clin.heat.kirc.var)))
dim(prob.kirc)
prob.kirc.roc=cbind(prob.kirc, clin.heat.kirc.var[rownames(prob.kirc),c(46,49)])
prob.kirc.roc$kaps.group.kirc..predict.numeric=as.numeric(paste(substr(prob.kirc.roc$kaps.group.kirc.predict,4,4)))
prob.kirc.roc=na.omit(prob.kirc.roc)
head(prob.kirc.roc)
dim(prob.kirc.roc)
#####
multiclass.roc(prob.kirc.roc$kaps.group.kirc~prob.kirc.roc$kaps.group.kirc..predict.numeric)
#
table(prob.kirc.roc$kaps.group.kirc,prob.kirc.roc$kaps.group.kirc.predict)
##################
#################
#################validation on CRC GSE39582
##########
dim(gse39582exp.agg.module)
gse39582exp.agg.module[1:3,1:4]
x.test.gse <- as.matrix(gse39582exp.agg.module)
prob.gse=predict(cvfit, newx = x.test.gse, s = "lambda.min", type = "class")
prob.gse=as.data.frame(prob.gse)
colnames(prob.gse)="kaps.group.gse.predict"
rownames(prob.gse)=rownames(gse39582exp.agg.module)
#####clin of GSE39582
clin.gse39582=read.csv("4.model/kaps/glmnet/GSE39582.csv",header = T,row.names = 1)
length(intersect(rownames(gse39582exp.agg.module), rownames(clin.gse39582)))
#
prob.gse.roc=cbind(prob.gse, clin.gse39582[rownames(prob.gse),])
#prob.gse.roc$kaps.group.gse.predict=as.numeric(paste(substr(prob.gse.roc$kaps.group.gse.predict,4,4)))
#prob.gse.roc=na.omit(prob.gse.roc)
head(prob.gse.roc)
dim(prob.gse.roc)
table(prob.gse.roc$kaps.group.gse.predict)
#####
#
table(prob.gse.roc$kaps.group.gse.predict,prob.gse.roc$msi)
#####
fit2<- survfit(Surv(dfs_days, dfs_status) ~ kaps.group.gse.predict, data = prob.gse.roc )
ggsurvplot(fit2, data =prob.gse.roc,
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
#
fit2<- survfit(Surv(dfs_days, dfs_status) ~ kaps.group.gse.predict, data = subset(prob.gse.roc, kaps.group.gse.predict=="set4"|kaps.group.gse.predict=="set3"))
ggsurvplot(fit2, data =subset(prob.gse.roc, kaps.group.gse.predict=="set4"|kaps.group.gse.predict=="set3"),
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TE.cluster.kaps",               # Change legend titles
           legend.labs = c("set3","set4"),  # Change legend labels
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
