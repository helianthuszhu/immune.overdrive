#
dim(traindata)
traindata[1:4,1:4]
class(traindata$kaps.group)
head(group_g.draw)
table(group_g.draw$group)
######
traindata.bo=traindata[,-1]
traindata.bo=as.data.frame(scale(traindata.bo))
traindata.bo=traindata.bo %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
rownames(traindata.bo)=rownames(traindata)
traindata.bo=cbind(traindata$kaps.group, traindata.bo)
colnames(traindata.bo)[1]="kaps.group"
traindata.bo[1:4,1:4]
traindata.bo$kaps.group=ifelse(traindata.bo$kaps.group=="set4",1,0)
traindata.bo$kaps.group=as.factor(traindata.bo$kaps.group)
#####KIRC
dim(glmnet.data.module.kirc)
head(clin.heat.kirc.var)
library(dplyr)
df= glmnet.data.module.kirc %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
rownames(df)=rownames(glmnet.data.module.kirc)
df=as.data.frame(scale(df))
class(df)
df[1:4,1:3]
testdata.kirc=cbind(clin.heat.kirc.var$kaps.group.kirc,df[rownames(clin.heat.kirc.var),])
colnames(testdata.kirc)[1]="kaps.group"
testdata.kirc$kaps.group=ifelse(testdata.kirc$kaps.group=="set4",1,0)
testdata.kirc$kaps.group=as.factor(testdata.kirc$kaps.group)
testdata.kirc[1:4,1:4]
table(testdata.kirc$kaps.group)

#
traindata.gy=traindata[, colnames(traindata) %in% rownames(group_g.draw[which(group_g.draw$group=="greenyellow"),])]
traindata.gy$kaps.group=traindata$kaps.group
traindata.gy$kaps.group=ifelse(traindata.gy$kaps.group=="set4",1,0)
traindata.gy$kaps.group=as.factor(traindata.gy$kaps.group)
class(traindata.gy$kaps.group)
dim(traindata.gy)
#
#
library(e1071)
model.svm = svm(formula = kaps.group ~ ., kernel = "linear", # 这里的待预测的变量Type必须是Factor
            data = traindata.gy)
summary(model.svm) 
train.pred = predict(model.svm, traindata.gy)
table(real=traindata.gy$kaps.group, predict=train.pred) 
confus.matrix = table(real=traindata.gy$kaps.group, predict=train.pred)
sum(diag(confus.matrix))/sum(confus.matrix)
#######
##train
model.svm = svm(formula = kaps.group ~ ., kernel = "linear", # 这里的待预测的变量Type必须是Factor
                data = traindata.bo)
summary(model.svm) 
train.pred = predict(model.svm, traindata.bo)
table(real=traindata.bo$kaps.group, predict=train.pred) 
confus.matrix = table(real=traindata.bo$kaps.group, predict=train.pred)
sum(diag(confus.matrix))/sum(confus.matrix)
##test kirc
test.pred = predict(model.svm,testdata.kirc)
confus.matrix = table(real=testdata.kirc$kaps.group, predict=test.pred)
table(real=testdata.kirc$kaps.group, predict=test.pred)
sum(diag(confus.matrix))/sum(confus.matrix)
#####test on GSE39582
#####
gse39582exp.agg.module[1:4,1:4]
dim(gse39582exp.agg.module)
test.gse=gse39582exp.agg.module %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
rownames(test.gse)=rownames(gse39582exp.agg.module)
dim(test.gse)
test.gse[1:4,1:4]
#####
test.pred.gse = predict(model.svm,test.gse)
test.pred.gse=as.data.frame(test.pred.gse)
head(test.pred.gse)
#####
head(clin.gse39582)
gse.svm.clin=cbind(test.pred.gse,clin.gse39582[rownames(test.pred.gse),])
gse.svm.clin$svm.group=ifelse(gse.svm.clin$test.pred.gse=="1","set4","set.other")
head(gse.svm.clin)
###
table(gse.svm.clin$test.pred.gse,gse.svm.clin$svm.group)
###
fit2<- survfit(Surv(dfs_days, dfs_status) ~ svm.group, data = gse.svm.clin)
ggsurvplot(fit2, data =gse.svm.clin,
               #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
               pval = TRUE, pval.method = TRUE, # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "TE.cluster.kaps",               # Change legend titles
               legend.labs = c("set.other","set4"),  # Change legend labels
               palette = c("#E69F00","#d73027"),  # Use JCO journal color palette,
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
##################################
##################################predict TE score
traindata.bo[1:4,1:4]
traindata.bo.score=traindata.bo
traindata.bo.score$kaps.group=kaps.td[rownames(traindata.bo.score),]$z.of.mean.exp
traindata.bo.score[1:3,1:4]
###
model.svm.score = svm(formula = kaps.group ~ ., # 这里的待预测的变量Type必须是Factor
                data = traindata.bo.score)
svr.pred = predict(model.svm.score,traindata.bo.score)
######
cor.test(svr.pred, traindata.bo.score$kaps.group)
#####
svr.pred.kirc = predict(model.svm.score,testdata.kirc[,-1])
cor.test(svr.pred.kirc, clin.heat.kirc.var$z.of.mean.exp)
#####
svr.pred.gse = predict(model.svm.score,test.gse)
gse.svm.clin.score=cbind(svr.pred.gse,clin.gse39582[names(svr.pred.gse),])
head(gse.svm.clin.score)
gse.svm.clin.score$group=cut(as.numeric(as.character(gse.svm.clin.score$svr.pred.gse)),
    breaks=quantile(as.numeric(as.character(gse.svm.clin.score$svr.pred.gse)), 
                    c(0, 0.25,0.5,0.75, 1), na.rm=T),labels=c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"), include.lowest = TRUE)
table(gse.svm.clin.score$group)
#
ss=Surv(gse.svm.clin.score$dfs_days, gse.svm.clin.score$dfs_status)
summary(coxph(ss~svr.pred.gse,data = gse.svm.clin.score))
#
fit2<- survfit(Surv(dfs_days, dfs_status) ~ group, data = gse.svm.clin.score)
ggsurvplot(fit2, data =gse.svm.clin.score,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TE.cluster.kaps",               # Change legend titles
           legend.labs = c("lmscore.set1","lmscore.set2","lmscore.set3","lmscore.set4"),  # Change legend labels
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
fit2<- survfit(Surv(dfs_days, dfs_status) ~ group, data = subset(gse.svm.clin.score, group=="lmscore.set4"|group=="lmscore.set3"))
ggsurvplot(fit2, data =subset(gse.svm.clin.score, group=="lmscore.set4"|group=="lmscore.set3"),
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TE.cluster.kaps",               # Change legend titles
           legend.labs = c("lmscore.set3","lmscore.set4"),  # Change legend labels
           palette = c("#E69F00","#d73027"),  # Use JCO journal color palette,
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
