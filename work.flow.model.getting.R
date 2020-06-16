#load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/tmp.sur.immu.overlap.RData")


#######correlate immune pathways
################
#load pathway scores
#####
load("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/genesets17s.score.RData")
######immuscore
immuscore=read.table("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/estimate/CRC_estimate_score.illumina.txt",header=T,row.names = 1)
rownames(immuscore)=gsub("[.]","-",rownames(immuscore))
#load CRC gene expression matrix
######
load("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/expMatrixCRC.RData")
####

immg=c("LAG3","TNFRSF14","BTLA","CD86","CD80","CTLA4","PDCD1LG2","CD274","PDCD1","CD8A","HAVCR2")
#####
immgM=as.data.frame(t(expM[immg,]))
######
score.data=cbind(immuscore$ImmuneScore, genesets17s.score[rownames(immuscore),], immgM[rownames(immuscore),] )
colnames(score.data)[1]="immuscore"
dim(score.data)
colnames(score.data)
###get the TE expression
load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.plot/tmp.RData")
load('~/nas/Xiaoqiang/opti.data//1.CGC.TE.CRC.REdiscoverTE/stats.20200311/1.TE.used/TE.label.used.RData')
###inclue TE used
#head(cdrep)
#table(cdrep$repClass.new2)
cdrep.sel=subset(cdrep, !(repClass.new2=="other.repeats"))
dim(cdrep.sel)
###TE matrix selection
te.tumorexp.sel=te.tumorexp[, colnames(te.tumorexp) %in% rownames(cdrep.sel)]
te.tumorexp.sel[1:4,1:4]
dim(te.tumorexp.sel)
######
ids=intersect(rownames(te.tumorexp.sel), rownames(score.data))
#score.data[which(rownames(score.data)=="TCGA-CM-4748-01"),1:4]

#####get the R and pvalue of spearman correlation
cor.score=score.data[ids,]
cor.exp=te.tumorexp.sel[ids,]

cor.score=cor.score[-which(rownames(cor.score)=="TCGA-CM-4748-01"),]
cor.exp=cor.exp[-which(rownames(cor.exp)=="TCGA-CM-4748-01"),]
cor.score=cor.score[rownames(cor.exp),]
#####
dim(cor.score)
dim(cor.exp)
cor.score[1:4,1:4]
cor.exp[1:4,1:4]
#
datalist=list()
oklist=list()
for (i in 1:ncol(cor.exp)) {
  for (j in 1:ncol(cor.score)) {
    res <- cor.test(cor.exp[,i], cor.score[,j],  method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    rownames(res)=colnames(cor.score)[j]
    rownames(res)=paste(colnames(cor.exp)[i], rownames(res),sep = "&")
    datalist[[j]] <- res
  }
  oklist[[i]]=do.call(rbind, datalist)
}
res.cor=do.call(rbind, oklist)
res.cor.29s=res.cor
save(res.cor.29s, 
     file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/3.immu.correlated/res.cor.29s.RData")
summary(res.cor$cor)
####################################

#use tumor purity as covariance
purdata=read.csv("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.plot/immune/tcga.tumor.purity.csv",header = T)
purdata=subset(purdata, Cancer.type=="COAD"|Cancer.type=="READ")
purdata$id=substr(purdata$Sample.ID,1,15)
purdata$idd=substr(purdata$Sample.ID,16,16)
purdata=subset(purdata,idd=="A")
rownames(purdata)=purdata$id
length(intersect(rownames(purdata),rownames(cor.score)))
cor.score.new=cbind(purdata[rownames(cor.score),]$IHC,cor.score)
colnames(cor.score.new)[1]="IHC"
#summary(cor.score.new$IHC)
cor.score.new=subset(cor.score.new, IHC>0)
dim(cor.score.new)
cor.exp.new=cor.exp[rownames(cor.score.new),]
cor.score.new[1:4,1:4]
cor.exp.new[1:4,1:4]
####adjusted by purity
datalist=list()
oklist=list()
library(ppcor)
for (i in 2:ncol(cor.score.new)) {
  subdata=cbind(cor.score.new[,c(1,i)],cor.exp.new)
  subdata=na.omit(subdata,cols=colnames(subdata)[2])
  #
  for (j in 3:ncol(subdata)) {
    #res <- pcor.test(cor.exp.new[,i], cor.score.new[,j], cor.score.new[,1],method = "spearman")
    res <- pcor.test(subdata[,2], subdata[,j], subdata[,1],method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    rownames(res)=colnames(subdata)[j]
    rownames(res)=paste(colnames(subdata)[2], rownames(res),sep = "&")
    datalist[[j]] <- res
  }
  oklist[[i]]=do.call(rbind, datalist)
}

res.Pcor.29s=do.call(rbind, oklist)
save(res.Pcor.29s, 
     file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/3.immu.correlated/res.pcor.29s.RData")
#################
#res.Pcor.29s$id=rownames(res.Pcor.29s)
library(tidyr)
#res.Pcor.29s=tidyr::separate(res.Pcor.29s, col = id, into = c("geneset","repName"),sep = "&")
#res.Pcor.29s.pos.sig=subset(res.Pcor.29s, cor >=0.3)
#head(res.Pcor.29s.pos.sig)
#res.Pcor.29s.pos.sig=subset(res.Pcor.29s, abs(cor)>=0.2)
##
###
#####
####select out sig tes
head(res.cor.29s)
res.cor.29s$id=rownames(res.cor.29s)
library(tidyr)
res.cor.29s=tidyr::separate(res.cor.29s, col = id, into = c("repName","geneset"),sep = "&")
res.Pcor.29s.pos.sig=subset(res.cor.29s, cor >=0.4)  ####good
#res.Pcor.29s.pos.sig=subset(res.cor.29s, cor >=0.3)
#head(res.Pcor.29s.pos.sig)
########
#####TE counts
te.counts=as.data.frame(table(res.Pcor.29s.pos.sig$repName))

colnames(te.counts)=c("repName","count")

rownames(te.counts)=te.counts$repName

te.counts=te.counts[order(te.counts$count,decreasing = T),]

te.counts=cbind(te.counts, cdrep[rownames(te.counts),])
head(te.counts)
dim(te.counts)
########bar plot of TE counts
library(ggplot2)
library(ggpubr)
cbPalette= c("#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e")
p1=ggbarplot(te.counts[,-1], x = "repName", y = "count",
             fill = "repClass.new2",               # change fill color by cyl
             color = "white",            # Set bar border colors to white
             #color = "#525252", 
             palette = cbPalette,            # jco journal color palett. see ?ggpar
             sort.val = "asc",          # Sort the value in dscending order
             sort.by.groups = TRUE,     # Don't sort inside each group
             x.text.angle = 90,           # Rotate vertically x axis texts
             ggtheme = theme_pubclean()
)+
  #font("x.text", size = 5, vjust = 0.5)+
  theme(axis.text.x = element_text(size=8))+
  ylab(paste0("Potentially immunogenic","\n", "with number of immue sets"))+ggtitle("TEs at least postively (cor>=0.4)")
generate.PDF <- function(pbmc) {
  #pdf("bar.top.200.of.781.tes.counts.geneset.pdf",width = 10,height = 5)
  pdf("3.immu.correlated/bar.total.14.tes.counts.geneset.2.pdf",width = 6,height = 5)
  print(p1)
  dev.off()
}
generate.PDF(pbmc)

#pie chart of TE 
library(ggplot2)
library(scales)
dfpie.pos <- as.data.frame(table(te.counts$repClass.new2))
colnames(dfpie.pos)[1]="TE.class"
head(dfpie.pos)
rownames(dfpie.pos)=dfpie.pos$TE.class
#############################################
colors <- c("#377eb8","#4daf4a","#e41a1c",'#27408B', '#FF0000', '#2E8B57', '#CD00CD', '#EE7942','#009ACD')
colors=c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")
slices <- dfpie.pos$Freq
lbls <- dfpie.pos$TE.class
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
library(ggplotify)
pdf("3.immu.correlated/pie.chart.613.tes.positive.2.pdf",width = 7,height = 7)
pie(slices,labels = lbls, col=colors,radius = 0.5,
    main=paste0("Proportion of 14", "\n","positive TEs at class level"))
dev.off()
###########draw individual cor of sig TE with immune sets
#########correlation
bbb=res.cor.29s
dim(te.counts)
head(te.counts)
dim(cor.exp)
cor.exp.draw=cor.exp[, colnames(cor.exp) %in% rownames(te.counts)]
library("ggpubr")
for (i in 1:ncol(cor.exp.draw)) {
  for (j in 1:ncol(cor.score)) {
    cd=data.frame(exp=cor.exp.draw[,i],score=cor.score[,j])
    colnames(cd)=c("exp","score")
    idxcor=paste0(colnames(cor.exp.draw)[i],"&",colnames(cor.score)[j])
    pvalue=bbb[rownames(bbb)==idxcor,]$pvalue
    coef=bbb[rownames(bbb)==idxcor,]$cor
    pdf(paste0("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/3.immu.correlated/indi.cors/",
               colnames(cor.exp.draw)[i],"-vs-",colnames(cor.score)[j],".pdf"))
    grob <- grobTree(textGrob(paste0("pcor =",coef, "\n", "pvalue =",pvalue), x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=13, fontface="italic")))
    print(ggscatter(cd, x = "score", y = "exp", 
                    add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
                    cor.coef = T, cor.method = "spearman",
                    xlab = colnames(cor.score)[j], ylab = colnames(cor.exp.draw)[i])+
            ggtitle(paste0(colnames(cor.exp.draw)[i],"-vs-",colnames(cor.score)[j]))+#geom_point(aes(color="#dd1c77"))
            geom_point(fill="#dd1c77",color="#dd1c77")+
            annotation_custom(grob)
          #annotate(geom="text", x=quantile(cd$score)[1], y=quantile(cd$exp)[4], label=paste0("pcor =",coef, "\n", "pvalue =",pvalue),color="red")
    )
    dev.off()
  }
}
#####load survival stats
#######
#load("2.survival.TE.screen/sur.res.combined.endpoints.RData")
load("2.survival.TE.screen/sur.res.combined.endpoints.2.RData")
write.csv(cb.data, "2.survival.TE.screen/cb.data.csv")
head(cb.data)
a1=subset(cb.data, endpoint=="dfi")
a2=subset(cb.data, endpoint=="dss")
a3=subset(cb.data, endpoint=="os")
a4=subset(cb.data, endpoint=="pfi")
table(cb.data$endpoint)
#b1=intersect(a3$repName, unique(subset(cb.data, !(endpoint=="os"))$repName))
#b2=intersect(a1$repName, a4$repName)
#b3=c(b1,b2)
#length(b3)
#b4=Reduce(intersect, list(a1$repName,a2$repName,a3$repName,a4$repName))
#length(b4)
#########
head(cb.data)
sur.total.id=unique(cb.data$repName)
#sur.total.id=unique(cb.data[which(cb.data$HR>=1.5),]$repName)
#sur.total.id=unique(cb.data[-which(cb.data$endpoint=="dfi"),]$repName)
length(intersect(rownames(te.counts), sur.total.id))
cndi.id=intersect(rownames(te.counts), sur.total.id)
#cndi.id=intersect(rownames(te.counts[which(te.counts$count>7),-1]), sur.total.id)
######
#sur.total.id=unique(b4)
#sur.total.id=unique(b3)
#length(intersect(rownames(te.counts), sur.total.id))
#cndi.id=intersect(rownames(te.counts), sur.total.id)
#cndi.id=intersect(rownames(te.counts[which(te.counts$count>4),-1]), sur.total.id)
cndi.rep=cdrep[rownames(cdrep) %in% cndi.id, ]
table(cndi.rep$repClass.new2)
#head(cndi.rep)
#rownames(cndi.rep)
dim(cndi.rep)
head(cndi.rep)
##############
head(te.counts)
dim(te.counts)
head(res.Pcor.29s.pos.sig)
class(res.Pcor.29s.pos.sig)
#
cds=res.Pcor.29s.pos.sig[res.Pcor.29s.pos.sig$repName %in% cndi.rep$repName,]
cds=cds[order(cds$cor,decreasing = T),]
head(cds)
summary(cds$cor)
#
#heatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
library(RColorBrewer)
library(pheatmap)
dim(cor.exp)
#da1=te.tumorexp.sel[, colnames(te.tumorexp.sel) %in% cndi.rep$repName]
da1=cor.exp[, colnames(cor.exp) %in% cndi.rep$repName]
dat <- as.data.frame(t(da1))
dim(dat)
dim(da1)
######get the clinical infor
ann.col=data.frame(msi=clin.CRC.full.sel$msi,
                   os_status=clin.CRC.full.sel$os_status,
                   dfs_status=clin.CRC.full.sel$dfs_status)
rownames(ann.col)=rownames(clin.CRC.full.sel)
ann.col=ann.col[rownames(da1),]

head(ann.col)

pp1=pheatmap(dat,show_rownames = T,border_color = NA,
             show_colnames = F,
             #cluster_rows = hc,
             #cluster_cols =hc,
             annotation_col = ann.col,
             scale="row",
             color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256), 
             annotation_legend = T
             #annotation_colors  = list(response = c("PD" = "#dd1c77","PRCR"="#74c476","SD"="#feb24c"),
             #                         treatment.status=c("On"="#74a9cf", "Pre"="#f768a1"))
)
######
clu=as.data.frame(sort(cutree(pp1$tree_col, k=3)))
colnames(clu)="pt.clu3"
clu$pt.clu3=paste("cluster",clu$pt.clu3,sep = "_")
clu$idss=rownames(clu)
table(clu$pt.clu3)
head(clu)
length(intersect(rownames(clu), rownames(ann.col)))
dim(clu)
dim(ann.col)
clu=cbind(clu, ann.col[rownames(clu),])
table(clu$pt.clu3, clu$msi)
#clu[which(clu$pt.clu3=="cluster_3"),]
head(clu)
#clu$pt.clu3.new=gsub("cluster_3","cluster_2",clu$pt.clu3)
#clu$pt.clu3.new=gsub("cluster_4","cluster_2",clu$pt.clu3.new)
#clu$pt.clu3.new=gsub("cluster_6","cluster_2",clu$pt.clu3.new)
#clu$pt.clu3.new=gsub("cluster_5","cluster_2",clu$pt.clu3.new)
#clu$pt.clu3.new=gsub("cluster_4","cluster_2",clu$pt.clu3.new)
#clu$pt.clu3.new=gsub("cluster_1","cluster_3",clu$pt.clu3)
######
clus=clu[,c(1,3)]
clus$msi=gsub("Indeterminate","MSS",clus$msi)
colnames(clus)=c("TE.cluster","MSI.status")
clus$TE.cluster.agg=ifelse(clus$TE.cluster=="cluster_1","TE.high","TE.low")
head(clus)
head(cndi.rep)
rownames(dat)
pp1=pheatmap(dat,show_rownames = T,border_color = NA,
             show_colnames = F,fontsize = 7,
             #cluster_rows = hc,
             #cluster_cols =hc,
             annotation_col = clus[,-2],
             annotation_row = cndi.rep[,c(2,3)],
             scale="row",
             #color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8), 
             color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             annotation_legend = T,
             annotation_colors  = list( TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                        TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                   #,MSI.status=c("MSI-H"="#a65628", "MSI-L"="#ffd92f","MSS"="#8da0cb"),
                                   repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                   repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
                                   )
)

save_pheatmap_pdf(pp1,filename = "4.model/heatmpa.9.TEs.full.pdf",height = 4)
########test the immune difference among clusters
dim(cor.exp)
head(cds)
drawbox=cbind(clu, cor.exp[rownames(clu),],cor.score[rownames(clu),])

ggplot(drawbox, aes_string(x="pt.clu3", y="LTR21B", group="pt.clu3")) + 
  geom_boxplot(aes(fill = pt.clu3), alpha = 0.5,outlier.colour = "black",outlier.size = 1,notch = F)+
  stat_compare_means(label = "p.signif")+
  geom_jitter(aes(color = pt.clu3),width = 0.15,size=0.3)#+
  #facet_grid(. ~ endpoint)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  #scale_color_manual(values = cbPalette)+ylab("Hazard ratio")
ggplot(drawbox, aes_string(x="pt.clu3", y="LAG3", group="pt.clu3")) + 
  geom_boxplot(aes(fill = pt.clu3), alpha = 0.5,outlier.colour = "black",outlier.size = 1,notch = F)+
  stat_compare_means(label = "p.signif")+
  geom_jitter(aes(color = pt.clu3),width = 0.15,size=0.3)

#cor.test(drawbox$LTR21B, drawbox$immuscore,method = "spearman")
#head(clu)
idclu=intersect(rownames(clu), rownames(CRC.survivaldata))
#####
surstats=cbind(clu[idclu,], CRC.survivaldata[idclu,])
surstats$TE.culster.agg=ifelse(surstats$pt.clu3=="cluster_1","TE.high","TE.low")
head(surstats)
dim(surstats)
table(surstats$pt.clu3,surstats$TE.culster.agg)
#head(surstats)
require("survival")
require("survminer")
######only cluster 1 vs cluster 2
surstats.sel=subset(surstats, pt.clu3=="cluster_1" | pt.clu3=="cluster_2")
########DSS
ss=Surv(surstats.sel$DSS.time, surstats.sel$DSS)
cox = summary(coxph(ss~factor(surstats.sel$pt.clu3,levels = c("cluster_2","cluster_1"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(DSS.time, DSS) ~ pt.clu3, data = surstats.sel)
dss.two=ggsurvplot(fit1, data =surstats.sel,
                title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                surv.median.line = "hv",            # Add median survival lines
                legend.title = "TE.cluster",               # Change legend titles
                legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                palette = c("#ca0020","#0571b0","#4daf4a"),
                risk.table = TRUE,                  # Add No at risk table
                cumevents = TRUE,                   # Add cumulative No of events table
                tables.height = 0.15,               # Specify tables height
                tables.theme = theme_cleantable(),  # Clean theme for tables
                tables.y.text = FALSE,             # Hide tables y axis text
                ylab="Disease specific survival probability",
                #ylab="Overall survival probability",
                #ylab="progression-free interval probability",
                #ylab="disease-free interval probability",
                xlab="Time(days)"
)
#########
fit2<- survfit(Surv(DSS.time, DSS) ~ pt.clu3, data = surstats )
dss.three=ggsurvplot(fit2, data =surstats,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TE.cluster",               # Change legend titles
           legend.labs = c("cluster_1","cluster_2","cluster_3"),  # Change legend labels
           #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
           palette = c("#ca0020","#0571b0","#4daf4a"),
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE,             # Hide tables y axis text
           ylab="Disease specific survival probability",
           #ylab="Overall survival probability",
           #ylab="progression-free interval probability",
           #ylab="disease-free interval probability",
           xlab="Time(days)"
)
#
ss=Surv(surstats.sel$DSS.time, surstats.sel$DSS)
cox = summary(coxph(ss~factor(surstats$TE.culster.agg,levels = c("TE.high","TE.low"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
fit2<- survfit(Surv(DSS.time, DSS) ~ TE.culster.agg, data = surstats )
dss.agg=ggsurvplot(fit2, data =surstats,
                     title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                     pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                     surv.median.line = "hv",            # Add median survival lines
                     legend.title = "TE.cluster.agg",               # Change legend titles
                     legend.labs = c("TE.high","TE.low"),  # Change legend labels
                     #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                     palette =c("#d73027","#E69F00","#00AFBB"),
                     #palette = c("#ca0020","#0571b0","#4daf4a"),
                     risk.table = TRUE,                  # Add No at risk table
                     cumevents = TRUE,                   # Add cumulative No of events table
                     tables.height = 0.15,               # Specify tables height
                     tables.theme = theme_cleantable(),  # Clean theme for tables
                     tables.y.text = FALSE,             # Hide tables y axis text
                     ylab="Disease specific survival probability",
                     #ylab="Overall survival probability",
                     #ylab="progression-free interval probability",
                     #ylab="disease-free interval probability",
                     xlab="Time(days)"
)
generate.PDF <- function(fig) {
  pdf("4.model/surC.dss.agg.pdf",width = 5,height = 8)
  print(dss.two)
  print(dss.three)
  print(dss.agg)
  dev.off()
}
generate.PDF(fig)
###OS
ss=Surv(surstats.sel$OS.time, surstats.sel$OS)
cox = summary(coxph(ss~factor(surstats.sel$pt.clu3,levels = c("cluster_2","cluster_1"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(OS.time, OS) ~ pt.clu3, data = surstats.sel)
os.two=ggsurvplot(fit1, data =surstats.sel,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster",               # Change legend titles
                   legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                   palette = c("#ca0020","#0571b0","#4daf4a"),
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
#########
fit2<- survfit(Surv(OS.time, OS) ~ pt.clu3, data = surstats )
os.three=ggsurvplot(fit2, data =surstats,
                     #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                     pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                     surv.median.line = "hv",            # Add median survival lines
                     legend.title = "TE.cluster",               # Change legend titles
                     legend.labs = c("cluster_1","cluster_2","cluster_3"),  # Change legend labels
                     #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                     palette = c("#ca0020","#0571b0","#4daf4a"),
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
ss=Surv(surstats.sel$OS.time, surstats.sel$OS)
cox = summary(coxph(ss~factor(surstats$TE.culster.agg,levels = c("TE.high","TE.low"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(OS.time, OS) ~ TE.culster.agg, data = surstats )
os.agg=ggsurvplot(fit2, data =surstats,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster.agg",               # Change legend titles
                   legend.labs = c("TE.high","TE.low"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                   palette =c("#d73027","#E69F00","#00AFBB"),
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
generate.PDF <- function(fig) {
  pdf("4.model/surC.os.agg.pdf",width = 5,height = 8)
  print(os.two)
  print(os.three)
  print(os.agg)
  dev.off()
}
generate.PDF(fig)
###DFI
ss=Surv(surstats.sel$DFI.time, surstats.sel$DFI)
cox = summary(coxph(ss~factor(surstats.sel$pt.clu3,levels = c("cluster_2","cluster_1"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(DFI.time, DFI) ~ pt.clu3, data = surstats.sel)
dfi.two=ggsurvplot(fit1, data =surstats.sel,
                  title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                  pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                  surv.median.line = "hv",            # Add median survival lines
                  legend.title = "TE.cluster",               # Change legend titles
                  legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                  #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                  palette = c("#ca0020","#0571b0","#4daf4a"),
                  risk.table = TRUE,                  # Add No at risk table
                  cumevents = TRUE,                   # Add cumulative No of events table
                  tables.height = 0.15,               # Specify tables height
                  tables.theme = theme_cleantable(),  # Clean theme for tables
                  tables.y.text = FALSE,             # Hide tables y axis text
                  #ylab="Disease specific survival probability",
                  #ylab="Overall survival probability",
                  #ylab="progression-free interval probability",
                  ylab="disease-free interval probability",
                  xlab="Time(days)"
)
#########
fit2<- survfit(Surv(DFI.time, DFI) ~ pt.clu3, data = surstats )
dfi.three=ggsurvplot(fit2, data =surstats,
                    #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                    pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                    surv.median.line = "hv",            # Add median survival lines
                    legend.title = "TE.cluster",               # Change legend titles
                    legend.labs = c("cluster_1","cluster_2","cluster_3"),  # Change legend labels
                    #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                    palette = c("#ca0020","#0571b0","#4daf4a"),
                    risk.table = TRUE,                  # Add No at risk table
                    cumevents = TRUE,                   # Add cumulative No of events table
                    tables.height = 0.15,               # Specify tables height
                    tables.theme = theme_cleantable(),  # Clean theme for tables
                    tables.y.text = FALSE,             # Hide tables y axis text
                    #ylab="Disease specific survival probability",
                    #ylab="Overall survival probability",
                    #ylab="progression-free interval probability",
                    ylab="disease-free interval probability",
                    xlab="Time(days)"
)
#
ss=Surv(surstats.sel$DFI.time, surstats.sel$DFI)
cox = summary(coxph(ss~factor(surstats$TE.culster.agg,levels = c("TE.high","TE.low"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)

fit2<- survfit(Surv(DFI.time, DFI) ~ TE.culster.agg, data = surstats )
dfi.agg=ggsurvplot(fit2, data =surstats,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster.agg",               # Change legend titles
                   legend.labs = c("TE.high","TE.low"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                   palette =c("#d73027","#E69F00","#00AFBB"),
                   #palette = c("#ca0020","#0571b0","#4daf4a"),
                   risk.table = TRUE,                  # Add No at risk table
                   cumevents = TRUE,                   # Add cumulative No of events table
                   tables.height = 0.15,               # Specify tables height
                   tables.theme = theme_cleantable(),  # Clean theme for tables
                   tables.y.text = FALSE,             # Hide tables y axis text
                   #ylab="Disease specific survival probability",
                   #ylab="Overall survival probability",
                   #ylab="progression-free interval probability",
                   ylab="disease-free interval probability",
                   xlab="Time(days)"
)
generate.PDF <- function(fig) {
  pdf("4.model/surC.dfi.agg.pdf",width = 5,height = 8)
  print(dfi.two)
  print(dfi.three)
  print(dfi.agg)
  dev.off()
}
generate.PDF(fig)
#####
###PFI
ss=Surv(surstats.sel$PFI.time, surstats.sel$PFI)
cox = summary(coxph(ss~factor(surstats.sel$pt.clu3,levels = c("cluster_2","cluster_1"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
#
fit1<- survfit(Surv(PFI.time, PFI) ~ pt.clu3, data = surstats.sel)
pfi.two=ggsurvplot(fit1, data =surstats.sel,
                  title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                  pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                  surv.median.line = "hv",            # Add median survival lines
                  legend.title = "TE.cluster",               # Change legend titles
                  legend.labs = c("cluster_1","cluster_2"),  # Change legend labels
                  #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                  palette = c("#ca0020","#0571b0","#4daf4a"),
                  risk.table = TRUE,                  # Add No at risk table
                  cumevents = TRUE,                   # Add cumulative No of events table
                  tables.height = 0.15,               # Specify tables height
                  tables.theme = theme_cleantable(),  # Clean theme for tables
                  tables.y.text = FALSE,             # Hide tables y axis text
                  #ylab="Disease specific survival probability",
                  #ylab="Overall survival probability",
                  ylab="progression-free interval probability",
                  #ylab="disease-free interval probability",
                  xlab="Time(days)"
)
#########
fit2<- survfit(Surv(PFI.time, PFI) ~ pt.clu3, data = surstats )
pfi.three=ggsurvplot(fit2, data =surstats,
                    #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                    pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                    surv.median.line = "hv",            # Add median survival lines
                    legend.title = "TE.cluster",               # Change legend titles
                    legend.labs = c("cluster_1","cluster_2","cluster_3"),  # Change legend labels
                    #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                    palette = c("#ca0020","#0571b0","#4daf4a"),
                    risk.table = TRUE,                  # Add No at risk table
                    cumevents = TRUE,                   # Add cumulative No of events table
                    tables.height = 0.15,               # Specify tables height
                    tables.theme = theme_cleantable(),  # Clean theme for tables
                    tables.y.text = FALSE,             # Hide tables y axis text
                    #ylab="Disease specific survival probability",
                    #ylab="Overall survival probability",
                    ylab="progression-free interval probability",
                    #ylab="disease-free interval probability",
                    xlab="Time(days)"
)
#
ss=Surv(surstats.sel$PFI.time, surstats.sel$PFI)
cox = summary(coxph(ss~factor(surstats$TE.culster.agg,levels = c("TE.high","TE.low"))))
pvalue=cox$sctest['pvalue']
hr = round(cox$conf.int[1],2)
hr_left = round(cox$conf.int[3],2)
hr_right = round(cox$conf.int[4],2)
conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
list(pvalue, hr, hr_left, hr_right)
fit2<- survfit(Surv(PFI.time, PFI) ~ TE.culster.agg, data = surstats )
pfi.agg=ggsurvplot(fit2, data =surstats,
                   title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster.agg",               # Change legend titles
                   legend.labs = c("TE.high","TE.low"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                   palette =c("#d73027","#E69F00","#00AFBB"),
                   #palette = c("#ca0020","#0571b0","#4daf4a"),
                   risk.table = TRUE,                  # Add No at risk table
                   cumevents = TRUE,                   # Add cumulative No of events table
                   tables.height = 0.15,               # Specify tables height
                   tables.theme = theme_cleantable(),  # Clean theme for tables
                   tables.y.text = FALSE,             # Hide tables y axis text
                   #ylab="Disease specific survival probability",
                   #ylab="Overall survival probability",
                   ylab="progression-free interval probability",
                   #ylab="disease-free interval probability",
                   xlab="Time(days)"
)
generate.PDF <- function(fig) {
  pdf("4.model/surC.pfi.agg.pdf",width = 5,height = 8)
  print(pfi.two)
  print(pfi.three)
  print(pfi.agg)
  dev.off()
}
generate.PDF(fig)
###########

#####################

#########immune pathway
####
score=score.data
score=as.data.frame(scale(score))
head(score)
ids=intersect(rownames(clu), rownames(score))

cbscoreshow=cbind(clu[ids,]$pt.clu3, score[ids,])
colnames(cbscoreshow)[1]="cluster"
cbscoreshow$cluster.agg=ifelse(cbscoreshow$cluster=="cluster_1","TE.high","TE.low")
head(cbscoreshow)
dim(cbscoreshow)
#cbscoreshow=subset(cbscoreshow, cluster=="cluster_1" | cluster=="cluster_2")
#######
library(tidyr)
drawdatashow <- cbscoreshow %>% pivot_longer(cols=colnames(cbscoreshow)[-c(1,31)],
                                             names_to= "signature",
                                             values_to = "score")
drawdatashow=as.data.frame(drawdatashow)

#drawdatashow$signature <- factor(drawdatashow$signature, levels = rownames(sigshow))
farb=c("#ca0020","#0571b0","#4daf4a","#27408B","#FF0000","#2E8B57","#CD00CD")
farb=c("#d73027", "#E69F00","#00AFBB")
#head(drawdatashow)
#table(drawdatashow$signature)
drawdatashow$signature <- factor(drawdatashow$signature, levels = colnames(cbscoreshow))
p1=ggplot(drawdatashow, aes(x=cluster, y=score, group=cluster)) + 
  geom_boxplot(aes(fill=cluster),alpha = 0.8,outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.4)+
  stat_compare_means(label = "p.signif")+
  facet_grid(signature~. )+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+rotate()+xlab("TE.cluster")+ylab("Z.score")
p2=ggplot(subset(drawdatashow, !(cluster=="cluster_3")), aes(x=cluster, y=score, group=cluster)) + 
  geom_boxplot(aes(fill=cluster),alpha = 0.8,outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.4)+
  stat_compare_means(label = "p.signif")+
  facet_grid(signature~. )+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+rotate()+xlab("TE.cluster")+ylab("Z.score")
p3=ggplot(drawdatashow, aes(x=cluster.agg, y=score, group=cluster.agg)) + 
  geom_boxplot(aes(fill=cluster.agg),alpha = 0.8,outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.4)+
  stat_compare_means(label = "p.signif")+
  facet_grid(signature~. )+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+rotate()+xlab("TE.cluster")+ylab("Z.score")
#+theme_classic()
#theme(strip.text.x = element_text(size = 8, colour = "orange", angle = 90))
generate.PDF <- function(fig) {
  pdf("4.model/immune/genesets.29s.three.clusters.agg.pdf",height  = 8)
  print(p1)
  print(p2)
  print(p3)
  dev.off()
}
generate.PDF(fig)
############cell types
load("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/cell.type.28s.gsva.RData")
score=genesets17s.score
#
score=as.data.frame(scale(score))
head(score)
ids=intersect(rownames(clu), rownames(score))

cbscoreshow=cbind(clu[ids,]$pt.clu3, score[ids,])
colnames(cbscoreshow)[1]="cluster"
cbscoreshow$cluster.agg=ifelse(cbscoreshow$cluster=="cluster_1","TE.high","TE.low")
head(cbscoreshow)
#cbscoreshow=subset(cbscoreshow, cluster=="cluster_1" | cluster=="cluster_2")
#######
library(tidyr)
drawdatashow <- cbscoreshow %>% pivot_longer(cols=colnames(cbscoreshow)[-c(1,30)],
                                             names_to= "signature",
                                             values_to = "score")
drawdatashow=as.data.frame(drawdatashow)

#drawdatashow$signature <- factor(drawdatashow$signature, levels = rownames(sigshow))
#farb=c("#ca0020","#0571b0","#4daf4a","#27408B","#FF0000","#2E8B57","#CD00CD")
#head(drawdatashow)
#table(drawdatashow$signature)
farb=c("#d73027","#E69F00","#00AFBB")
p1=ggplot(drawdatashow, aes(x=cluster, y=score, group=cluster)) + 
  geom_boxplot(aes(fill=cluster),alpha = 0.8,outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.4)+
  stat_compare_means(label = "p.signif")+
  facet_grid(signature~. )+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+rotate()+xlab("TE.cluster")+ylab("Z.score")
p2=ggplot(drawdatashow, aes(x=cluster.agg, y=score, group=cluster.agg)) + 
  geom_boxplot(aes(fill=cluster.agg),alpha = 0.8,outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.4)+
  stat_compare_means(label = "p.signif")+
  facet_grid(signature~. )+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+rotate()+xlab("TE.cluster")+ylab("Z.score")

#+theme_classic()
#theme(strip.text.x = element_text(size = 8, colour = "orange", angle = 90))
generate.PDF <- function(fig) {
  pdf("4.model/immune/genesets.cell.type.agg.pdf",height  = 8)
  print(p1)
  print(p2)
  dev.off()
}
generate.PDF(fig)
####
#######
#tclu=subset(clu, !(pt.clu3=="cluster_3")& !(msi=="Indeterminate") )
tclu=subset(clu, !(pt.clu3=="cluster_3")& !(msi=="Indeterminate") )
table(tclu$pt.clu3,tclu$msi)
chisq.test(tclu$pt.clu3,tclu$msi)
######
save(res.Pcor.29s.pos.sig,
     te.counts,
     cndi.rep,
     clu,
     dat,
     surstats,
     cbscoreshow,
     sur.kirc,
     stat.exp.KIRC,
     score.data,file = "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/model.cor.>=0.4.totoal.surTEs.RData"
     )
save.image(file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/tmp.sur.immu.overlap.RData")
#load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/tmp.sur.immu.overlap.RData")
#setwd("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/")
#######
save(cor.exp, expM, file = "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/5.TE.cor.with.gene/TE.gene.CRC.matrix.RData")
########
######
#################compare the averaage expression of 9 tes
#######
mean.stats=cbind(clus, as.data.frame(t(dat))[rownames(clus),])
head(mean.stats)
mean.stats$mean.exp=rowSums(mean.stats[,3:11])/9
mean.stats$TE.cluster.agg=ifelse(mean.stats$TE.cluster=="cluster_1","TE.high","TE.low")
save(mean.stats, file="4.model/Mean.expression.of.9.TE.model.RData")
#####
my_comparisons <- list( c("cluster_1", "cluster_2"), c("cluster_1", "cluster_3"), c("cluster_2", "cluster_3") )
gepm1=ggviolin(mean.stats,x = "TE.cluster", y ="mean.exp" , fill = "TE.cluster",alpha = 1,size = 0.3,
              palette = c("#d73027","#00AFBB", "#E69F00"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+theme_bw()+xlab("TE.cluster")+ylab("mean.expression.of.9.TEs")

gepm2=ggviolin(mean.stats,x = "TE.cluster.agg", y ="mean.exp" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
              palette = c("#d73027","#E69F00","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster.agg")+
  stat_compare_means(label = "p.signif")+theme_bw()+xlab("TE.cluster")+ylab("mean.expression.of.9.TEs")+ggtitle("CRC")
########
generate.PDF <- function(fig) {
  pdf("4.model/Mean.expression.of.9.TE.model.pdf",height  = 4,width = 3)
  print(gepm1)
  print(gepm2)
  dev.off()
}
generate.PDF(fig)

