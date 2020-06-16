########validation on other cancers
########
#pan cancer data
load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/pancancer.TEM.clin.RData")
dim(stat.rexp)
dim(stat.surclin)
#
head(stat.surclin)
table(stat.surclin$type)
###########calculate for each cancer type
###########
exp9TEs.pancacner=stat.rexp[, colnames(stat.rexp) %in% rownames(cndi.rep) ]
head(exp9TEs.pancacner)
#########delete the CRC data
stat.surclin.noCRC=subset(stat.surclin, !(type=="COAD"|type=="READ"))
table(stat.surclin.noCRC$type)
exp9TEs.pancacner.noCRC=exp9TEs.pancacner[rownames(stat.surclin.noCRC),]
#########add 590 CRC back clin
#########
stat.surclin.pan.fullCRC=data.frame(type=stat.surclin.noCRC$type,
                                    OS=stat.surclin.noCRC$OS,
                                    OS.time=stat.surclin.noCRC$OS.time,
                                    DSS=stat.surclin.noCRC$DSS,
                                    DSS.time=stat.surclin.noCRC$DSS.time,
                                    DFI=stat.surclin.noCRC$DFI,
                                    DFI.time=stat.surclin.noCRC$DFI.time,
                                    PFI=stat.surclin.noCRC$PFI,
                                    PFI.time=stat.surclin.noCRC$PFI.time
                                    )
rownames(stat.surclin.pan.fullCRC)=rownames(stat.surclin.noCRC)
#####deal with NA
stat.surclin.pan.fullCRC[,-1] <- lapply(stat.surclin.pan.fullCRC[,-1], function(x) as.numeric(as.character(x)))
head(stat.surclin.pan.fullCRC)
#
ttp=kaps.td[,c(30:37)]
ttp$type=rep("CRC",times=nrow(ttp))
head(ttp)
#####
stat.surclin.pan.fullCRC=rbind(stat.surclin.pan.fullCRC, ttp)
dim(stat.surclin.pan.fullCRC)
#####full 9 TE expression
exp9TEs.pancacner.fullCRC=rbind(exp9TEs.pancacner.noCRC, kaps.td[,c(3:11)])
dim(exp9TEs.pancacner.fullCRC)
all(rownames(exp9TEs.pancacner.fullCRC) %in% rownames(stat.surclin.pan.fullCRC))
####################
#######draw the correlation matrix
#######
head(exp9TEs.pancacner.fullCRC)
head(stat.surclin.pan.fullCRC)
save(exp9TEs.pancacner.fullCRC,stat.surclin.pan.fullCRC,file="4.model/kaps/other.cancers/TE.9s.exp.and.clin.pancancer.RData" )
#
cancerT.9TE.matrix=cbind(stat.surclin.pan.fullCRC,exp9TEs.pancacner.fullCRC[rownames(stat.surclin.pan.fullCRC),])
cancerT.9TE.matrix=cancerT.9TE.matrix[,-c(2:9)]
#
library(dplyr)
cancerT.9TE.matrix.agg= cancerT.9TE.matrix %>% group_by(type) %>% summarise_all(mean)
cancerT.9TE.matrix.agg=as.data.frame(cancerT.9TE.matrix.agg)
rownames(cancerT.9TE.matrix.agg)=cancerT.9TE.matrix.agg$type
cancerT.9TE.matrix.agg=cancerT.9TE.matrix.agg[,-1]
head(cancerT.9TE.matrix.agg)
##
dat.pan <- as.data.frame(t(cancerT.9TE.matrix.agg))
dist_dat.pan  <- cor(dat.pan,method="spearman")
hc.pan <- hclust(as.dist(1-dist_dat.pan ),method="ward.D2")
library(viridis)
library(pheatmap)
mycolor=col = rev(viridis(10))
#col_row$ <- factor(col_row$tumor)
pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/other.cancers/corrlation.Matrix.pdf")
pheatmap(dist_dat.pan,show_rownames = T,border_color = NA,
         show_colnames = F,cluster_rows = hc.pan,
         cluster_cols =hc.pan,#annotation_col = rclin[,c(1,6,7)],
         color=mycolor, annotation_legend = T
         #annotation_colors  = list(response = c("PD" = "#dd1c77","PRCR"="#74c476","SD"="#feb24c"),
         #                         treatment.status=c("On"="#74a9cf", "Pre"="#f768a1"))
)
dev.off()
pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/other.cancers/ht.9tes.across.pancancer.pdf")
pheatmap(cancerT.9TE.matrix.agg,show_rownames = T,border_color = NA,
         clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",
         show_colnames = T,#cluster_rows = hc.pan,
         #cluster_cols =hc.pan,#annotation_col = rclin[,c(1,6,7)],
         color=mycolor, annotation_legend = T
         #annotation_colors  = list(response = c("PD" = "#dd1c77","PRCR"="#74c476","SD"="#feb24c"),
         #                         treatment.status=c("On"="#74a9cf", "Pre"="#f768a1"))
)
pheatmap(scale(cancerT.9TE.matrix.agg),show_rownames = T,border_color = NA,
         clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",
         show_colnames = T,#cluster_rows = hc.pan,
         #cluster_cols =hc.pan,#annotation_col = rclin[,c(1,6,7)],
         color=mycolor, annotation_legend = T
         #annotation_colors  = list(response = c("PD" = "#dd1c77","PRCR"="#74c476","SD"="#feb24c"),
         #                         treatment.status=c("On"="#74a9cf", "Pre"="#f768a1"))
)
dev.off()

########################
#################mean exp
#########
head(exp9TEs.pancacner.fullCRC)
head(stat.surclin.pan.fullCRC)
#
pan9TE.stat=exp9TEs.pancacner.fullCRC
pan9TE.stat$mean.exp=rowMeans(pan9TE.stat,na.rm = T)
pan9TE.stat=cbind(pan9TE.stat, stat.surclin.pan.fullCRC[rownames(pan9TE.stat),])
pan9TE.stat$z.of.mean.exp=(pan9TE.stat$mean.exp - mean(pan9TE.stat$mean.exp))/sd(pan9TE.stat$mean.exp)
head(pan9TE.stat)
#
datalist=list()
panidx=unique(pan9TE.stat$type)
for (i in 1:length(panidx)) {
  aa=subset(pan9TE.stat, type==panidx[i])
  ss=Surv(aa$OS.time, aa$OS)
  cox = summary(coxph(ss~mean.exp,data = aa))
  hr_left=round(cox$conf.int[3],2)
  hr_right=round(cox$conf.int[4],2)
  bb=data.frame(pvalue=cox$sctest['pvalue'],
                hr = round(cox$conf.int[1],2),
                hr_left = round(cox$conf.int[3],2),
                hr_right = round(cox$conf.int[4],2),
                conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep="")
                )
  rownames(bb)=panidx[i]
  datalist[[i]]=bb
}
cox.meanexp.pancancer=do.call(rbind, datalist)
subset(cox.meanexp.pancancer, pvalue<0.05)
subset(cox.meanexp.pancancer, hr>1)
save(pan9TE.stat, cox.meanexp.pancancer,file="4.model/kaps/other.cancers/cox.meanexp.pancancer.RData")
################
############9TE mean expression among cancer types
############
pdf("4.model/kaps/other.cancers/box.meanexp.cancer.type.pdf",width = 10,height = 5)
ggboxplot(pan9TE.stat, x = "type", y = "mean.exp",
          color = "type"#, #palette = c("Yes"="#d01c8b","No"="#4dac26"),
          #add = "jitter")+stat_compare_means(aes(group = DFS.status)
          )+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")

aa=aggregate(mean.exp~type,pan9TE.stat,median)
aa=aa[order(aa$mean.exp,decreasing = T),]
ggplot(pan9TE.stat, aes(x=type, y=mean.exp, fill=type)) + 
  geom_boxplot()+scale_x_discrete(limits=aa$type)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
        )#+
  #scale_fill_manual(values = brewer.pal(24, "gist_ncar"), guide = guide_legend(title = "cancer.type"))
ggplot(pan9TE.stat, aes(x=type, y=z.of.mean.exp, fill=type)) + 
  geom_boxplot()+scale_x_discrete(limits=aa$type)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())
#colorRampPalette(rev(brewer.pal(n = 24, name ="gist_ncar")))(8)
dev.off()
########plot the forest
library(magrittr)
library(ggforestplot)
########
cox.meanexp.pancancer$type=rownames(cox.meanexp.pancancer)
cox.meanexp.pancancer=cox.meanexp.pancancer[order(cox.meanexp.pancancer$type),]
cox.meanexp.pancancer$pvalue.round=round(cox.meanexp.pancancer$pvalue,digits = 4)
cox.meanexp.pancancer$pos=seq(1:nrow(cox.meanexp.pancancer))
pdf("4.model/kaps/other.cancers/Univariate.Cox.HR.distribution.new.pdf",width = 14,height = 5)
ggplot(cox.meanexp.pancancer, aes( type,hr, colour = type)) + 
  geom_pointrange(aes(ymin = hr_left, ymax = hr_right),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  #scale_color_manual(values = cbPalette)+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Hazard ratio")+ggtitle("univariate Cox regression of z.TE.score across pan cancer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
        )+
        # Add striped background
        #geom_stripes(odd = "#33333333", even = "#00000000")+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+scale_y_continuous(trans='log2',breaks=c(0,1,2,3, 25))+
  annotate("text", x=1, y=26, label=paste0("ns","\n","n=79"),size=3)+
  annotate("text", x=2, y=26, label=paste0("ns","\n","n=211"),size=3)+
  annotate("text", x=3, y=26, label=paste0("***","\n","n=982"),size=3)+
  annotate("text", x=4, y=26, label=paste0("*","\n","n=144"),size=3)+
  annotate("text", x=5, y=26, label=paste0("**","\n","n=590"),size=3)+
  annotate("text", x=6, y=26, label=paste0("ns","\n","n=27"),size=3)+
  annotate("text", x=7, y=26, label=paste0("ns","\n","n=154"),size=3)+
  annotate("text", x=8, y=26, label=paste0("***","\n","n=403"),size=3)+
  annotate("text", x=9, y=26, label=paste0("ns","\n","n=66"),size=3)+
  annotate("text", x=10, y=26, label=paste0("****","\n","n=494"),size=3)+
  annotate("text", x=11, y=26, label=paste0("ns","\n","n=161"),size=3)+
  annotate("text", x=12, y=26, label=paste0("ns","\n","n=300"),size=3)+
  annotate("text", x=13, y=26, label=paste0("*","\n","n=147"),size=3)+
  annotate("text", x=14, y=26, label=paste0("ns","\n","n=486"),size=3)+
  annotate("text", x=15, y=26, label=paste0("ns","\n","n=481"),size=3)+
  annotate("text", x=16, y=26, label=paste0("**","\n","n=404"),size=3)+
  annotate("text", x=17, y=26, label=paste0("ns","\n","n=56"),size=3)+
  annotate("text", x=18, y=26, label=paste0("ns","\n","n=246"),size=3)+
  annotate("text", x=19, y=26, label=paste0("ns","\n","n=103"),size=3)+
  annotate("text", x=20, y=26, label=paste0("ns","\n","n=73"),size=3)+
  annotate("text", x=21, y=26, label=paste0("ns","\n","n=282"),size=3)+
  annotate("text", x=22, y=26, label=paste0("ns","\n","n=491"),size=3)+
  annotate("text", x=23, y=26, label=paste0("ns","\n","n=118"),size=3)+
  annotate("text", x=24, y=26, label=paste0("ns","\n","n=56"),size=3)
  
dev.off()
#ns: p > 0.05
#*: p <= 0.05
#**: p <= 0.01
#***: p <= 0.001
#****: p <= 0.0001
##############
##############
write.csv(cox.meanexp.pancancer,"4.model/kaps/other.cancers/Univariate.Cox.HR.distribution.csv")
