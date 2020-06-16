######################
######################TE score, TE five clss and 9 individual TEs correlated with GEP at pan cancer
######################
#########load GEP score at pan cancer level
load("~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/GEP.ICB.score.pancancer.RData")
head(preexp.cal.pancancer)
pancancer.GEP=preexp.cal.pancancer
pancancer.GEP$id=rownames(pancancer.GEP)
head(pancancer.GEP)
######### load TE score pancancer
load("4.model/kaps/other.cancers/cox.meanexp.pancancer.RData")
head(pan9TE.stat)
pancancer.TEscore=pan9TE.stat
pancancer.TEscore$id=substr(rownames(pancancer.TEscore),1,12)
head(pancancer.TEscore)
dim(pancancer.TEscore)
######## load the five class level exp pancancer
load("8.KIRC/TE.class.compate.kaps.group/TEclass.level.expM.pancancer.RData")
head(TE.classexp.kirc.agg)
pancancer.TE.classEXP=TE.classexp.kirc.agg[,c(2,4,6,10,14)]
head(pancancer.TE.classEXP)
######## add more CRC samples for five class expression
load("4.model/kaps/TEclass.level.compare.kaps.groups.RData")
head(kaps.te.class.compare)
class.TEexp.CRC= kaps.te.class.compare[,colnames(kaps.te.class.compare) %in% colnames(pancancer.TE.classEXP)]
head(class.TEexp.CRC)
#
########combine TE expression firstly
TE.score.class.exp=rbind(class.TEexp.CRC, pancancer.TE.classEXP)
TE.score.class.exp$id=rownames(TE.score.class.exp)
TE.score.class.exp.full.CRC=TE.score.class.exp[!(duplicated(TE.score.class.exp$id)),]
head(TE.score.class.exp.full.CRC)
dim(TE.score.class.exp.full.CRC)
####
dim(pancancer.TEscore)
head(pancancer.TEscore)
length(intersect(rownames(TE.score.class.exp.full.CRC),rownames(pancancer.TEscore)))   ###6554
all(rownames(pancancer.TEscore) %in% rownames(TE.score.class.exp.full.CRC))
TE.score.class.exp.full.CRC.ok=cbind(TE.score.class.exp.full.CRC[rownames(pancancer.TEscore),-6], pancancer.TEscore)
head(TE.score.class.exp.full.CRC.ok)
write.csv(pancancer.TEscore,"4.model/kaps/other.cancers/GEP/pancancer.TEscore.csv")
#########combine GEP
length(intersect(pancancer.GEP$id, TE.score.class.exp.full.CRC.ok$id))
stat.TEscore.vs.GEP=merge(pancancer.GEP[,c(19,20)], TE.score.class.exp.full.CRC.ok, by="id",all=FALSE)
dim(stat.TEscore.vs.GEP)
#########
#########
head(stat.TEscore.vs.GEP)
table(stat.TEscore.vs.GEP$type)
class(stat.TEscore.vs.GEP$DNA)
#########
library(ggpubr)
pdf("4.model/kaps/other.cancers/GEP/GEP.vs.TE.score.label.corrected.pdf",width = 7,height = 7)
ggscatter(stat.TEscore.vs.GEP, x = "z.of.mean.exp", y = "GEP",
        add = "reg.line",                                 # Add regression line
        conf.int = TRUE,                                  # Add confidence interval
        add.params = list(color = "black",fill = "lightgray"),
        color = "type"#, palette =  c("#27408B","#FF0000","#2E8B57","#CD00CD","#EE7942","#009ACD")
        )+xlab("z of TE score")+ylab("GEP")+
 stat_cor(method = "spearman", label.x = 2, label.y = 2)
dev.off()
#########calculate correlation for TE score and five class vs GEP in each cancer type
#########
datalist1=list()
datalist2=list()
canceridx=unique(stat.TEscore.vs.GEP$type)
for (i in 1:length(canceridx)) {
  for (j in c(3:16,27)) {
    aa=subset(stat.TEscore.vs.GEP, type==canceridx[i])
    res <- cor.test(aa[,2], aa[,j],  method = "spearman")
    pvalue=res$p.value
    cor=res$estimate
    res=data.frame(pvalue, cor)
    rownames(res)=paste(colnames(stat.TEscore.vs.GEP)[j],sep = "_")
    colnames(res)=paste(colnames(res),canceridx[i],sep = "_")
    datalist1[[j]] <- res
  }
  res1=do.call(rbind, datalist1)
  datalist2[[i]]=res1
}
res.cor.GEP.vs.TE.pancacner=do.call(cbind, datalist2)
################draw balloon plot
################
aa=res.cor.GEP.vs.TE.pancacner[,seq(2,ncol(res.cor.GEP.vs.TE.pancacner),2)]
colnames(aa)=substr(colnames(aa),5, nchar(as.character(colnames(aa))))
aa
pdf("4.model/kaps/other.cancers/GEP/GEP.vs.TE.score.class.pdf",width = 10,height = 7)
ggballoonplot(aa, fill = "value",size.range = c(1, 5))+
  scale_fill_gradientn(colors = c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
                                  "#E16462FF", "#FCA636FF", "#F0F921FF"))+
  ggtitle(paste0("correlation of TE ","\n","expression vs GEP"))+theme_bw()+xlab("cohort")+ylab("TE variable")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
#
ggballoonplot(aa[-c(1:5),], fill = "value",size.range = c(1, 5))+
  scale_fill_gradientn(colors = c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
                                  "#E16462FF", "#FCA636FF", "#F0F921FF"))+
  ggtitle(paste0("correlation of TE ","expression vs GEP"))+theme_bw()+xlab("cohort")+ylab("TE variable")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#
pdf("4.model/kaps/other.cancers/GEP/GEP.vs.TE.score.only.pdf",width = 10,height = 3)
ggballoonplot(aa[-c(1:14),], fill = "value",size.range = c(5, 10))+
  scale_fill_gradientn(colors = c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
                                  "#E16462FF", "#FCA636FF", "#F0F921FF"))+
  ggtitle(paste0("correlation of TE ","expression vs GEP"))+theme_bw()+xlab("cohort")+ylab("TE variable")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
#
pdf("4.model/kaps/other.cancers/GEP/GEP.vs.TE.indivi.pdf",width = 10,height = 5)
ggballoonplot(aa[-c(1:5,15),], fill = "value",size.range = c(5, 10))+
  scale_fill_gradientn(colors = c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
                                  "#E16462FF", "#FCA636FF", "#F0F921FF"))+
  ggtitle(paste0("correlation of TE ","expression vs GEP"))+theme_bw()+xlab("cohort")+ylab("TE variable")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
###############
save(stat.TEscore.vs.GEP, TE.score.class.exp.full.CRC.ok, file="4.model/kaps/other.cancers/GEP/TE.score.class.exp.full.CRC.RData")
###############
###############
###############correlation volcano plot for TE score just
###############
corvol.pancaner=res.cor.GEP.vs.TE.pancacner[-c(1:14),]
corvol.pancaner=as.data.frame(t(corvol.pancaner))
corvol.pancaner$id=rownames(corvol.pancaner)
corvol.pancaner1=corvol.pancaner[seq(2,nrow(corvol.pancaner),2),]
colnames(corvol.pancaner1)[1]="cor"
corvol.pancaner2=corvol.pancaner[seq(1,nrow(corvol.pancaner),2),]
colnames(corvol.pancaner2)[1]="p.value"
####
corvol.pancaner.drawdata=cbind(corvol.pancaner1,corvol.pancaner2)
corvol.pancaner.drawdata$type=substr(corvol.pancaner.drawdata$id,5, nchar(as.character(corvol.pancaner.drawdata$id)))
corvol.pancaner.drawdata=corvol.pancaner.drawdata[,c(1,3,5)]
corvol.pancaner.drawdata$pvalue.2=-log10(corvol.pancaner.drawdata$p.value)
corvol.pancaner.drawdata$pvalue.2=gsub("Inf",1000,corvol.pancaner.drawdata$pvalue.2)
corvol.pancaner.drawdata$pvalue.2=as.numeric(paste(corvol.pancaner.drawdata$pvalue.2))
head(corvol.pancaner.drawdata)
summary((corvol.pancaner.drawdata$cor))
class(corvol.pancaner.drawdata$cor)
write.csv(corvol.pancaner.drawdata,"4.model/kaps/other.cancers/GEP/corvol.pancaner.drawdata.csv")
#####
ggplot(corvol.pancaner.drawdata, aes(x=p.value, y=cor)) + geom_point()
library(easyGgplot2)
require("ggrepel")
pdf("4.model/kaps/other.cancers/GEP/cor.volcano.TE.score.vs.GEP.pancancer.new2.pdf",width = 5,height = 5)
ggplot2.scatterplot(data=corvol.pancaner.drawdata, xName='p.value',yName='cor', size=9, 
                    mapping=aes(alpha = cor)#, color="purple"
                    )+
  geom_point(size=9,fill="purple",shape = 21,colour="black",stroke = 0.5)+
  theme(axis.title=element_text(size=20),
        axis.text = element_text(face = "bold",size = 16),
        axis.ticks.length=unit(.4, "cm"),
        axis.ticks = element_line(colour = "black", size = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        plot.margin = margin(1, 1, 1, 1, "cm")
  )+
  geom_text_repel(aes(label = corvol.pancaner.drawdata$type),size = 3.5)+
  labs(x="P.value",y="Spearman correlation (r)")+ggtitle("TE score vs GEP across pan cancer")+
  geom_vline(xintercept=0.05,colour="gray50",linetype=2,size=0.6)+
  scale_x_continuous(breaks=c(0,0.05,0.2,0.4,0.6))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.5,0.55))
  #labs(x=bquote(-log[10]~italic("P")),y="Spearman correlation (r)")+
  #geom_text(label=corvol.pancaner.drawdata$type,vjust=0,hjust=-0.15 )+
dev.off()


  


