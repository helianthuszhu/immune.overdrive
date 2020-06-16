#######
my_comparisons <- list( c("PD", "PRCR"), c("PD", "SD"))
pdf("6.immnue.treat.data/stat.gse91061.20200523/Pre.samples.TE.score.comparison.pdf",width = 4,height = 5)
ggboxplot(subset(stat.reponse.data.gse91061,treat.status=="Pre" ), x = "reponse2", y = "z.of.mean.exp",
          color = "reponse2", palette = c("PD"="#d01c8b","PRCR"="#4dac26","SD"="#00AFBB"),
        add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
         stat_compare_means(label.y = 7)+ggtitle("Pre samples")

ggboxplot(subset(stat.reponse.data.gse91061,treat.status=="Pre" ), x = "reponse2", y = "MER57F",
          color = "reponse2", palette = c("PD"="#d01c8b","PRCR"="#4dac26","SD"="#00AFBB"),
          add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
          stat_compare_means(label.y = 7)+ggtitle("Pre samples")

ggboxplot(subset(stat.reponse.data.gse91061,treat.status=="Pre" ), x = "reponse2", y = "LTR21B",
          color = "reponse2", palette = c("PD"="#d01c8b","PRCR"="#4dac26","SD"="#00AFBB"),
          add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 7)+ggtitle("Pre samples")
dev.off()
#
pdf("6.immnue.treat.data/stat.gse91061.20200523/On.samples.TE.score.comparison.pdf",width = 4,height = 5)
ggboxplot(subset(stat.reponse.data.gse91061,treat.status=="On" ), x = "reponse2", y = "z.of.mean.exp",
          color = "reponse2", palette = c("PD"="#d01c8b","PRCR"="#4dac26","SD"="#00AFBB"),
          add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 7)+ggtitle("Pre samples")

ggboxplot(subset(stat.reponse.data.gse91061,treat.status=="On" ), x = "reponse2", y = "MER57F",
          color = "reponse2", palette = c("PD"="#d01c8b","PRCR"="#4dac26","SD"="#00AFBB"),
          add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 7)+ggtitle("Pre samples")

ggboxplot(subset(stat.reponse.data.gse91061,treat.status=="On" ), x = "reponse2", y = "LTR21B",
          color = "reponse2", palette = c("PD"="#d01c8b","PRCR"="#4dac26","SD"="#00AFBB"),
          add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 7)+ggtitle("Pre samples")
dev.off()

#########################
#########################
head(clin.GSE91061)
clin.GSE91061.full=read.table("6.immnue.treat.data/stat.gse91061.20200523/clin.gse91061.txt",header = T,sep = "\t")
colnames(clin.GSE91061.full)[1]="pt"
head(clin.GSE91061.full)
dim(clin.GSE91061.full)
length(unique(clin.GSE91061.full$Patient))
length(intersect(unique(stat.reponse.data.gse91061$pt), clin.GSE91061.full$Patient))
stat.reponse.data.gse91061.Pre=subset(stat.reponse.data.gse91061, treat.status=="Pre")
##############
stat.clin.gse91061=merge(stat.reponse.data.gse91061.Pre, clin.GSE91061.full,by="pt",all=FALSE)
head(stat.clin.gse91061)
##############tide score
tidescore.gse91061=read.table("6.immnue.treat.data/stat.gse91061.20200523/TIDE.score.GSE91061.txt",header = T,sep = "\t")
colnames(tidescore.gse91061)[1]="pt"
head(tidescore.gse91061)
dim(tidescore.gse91061)
stat.clin.gse91061=merge(stat.clin.gse91061, tidescore.gse91061, by="pt",all.x=TRUE)
##########
#########expMatrix
#
geneexpM.gse91061=read.csv("~/nas/Xiaoqiang/opti.data/data.melanoma/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",header = T)
geneexpM.gse91061[1:4,1:6]
#
library(clusterProfiler)
gidgse91061 = bitr(geneexpM.gse91061$X, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
head(gidgse91061)
#
colnames(geneexpM.gse91061)[1]="ENTREZID"
geneexpM.gse91061=merge(gidgse91061, geneexpM.gse91061, by="ENTREZID",all.x=TRUE)
dim(geneexpM.gse91061)
save(geneexpM.gse91061, file="6.immnue.treat.data/stat.gse91061.20200523/geneexpMatrix.gse91061.RData")
########PDL1 exp
rownames(geneexpM.gse91061)=geneexpM.gse91061$SYMBOL
geneexpM.gse91061=geneexpM.gse91061[,-c(1,2)]
#
pdl1.gse91061=as.data.frame(t(geneexpM.gse91061[which(rownames(geneexpM.gse91061)=="CD274"),]))
pdl1.gse91061$id=rownames(pdl1.gse91061)
pdl1.gse91061=separate(pdl1.gse91061,"id", into=c("pt","ss","idn"),sep = "_")
head(pdl1.gse91061)
dim(pdl1.gse91061)
pdl1.gse91061.pre=subset(pdl1.gse91061, ss=="Pre")
######
######
stat.clin.gse91061.ok=merge(stat.clin.gse91061, pdl1.gse91061.pre, by="pt",all.x=TRUE)
head(stat.clin.gse91061.ok)
table(stat.clin.gse91061.ok$response)
###########
#pROC
library(ggplot2)
library(pROC)
head(stat.clin.gse91061)
stat.clin.gse91061.ok$response.status=ifelse(stat.clin.gse91061.ok$response=="PD",0,1)
roc.list.gse91061.Pre <- roc(response.status ~z.of.mean.exp+M.stage.new+Mutation.Load+Neo.antigen.Load+Neo.peptide.Load+ Dysfunction+CD274
                               , 
                         auc=TRUE, ci=TRUE,data = stat.clin.gse91061.ok)


var.roc.gse91061=glm(response.status~z.of.mean.exp+M.stage.new+Mutation.Load+Neo.antigen.Load+Neo.peptide.Load+ Dysfunction+CD274,family=binomial(link="logit"),
                     data=stat.clin.gse91061.ok)

stat.clin.gse91061.ok$combined.all=1.7538564+ stat.clin.gse91061.ok$z.of.mean.exp*0.6688903+ stat.clin.gse91061.ok$M.stage.new*(-0.4239713)+
                     stat.clin.gse91061.ok$Mutation.Load*(-0.0120179)+ stat.clin.gse91061.ok$Neo.antigen.Load*0.0179395+ 
                     stat.clin.gse91061.ok$Neo.peptide.Load*0.0009811+ stat.clin.gse91061.ok$Dysfunction*(-0.8245884)+ stat.clin.gse91061.ok$CD274*0.0271402

roc.list.gse91061.Pre <- roc(response.status ~z.of.mean.exp+M.stage.new+Mutation.Load+Neo.antigen.Load+Neo.peptide.Load+ Dysfunction+CD274+combined.all
                             , 
                             auc=TRUE, ci=TRUE,data = stat.clin.gse91061.ok)
pdf("6.immnue.treat.data/stat.gse91061.20200523/ROC.GSE91061.pdf",width = 7,height = 7)
ggroc(roc.list.gse91061.Pre,linetype=7, size=1, alpha=0.8)+ geom_abline(slope=1, intercept=1,color="grey") + 
  theme_classic()+
  scale_colour_manual(values = c("z.of.mean.exp"="#e31a1c","M.stage.new"="#d95f02","Mutation.Load"="#7570b3",
                                 "Neo.antigen.Load"="#e7298a","Neo.peptide.Load"="#66a61e","Dysfunction"="#e6ab02","CD274"="#00AFBB", "combined.all"="black"), 
                      name = "risk factors") +
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$z.of.mean.exp$auc,2)," 95%CI: 0.49-0.73"), color="#e31a1c", x= 0.22, y=0.40, size=4) +
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$M.stage.new$auc,2), " 95%CI: 0.52-0.76"), color="#d95f02", x= 0.22, y=0.35, size =4) +
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$Mutation.Load$auc,2), " 95%CI: 0.39-0.65"), color="#7570b3", x= 0.22, y=0.30, size =4)+
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$Neo.antigen.Load$auc,2), " 95%CI: 0.42-0.67"), color="#e7298a", x= 0.22, y=0.25, size =4)+
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$Neo.peptide.Load$auc,2), " 95%CI: 0.43-0.68"), color="#66a61e", x= 0.22, y=0.20, size =4)+
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$Dysfunction$auc,2), " 95%CI: 0.45-0.69"), color="#e6ab02", x= 0.22, y=0.15, size =4)+
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$CD274$auc,2), " 95%CI: 0.46-0.70"), color="#00AFBB", x= 0.22, y=0.10, size =4)+
  annotate("text", label=paste0("AUC=",signif(roc.list.gse91061.Pre$combined.all$auc,2), " 95%CI: 0.57-0.91"), color="black", x= 0.22, y=0.05, size =4)+
  labs( x = "Specificity", y = "Sensitivity", size =5) + theme(axis.text=element_text(size=12),
                                                               axis.title=element_text(size=14,face="bold"), 
                                                               legend.text=element_text(size=14),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()
#######
dim(stat.clin.gse91061.ok)
#######
#######
table(stat.clin.gse91061.ok$response.status)
#######
oddsmodel.gse91061 <- glm( response.status ~ z.of.mean.exp+M.stage.new+Mutation.Load+Neo.antigen.Load+Neo.peptide.Load+ Dysfunction+CD274, 
                           data = stat.clin.gse91061.ok, family = binomial)
oddsmodel.gse91061 <- glm( response.status ~ z.of.mean.exp+M.stage.new+Mutation.Load+Neo.antigen.Load, 
                           data = stat.clin.gse91061.ok, family = binomial)
oddres1gse91061=as.data.frame(summary(oddsmodel.gse91061)$coef)
oddres1gse91061$odd.ratio=exp(oddres1gse91061$Estimate)
oddres2gse91061=as.data.frame(exp(confint(oddsmodel.gse91061,level = 0.9)))
oddres.ok.gse91061=cbind(oddres1gse91061,oddres2gse91061)
oddres.ok.gse91061=oddres.ok.gse91061[-1,]
oddres.ok.gse91061$variable=rownames(oddres.ok.gse91061)
#
oddres.ok.gse91061
#
pdf("6.immnue.treat.data/stat.gse91061.20200523/odds.gse91061.pdf",width = 6,height = 3)
ggplot(oddres.ok.gse91061, aes( variable, odd.ratio, colour = variable)) + 
  geom_pointrange(aes(ymin = `5 %`, ymax = `95 %`),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#a6761d","#e7298a","#b10026","#fc4e2a","#feb24c"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Odds ratio")+ggtitle("Odds ratio on predicting immune response (GSE91061)")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="right",
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
  geom_hline(yintercept=1, linetype="dashed", color = "black")+scale_y_continuous(breaks=c(0,1,2,4,6,8,10))+coord_flip()
dev.off()
