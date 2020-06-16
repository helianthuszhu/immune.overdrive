###########odd ratio
###########
dim(clin.heat)
#######
head(GEP.stats.kaps)
###########
oddratiodata=cbind(GEP=GEP.stats.kaps$GEP,clin.heat[rownames(GEP.stats.kaps),])
oddratiodata$GEP.group=ifelse(oddratiodata$GEP> median(oddratiodata$GEP), 1,0)
oddratiodata$msi.bin.bin=ifelse(oddratiodata$MSI.status.bin=="MSI","2MSI","1MSS")
#
#oddratiodata$AJCC.stage.new=gsub("not_available","NA",oddratiodata$AJCC.stage)
oddratiodata= oddratiodata %>%  mutate_at(vars(AJCC.stage), na_if, "not_available")
#oddratiodata$AJCC.stage.new=substr(oddratiodata$AJCC.stage.new, 6, nchar(as.character(oddratiodata$AJCC.stage.new)))

############
############
oddsmodel <- glm( GEP.group ~ kaps.group +msi.bin.bin+AJCC.stage+age.bi+gender, data = oddratiodata, family = binomial)
oddres1=as.data.frame(summary(oddsmodel)$coef)
oddres1$odd.ratio=exp(oddres1$Estimate)
oddres2=as.data.frame(exp(confint(oddsmodel,level = 0.9)))
oddres.ok=cbind(oddres1,oddres2)
oddres.ok=oddres.ok[-1,]
oddres.ok$variable=rownames(oddres.ok)
###
oddres.ok$variable.name=c("TE.cluster2","TE.cluster3","TE.cluster4","MSI.status (MSI)","stageII","stageIII","stageIV","age(high)","gender(male)")
oddres.ok$variable.name=factor(oddres.ok$variable.name,levels = rev(c("TE.cluster2","TE.cluster3","TE.cluster4","MSI.status (MSI)","stageII","stageIII","stageIV","age(high)","gender(male)")))
oddres.ok.msi=oddres.ok
####draw
pdf("4.model/kaps/odds.ratio/oddsratio.GEP.against.TEscore.MSI.pdf",width = 6,height = 3)
ggplot(oddres.ok.msi, aes( variable.name,odd.ratio, colour = variable.name)) + 
  geom_pointrange(aes(ymin = `5 %`, ymax = `95 %`),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#a6761d","#e7298a","#b10026","#fc4e2a","#feb24c"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Odds ratio")+ggtitle("Odds ratio on predicting immune infiltration (GEP)")+
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
#
############
oddsmodel <- glm( GEP.group ~ z.of.mean.exp + TMB+AJCC.stage +age.bi+gender, data = oddratiodata, family = binomial)
oddres1=as.data.frame(summary(oddsmodel)$coef)
oddres1$odd.ratio=exp(oddres1$Estimate)
oddres2=as.data.frame(exp(confint(oddsmodel,level = 0.9)))
oddres.ok=cbind(oddres1,oddres2)
oddres.ok=oddres.ok[-1,]
oddres.ok$variable=rownames(oddres.ok)
oddres.ok$variable.name=c("TE.score","TMB","stageII","stageIII","stageIV","age(high)","gender(male)")
oddres.ok$variable.name=factor(oddres.ok$variable.name,levels = rev(c("TE.score","TMB","stageII","stageIII","stageIV","age(high)","gender(male)")))
oddres.ok.TMB=oddres.ok
############
#cbPalette=c("#1f78b4" "#5ab4ac" "#c51b7d" "#d95f02" "#7570b3" "#e7298a" "#e6ab02" "#66a61e")
#"#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#e7298a","#a65628"
pdf("4.model/kaps/odds.ratio/oddsratio.GEP.against.TEscore.TMB.pdf",width = 6,height = 3)
ggplot(oddres.ok.TMB, aes( variable.name,odd.ratio, colour = variable.name)) + 
  geom_pointrange(aes(ymin = `5 %`, ymax = `95 %`),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#e7298a","#b10026","#a65628"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Odds ratio")+ggtitle("Odds ratio on predicting immune infiltration (GEP)")+
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
  geom_hline(yintercept=1, linetype="dashed", color = "black")+scale_y_continuous(breaks=c(0,1,2,3,4))+coord_flip()
dev.off()
##########
save(oddres.ok.TMB,oddres.ok.msi,file="4.model/kaps/odds.ratio/oddsratio.GEP.against.TEscore.RData")
###############################
#########################
###################
###########predict CD8A
###########
oddratiodata$CD8A=as.numeric(paste(oddratiodata$CD8A))
oddratiodata$CD8A.group=ifelse(oddratiodata$CD8A> median(oddratiodata$CD8A),1,0)
###########
oddsmodel.CD8A <- glm( CD8A.group ~ kaps.group +msi.bin.bin+AJCC.stage+age.bi+gender, data = oddratiodata, family = binomial)
oddres1CD8A=as.data.frame(summary(oddsmodel.CD8A)$coef)
oddres1CD8A$odd.ratio=exp(oddres1CD8A$Estimate)
oddres2CD8A=as.data.frame(exp(confint(oddsmodel.CD8A,level = 0.9)))
oddres.ok.CD8A=cbind(oddres1CD8A,oddres2CD8A)
oddres.ok.CD8A=oddres.ok.CD8A[-1,]
oddres.ok.CD8A$variable=rownames(oddres.ok.CD8A)
###
oddres.ok.CD8A$variable.name=c("TE.cluster2","TE.cluster3","TE.cluster4","MSI.status (MSI)","stageII","stageIII","stageIV","age(high)","gender(male)")
oddres.ok.CD8A$variable.name=factor(oddres.ok.CD8A$variable.name,levels = rev(c("TE.cluster2","TE.cluster3","TE.cluster4","MSI.status (MSI)","stageII","stageIII","stageIV","age(high)","gender(male)")))
oddres.ok.CD8A.msi=oddres.ok.CD8A
####draw
pdf("4.model/kaps/odds.ratio/oddsratio.CD8A.against.TEscore.MSI.pdf",width = 6,height = 3)
ggplot(oddres.ok.CD8A.msi, aes( variable.name,odd.ratio, colour = variable.name)) + 
  geom_pointrange(aes(ymin = `5 %`, ymax = `95 %`),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#a6761d","#e7298a","#b10026","#fc4e2a","#feb24c"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Odds ratio")+ggtitle("Odds ratio on predicting immune infiltration (CD8A)")+
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
  geom_hline(yintercept=1, linetype="dashed", color = "black")+scale_y_continuous(breaks=c(0,1,2,4,6,8,10,12))+coord_flip()
dev.off()
##############tmb
oddsmodel.CD8A <- glm( CD8A.group ~ z.of.mean.exp +TMB+AJCC.stage+age.bi+gender, data = oddratiodata, family = binomial)
oddres1CD8A=as.data.frame(summary(oddsmodel.CD8A)$coef)
oddres1CD8A$odd.ratio=exp(oddres1CD8A$Estimate)
oddres2CD8A=as.data.frame(exp(confint(oddsmodel.CD8A,level = 0.9)))
oddres.ok.CD8A=cbind(oddres1CD8A,oddres2CD8A)
oddres.ok.CD8A=oddres.ok.CD8A[-1,]
oddres.ok.CD8A$variable=rownames(oddres.ok.CD8A)
###
oddres.ok.CD8A$variable.name=c("TE.score","TMB","stageII","stageIII","stageIV","age(high)","gender(male)")
oddres.ok.CD8A$variable.name=factor(oddres.ok.CD8A$variable.name,levels = rev(c("TE.score","TMB","stageII","stageIII","stageIV","age(high)","gender(male)")))
oddres.ok.CD8A.TMB=oddres.ok.CD8A
####draw
pdf("4.model/kaps/odds.ratio/oddsratio.CD8A.against.TEscore.TMB.pdf",width = 6,height = 3)
ggplot(oddres.ok.CD8A.TMB, aes( variable.name,odd.ratio, colour = variable.name)) + 
  geom_pointrange(aes(ymin = `5 %`, ymax = `95 %`),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#a6761d","#e7298a","#b10026","#fc4e2a","#feb24c"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Odds ratio")+ggtitle("Odds ratio on predicting immune infiltration (CD8A)")+
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
  geom_hline(yintercept=1, linetype="dashed", color = "black")+scale_y_continuous(breaks=c(0,1,2,2.5))+coord_flip()
dev.off()
#
save(oddres.ok.CD8A.TMB,oddres.ok.CD8A.msi,file="4.model/kaps/odds.ratio/oddsratio.CD8A.against.TEscore.RData")
