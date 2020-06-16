########response data of GSE91061 and GSE78220
####load gse91061 data
#
load("6.immnue.treat.data/repM.clin.GSE91061.RData")
dim(repM.GSE91061)
#####
response.data.gse91061=repM.GSE91061[, colnames(repM.GSE91061) %in% rownames(cndi.rep)]
response.data.gse91061$mean.exp=rowMeans(response.data.gse91061,na.rm = T)
response.data.gse91061$z.of.mean.exp=(response.data.gse91061$mean.exp - mean(response.data.gse91061$mean.exp))/sd(response.data.gse91061$mean.exp)
head(response.data.gse91061)
#
head(clin.GSE91061)
all(rownames(clin.GSE91061) %in% rownames(response.data.gse91061))
stat.reponse.data.gse91061=cbind(response.data.gse91061, clin.GSE91061[rownames(response.data.gse91061),])
colnames(stat.reponse.data.gse91061)=gsub("-","_",colnames(stat.reponse.data.gse91061))
stat.reponse.data.gse91061$response=ifelse(stat.reponse.data.gse91061$reponse2=="PD","PD","PRCR.SD")
head(stat.reponse.data.gse91061)
table(stat.reponse.data.gse91061$response)
#######
#######
pdf("6.immnue.treat.data/boxplot.TEscore.vs.response.GSE91061.pdf",width = 5)
ggboxplot(stat.reponse.data.gse91061, x = "response", y = "z.of.mean.exp",
          color = "response", palette = c("PD"="#d01c8b","PRCR.SD"="#4dac26")
          #,add = "jitter"
)+stat_compare_means(aes(group = response))+
  geom_point(data=stat.reponse.data.gse91061,aes(x=as.factor(response),y=z.of.mean.exp,colour=response),
             position=position_jitter(width=0.2),alpha=0.75,size=2)+
  theme_classic() + 
  labs(x="response",y="TE score",title=paste0("TE score"," comparison"))+
  theme(legend.position = "right",
        #legend.justification = c(0,1),
        #legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        aspect.ratio = 1,
        plot.margin = margin(t=5,r=40,b=5,l=40,unit="pt"),
        axis.title.y = element_text(size=8,colour="black"),
        axis.title.x = element_text(size=8,colour="black"),
        #axis.title.x = element_blank(),
        axis.text = element_text(size=8,colour="black",angle=0),
        axis.text.x = element_text(size=8,colour="black",angle=0,hjust=1,vjust=0.5),
        #axis.text.x = element_blank(),
        plot.title = element_text(size=8,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.1, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5)
  )
dev.off()
#####
stat.reponse.data.gse91061.draw=  stat.reponse.data.gse91061[,c(1:9,11,20)] %>% pivot_longer(cols=colnames(stat.reponse.data.gse91061)[c(1:9,11)],
                                                                                                   names_to= "TE",
                                                                                                   values_to = "expression")
stat.reponse.data.gse91061.draw=as.data.frame(stat.reponse.data.gse91061.draw)
head(stat.reponse.data.gse91061.draw)
pdf("6.immnue.treat.data/boxplot.TEscore.individualTEs.vs.response.GSE91061.pdf",width = 12,height = 5)
ggboxplot(stat.reponse.data.gse91061.draw, x = "response", y = "expression",
          color = "response", palette = c("PD"="#d01c8b","PRCR.SD"="#4dac26"),
          add = "jitter")+
  #stat_compare_means(aes(group = response))+
          facet_grid(. ~ TE)+stat_compare_means(label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()          
######################
######################
#############gse78220
###
load("6.immnue.treat.data/validata.gse78220.RData")
dim(validata.gse78220)
############
response.data.gse78220=validata.gse78220[, colnames(validata.gse78220) %in% rownames(cndi.rep)]
response.data.gse78220$mean.exp=rowMeans(response.data.gse78220,na.rm = T)
response.data.gse78220$z.of.mean.exp=(response.data.gse78220$mean.exp - mean(response.data.gse78220$mean.exp))/sd(response.data.gse78220$mean.exp)
head(response.data.gse78220)
#
stat.reponse.data.gse78220=cbind(response.data.gse78220, validata.gse78220[,1])
colnames(stat.reponse.data.gse78220)=gsub("-","_",colnames(stat.reponse.data.gse78220))
colnames(stat.reponse.data.gse78220)[12]="response"
head(stat.reponse.data.gse78220)
table(stat.reponse.data.gse78220$response)
############
pdf("6.immnue.treat.data/boxplot.TEscore.vs.response.GSE78220.pdf",width = 5)
ggboxplot(stat.reponse.data.gse78220, x = "response", y = "z.of.mean.exp",
          color = "response", palette = c("res"="#d01c8b","nonres"="#4dac26")
          #,add = "jitter"
)+stat_compare_means(aes(group = response))+
  geom_point(data=stat.reponse.data.gse78220,aes(x=as.factor(response),y=z.of.mean.exp,colour=response),
             position=position_jitter(width=0.2),alpha=0.75,size=2)+
  theme_classic() + 
  labs(x="response",y="TE score",title=paste0("TE score"," comparison"))+
  theme(legend.position = "right",
        #legend.justification = c(0,1),
        #legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        aspect.ratio = 1,
        plot.margin = margin(t=5,r=40,b=5,l=40,unit="pt"),
        axis.title.y = element_text(size=8,colour="black"),
        axis.title.x = element_text(size=8,colour="black"),
        #axis.title.x = element_blank(),
        axis.text = element_text(size=8,colour="black",angle=0),
        axis.text.x = element_text(size=8,colour="black",angle=0,hjust=1,vjust=0.5),
        #axis.text.x = element_blank(),
        plot.title = element_text(size=8,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.1, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5)
  )
dev.off()
#########
#########
stat.reponse.data.gse78220.draw=  stat.reponse.data.gse78220[,c(1:9,11,12)] %>% pivot_longer(cols=colnames(stat.reponse.data.gse78220)[c(1:9,11)],
                                                                                             names_to= "TE",
                                                                                             values_to = "expression")
stat.reponse.data.gse78220.draw=as.data.frame(stat.reponse.data.gse78220.draw)
head(stat.reponse.data.gse78220.draw)
pdf("6.immnue.treat.data/boxplot.TEscore.individualTEs.vs.response.GSE78220.pdf",width = 12,height = 5)
ggboxplot(stat.reponse.data.gse78220.draw, x = "response", y = "expression",
          color = "response", palette = c("res"="#d01c8b","nonres"="#4dac26"),
          add = "jitter")+
  #stat_compare_means(aes(group = response))+
  facet_grid(. ~ TE)+stat_compare_means(label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off() 
#
save(stat.reponse.data.gse78220,stat.reponse.data.gse91061, stat.reponse.data.gse78220.draw, stat.reponse.data.gse91061.draw,file="6.immnue.treat.data/stat.TE.score.vs.response.two.gsedatasets.RData")
