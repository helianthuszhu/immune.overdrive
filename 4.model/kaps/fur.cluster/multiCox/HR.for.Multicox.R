#######draw multiCox for TE group and MSI clinical vraibles
#######PFI
aa=read.csv("4.model/kaps/fur.cluster/multiCox/MulCox.PFI.res.csv", row.names = 1, header = T)
aa$variable.name=rownames(aa)
aa=aa %>% separate(col="CI95.PFI",into = c("lower","upper"),sep = "[-]")
aa$lower=as.numeric(paste(aa$lower))
aa$upper=as.numeric(paste(aa$upper))
aa$variable.name=c("TE.cluster1 (ref=cluster3)","TE.cluster2 (ref=cluster3)","TE.cluster4 (ref=cluster4)",
                   "gender(male vs female)","age(high vs low)","MSI.status (MSS vs MSI)","lymphhatic.invasion (yes vs no)",
                   "location (proximal vs distal)","AJCC stage (III&IV vs I&II)")
aa$variable.name=factor(aa$variable.name,levels = rev(c("TE.cluster1 (ref=cluster3)","TE.cluster2 (ref=cluster3)","TE.cluster4 (ref=cluster4)",
                                                        "gender(male vs female)","age(high vs low)","MSI.status (MSS vs MSI)","lymphhatic.invasion (yes vs no)",
                                                        "location (proximal vs distal)","AJCC stage (III&IV vs I&II)")))
aa
#######
pdf("4.model/kaps/fur.cluster/multiCox/HR.PFI.pdf",width = 8,height = 3)
ggplot(aa, aes( variable.name,Hazard_Ratio.PFI, colour = variable.name)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#a6761d","#e7298a","#b10026","#fc4e2a","#feb24c"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Hazard ratio")+ggtitle("Hazard ratio of Multicox for PFI")+
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
  geom_hline(yintercept=1, linetype="dashed", color = "black")+#scale_y_continuous(breaks=c(0,1,2,3,5,7))+
  coord_flip()+
  scale_y_continuous(trans='log2',breaks = c(0,1,2,3,5,7))+
  annotate("text", x=9, y=8, label=paste0("ns"),size=5)+
  annotate("text", x=8, y=8, label=paste0("ns"),size=5)+
  annotate("text", x=7, y=8, label=paste0("*"),size=5)+
  annotate("text", x=6, y=8, label=paste0("ns"),size=5)+
  annotate("text", x=5, y=8, label=paste0("ns"),size=5)+
  annotate("text", x=4, y=8, label=paste0("ns"),size=5)+
  annotate("text", x=3, y=8, label=paste0("*"),size=5)+
  annotate("text", x=2, y=8, label=paste0("ns"),size=5)+
  annotate("text", x=1, y=8, label=paste0("****"),size=5)
dev.off()
  #ns: p > 0.05
  #*: p <= 0.05
  #**: p <= 0.01
  #***: p <= 0.001
  #****: p <= 0.0001
###OS
aa=read.csv("4.model/kaps/fur.cluster/multiCox/MulCox.OS.res.csv", row.names = 1, header = T)
aa$variable.name=rownames(aa)
aa=aa %>% separate(col="CI95.OS",into = c("lower","upper"),sep = "[-]")
aa$lower=as.numeric(paste(aa$lower))
aa$upper=as.numeric(paste(aa$upper))
aa$variable.name=c("TE.cluster1 (ref=cluster3)","TE.cluster2 (ref=cluster3)","TE.cluster4 (ref=cluster4)",
                   "gender(male vs female)","age(high vs low)","MSI.status (MSS vs MSI)","lymphhatic.invasion (yes vs no)",
                   "location (proximal vs distal)","AJCC stage (III&IV vs I&II)")
aa$variable.name=factor(aa$variable.name,levels = rev(c("TE.cluster1 (ref=cluster3)","TE.cluster2 (ref=cluster3)","TE.cluster4 (ref=cluster4)",
                                                        "gender(male vs female)","age(high vs low)","MSI.status (MSS vs MSI)","lymphhatic.invasion (yes vs no)",
                                                        "location (proximal vs distal)","AJCC stage (III&IV vs I&II)")))
aa
#######
pdf("4.model/kaps/fur.cluster/multiCox/HR.OS.pdf",width = 8,height = 3)
ggplot(aa, aes( variable.name,Hazard_Ratio.OS, colour = variable.name)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#a6761d","#e7298a","#b10026","#fc4e2a","#feb24c"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Hazard ratio")+ggtitle("Hazard ratio of Multicox for OS")+
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
  geom_hline(yintercept=1, linetype="dashed", color = "black")+#scale_y_continuous(breaks=c(0,1,2,3,5,7,10,14))+
  coord_flip()+
  scale_y_continuous(trans='log2',breaks = c(0,1,2,3,5,7,10,14))+
  annotate("text", x=9, y=15, label=paste0("ns"),size=5)+
  annotate("text", x=8, y=15, label=paste0("*"),size=5)+
  annotate("text", x=7, y=15, label=paste0("*"),size=5)+
  annotate("text", x=6, y=15, label=paste0("ns"),size=5)+
  annotate("text", x=5, y=15, label=paste0("***"),size=5)+
  annotate("text", x=4, y=15, label=paste0("ns"),size=5)+
  annotate("text", x=3, y=15, label=paste0("*"),size=5)+
  annotate("text", x=2, y=15, label=paste0("ns"),size=5)+
  annotate("text", x=1, y=15, label=paste0("****"),size=5)
dev.off()
summary(aa$upper)
#DSS
aa=read.csv("4.model/kaps/fur.cluster/multiCox/MulCox.DSS.res.csv", row.names = 1, header = T)
aa$variable.name=rownames(aa)
aa=aa %>% separate(col="CI95.DSS",into = c("lower","upper"),sep = "[-]")
aa$lower=as.numeric(paste(aa$lower))
aa$upper=as.numeric(paste(aa$upper))
aa$variable.name=c("TE.cluster1 (ref=cluster3)","TE.cluster2 (ref=cluster3)","TE.cluster4 (ref=cluster4)",
                   "gender(male vs female)","age(high vs low)","MSI.status (MSS vs MSI)","lymphhatic.invasion (yes vs no)",
                   "location (proximal vs distal)","AJCC stage (III&IV vs I&II)")
aa$variable.name=factor(aa$variable.name,levels = rev(c("TE.cluster1 (ref=cluster3)","TE.cluster2 (ref=cluster3)","TE.cluster4 (ref=cluster4)",
                                                        "gender(male vs female)","age(high vs low)","MSI.status (MSS vs MSI)","lymphhatic.invasion (yes vs no)",
                                                        "location (proximal vs distal)","AJCC stage (III&IV vs I&II)")))
aa

pdf("4.model/kaps/fur.cluster/multiCox/HR.DSS.pdf",width = 8,height = 3)
ggplot(aa, aes( variable.name,Hazard_Ratio.DSS, colour = variable.name)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper),fatten = 2)+
  #facet_grid(. ~ endpoint)+
  #scale_fill_manual(values=cbPalette)+ # Boxplot fill color
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#a6761d","#e7298a","#b10026","#fc4e2a","#feb24c"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Hazard ratio")+ggtitle("Hazard ratio of Multicox for DSS")+
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
  geom_hline(yintercept=1, linetype="dashed", color = "black")+#scale_y_continuous(breaks=c(0,1,2,3,5,7,10))+
  coord_flip()+
  #coord_cartesian(ylim = c(0, 10))+
  scale_y_continuous(trans='log2',breaks = c(0,1,2,4,5,7,9,16,64))+
  annotate("text", x=9, y=75, label=paste0("ns"),size=5)+
  annotate("text", x=8, y=75, label=paste0("*"),size=5)+
  annotate("text", x=7, y=75, label=paste0("*"),size=5)+
  annotate("text", x=6, y=75, label=paste0("ns"),size=5)+
  annotate("text", x=5, y=75, label=paste0("*"),size=5)+
  annotate("text", x=4, y=75, label=paste0("ns"),size=5)+
  annotate("text", x=3, y=75, label=paste0("*"),size=5)+
  annotate("text", x=2, y=75, label=paste0("ns"),size=5)+
  annotate("text", x=1, y=75, label=paste0("****"),size=5)
dev.off()
