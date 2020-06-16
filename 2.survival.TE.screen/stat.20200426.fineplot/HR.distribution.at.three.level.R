########draw at subfamily
head(cb.data)
HRdata.subfamily=cb.data
pdf("2.survival.TE.screen/stat.20200426.fineplot/HR.subfamily.pdf",width = 8,height = 5)
ggboxplot(HRdata.subfamily, x = "repClass", y = "HR",size=0.5,
          color = "repClass", palette = c("DNA"="#1f78b4","LINE"="#d95f02","LTR"="#7570b3","Retroposon"="#e7298a","Satellite"="#66a61e","SINE"="#e6ab02"),
          #add = "jitter"
          )+
  #stat_compare_means()+
  geom_point(data=cb.data,aes(x=as.factor(repClass),
                              y=HR,colour=repClass),
             position=position_jitter(width=0.15),alpha=0.75,size=0.5)+
  theme(legend.position = "right",axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1))+
  ggtitle(paste0("HR distribution at subfamily"))+
  xlab("endpoint")+ylab("HR")+
  facet_grid(. ~ endpoint)+
  geom_hline(yintercept=1, linetype="dashed", color = "black")
dev.off()
####### draw family level
######
"DNA"="#1f78b4","LINE"="#d95f02","LTR"="#7570b3","Retroposon"="#e7298a","Satellite"="#66a61e","SINE"="#e6ab02"
#cbPalette= c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e", "#e6ab02")
head(cb.data.family)
pdf("2.survival.TE.screen/stat.20200426.fineplot/HR.family.pdf",width = 7,height = 4)
ggplot(cb.data.family, aes( repFamily,HR, colour = repClass)) + 
  geom_pointrange(aes(ymin = HRCILL, ymax = HRCIUL),fatten = 2)+
  facet_grid(. ~ endpoint)+
  scale_fill_manual(values=c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  scale_color_manual(values = c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Hazard ratio")+ggtitle("HR distribution significant at least in one endpoint")+
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6))+
    geom_hline(yintercept = 1, linetype="dashed", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
  )+coord_flip()
dev.off()
################
################draw at class level
#
head(cb.data.class)
pdf("2.survival.TE.screen/stat.20200426.fineplot/HR.class.pdf",width = 7,height = 4)
ggplot(cb.data.class, aes( repClass,HR, colour = repClass)) + 
  geom_pointrange(aes(ymin = HRCILL, ymax = HRCIUL),fatten = 2)+
  facet_grid(. ~ endpoint)+
  scale_fill_manual(values=c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  scale_color_manual(values = c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Hazard ratio")+ggtitle("HR distribution significant at least in one endpoint")+
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6))+
  geom_hline(yintercept = 1, linetype="dashed", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
  )+coord_flip()
dev.off()
############
save(cb.data, cb.data.class, cb.data.family, file="2.survival.TE.screen/stat.20200426.fineplot/HR.distribution.at.three.level.RData")
