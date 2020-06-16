###################GEP
###################
head(stat.kirc.kaps.vali)
head(preexp.cal.kircfh)
#
length(intersect(rownames(stat.kirc.kaps.vali), rownames(preexp.cal.kircfh)))
#####
stat.kirc.kaps.vali.GEP=cbind(preexp.cal.kircfh[rownames(stat.kirc.kaps.vali),]$GEP,stat.kirc.kaps.vali)
colnames(stat.kirc.kaps.vali.GEP)[1]="GEP"
head(stat.kirc.kaps.vali.GEP)

#
stat.kirc.kaps.vali.GEP$kaps.group.kirc=factor(stat.kirc.kaps.vali.GEP$kaps.group.kirc,levels = c("set4","set3","set2","set1"))

my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
######
pdf("8.KIRC/2.GEP/GEP.kaps.kirc.own.cutoff.new.pdf",height  = 5,width = 4)
ggviolin(stat.kirc.kaps.vali.GEP,x = "kaps.group.kirc", y ="GEP" , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
         palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("kaps.group")+ylab("GEP") +
  ggtitle(paste0("GEP",".KIRC.cutoff"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
dev.off()
save(stat.kirc.kaps.vali,stat.kirc.kaps.vali.GEP,file="8.KIRC/2.GEP/GEP.kaps.kirc.own.cutoff.RData")
