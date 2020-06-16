head(GEP.stats)
head(kaps.td)
GEP.stats.kaps=cbind(GEP.stats[rownames(kaps.td),]$GEP, kaps.td)
colnames(GEP.stats.kaps)[1]="GEP"
head(GEP.stats.kaps)
#
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
GEP.stats.kaps$kaps.group=factor(GEP.stats.kaps$kaps.group,levels = c("set4","set3","set2","set1"))
ge1=ggviolin(GEP.stats.kaps,x = "kaps.group", y ="GEP" , fill = "kaps.group",alpha = 1,size = 0.3,
             #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
             palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
             add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.group.kaps")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.group.kaps")+ylab("GEP")+ggtitle("GEP")

generate.PDF <- function(fig) {
  pdf("4.model/kaps/GEP/GEPvs.kaps.group.pdf",width = 4,height = 6)
  print(ge1)
  dev.off()
}
generate.PDF(fig)

pdf(file = paste0("4.model/kaps/GEP/GEPvs.kaps.group.new.pdf"),height  = 5,width = 4)
print(ggviolin(GEP.stats.kaps,x = "kaps.group", y ="GEP" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
               palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
               add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
        stat_compare_means(comparisons = my_comparisons,label = "p.signif")+xlab("TE.cluster")+ylab("GEP") +
        ggtitle(paste0("GEP"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
)
print(ggviolin(GEP.stats.kaps,x = "kaps.group", y ="GEP" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
               palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
               add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.new")+
        stat_compare_means(label = "p.signif",method = "kruskal.test")+xlab("TE.cluster")+ylab("GEP") +
        ggtitle(paste0("GEP"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
)
dev.off()
#######
save(GEP.stats.kaps,file="4.model/kaps/GEP/GEPvs.kaps.group.new.RData")
##########
##########
##########coreelation plot
####GEP
pdf("4.model/kaps/GEP/cor.GEP.vs.TEscore.pdf",width = 5,height = 5)
ggscatter(GEP.stats.kaps, x = "z.of.mean.exp", y = "GEP", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "TE score", ylab = "GEP")+
  ggtitle("TE score vs GEP in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
####CD8A
pdf("4.model/kaps/GEP/cor.CD8A.vs.TEscore.pdf",width = 5,height = 5)
ggscatter(GEP.stats.kaps, x = "z.of.mean.exp", y = "CD8A", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "TE score", ylab = "CD8A")+
  ggtitle("TE score vs CD8A in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
