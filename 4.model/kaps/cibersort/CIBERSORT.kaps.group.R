#####cibersort
##########
#ciberdata=read.table("4.model/kaps/cibersort/CIBERSORT.Output_.no.missing.txt",header = T,row.names = 1,sep = "\t")
ciberdata=read.table("4.model/kaps/cibersort/CIBERSORT.Output_imputed.txt",header = T,row.names = 1,sep = "\t")
head(ciberdata)
summary(ciberdata$P.value)
#####
length(intersect(rownames(kaps.td), rownames(ciberdata)))
stat.ciber=cbind(kaps.td$kaps.group,ciberdata[rownames(kaps.td),])
colnames(stat.ciber)[1]="kaps.group"
head(stat.ciber)
colnames(stat.ciber)
#####
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
index.ciber=colnames(stat.ciber)[c(2:23)]
for (i in 1:length(index.ciber)) {
  pdf(file = paste0("4.model/kaps/cibersort/indi.imputed/",index.ciber[i],".cibersort.imputed.four.pdf"),height  = 4,width = 3)
  print(ggviolin(stat.ciber,x = "kaps.group", y =index.ciber[i] , fill = "kaps.group",alpha = 1,size = 0.3,
                 palette = c("#d73027","#E69F00","#00AFBB","#9970ab"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons =my_comparisons ,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab(index.ciber[i]) +
          ggtitle(paste0(index.ciber[i],"CRC.cibersort.no.filtering"))
  )
  dev.off()
}
#####draw together
#########
#######
library(tidyr)
drawdatashow.ciber <- stat.ciber[,c(1:23)] %>% pivot_longer(cols=colnames(stat.ciber)[c(2:23)],
                                             names_to= "signature",
                                             values_to = "score")
drawdatashow.ciber=as.data.frame(drawdatashow.ciber)

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