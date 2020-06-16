############MCP counter
load("4.model/kaps/MCP.counter/socre.Mcpcounter.CRC.RData")
score=as.data.frame(t(score.mcp.CRC))
#
score=as.data.frame(scale(score))
head(score)
ids=intersect(rownames(kaps.td), rownames(score))

cbscoreshow=cbind(kaps.td$kaps.group, score[rownames(kaps.td),])
colnames(cbscoreshow)[1]="kaps.group"

head(cbscoreshow)
#cbscoreshow=subset(cbscoreshow, cluster=="cluster_1" | cluster=="cluster_2")
#######
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
index.sets=colnames(cbscoreshow)[-1]
for (i in 1:length(index.sets)) {
  pdf(file = paste0("4.model/kaps/MCP.counter/indi/",index.sets[i],".MCPcounter.four.pdf"),height  = 4,width = 3)
  print(ggviolin(cbscoreshow,x = "kaps.group", y =index.sets[i] , fill = "kaps.group",alpha = 1,size = 0.3,
                 palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
          stat_compare_means(comparisons =my_comparisons ,label = "p.signif")+theme_bw()+xlab("kaps.group")+ylab(index.sets[i]) +
          ggtitle(paste0(index.sets[i],"CRC.MCPcounter"))
  )
  dev.off()
}
#
library(tidyr)
drawdatashow <- cbscoreshow %>% pivot_longer(cols=colnames(cbscoreshow)[-1],
                                             names_to= "signature",
                                             values_to = "score")
drawdatashow=as.data.frame(drawdatashow)

#drawdatashow$signature <- factor(drawdatashow$signature, levels = rownames(sigshow))
#farb=c("#ca0020","#0571b0","#4daf4a","#27408B","#FF0000","#2E8B57","#CD00CD")
#head(drawdatashow)
#table(drawdatashow$signature)
farb=c("#00AFBB","#756bb1","#E69F00","#d73027")
p1=ggplot(drawdatashow, aes(x=kaps.group, y=score, group=kaps.group)) + 
  geom_boxplot(aes(fill=kaps.group),alpha = 0.8,outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.4)+
  stat_compare_means(label = "p.signif")+
  facet_grid(signature~. )+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+rotate()+xlab("kaps.group")+ylab("Z.score")


p2=ggplot(drawdatashow, aes(x=kaps.group, y=score, group=kaps.group)) + 
  geom_boxplot(aes(fill=kaps.group),outlier.colour = "black",outlier.size = 0.1,width=0.7,size=0.3)+
  stat_compare_means(label = "p.signif")+
  facet_grid(signature~. )+
  theme(strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= farb)+rotate()+xlab("kaps.group")+ylab("Z.score")

#+theme_classic()
#theme(strip.text.x = element_text(size = 8, colour = "orange", angle = 90))
generate.PDF <- function(fig) {
  pdf("4.model/kaps/MCP.counter/MCPcounter.10.celltypes.pdf",height  = 8)
  print(p1)
  print(p2)
  dev.off()
}
generate.PDF(fig)
