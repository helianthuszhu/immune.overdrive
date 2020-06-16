##########genes correlated with mean TE expression
#############
load("5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/TEmean.cor.with.gene.RData")
head(TE.mean.corM.with.genes)
TE.mean.corM.with.genes$id=rownames(TE.mean.corM.with.genes)
TE.mean.corM.with.genes=separate(TE.mean.corM.with.genes, id, into = c("TE","symbol"),sep = "&")
TE.mean.corM.with.genes=TE.mean.corM.with.genes[order(TE.mean.corM.with.genes$cor,decreasing = T),]
head(TE.mean.corM.with.genes)
dim(TE.mean.corM.with.genes)
dim(subset(TE.mean.corM.with.genes, cor >= 0.6))
TE.mean.corM.with.genes[grep("APOBEC",TE.mean.corM.with.genes$symbol),]
TE.mean.corM.with.genes$pos=rev(seq(1:nrow(TE.mean.corM.with.genes)))
#TE.mean.corM.with.genes$set= Hmisc::cut2(TE.mean.corM.with.genes$cor, c(-0.2,0.2,0.4),labels=c("set1","set2","set3","set4"))
TE.mean.corM.with.genes$set <- cut(TE.mean.corM.with.genes$cor, 
                                   breaks = c(min(TE.mean.corM.with.genes$cor),-0.2,0.3,0.4,max(TE.mean.corM.with.genes$cor)),
                                   labels = c("set4","set3","set2","set1"), include.lowest = TRUE)
head(TE.mean.corM.with.genes)
table(TE.mean.corM.with.genes$set)
write.csv(TE.mean.corM.with.genes, "5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/TE.mean.cor.with.gene.csv")
#
farb=c("#E69F00","gray","#31a354","#d73027")
gcf=ggplot(TE.mean.corM.with.genes, aes(x=cor, y=pos, color=set)) +
  scale_color_manual(values=farb)+
  geom_point(size=0.5)+geom_vline(xintercept = 0.3, color = "#31a354", size=0.2)+theme()+
  geom_vline(xintercept = 0.4, color = "#d73027", size=0.2)+
geom_vline(xintercept = -0.2, color = "#E69F00", size=0.2)+theme_classic()+
annotate("text", x=0.5, y=16900, label=paste0("n=","261"),size=4)+
  annotate("text", x=0.35, y=16000, label=paste0("n=","885"),size=4)+
  annotate("text", x=0.2, y=10000, label=paste0("n=","15804"),size=4)+
  annotate("text", x=-0.3, y=350, label=paste0("n=","406"),size=4)+
  xlab("correlation with TE mean exp")+ylab("genes order by correlation coefs")+ggtitle("correlation coefs of TE mean exp with genes")
  
generate.PDF <- function(fig) {
    pdf("5.TE.cor.with.gene/TE.mean.exp.cor.with.gene/TE.mean.exp.cor.gene.distribution.pdf",height  = 4,width = 5)
    print(gcf)
    dev.off()
  }
generate.PDF(fig)
  