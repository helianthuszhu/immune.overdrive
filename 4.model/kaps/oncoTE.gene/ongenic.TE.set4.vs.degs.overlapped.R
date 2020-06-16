
####
####
head(kaps.)
####
####TE oncogenuc gene list from NGs wang ting
oncotegene=read.table("4.model/kaps/oncoTE.gene/ongenic.TEgenelist.txt",header = T,sep = "\t")
head(oncotegene)
length(unique(oncotegene$Oncogene))
length(intersect(rownames(kaps.four.degs), unique(oncotegene$Oncogene)))
##filter in CRC
oncotegene=subset(oncotegene, COAD>0)
####
degs.oncotegene=kaps.four.degs[intersect(rownames(kaps.four.degs), unique(oncotegene$Oncogene)),]
head(degs.oncotegene)
degs.oncotegene=degs.oncotegene[order(degs.oncotegene$set4logFC,decreasing = T),]
degs.oncotegene$pos=seq(1:nrow(degs.oncotegene))
degs.oncotegene$sig.padj.0.05=ifelse(degs.oncotegene$set4adj.P.Val<0.05,"TRUE","FALSE")
degs.oncotegene$geneName=rownames(degs.oncotegene)
table(degs.oncotegene$sig.padj.0.05)
dim(degs.oncotegene)
#####
pdf("4.model/kaps/oncoTE.gene/ongenic.TE.set4.vs.allother.CRC.only.pdf",width = 7,height = 5)
ggplot(degs.oncotegene, aes(x=pos, y=set4logFC, color=sig.padj.0.05)) +
  scale_color_manual(values=c("#00AFBB","#d73027","#756bb1","#E69F00"))+
  geom_point(size=2) +theme_classic()+
  geom_text(data=subset(degs.oncotegene, sig.padj.0.05=="TRUE"),
            aes(label=geneName))
dev.off()
####
degs.oncotegene=deg.res.set34[intersect(rownames(deg.res.set34), unique(oncotegene$Oncogene)),]
  dim(degs.oncotegene)
  head(degs.oncotegene)
  degs.oncotegene=degs.oncotegene[order(degs.oncotegene$logFC,decreasing = T),]
  degs.oncotegene$pos=seq(1:nrow(degs.oncotegene))
  degs.oncotegene$sig.padj.0.05=ifelse(degs.oncotegene$adj.P.Val<0.05,"TRUE","FALSE")
  degs.oncotegene$geneName=rownames(degs.oncotegene)
  summary(degs.oncotegene$adj.P.Val)
  #write.csv(degs.oncotegene, "4.model/kaps/oncoTE.gene/logfc.set3vs4.ontoTEgene.csv")
  pdf("4.model/kaps/oncoTE.gene/ongenic.TE.set4.vs.set3.CRC.only.pdf",width = 7,height = 5)
  ggplot(degs.oncotegene, aes(x=pos, y=logFC, color=sig.padj.0.05)) +
    scale_color_manual(values=c("#00AFBB","#d73027","#756bb1","#E69F00"))+
    geom_point(size=2) +theme_classic()+
    geom_text(data=subset(degs.oncotegene, sig.padj.0.05=="TRUE"),
              aes(label=geneName))
  dev.off()
###########
  #######
#
