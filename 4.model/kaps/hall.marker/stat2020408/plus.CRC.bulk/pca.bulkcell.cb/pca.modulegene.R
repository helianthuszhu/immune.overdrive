#######combine bulk CRC and cell line data
#######bulk
length(intersect(colnames(expM),rownames(kaps.td)))
kaps.td.exp=expM[, colnames(expM) %in% rownames(kaps.td)]
#######
dim(CB.cell.line.expMatrix.gse22250)
dim(CB.cell.line.expMatrix.gse5816)
#
dim((CB.cell.line.expMatrix.gse80137))
dim(aza.expmatrix.normd.symbol)####use this one
######
ccids=Reduce(intersect, list(rownames(kaps.td.exp),
                             rownames(CB.cell.line.expMatrix.gse5816),
                             rownames(aza.expmatrix.normd.symbol),
                             rownames(CB.cell.line.expMatrix.gse22250)
                             )
             )
length(ccids)
#####
library(dplyr)
#
aa1=as.data.frame(t(kaps.td.exp[ccids,]))
aa2=as.data.frame(t(CB.cell.line.expMatrix.gse5816[ccids,]))
aa3=as.data.frame(t(aza.expmatrix.normd.symbol[ccids,]))
aa4=as.data.frame(t(CB.cell.line.expMatrix.gse22250[ccids,]))
#change NA into median
aa11=aa1 %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
aa22=aa2 %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
aa33=aa3 %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
aa44=aa4 %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
class(aa11)
rownames(aa11)=rownames(aa1)
rownames(aa22)=rownames(aa2)
rownames(aa33)=rownames(aa3)
rownames(aa44)=rownames(aa4)
#
aa22[1:3,1:4]

bulkccle.cb.expM=cbind(kaps.td.exp[ccids,], CB.cell.line.expMatrix.gse5816[ccids,],
                       aza.expmatrix.normd.symbol[ccids,],CB.cell.line.expMatrix.gse22250[ccids,])

bulkccle.cb.expM=cbind(as.data.frame(t(aa11)),
                       as.data.frame(t(aa22)),
                       as.data.frame(t(aa33)),
                       as.data.frame(t(aa44)))
bulkccle.cb.expM=as.data.frame(bulkccle.cb.expM)
dim(bulkccle.cb.expM)
bulkccle.cb.expM[1:4,1:4]
#
bulkccle.cb.expM.tmp=na.omit(bulkccle.cb.expM)
dim(bulkccle.cb.expM.tmp)
######combine cell and samples infomation
anno1=data.frame(group=kaps.td$kaps.group)
rownames(anno1)=rownames(kaps.td)
anno1$batch=rep("bulk.CRC",times=nrow(anno1))
head(anno1)
#
ascore.five.CB[1:4,1:6]
anno2=data.frame(group=ascore.five.CB$treatment,
                 batch=ascore.five.CB$GEOid
                 )
rownames(anno2)=rownames(ascore.five.CB)
anno2=subset(anno2, batch=="GSE80137" | batch=="GSE5816" | batch=="GSE22250")
anno2$group=ifelse(anno2$group=="1Control" | anno2$group=="untreatment" | anno2$group== "WT", "untreated.cell","treated.cell")
head(anno2)
#
bulkccle.cb.annotation=rbind(anno1, anno2)
head(bulkccle.cb.annotation)
###############
###############
#PCA
#
library("FactoMineR")
library("factoextra")
#
bulkccle.pca <- PCA(as.matrix(t(bulkccle.cb.expM.tmp)), graph = FALSE)
#
pdf("4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/pca.bulkcell.cb/pca.before.batch.bulkcell.CB.pdf",width = 5,height = 4)
fviz_pca_ind(bulkccle.pca,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = bulkccle.cb.annotation$batch, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.bulkCRC.aza.cellline")+theme_bw()
fviz_pca_ind(bulkccle.pca,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = bulkccle.cb.annotation$group, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.bulkCRC.aza.cellline")+theme_bw()
dev.off()
###########
###########
#remove batch effect
#####
library(sva)
library(pamr)
library(limma)
bulkccle.batch = bulkccle.cb.annotation$batch
bulkccle.modcombat = model.matrix(~1, data=bulkccle.cb.annotation)
combat_edata.bulkccle.CB = ComBat(dat=as.matrix(bulkccle.cb.expM.tmp), batch=bulkccle.batch, mod=bulkccle.modcombat, par.prior=TRUE, prior.plots=FALSE)
dim(combat_edata.bulkccle.CB)
######
#Pca again
bulkccle.pca2 <- PCA(as.matrix(t(combat_edata.bulkccle.CB)), graph = FALSE)
#
pdf("4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/pca.bulkcell.cb/pca.before.batch.bulkcell.CB.batch.rm.pdf",width = 5,height = 4)
fviz_pca_ind(bulkccle.pca2,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = bulkccle.cb.annotation$batch, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.bulkCRC.aza.cellline")+theme_bw()
fviz_pca_ind(bulkccle.pca2,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = bulkccle.cb.annotation$group, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.bulkCRC.aza.cellline")+theme_bw()
dev.off()
##################
##################
##########use module genes 
##########
head(group_g.draw)
######
combat_edata.bulkccle.CB[1:4,1:4]
length(intersect(rownames(group_g.draw), rownames(combat_edata.bulkccle.CB)))
combat_edata.bulkccle.CB.modulegeneM=combat_edata.bulkccle.CB[rownames(combat_edata.bulkccle.CB) %in% rownames(group_g.draw), ]
#combat_edata.bulkccle.CB.modulegeneM=bulkccle.cb.expM.tmp[rownames(bulkccle.cb.expM.tmp) %in% rownames(group_g.draw), ]

dim(combat_edata.bulkccle.CB.modulegeneM)
#
bulkccle.pcamodulegene <- PCA(as.matrix(t(combat_edata.bulkccle.CB.modulegeneM)), graph = FALSE)
#
pdf("4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/pca.bulkcell.cb/pca.before.batch.bulkcell.CB.batch.rm.modulegene.pdf",width = 5,height = 4)
fviz_pca_ind(bulkccle.pcamodulegene,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = bulkccle.cb.annotation$batch, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.bulkCRC.aza.cellline")+theme_bw()
fviz_pca_ind(bulkccle.pcamodulegene,pointshape = 19,geom = "point",invisible="quali",
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = bulkccle.cb.annotation$group, # color by groups
             palette = c( "#E7B800", "#FC4E07","#00AFBB","#377eb8","#4daf4a","#984ea3","#ff7f00","#e41a1c"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",title="PCA.bulkCRC.aza.cellline")+theme_bw()
dev.off()
#####tSNE
tsne_bulkccle.cb <- Rtsne(t(combat_edata.bulkccle.CB.modulegeneM),initial_dims=40,theta=0.0,perplexity = 40)
tsne_output.bulkccle.cb <- data.frame(x=tsne_bulkccle.cb$Y[,1],y=tsne_bulkccle.cb$Y[,2],
                                      group = bulkccle.cb.annotation[colnames(combat_edata.bulkccle.CB.modulegeneM),])
write.csv(tsne_output.bulkccle.cb,"4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/pca.bulkcell.cb/tsne.modulegene.csv")
head(tsne_output.bulkccle.cb)
pdf("4.model/kaps/hall.marker/stat2020408/plus.CRC.bulk/pca.bulkcell.cb/tsne.modulegene.pdf",width = 5,height = 5)
ggplot(tsne_output.bulkccle.cb, aes(x, y, label=group.group)) + 
  #geom_point(aes(colour=factor(K8.detailed),shape=factor(K8.detailed)), size=2) +
  geom_point(aes(colour=factor(group.group)), size=1) +
  #scale_shape_manual(values=c(21,22,23,25,24))+
  #scale_shape_manual(values=c(10,4,8,3,6))+
  #scale_color_manual(values=wes_palette(n=5, name="GrandBudapest"))+
  #geom_point(aes(shape=factor(site)))+
  #geom_text(aes(colour=factor(control)), check_overlap = TRUE, size=2.2, hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.005) + 
  #scale_colour_discrete(name = "Patient.ID") +
  scale_color_manual(values=c("set4"="#e31a1c", "set3"="#e31a1c","set2"="#e31a1c","set1"="#377eb8", 
                              "treated.cell"="#4daf4a", "untreated.cell"="#984ea3"))+
  labs(x="tSNE1", y="tSNE2", title="Rtsne.CRC.cell.CB") + theme_bw()+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(0.3, "cm"))

ggplot(tsne_output.bulkccle.cb, aes(x, y, label=group.group)) + 
  #geom_point(aes(colour=factor(K8.detailed),shape=factor(K8.detailed)), size=2) +
  geom_point(aes(colour=factor(group.group)), size=1) +
  #scale_shape_manual(values=c(21,22,23,25,24))+
  #scale_shape_manual(values=c(10,4,8,3,6))+
  #scale_color_manual(values=wes_palette(n=5, name="GrandBudapest"))+
  #geom_point(aes(shape=factor(site)))+
  #geom_text(aes(colour=factor(control)), check_overlap = TRUE, size=2.2, hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.005) + 
  #scale_colour_discrete(name = "Patient.ID") +
  scale_color_manual(values=c("set4"="#e31a1c", "set3"="#984ea3","set2"="#4daf4a","set1"="#377eb8", 
                              "treated.cell"="#f781bf", "untreated.cell"="#ff7f00"))+
  labs(x="tSNE1", y="tSNE2", title="Rtsne.CRC.cell.CB") + theme_bw()+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(0.3, "cm"))


ggplot(tsne_output.bulkccle.cb) + geom_point(aes(x, y, colour = group.batch), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() 
dev.off()
