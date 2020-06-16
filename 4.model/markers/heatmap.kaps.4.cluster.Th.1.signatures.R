#save(msi.com, immgM,file = "4.model/markers/marker.expM.msi.TE.cluster.RData")
#load("4.model/markers/marker.expM.msi.TE.cluster.RData")
cor.immune.drive=cbind(msi.com, immgM[rownames(msi.com),])
head(cor.immune.drive)
tmp.idmsi=intersect(rownames(surstats), rownames(cor.immune.drive))
cor.immune.drive.sur=cbind(cor.immune.drive[tmp.idmsi,], surstats[tmp.idmsi,7:14])
head(cor.immune.drive.sur)
#####
load("4.model/kaps/fur.cluster/kaps.td.RData")
############################
#genestats=cbind(cor.immune.drive[,1:18], as.data.frame(t(expM))[rownames(cor.immune.drive),])
genestats=cbind(kaps.td[,c(1:18,30:40)], as.data.frame(t(expM))[rownames(kaps.td),])
#genestats=genestats[order(genestats$z.of.mean.exp,decreasing = T),]
genestats[1:4,1:40]
head(kaps.td)
table(kaps.td$kaps.group)
genestats$MSI.status.bin=ifelse(genestats$MSI.status.bin=="MSI","2MSI","1MSS")
#genestats.orderd =genestats[order(genestats$kaps.group,genestats$MSI.status.bin,decreasing = T),]
genestats.orderd =genestats[order(genestats$z.of.mean.exp,decreasing = T),]
library(plyr)
library("data.table")
#score_dt <- setDT(genestats)
#genestats= as.data.frame(setkey(score_dt,z.of.mean.exp,MSI.status.bin))
####input markers panel
markerpanel=read.table("4.model/markers/immue.gene.to.show.for.heatmap.txt",header = T,sep = "\t")
head(markerpanel)
####
idx.panel=unique(markerpanel$category)
length(subset(markerpanel, category=="NanoString.panel")$symbol)
#
panelM=genestats.orderd[, colnames(genestats.orderd) %in% subset(markerpanel, category=="Th.1.signatures")$symbol]
#panelM=genestats.orderd[, colnames(genestats.orderd) %in% subset(markerpanel, category=="CD.8.T.exhuasted")$symbol]
colnames(panelM)
dim(panelM)
#
col_ha.top=genestats.orderd[rownames(panelM),c(1,13,14,17,27)]
class(col_ha.top)
head(col_ha.top)
#
pk1=pheatmap(t(panelM),show_rownames = T,border_color = NA,main = "Th-1.signatures",
         show_colnames = F,fontsize = 4,#na_col = "black",
         #cluster_rows = hc,
         cluster_cols =F,
         annotation_col = col_ha.top,
         #annotation_row = cndi.rep[rownames(ppd),c(2,3)],
         scale="row",
         #color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8), 
         color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
         #annotation_legend = T,
         annotation_colors  = list( TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                    TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                    #MSI.status=c("MSI"="#a65628", "MSI"="#ffd92f","MSS"="#8da0cb"),
                                    MSI.status.bin=c("2MSI"="#e7298a","1MSS"="#66a61e"),
                                    kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")#,
                                    #repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                    #repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
         )
)
save_pheatmap_pdf(pk1, "4.model/markers/heatmap.kaps.4.cluster.Th.1.signatures.new.pdf",width = 7,height = 3)


col_ha.top1=columnAnnotation(z.of.mean.exp.score=anno_lines(col_ha.top$z.of.mean.exp),
                             TE.cluster = col_ha.top$TE.cluster,
                             TE.cluster.agg=col_ha.top$TE.cluster.agg,
                             MSI.status.bin=col_ha.top$MSI.status.bin,
                             z.of.mean.exp=col_ha.top$z.of.mean.exp,
                             kaps.group=col_ha.top$kaps.group,
                             col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                      TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                      MSI.status.bin=c("2MSI"="#e7298a","1MSS"="#66a61e"),
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                                      ),show_annotation_name = TRUE)

comp=Heatmap(t(scale(panelM)), name = "Th-1.signatures", 
        col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
        #col = rev(viridis(10)),border = F,
        show_column_names = F,show_row_names = T,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha.top1#, 
        #right_annotation = row_ha.right,
        #show_row_dend = T,show_column_dend = T,
        #row_names_side = "left",
        #left_annotation = row_ha.left
)
generate.PDF <- function(fig) {
  pdf("4.model/markers/complexheatmap.kaps.4.cluster.Th.1.signatures.new.pdf",height = 5,width = 10)
  print(comp)
  dev.off()
}
generate.PDF(fig)

save.image(file = "tmp.RData")
####################

#my_comparisons <- list( c("set1", "set2"), c("set3", "cluster_3"), c("cluster_2", "cluster_3") )
gidx=unique(subset(markerpanel, !(category=="NanoString.panel"))$symbol)

marid=cbind(col_ha.top, genestats[rownames(col_ha.top), colnames(genestats) %in% gidx])
marid$kaps.group.agg=gsub("set1","set_1&2",msisubsur$kaps.group)
marid$kaps.group.agg=gsub("set2","set_1&2",msisubsur$kaps.group.agg)

head(marid)
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
my_comparisons <- list( c("set_1&2", "set3"), c("set_1&2", "set4"),
                        c("set3", "set4"))
for (i in 6:ncol(marid)-1) {
  pdf(paste0("4.model/markers/indiMarker.plot.agg/",colnames(marid)[i],".kaps.group.agg.pdf"))
  cateN=markerpanel[which(markerpanel$symbol==colnames(marid)[i]),]$category[1]
  print(ggviolin(marid,x = "kaps.group.agg", y = colnames(marid)[i], fill = "kaps.group.agg",alpha = 1,size = 0.3,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group.agg")+
          stat_compare_means(comparisons =my_comparisons, label = "p.signif")+
          theme_bw()+xlab("kaps.group.agg")+ylab(colnames(marid)[i])+ggtitle(paste0(colnames(marid)[i],"_",cateN))
        )
  dev.off()
}

###################
############
################calculate the mean exp of Th-1 signautre
##################
colnames(panelM)

head(panelM)
Th1.sigM=panelM
Th1.sigM$mean.Th1.sig=rowSums(Th1.sigM,na.rm = T)/ncol(Th1.sigM)
dim(Th1.sigM)
stat.th1sig=cbind(Th1.sigM[rownames(kaps.td),]$mean.Th1.sig,kaps.td)

stat.th1sig$z.of.mean.Th1.sig=(stat.th1sig$mean.Th1.sig-mean(stat.th1sig$mean.Th1.sig))/sd(stat.th1sig$mean.Th1.sig)
head(stat.th1sig)
class(stat.th1sig$z.of.mean.Th1.sig)

colnames(stat.th1sig)[1]="mean.Th1.sig"
#####
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
stat.th1sig$kaps.group=factor(stat.th1sig$kaps.group,levels = c("set4", "set3","set2","set1"))
th1g=ggviolin(stat.th1sig,x = "kaps.group", y ="mean.Th1.sig" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("mean.Th1.sig")+ggtitle("mean.Th1.sig")

th2g=ggviolin(stat.th1sig,x = "kaps.group", y ="z.of.mean.Th1.sig" , fill = "kaps.group",alpha = 1,size = 0.3,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#E69F00","#00AFBB","#756bb1"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("z.of.mean.Th1.sig")+ggtitle("z.of.mean.Th1.sig")

th2g=ggviolin(stat.th1sig,x = "kaps.group", y ="z.of.mean.Th1.sig" , fill = "kaps.group",alpha = 1,size = 0.3,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("z.of.mean.Th1.sig")+ggtitle("z.of.mean.Th1.sig")


th3g=ggviolin(stat.th1sig,x = "kaps.group", y ="z.of.mean.Th1.sig" , fill = "kaps.group",alpha = 1,size = 0.01,width = 1,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
              add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+xlab("kaps.group")+ylab("z.of.mean.Th1.sig")+ggtitle("z.of.mean.Th1.sig")

generate.PDF <- function(fig) {
  pdf("4.model/markers/mean.Th1.sig.pdf",width = 4,height = 6)
  print(th1g)
  print(th2g)
  print(th3g)
  dev.off()
}
generate.PDF(fig)

#########################################################################
####################################################################
######immune regulators
########
regulator=read.table("4.model/markers/regulators/immuregulators.summary.for.kaps.txt",header = T,sep = "\t")
rownames(regulator)=regulator$Symbol
regulator=regulator[order(regulator$category,decreasing = T),]
dim(regulator)
head(regulator)
#
#write.csv(as.data.frame(t(genestats.orderd[1:4,])),"4.model/markers/regulators/tmp.geneName.mapping.csv")
idre=intersect(rownames(regulator), colnames(genestats.orderd))
panelM.regu=genestats.orderd[, idre]
#panelM=genestats.orderd[, colnames(genestats.orderd) %in% subset(markerpanel, category=="CD.8.T.exhuasted")$symbol]
dim(panelM.regu)
panelM.regu[1:4,1:4]
########count NA in these genes
############
count.regu=as.data.frame(colSums(is.na(panelM.regu)))
colnames(count.regu)="count"
count.regu.sel=subset(count.regu,count<300)
head(count.regu.sel)
#####
panelM.regu.sel=as.data.frame(t(panelM.regu))[rownames(count.regu.sel),]
panelM.regu.sel[1:4,1:4]

#
col_ha.top=genestats.orderd[rownames(panelM.regu),c(1,13,14,17,27)]
class(col_ha.top)
head(col_ha.top)
ann.row=regulator[colnames(panelM.regu),]
ann.row$category=gsub("immunoinhibitor","Immunoinhibitor",ann.row$category)
ann.row$category=gsub("immunostimulator","Immunostimulator",ann.row$category)
ann.row.sel=ann.row[rownames(count.regu.sel),]
unique(ann.row.sel$category)
head(ann.row.sel)
#
col_ha.top1=columnAnnotation(z.of.mean.exp.score=anno_lines(col_ha.top$z.of.mean.exp),
                             TE.cluster = col_ha.top$TE.cluster,
                             TE.cluster.agg=col_ha.top$TE.cluster.agg,
                             MSI.status.bin=col_ha.top$MSI.status.bin,
                             z.of.mean.exp=col_ha.top$z.of.mean.exp,
                             kaps.group=col_ha.top$kaps.group,
                             col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                      TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                      MSI.status.bin=c("2MSI"="#e7298a","1MSS"="#66a61e"),
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)
row_ha_right=rowAnnotation(category=ann.row.sel$category,
                           col=list(category=c("receptor"="#a6cee3",
                             "Immunostimulator"="#1f78b4",
                             "Immunoinhibitor"="#b2df8a",
                             "chemokine"="#33a02c",
                             "cancer.germline.antigen"="#fb9a99",
                             "MHC.class.II"="#e31a1c",
                             "MHC.class.I"="#fdbf6f",
                             "MHC.non.class.II"="#ff7f00",
                             "MHC.non.class.I"="#cab2d6"
                             ))
                           )

comp1=Heatmap(t(scale(as.data.frame(t(panelM.regu.sel)))), name = "Th-1.signatures", na_col = "black",
             col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             #col = rev(viridis(10)),border = F,
             show_column_names = F,show_row_names = T,
             cluster_columns = F,
             cluster_rows = F,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 8),
             top_annotation = col_ha.top1, 
             right_annotation = row_ha_right#,
             #show_row_dend = T,show_column_dend = T,
             #row_names_side = "left",
             #left_annotation = row_ha.left
)
######HLA
mhcgene=regulator[grep("MHC",regulator$category),]
#length(intersect(rownames(mhcgene), rownames(panelM.regu.sel)))
mchexp=panelM.regu.sel[rownames(mhcgene),]
mchexp[1:4,1:4]
row_ha_left.mhc=rowAnnotation(category=mhcgene$category,
                           col=list(category=c(
                                               "MHC.class.II"="#7fc97f",
                                               "MHC.class.I"="#beaed4",
                                               "MHC.non.class.II"="#fdc086",
                                               "MHC.non.class.I"="#ffff99"
                           ))
)
#
comp2=Heatmap(t(scale(t(mchexp))), name = "MHC class gene", na_col = "gray",
        col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
        #col = rev(viridis(10)),border = F,
        show_column_names = F,show_row_names = T,
        cluster_columns = F,
        cluster_rows = T,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha.top1, 
        #right_annotation = row_ha_right#,
        #show_row_dend = T,show_column_dend = T,
        #row_names_side = "left",
        left_annotation = row_ha_left.mhc
)
#
generate.PDF <- function(fig) {
  pdf("4.model/markers/regulators/heatmp.2s.regulator.MHC.pdf",height = 5,width = 10)
  print(comp1)
  print(comp2)
  dev.off()
}
generate.PDF(fig)
######################
########find the regulators overlapped with top genes correlated with TE mean expression
###########
length(intersect(rownames(regulator), subset(TE.mean.corM.with.genes,set=="set1" |set=="set2" )$symbol))

posiggene=regulator[intersect(rownames(regulator), subset(TE.mean.corM.with.genes,set=="set1" |set=="set2" )$symbol),]
posiggene$category=gsub("immunoinhibitor","Immunoinhibitor",posiggene$category)
posiggene$category=gsub("immunostimulator","Immunostimulator",posiggene$category)
#
dim(posiggene)
length(intersect(rownames(posiggene), rownames(panelM.regu.sel)))
posigexp=panelM.regu.sel[rownames(posiggene),]
posigexp[1:4,1:4]
posigexp=as.data.frame(t(posigexp))

row_ha_left.possig=rowAnnotation(category=posiggene$category,
                              col=list(category=c("receptor"="#a6cee3",
                                                  "Immunostimulator"="#1f78b4",
                                                  "Immunoinhibitor"="#b2df8a",
                                                  "chemokine"="#33a02c",
                                                  "cancer.germline.antigen"="#fb9a99",
                                                  "MHC.class.II"="#e31a1c",
                                                  "MHC.class.I"="#fdbf6f",
                                                  "MHC.non.class.II"="#ff7f00",
                                                  "MHC.non.class.I"="#cab2d6"
                              ))
)
#
comp3=Heatmap(t(scale(posigexp)), name = paste0("positively correlated with","\n","mean.TE.exp (Rho > 0.3)" ), na_col = "gray",
              col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
              #col = rev(viridis(10)),border = F,
              show_column_names = F,show_row_names = T,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 8),
              top_annotation = col_ha.top1, 
              #right_annotation = row_ha_right#,
              #show_row_dend = T,show_column_dend = T,
              #row_names_side = "left",
              left_annotation = row_ha_left.possig
)

generate.PDF <- function(fig) {
  pdf("4.model/markers/regulators/heatmp.overlapped.with.sig.cor.genes.pdf",height = 7,width = 10)
  print(comp3)
  dev.off()
}
generate.PDF(fig)
#
#
#



 ggviolin(genestats,x = "kaps.group", y ="AIM2" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("CD8A")+ggtitle("CD8A")

ggviolin(genestats,x = "kaps.group", y ="HAVCR2" , fill = "kaps.group",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("kaps.group")+ylab("KRT20")+ggtitle("KRT20")


#
load("4.model/kaps/for.kaps.classification.output.k3.RData")
fitkaps
summary(cor.immune.drive.sur$z.of.mean.exp)
length(intersect(rownames(cor.immune.drive.sur.stage),rownames(genestats)))

dim(genestats)
dim(cor.immune.drive)
dim(cor.immune.drive.sur)
####
head(msi.com)
dim(msi.com)
head(surstats)
#####
dim(clin.CRC.full.sel)
head(clin.CRC.full.sel)
########
library(ggpubr)
ggscatter(genestats, x = "CD274", y = "CMTM4",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black",fill = "lightgray"),#shape = "MSI.status.bin",
          color = "msi.TE.group.bin", palette =  c("#d73027","#00AFBB","#E69F00","#756bb1","#27408B",
                                                   "#FF0000","#2E8B57","#CD00CD","#EE7942","#009ACD")
)+
  stat_cor(method = "spearman")
  #geom_vline(xintercept = median(cor.immune.drive$CD274,na.rm = T), color = "#d73027", size=0.2)+
  #geom_hline(yintercept = median(cor.immune.drive$CD8A,na.rm = T), color = "#d73027", size=0.2)
  #annotate("text", x=0.2, y=10000, label=paste0("n=","15804"),size=4)+
######CRC cell difference
g1=ggviolin(genestats,x = "msi.TE.group.bin", y ="MS4A1" , fill = "msi.TE.group.bin",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.MSI.group")+ylab("MS4A1")+ggtitle("MS4A1")


g2=ggviolin(genestats,x = "TE.cluster.agg", y ="MS4A1" , fill = "TE.cluster.agg",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#E69F00","#756bb1","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(label = "p.signif")+
  theme_bw()+xlab("TE.MSI.group")+ylab("MS4A1")+ggtitle("MS4A1")

generate.PDF <- function(fig) {
  pdf("4.model/markers/KRT20.pdf",width = 5,height = 8)
  print(g1)
  print(g2)
  dev.off()
}
generate.PDF(fig)

celltypeMarkers=c("MS4A1","CD3E","CD68","ITGAX","FCGR3A","FOXP3")







library(plot3D)
x=genestats$z.of.mean.exp
y=genestats$CD274
z=genestats$CD8A
pdf("tmp.pdf")
scatter3D(x, y, z, pch = 18,  theta = 20, phi = 20,
          col.var = genestats$z.of.mean.exp, 
          col = c("#1B9E77", "#D95F02"),
          main = "TCGA.CRC", xlab = "Z.of.TE.mean",
          ylab ="CD274", zlab = "CD8A")
dev.off()


#save.image(file = "tmp.RData")


ggscatter(kaps.td, x = "z.of.mean.exp", y = "LAG3",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black",fill = "lightgray"),#shape = "MSI.status.bin",
          color = "msi.TE.group.bin", palette =  c("#d73027","#00AFBB","#E69F00","#756bb1","#27408B","#FF0000","#2E8B57","#CD00CD","#EE7942","#009ACD")
)+
  stat_cor(method = "spearman", label.x = 2, label.y = 2) +
  geom_vline(xintercept = -0.41, color = "#d73027", size=0.2)+
  geom_hline(yintercept = 6, color = "#d73027", size=0.2)
#annotate("text", x=0.2, y=10000, label=paste0("n=","15804"),size=4)+

ggviolin(cor.immune.drive,x = "msi.TE.group.bin", y ="CD274" , fill = "msi.TE.group.bin",alpha = 1,size = 0.3,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
         add = "boxplot", add.params = list(fill = "white"))+labs(fill = "TE.cluster")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_bw()+xlab("TE.cluster")+ylab("CD174")+ggtitle("CD274")



dim(cor.immune.drive)
#######

####
fit1<- survfit(Surv(OS.time, OS) ~ msi.TE.group.bin, data = cor.immune.drive.sur)
ggsurvplot(fit1, data =cor.immune.drive.sur,
                   #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "TE.cluster",               # Change legend titles
                   legend.labs = c("TE.high_MSI", "TE.high_MSS","TE.low_MSI", "TE.low_MSS"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
                   palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
                   risk.table = TRUE,                  # Add No at risk table
                   cumevents = TRUE,                   # Add cumulative No of events table
                   tables.height = 0.15,               # Specify tables height
                   tables.theme = theme_cleantable(),  # Clean theme for tables
                   tables.y.text = FALSE,             # Hide tables y axis text
                   #ylab="Disease specific survival probability",
                   ylab="Overall survival probability",
                   #ylab="progression-free interval probability",
                   #ylab="disease-free interval probability",
                   xlab="Time(days)"
)
############
##########
#########kaps
########
save(cor.immune.drive.sur, file="4.model/kaps/for.kaps.classification.RData")
library(kaps)
fitkaps <- kaps(survival::Surv(OS.time, OS) ~ z.of.mean.exp, data = cor.immune.drive.sur, K = 2)

res.cut <- surv_cutpoint(cor.immune.drive.sur, time = "OS.time", event = "OS",
                         variables = c("z.of.mean.exp","CD274","CD8A"))
summary(res.cut)
plot(res.cut, "CD274", palette = "npg")


order(res.cut$CD8A$stats, decreasing = T)

median(aaaa$z.of.mean.exp)

aaaa=cor.immune.drive.sur
aaaa$os.group.bi=ifelse(aaaa$CD274> 6,"2high","1low")
fit1<- survfit(Surv(OS.time, OS) ~ os.group.bi, data =aaaa)
ggsurvplot(fit1, data =aaaa,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TE.cluster",               # Change legend titles
           legend.labs = c("2high", "1low"),  # Change legend labels
           #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette
           palette = c("#d73027","#00AFBB","#E69F00","#756bb1"),
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE,             # Hide tables y axis text
           #ylab="Disease specific survival probability",
           ylab="Overall survival probability",
           #ylab="progression-free interval probability",
           #ylab="disease-free interval probability",
           xlab="Time(days)"
)
