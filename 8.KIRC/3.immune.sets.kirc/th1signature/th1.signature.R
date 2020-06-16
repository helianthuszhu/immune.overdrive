##########Th1 signature for kirc
#####
dim(kirc.expfireh)
dim(stat.kirc.kaps.vali)
####
kirc.expfireh[1:4,1:4]
stat.kirc.kaps.vali[1:4,1:4]
length(intersect(rownames(stat.kirc.kaps.vali), colnames(kirc.expfireh)))
#
kirc.expfireh.sel=as.data.frame(t(kirc.expfireh))[rownames(stat.kirc.kaps.vali),]
#########
markerpanel=read.table("4.model/markers/immue.gene.to.show.for.heatmap.txt",header = T,sep = "\t")
panelM.th1.kirc=kirc.expfireh.sel[, colnames(kirc.expfireh.sel) %in% subset(markerpanel, category=="Th.1.signatures")$symbol]
#panelM.th1.kirc=as.data.frame(t(panelM.th1.kirc))
#z_panelM.th1.kirc=(panelM.th1.kirc - rowMeans(panelM.th1.kirc))/apply(panelM.th1.kirc,1,sd)
z_panelM.th1.kirc=as.data.frame(t(scale(panelM.th1.kirc)))
z_panelM.th1.kirc[z_panelM.th1.kirc< -4] <- -4
z_panelM.th1.kirc[z_panelM.th1.kirc>4] <- 4
dim(z_panelM.th1.kirc)
#
col_ha.top.th1.kirc=columnAnnotation(z.of.mean.exp.score=anno_lines(stat.kirc.kaps.vali$z.of.mean.exp),
                             z.of.mean.exp=stat.kirc.kaps.vali$z.of.mean.exp,
                             kaps.group=stat.kirc.kaps.vali$kaps.group.kirc,
                             
                             col=list(z.of.mean.exp=colorRamp2(c(-4,-2,0,2,4), 
                                                                c("#ffffcc","#d9f0a3","#addd8e","#78c679","#006837")),
                                      kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)

th1kirc=Heatmap(z_panelM.th1.kirc, name = "Th-1.signatures", 
             col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
             #col = rev(viridis(10)),border = F,
             show_column_names = F,show_row_names = T,
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 8),
             top_annotation = col_ha.top.th1.kirc#, 
             #right_annotation = row_ha.right,
             #show_row_dend = T,show_column_dend = T,
             #row_names_side = "left",
             #left_annotation = row_ha.left
)
generate.PDF <- function(fig) {
  pdf("8.KIRC/3.immune.sets.kirc/th1signature/ht.th1.signature.pdf",height = 5,width = 10)
  print(th1kirc)
  dev.off()
}
generate.PDF(fig)

pdf("8.KIRC/3.immune.sets.kirc/th1signature/ht.th1.signature.2.pdf",width = 8,height = 6)
draw(th1kirc, padding = unit(c(45, 5, 45,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     ,heatmap_legend_side = "right"
)
dev.off()

########################
############calculate the mean of TH1 signautre
############
colnames(panelM.th1.kirc)

head(panelM.th1.kirc)
Th1.sigM.kirc=panelM.th1.kirc
Th1.sigM.kirc$mean.Th1.sig=rowMeans(Th1.sigM.kirc,na.rm = T)
dim(Th1.sigM.kirc)
stat.th1sig.kirc=cbind(Th1.sigM.kirc[rownames(stat.kirc.kaps.vali),]$mean.Th1.sig,stat.kirc.kaps.vali)
colnames(stat.th1sig.kirc)[1]="mean.Th1.sig"

stat.th1sig.kirc$z.of.mean.Th1.sig=(stat.th1sig.kirc$mean.Th1.sig-mean(stat.th1sig.kirc$mean.Th1.sig))/sd(stat.th1sig.kirc$mean.Th1.sig)
head(stat.th1sig.kirc)

#####
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
stat.th1sig.kirc$kaps.group.kirc=factor(stat.th1sig.kirc$kaps.group,levels = c("set4", "set3","set2","set1"))

pdf("8.KIRC/3.immune.sets.kirc/th1signature/violin.th1.signature.pdf",height  = 5,width = 4)
ggviolin(stat.th1sig.kirc,x = "kaps.group.kirc", y ="z.of.mean.Th1.sig" , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
              #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
              palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
              add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+xlab("kaps.group")+ylab("z.of.mean.Th1.sig")+ggtitle("z.of.mean.Th1.sig.KIRC")

ggviolin(stat.th1sig.kirc,x = "kaps.group.kirc", y ="mean.Th1.sig" , fill = "kaps.group.kirc",alpha = 1,size = 0.01,width = 1,
         #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
         palette = c("set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" ),
         add = "boxplot", add.params = list(fill = "white",color="black",size=0.5,width=0.2))+labs(fill = "kaps.group.kirc")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")+xlab("kaps.group")+ylab("mean.Th1.sig")+ggtitle("mean.Th1.sig.KIRC")
dev.off()
#
save(stat.th1sig.kirc,panelM.th1.kirc,stat.kirc.kaps.vali,markerpanel,file="8.KIRC/3.immune.sets.kirc/th1signature/violin.th1.signature.RData")
