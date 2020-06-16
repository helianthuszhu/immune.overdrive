#########heatmap of IPRES in KIRC
#########
head(ipres.score.kirc)
head(stat.kirc.kaps.vali)
dim(ipres.score.kirc)
dim(stat.kirc.kaps.vali)

ht.ipres.kirc=cbind(stat.kirc.kaps.vali[,c("z.of.mean.exp","kaps.group.kirc")],ipres.score.kirc[rownames(stat.kirc.kaps.vali),] )
ht.ipres.kirc=ht.ipres.kirc[order(ht.ipres.kirc$z.of.mean.exp,decreasing = T),]
head(ht.ipres.kirc)
#
col_ha.ipres.kirc=columnAnnotation(z.mean.IPRES=ht.ipres.kirc$z.mean.IPRES,
                              kaps.group=ht.ipres.kirc$kaps.group.kirc,
                              col=list(z.mean.IPRES=colorRamp2(c(-4,-2,0,2,4), 
                                                                c("#ffffcc","#d9f0a3","#addd8e","#78c679","#006837")),
                                kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                              ),show_annotation_name = TRUE)
#
ada.kirc=as.data.frame(t(ht.ipres.kirc[,-c(1:4)]))
z_ada.kirc=(ada.kirc - rowMeans(ada.kirc))/apply(ada.kirc,1,sd)
z_ada.kirc[z_ada.kirc< -3] <- -3
z_ada.kirc[z_ada.kirc>3] <- 3

min_cor = min(as.vector(z_ada.kirc))
max_cor = max(as.vector(z_ada.kirc))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=100)
#col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "Spectral")))(100))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#7fcdbb", "white","#ce1256"))(100))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#4393c3", "white","#ce1256"))(100))


pdf("8.KIRC/3.immune.sets.kirc/ipres/htmap.ipres.KIRC.pdf",width = 8,height = 3)
Heatmap(z_ada.kirc, name = "ipres", 
        #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
        #col = rev(viridis(10)),border = F,
        col=col.pal_cor,
        show_column_names = F,show_row_names = T,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha.ipres.kirc#, 
        #right_annotation = row_ha.right,
        #show_row_dend = T,show_column_dend = T,
        #row_names_side = "left",
        #left_annotation = row_ha.left
)
dev.off()
save(ht.ipres.kirc, file="8.KIRC/3.immune.sets.kirc/ipres/htmap.ipres.KIRC.RData")
