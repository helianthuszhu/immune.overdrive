########### compare pan cancer pathways activities
##############
pathypancancer.19k=read.table("4.model/immune/pancancer.pathways.xena/19k.pathways/merge_merged_reals_sample_level.txt.gz",header = T,row.names = 1,sep = "\t")
pathypancancer.19k=as.data.frame(t(pathypancancer.19k))
rownames(pathypancancer.19k)=gsub("[.]","-",rownames(pathypancancer.19k))
pathypancancer.19k[1:4,1:4]
#####
id19k=intersect(rownames(kaps.td), rownames(pathypancancer.19k))
#
eset.kaps.19k=pathypancancer.19k[id19k,]
eset.kaps.19k[1:3,1:3]
dim(eset.kaps.19k)
kaps.td.19k=kaps.td[id19k,]
table(kaps.td.19k$kaps.group)
#######
####
####构建对应的group_list
subtype <- as.factor(unique(kaps.td.19k$kaps.group))
group_list <- c()
for( i in 1:length(subtype)){
  tmp<- ifelse(kaps.td.19k$kaps.group==subtype[i],paste0(subtype[i]),'other')
  group_list<- as.data.frame(cbind(group_list,tmp))
  colnames(group_list)[i] <- paste0(subtype[i])
}
head(group_list)
dim(group_list)
deglist=list()
for( i in 1:ncol(group_list)){
  design2=model.matrix(~0+ factor(group_list[,i]))
  colnames(design2)=levels(factor(group_list[,i]))
  rownames(design2)=rownames(kaps.td.19k)
  fit=lmFit(t(eset.kaps.19k),design2)
  constrsts<- paste0(colnames(group_list)[i],"-other")
  contrast.matrix<-makeContrasts(contrasts=constrsts,levels = design2)
  fit=contrasts.fit(fit,contrast.matrix)
  fit=eBayes(fit)
  ####上面得到了p值等统计的结果，topTable对p值校验，对基因排序
  tT_tmp <- topTable(fit, adjust="fdr", number=nrow(fit))
  tT_tmp<- subset(tT_tmp, select=c('adj.P.Val',"P.Value","logFC"))
  colnames(tT_tmp)<- paste0(colnames(group_list)[i],colnames(tT_tmp))
  deglist[[i]] <- tT_tmp
}
######
degs.stat=list()
siglist=list()
for(i in 1:length(subtype)){
  tT<- deglist[[i]]
  #colnames(tT)<- gsub(subtype[1],'',colnames(tT))
  #tT$change = ifelse(tT[,1] < 0.05 & tT[,3]>1,TRUE,FALSE)
  tT$change = ifelse(tT[,1] < 0.001 ,TRUE,FALSE)   #change
  colnames(tT)[4]=paste0(subtype[i],"change")
  #siggene=rownames(subset(tT,change=="TURE"))
  sub.tt=subset(tT,tT[,1] < 0.001 )
  degs.stat[[i]]=tT
  siglist[[i]]=sub.tt
}
kaps.four.degs.19k=cbind(degs.stat[[1]][colnames(eset.kaps.19k),],degs.stat[[2]][colnames(eset.kaps.19k),],
                         degs.stat[[3]][colnames(eset.kaps.19k),],degs.stat[[4]][colnames(eset.kaps.19k),])
head(kaps.four.degs.19k)
dim(kaps.four.degs.19k)
table(kaps.four.degs.19k$set3change)
#
genesig.19k=c(rownames(siglist[[1]]), rownames(siglist[[2]]),rownames(siglist[[3]]),rownames(siglist[[4]]))
genesigexp.19k=as.data.frame(t(eset.kaps.19k))[unique(genesig.19k),]
dim(genesigexp.19k)

#
degp1=pheatmap(as.matrix(genesig.19k),show_rownames = F,border_color = NA,
               show_colnames = F,fontsize = 6,
               #cluster_rows = hc,
               cluster_cols =F,
               annotation_col = kaps.td.19k[,c(1,13,14,17,38)],
               #annotation_row = cndi.rep[rownames(ppd),c(2,3)],
               scale="row",
               #color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8), 
               color=colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
               annotation_legend = T,
               annotation_colors  = list( TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                          TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                          #MSI.status=c("MSI"="#a65628", "MSI"="#ffd92f","MSS"="#8da0cb"),
                                          MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                          kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027"),
                                          repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                          repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
               )
)
generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/19k.pathways/heatmp.degs.sig.among.kaps.group.19k.pdf",height = 5,width = 10)
  print(degp1)
  dev.off()
}
generate.PDF(fig)
###
temexp.19k=as.matrix(genesigexp.19k)

temexp.19k=(temexp.19k - rowMeans(temexp.19k))/apply(temexp.19k,1,sd)
temexp.19k[temexp.19k< -3] <- -3
temexp.19k[temexp.19k>3] <- 3
#temexp[is.na(temexp)] <- 0
temexp.19k[1:4,1:4]
max(temexp.19k)
col_ha.top1.pathyp.19k=columnAnnotation(z.of.mean.exp.score=anno_lines(kaps.td.19k$z.of.mean.exp),
                                         TE.cluster = kaps.td.19k$TE.cluster,
                                         TE.cluster.agg=kaps.td.19k$TE.cluster.agg,
                                         MSI.status.bin=kaps.td.19k$MSI.status.bin,
                                         z.of.mean.exp=kaps.td.19k$z.of.mean.exp,
                                         kaps.group=kaps.td.19k$kaps.group,
                                         col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                                  TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                                  MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                                  kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                                         ),show_annotation_name = TRUE)
pathpp2=Heatmap(temexp.19k, 
                #name = paste0("significant among kaps group","\n","(anova < 0.05)","\n",
                 #             "24 out of 40 pathways"), 
                na_col = "gray",
                col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                #col = rev(viridis(10)),border = F,
                show_column_names = F,show_row_names = F,
                cluster_columns = F,
                cluster_rows = T,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                top_annotation = col_ha.top1.pathyp.19k#, 
                #right_annotation = row_ha_right#,
                #show_row_dend = T,show_column_dend = T,
                #row_names_side = "left",
                #left_annotation = row_ha_left.pathyp.drug
)
generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/19k.pathways/heatmp.degs.sig.among.kaps.group.19k.scale.pdf",height = 5,width = 10)
  print(pathpp2)
  dev.off()
}
generate.PDF(fig)
save(pathypancancer.19k,eset.kaps.19k,kaps.td.19k, kaps.four.degs.19k,file = "4.model/immune/pancancer.pathways.xena/19k.pathways/heatmp.degs.sig.among.kaps.group.19k.scale.RData")
########################
##################show top 10 paths
########
topset4.19k=rownames(subset(kaps.four.degs.19k[order(kaps.four.degs.19k$set4logFC,decreasing = T),], set4adj.P.Val < 0.01))[1:10]

topset3.19k=rownames(subset(kaps.four.degs.19k[order(kaps.four.degs.19k$set3logFC,decreasing = T),], set3adj.P.Val < 0.01))[1:10]

topset2.19k=rownames(subset(kaps.four.degs.19k[order(kaps.four.degs.19k$set2logFC,decreasing = T),], set2adj.P.Val < 0.01))[1:10]

topset1.19k=rownames(subset(kaps.four.degs.19k[order(kaps.four.degs.19k$set1logFC,decreasing = T),], set1adj.P.Val < 0.01))[1:10]

genesig.19k.top=c(topset1.19k,topset2.19k,topset3.19k,topset4.19k)
genesigexp.19k.top=as.data.frame(t(eset.kaps.19k))[unique(genesig.19k.top),]
#
temexp.19k.top=as.matrix(genesigexp.19k.top)

temexp.19k.top=(temexp.19k.top - rowMeans(temexp.19k.top))/apply(temexp.19k.top,1,sd)
temexp.19k.top[temexp.19k.top< -3] <- -3
temexp.19k.top[temexp.19k.top>3] <- 3
#temexp[is.na(temexp)] <- 0
temexp.19k.top[1:4,1:4]

col_ha.top1.pathyp.19k=columnAnnotation(z.of.mean.exp.score=anno_lines(kaps.td.19k$z.of.mean.exp),
                                        TE.cluster = kaps.td.19k$TE.cluster,
                                        TE.cluster.agg=kaps.td.19k$TE.cluster.agg,
                                        MSI.status.bin=kaps.td.19k$MSI.status.bin,
                                        z.of.mean.exp=kaps.td.19k$z.of.mean.exp,
                                        kaps.group=kaps.td.19k$kaps.group,
                                        col=list(TE.cluster= c("cluster_1" = "#d73027","cluster_2"="#00AFBB","cluster_3"="#E69F00"),
                                                 TE.cluster.agg=c("TE.high"="#d73027","TE.low"="#E69F00"),
                                                 MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e"),
                                                 kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                                        ),show_annotation_name = TRUE)
pathpp3=Heatmap(temexp.19k.top, 
                #name = paste0("significant among kaps group","\n","(anova < 0.05)","\n",
                #             "24 out of 40 pathways"), 
                na_col = "gray",
                col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                #col = rev(viridis(10)),border = F,
                show_column_names = F,show_row_names = T,
                cluster_columns = F,
                cluster_rows = T,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                top_annotation = col_ha.top1.pathyp.19k#, 
                #right_annotation = row_ha_right#,
                #show_row_dend = T,show_column_dend = T,
                #row_names_side = "left",
                #left_annotation = row_ha_left.pathyp.drug
)
generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/19k.pathways/heatmp.degs.sig.among.kaps.group.19k.top10.pdf",height = 5,width = 10)
  print(pathpp3)
  dev.off()
}
generate.PDF(fig)
#####################show top 10 paths
##########
showingp.19k=genesig.19k.top
idxpath.19k=match(showingp.19k,rownames(temexp.19k))
#
row_ha.right.19k = rowAnnotation(foo2 = anno_mark(at = idxpath.19k, labels = showingp.19k,padding = unit(4, "mm"),link_height=unit(7, "mm"),
                                              link_width = unit(20, "mm")), show_annotation_name = FALSE)


pathpp4=Heatmap(temexp.19k, 
                name = paste0("significant among kaps group","\n",
                              "(limma adj.p.value < 0.01)","\n",
                            "3915 out of 19k pathways","\n",
                            "(top 10 logFC in each group)"), 
                na_col = "gray",
                col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                #col = rev(viridis(10)),border = F,
                show_column_names = F,show_row_names = F,
                cluster_columns = F,
                cluster_rows = T,
                row_names_gp = gpar(fontsize = 4),
                column_names_gp = gpar(fontsize = 8),
                top_annotation = col_ha.top1.pathyp.19k, 
                right_annotation = row_ha.right.19k,
                #show_row_dend = T,show_column_dend = T,
                #row_names_side = "left",
                #left_annotation = row_ha_left.pathyp.drug
)
generate.PDF <- function(fig) {
  pdf("4.model/immune/pancancer.pathways.xena/19k.pathways/heatmp.degs.sig.among.kaps.group.19k.scale.anno.top10.pdf",height = 10,width = 18)
  print(pathpp4)
  dev.off()
}
generate.PDF(fig)

#
