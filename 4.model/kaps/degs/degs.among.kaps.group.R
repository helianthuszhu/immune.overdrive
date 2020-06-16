######DEGs among kaps groups
########load gene expM
#
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/5.TE.cor.with.gene/TE.gene.CRC.matrix.RData")
gexp=as.data.frame(t(expM))
gexp[1:4,1:4]
dim(gexp)
###
############
count.gene=as.data.frame(colSums(is.na(gexp)))
colnames(count.gene)="count"
count.gene.sel=subset(count.gene,count<339)
head(count.gene.sel)
dim(count.gene.sel)
###
dim(kaps.td)
colnames(kaps.td)
degkapsid=intersect(rownames(kaps.td), rownames(gexp))
#####
eset.kaps=gexp[rownames(kaps.td), colnames(gexp) %in% rownames(count.gene.sel)]
eset.kaps[1:3,1:3]
dim(eset.kaps)
#######
####
####构建对应的group_list
subtype <- as.factor(unique(kaps.td$kaps.group))
group_list <- c()
for( i in 1:length(subtype)){
  tmp<- ifelse(kaps.td$kaps.group==subtype[i],paste0(subtype[i]),'other')
  group_list<- as.data.frame(cbind(group_list,tmp))
  colnames(group_list)[i] <- paste0(subtype[i])
}
head(group_list)
dim(group_list)
deglist=list()
for( i in 1:ncol(group_list)){
  design2=model.matrix(~0+ factor(group_list[,i]))
  colnames(design2)=levels(factor(group_list[,i]))
  rownames(design2)=rownames(kaps.td)
  fit=lmFit(t(eset.kaps),design2)
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
kaps.four.degs=cbind(degs.stat[[1]][colnames(eset.kaps),],degs.stat[[2]][colnames(eset.kaps),],degs.stat[[3]][colnames(eset.kaps),],degs.stat[[4]][colnames(eset.kaps),])
head(kaps.four.degs)
#write.csv(kaps.four.degs,"4.model/kaps/degs/kaps.four.degs.csv")
#dim(kaps.four.degs)
####
#genesig=c(rownames(siglist[[1]]), rownames(siglist[[2]]),rownames(siglist[[3]]),rownames(siglist[[4]]))
#genesig=Reduce(intersect, list(rownames(siglist[[1]]), rownames(siglist[[2]]),rownames(siglist[[3]]),rownames(siglist[[4]])))
#get the top 10 genes in each group
#####
topset4=rownames(siglist[[1]][order(siglist[[1]]$set4logFC,decreasing = T),])[1:20]

topset3=rownames(siglist[[2]][order(siglist[[2]]$set3logFC,decreasing = T),])[1:20]

topset2=rownames(siglist[[3]][order(siglist[[3]]$set2logFC,decreasing = T),])[1:20]

topset1=rownames(siglist[[4]][order(siglist[[4]]$set1logFC,decreasing = T),])[1:20]

genesig=c(topset1,topset2,topset3,topset4)

genesigexp=as.data.frame(t(eset.kaps))[unique(genesig),]
dim(genesigexp)
########
######
####check the overlapped te among these four endpoints
####
library(VennDiagram)
deg.set4=rownames(siglist[[1]])
deg.set3=rownames(siglist[[2]])
deg.set2=rownames(siglist[[3]])
deg.set1=rownames(siglist[[4]])
#
venn.diagram(x= list(deg.set4 = deg.set4, deg.set3 = deg.set3,deg.set2 = deg.set2, deg.set1=deg.set1), 
             compression ="lzw",main="sig TEs overlapped among four kaps groups",main.cex = 0.45,
             filename = "4.model/kaps/degs/venn.degs.tif", height = 450,
             width = 450,resolution =300,  col ="transparent", 
             fill =c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"),alpha = 0.5, 
             cex = 0.45,fontfamily = "serif", fontface = "bold",
             #cat.col =c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"), 
             cat.cex = 0.45,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif", 
             rotation.degree = 0)

######
#genesigexp[is.na(genesigexp)] <- 0
#colnames(kaps.td)
temexp=genesigexp
temexp[temexp< -2] <- -2
temexp[temexp>2] <- 2
temexp[is.na(temexp)] <- 0
temexp[1:4,1:4]

degp1=pheatmap(as.matrix(genesigexp),show_rownames = T,border_color = NA,
         show_colnames = F,fontsize = 6,
         #cluster_rows = hc,
         cluster_cols =F,
         annotation_col = kaps.td[,c(1,13,14,17,38)],
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
  pdf("4.model/kaps/degs/heatmp.degs.sig.among.kaps.group.top.20.pdf",height = 6,width = 10)
  print(degp1)
  dev.off()
}
generate.PDF(fig)

#
#
#
############compare between set4 and set3
#####
set34.ann=subset(kaps.td, kaps.group=="set3"| kaps.group=="set4")
set34.exp=eset.kaps[rownames(set34.ann),]
#
design2=model.matrix(~0+ factor(set34.ann$kaps.group))
colnames(design2)=levels(factor(set34.ann$kaps.group))
rownames(design2)=rownames(set34.ann)
fit=lmFit(t(set34.exp),design2)
contrast.matrix<-makeContrasts(contrasts='set4-set3',levels = design2)
fit=contrasts.fit(fit,contrast.matrix)
fit=eBayes(fit)
####上面得到了p值等统计的结果，topTable对p值校验，对基因排序
tT_tmp <- topTable(fit, adjust="fdr", number=nrow(fit))
tT_tmp$change = ifelse(tT_tmp[,4] < 0.01 & tT_tmp[,1]>1,TRUE,FALSE)
head(tT_tmp)
#tT_tmp<- subset(tT_tmp, select=c('adj.P.Val',"P.Value","logFC"))
deg.res.set34=tT_tmp
write.csv(deg.res.set34, "4.model/kaps/degs/kaps.set34.degs.csv")
###################################################################################################
#####GSEA analysis
###################################################################################################
###################################################################################################
###################################################################################################
library(topGO)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)#人类基因组注释相关的包
library(DO.db)
library(clusterProfiler)
####
head(kaps.four.degs)

#setres=kaps.four.degs[,1:4]
setres=kaps.four.degs[,5:8]
colnames(setres)=substr(colnames(setres),5, nchar(colnames(setres)))
head(setres)

#####
###########差异分析结果
AML_genelist_up<-setres[(which(setres$logFC>0 & setres$adj.P.Val < 0.05)),]
head(AML_genelist_up)
dim(AML_genelist_up)
###########ID转换
AML_genelist_up_l<-bitr(unique(row.names(AML_genelist_up)),
                        fromType="SYMBOL",toType="ENTREZID",
                        OrgDb="org.Hs.eg.db",drop = TRUE)#转换ID
AML_genelist_up_l<- dplyr::distinct(AML_genelist_up_l,SYMBOL,.keep_all=TRUE)
rownames(AML_genelist_up_l)=AML_genelist_up_l$SYMBOL
head(AML_genelist_up_l)
dim(AML_genelist_up_l)
#length(intersect(rownames(AML_genelist_up_l), rownames(AML_genelist_up)))

###########信息合并
AML_genelist_up<-cbind(AML_genelist_up_l, AML_genelist_up[rownames(AML_genelist_up_l),])
head(AML_genelist_up)
#
AML_up<-AML_genelist_up[,c(2,5)]
colnames(AML_up)<-c("ENSEMBL","logFC")
head(AML_up)
###########排序
geneList<-as.numeric(as.character(AML_up[,2]))
names(geneList) = as.character(AML_up[,1])
geneList= sort(geneList, decreasing = TRUE)#构建geneList,并根据logFC由高到低排列
#
AML_GSEA_KEGG_up<-gseKEGG(
  geneList =geneList,
  nPerm = 1000,#
  keyType = 'kegg',#可以选择"kegg",'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
  organism = 'hsa'#定义物种,
  #pvalueCutoff = 0.05, #自定义pvalue的范围
  #pAdjustMethod     = "BH" #校正p值的方法
)
AML_GSEA_KEGG_up<-gseGO(geneList     = geneList,
                        OrgDb        = org.Hs.eg.db,
                        ont          = "ALL",
                        nPerm        = 1000,
                        minGSSize    = 100,
                        maxGSSize    = 500,
                        pvalueCutoff = 1,
                        verbose      = FALSE)



sortAML_GSEA_KEGG_up<-AML_GSEA_KEGG_up[order(AML_GSEA_KEGG_up$enrichmentScore,decreasing=T)]

gseaplot2(AML_GSEA_KEGG_up,row.names(sortAML_GSEA_KEGG_up))

gseaplot2(AML_GSEA_KEGG_up,row.names(sortAML_GSEA_KEGG_up)[1:10],#只显示前三个GSEA的结果
          title="AML_GSEA_KEGG_up",#标题
          color = c("#626262","#8989FF","#FF0404"),#颜色
          pvalue_table = FALSE,
          ES_geom = "line"#enrichment scored的展现方式 'line' or 'dot'
          )
class(AML_GSEA_KEGG_up)
########
library(clusterProfiler)
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("4.model/kaps/degs/gsea.based.on.degs/h.all.v7.0.entrez.gmt")
# 需要网络
y <- GSEA(geneList,TERM2GENE =hallmarks)
y
class(y)
sort.y<-y[order(y$enrichmentScore,decreasing=T)]

gseaplot2(y,row.names(sort.y)[1],#只显示前三个GSEA的结果
          title="AML_GSEA_KEGG_up",#标题
          color = c("#626262","#8989FF","#FF0404"),#颜色
          pvalue_table = FALSE,
          ES_geom = "line"#enrichment scored的展现方式 'line' or 'dot'
)
