################
################run WGCNA find module gene
################
dim(expM)
library(WGCNA)
#group information
datTraits=kaps.td[,c(14,17,38)]
#datTraits$kaps.group=substr(datTraits$kaps.group,4,4)
#datTraits$MSI.status.bin=ifelse(datTraits$MSI.status.bin=="MSI",1,2)
#datTraits$kaps.group=as.numeric(paste(datTraits$kaps.group))
datTraits[1:5,]
#
#gene exp Matrix
datExpr0=as.data.frame(t(expM))
datExpr0=datExpr0[rownames(datTraits),]
datExpr0[1:4,1:4]
dim(datExpr0)
#
#filter genes
gsg = goodSamplesGenes(datExpr = datExpr0, minNSamples = 295,verbose = 3)
gsg$allOK
#
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  #if (sum(!gsg$goodSamples)>0)
   # printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  #datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  datExpr0 = datExpr0[, gsg$goodGenes]
}
#
datExpr0=as.data.frame(t(datExpr0))
datExpr0 = t(datExpr0[order(apply(datExpr0,1,mad), decreasing = T)[1:5000],])
###cluster
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
pdf("4.model/kaps/wgcna/cluster.tree.pdf",width = 20,height = 9)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

####filter outliers
# Plot a line to show the cut
abline(h = 180, col = "red")

#traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
#plotDendroAndColors(sampleTree, traitColors,
                    #groupLabels = names(datTraits),
 #                   main = "Sample dendrogram and trait heatmap")
dev.off()
############
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datExpr[1:4,1:4]
####
#datTraits=datTraits[rownames(datExpr),]
#table(datTraits$kaps.group)
#datTraits$kaps.group=paste0("set",datTraits$kaps.group)
datTraits=datTraits[rownames(datExpr),]
datTraits.score=as.data.frame(datTraits$z.of.mean.exp)
rownames(datTraits.score)=rownames(datTraits)
colnames(datTraits.score)="TE.score"
head(datTraits.score)
dim(datExpr)
dim(datTraits)
#
save(datExpr, datTraits, file = "4.model/kaps/wgcna/wgcna.CR.dataInput.RData")
######
######
######
##########################
##################step 2 ----beta value selection
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
# Plot the results:
##sizeGrWindow(9, 5)
pdf("4.model/kaps/wgcna/beta.value.pdf")
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#
#################
########
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "4.model/kaps/wgcna/CRC.kaps.wgcna",
  verbose = 3
)
table(net$colors) 
#
# open a graphics window
pdf("4.model/kaps/wgcna/Cluster.Dendrogram.net.pdf",width = 12,height = 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,net,
     file = "4.model/kaps/wgcna/CRC.kaps-02-networkConstruction-auto.RData")
##################
load("4.model/kaps/wgcna/CRC.kaps-02-networkConstruction-auto.RData")
load("4.model/kaps/wgcna/wgcna.CR.dataInput.RData")
##################Relating modules to external clinical traits
###########
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits.score, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


pdf("4.model/kaps/wgcna/module.trait.2.pdf")
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits.score),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
########################
################3
############
# Define variable weight containing the weight column of datTrait
TE.score = as.data.frame(datTraits.score$TE.score);
names(TE.score) = "TE.score"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, TE.score, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(TE.score), sep="");
names(GSPvalue) = paste("p.GS.", names(TE.score), sep="");
#################
################module1
module = "greenyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
pdf("4.model/kaps/wgcna/module.memnership.pdf",width = 5,height = 7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
##################
#################
#
module = "greenyellow";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
modProbes=as.data.frame(modProbes)
write.csv(modProbes,"4.model/kaps/wgcna/gene.module.1.csv" )
#########
#########
##########module2
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
pdf("4.model/kaps/wgcna/module.memnership.2.pdf",width = 5,height = 7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TE score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
##################
#################
module = "brown";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
modProbes=as.data.frame(modProbes)
write.csv(modProbes,"4.model/kaps/wgcna/gene.module.2.csv" )
##############################
##############################annotation
#########
group_g <- data.frame(gene=colnames(datExpr),
                      group=moduleColors)
write.csv(group_g,file = "4.model/kaps/wgcna/group_g.module.vs.genes.csv")
######
library(clusterProfiler)
library(org.Hs.eg.db)
tmp <- bitr(group_g$gene,fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
colnames(tmp)[1]="gene"
head(tmp)
de_gene_cluster <- merge(tmp, group_g, by="gene",all.x=TRUE)
table(de_gene_cluster$group)
head(de_gene_cluster)
######only brown and greenyellow modules
####3
de_gene_cluster.sel=subset(de_gene_cluster, group=="brown"| group=="greenyellow")
table(de_gene_cluster.sel$group)
###run go analysis
# Run full GO enrichment test
formula_res <- compareCluster(
  ENTREZID~group, 
  data=de_gene_cluster.sel, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)


pdf('4.model/kaps/wgcna/GOenrich.all.pdf',width = 8,height = 8)
dotplot(formula_res, showCategory=8)
dev.off()
######
go.res=as.data.frame(formula_res)
go.res=separate(go.res, col="GeneRatio",into = c("counts","total"),sep = "[/]")
go.res$GeneRatio=as.numeric(paste(go.res$counts))/as.numeric(paste(go.res$total))
head(go.res)
save(group_g,de_gene_cluster,de_gene_cluster.sel,formula_res,go.res,file="4.model/kaps/wgcna/step5_GOananlysis.enrichGOall.Rdata")
write.csv(go.res, "4.model/kaps/wgcna/step5_GOananlysis.enrichGOall.res.csv")
######individually draw
go.res.sel=read.table("4.model/kaps/wgcna/GOenrich.res.hand.selected.txt",header = T,row.names = 1,sep = "\t")
head(go.res.sel)
#
pdf('4.model/kaps/wgcna/GOenrich.defined.pdf',width = 10,height = 8)
ggplot(go.res.sel, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO enrichment")+ facet_grid(.~Cluster)
ggplot(go.res.sel, aes(x = group, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO enrichment")
dev.off()

#########KEGG
formula_res.kegg <- compareCluster(
  ENTREZID~group, 
  data=de_gene_cluster.sel, 
  fun="enrichKEGG", 
  #OrgDb="org.Hs.eg.db",
  #ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

pdf('4.model/kaps/wgcna/KEGGenrich.pdf',width = 8,height = 8)
dotplot(formula_res.kegg,showCategory=11)
dev.off()
#########################################
kegg.res=as.data.frame(formula_res.kegg)
kegg.res=separate(kegg.res, col="GeneRatio",into = c("counts","total"),sep = "[/]")
kegg.res$GeneRatio=as.numeric(paste(kegg.res$counts))/as.numeric(paste(kegg.res$total))
save(group_g,de_gene_cluster,de_gene_cluster.sel,formula_res.kegg,kegg.res,file="4.model/kaps/wgcna/step5_GOananlysis.enrichKEGG.Rdata")
write.csv(kegg.res, "4.model/kaps/wgcna/step5_GOananlysis.enrichKEGG.res.csv")
#individually draw
dot_df=kegg.res[c(4:6,11,15,21,25,26,28,32,33,36:38,41,43,45,47:57),]
## from Tommy's code
pdf('4.model/kaps/wgcna/KEGGenrich.defined.pdf',width = 10,height = 8)
ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("KEGG pathway enrichment")+ facet_grid(.~Cluster)
ggplot(dot_df, aes(x = group, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("KEGG pathway enrichment")
dev.off()

#formula_res.kegg@compareClusterResult$Description
###################
####################
####################
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 10); 
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
#
png("4.model/kaps/wgcna/step6_all_Network-heatmap.1.png")
TOMplot(plotTOM,geneTree,moduleColors,main="Network heapmap plot of all genes")
dev.off()
#
png("4.model/kaps/wgcna/step6_all_Network-heatmap.2.png")
TOMplot(plotTOM,geneTree,main="Network heapmap plot of all genes",
        Colors =ifelse(moduleColors=="brown"|moduleColors=="greenyellow"|moduleColors=="grey",moduleColors,"#3288bd")
)
dev.off()
### 我这里只取了1000个基因哈，我试了一下全部基因，结果跑了半个小时没跑完，被我强行退出！
  nSelect =1000
  set.seed(20)
  #select=sample(nGenes,size = nSelect)
  select = (moduleColors=="brown"|moduleColors=="greenyellow");
  table(select)
  selectTOM = dissTOM[select,select]
  selectTree = hclust(as.dist(selectTOM),method = "average")
  selectColors = moduleColors[select]
  plotDiss=selectTOM^7
  diag(plotDiss)=NA
  png("4.model/kaps/wgcna/step6_select_Network-heatmap.png",width = 300,height=300)
  TOMplot(plotDiss,selectTree,selectColors,main="Network heapmap of brown and greenyellow module")
          #col=colorRampPalette(c("#FF3E96","white","#00F5FF"))(256))
          #col=colorRampPalette(c("#e7298a","#ffffcc"))(256)
  
  dev.off()
  pdf("4.model/kaps/wgcna/step6_select_Network-heatmap.pdf",width =8,height=8)
  TOMplot(plotDiss,selectTree,selectColors,main="Network heapmap of brown and greenyellow module")
  #col=colorRampPalette(c("#FF3E96","white","#00F5FF"))(256))
  #col=colorRampPalette(c("#e7298a","#ffffcc"))(256)
  
  dev.off()
  ###############3
  ##############
  ###############heatmap of module gene
  ########
head(group_g)
group_g.draw=subset(group_g, group=="brown"|group=="greenyellow")
rownames(group_g.draw)=group_g.draw$gene
htmodule=datExpr[,colnames(datExpr) %in% group_g.draw$gene]
z.htmodule=as.data.frame(t(scale(htmodule))) # 'scale'可以对log-ratio数值进行归一化
z.htmodule[z.htmodule>4]=4 
z.htmodule[z.htmodule< -4]= -4
#
min_cor = min(as.vector(z.htmodule))
max_cor = max(as.vector(z.htmodule))
range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=50)
col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(50))
col.pal_cor = colorRamp2(range_cor,colorRampPalette(c("#3288bd", "white","#ae017e"))(50))

head(group_g.draw)

col_ha.top.module=columnAnnotation(z.of.mean.exp.score=anno_lines(datTraits$z.of.mean.exp),
                                     z.of.mean.exp=datTraits$z.of.mean.exp,
                                     kaps.group=datTraits$kaps.group,
                                     col=list(z.of.mean.exp=colorRamp2(c(-4,-2,0,2,4), 
                                                                       c("#ffffcc","#d9f0a3","#addd8e","#78c679","#006837")),
                                              kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                                     ),show_annotation_name = TRUE)
row_ha.module.gene=rowAnnotation(module=group_g.draw$group,
                                 col=list(module=c("brown"="brown","greenyellow"="greenyellow")))


htmodule=Heatmap(z.htmodule, name = "module.gene expression", 
                #col = colorRampPalette(c("#00F5FF", "white","#FF3E96"))(256),
                #col = rev(viridis(10)),border = F,
                col=col.pal_cor,show_row_dend = F,
                show_column_names = F,show_row_names = F,
                cluster_columns = F,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 8),
                top_annotation = col_ha.top.module, 
                #right_annotation = row_ha.right,
                #show_row_dend = T,show_column_dend = T,
                #row_names_side = "left",
                left_annotation = row_ha.module.gene
)
pdf("4.model/kaps/wgcna/heatmap.of.module.gene.pdf",width = 7,height = 6)
draw(htmodule, padding = unit(c(20, 5, 20,5), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "bottom"
     ,heatmap_legend_side = "bottom"
)
dev.off()
#####
save.image(file="4.model/kaps/wgcna/containing.WGCNA.result.RData")
##########################################

