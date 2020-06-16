
#######
library(Biobase)
library(dplyr)
library(tidyr)
library(pheatmap)
aa=readRDS("~/nas/Xiaoqiang/R.pacakages/TE/REdiscoverTEdata/inst/Fig4_data/eset_TCGA_TE_intergenic_logCPM.RDS")
rexp=as.data.frame(t(exprs(aa)))
rexp[1:4,1:4]
rclin=pData(aa)
rclin$control=substr(rownames(rclin),14,15)
rclin$control=ifelse(rclin$control=="01","tumor","normal")
head(rclin)

head(rclin.agg.tumor)



load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.plot/repinfor.sel.RData")
dim(rclin)
dim(rexp)
############
table(repinfo.sel$repFamily,repinfo.sel$repClass)
#LTR=subset(repinfo.sel, repClass=="LTR")
#LTR=subset(repinfo.sel, repFamily=="Alu")
#LTR=subset(repinfo.sel, repFamily=="L1")
LTR=subset(repinfo.sel, repFamily=="SVA")
head(LTR)
#####
da1=rexp[, colnames(rexp) %in% rownames(LTR)]
dim(da1)
#######
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
png("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/pan.ht.SVA.png")
pheatmap(as.matrix(t(da1)),color = colorRampPalette(c("firebrick3", "white","navy"))(50),
             show_rownames = F,show_colnames = F,
             #color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(256),
             annotation_row = LTR[,2:3],annotation_col = rclin[,c(1,6)],
             scale = "row",
             #clustering_method = "median",
             #method = c("pearson"),
             fontsize = 7
             #clustering_distance_rows = "euclidean",clustering_distance_cols = "minkowski",
            
)
dev.off()
#save_pheatmap_pdf(ptt, "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/pan.ht.pdf")
######

#correlation matrix
library(RColorBrewer)
library(pheatmap)
dat <- as.data.frame(t(da1))
dist_dat  <- cor(dat,method="spearman")
hc <- hclust(as.dist(1-dist_dat),method="ward.D2")
#mycolor = colorRampPalette(c("white","red"))(10)
#mycolor=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8)
#mycolor=colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)
library(viridis)
library(pheatmap)
mycolor=col = rev(viridis(10))
#col_row$ <- factor(col_row$tumor)
png("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/2.PCA/ht.cor.LTR.2.png")
pheatmap(dist_dat,show_rownames = F,border_color = NA,
              show_colnames = T,cluster_rows = hc,
              cluster_cols =hc,annotation_col = rclin[,c(1,6,7)],
              color=mycolor, annotation_legend = T
              #annotation_colors  = list(response = c("PD" = "#dd1c77","PRCR"="#74c476","SD"="#feb24c"),
               #                         treatment.status=c("On"="#74a9cf", "Pre"="#f768a1"))
              )
dev.off()
#save_pheatmap_pdf(pcor, "/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/ht.cor.pdf",height = 6)
