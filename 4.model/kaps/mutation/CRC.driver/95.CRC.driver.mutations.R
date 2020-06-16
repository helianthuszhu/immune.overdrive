###########oncogenic mutation in CRC
#################
muted.anno.ready=read.csv("4.model/kaps/mutation/CRC.driver/process.muted.data.ready.csv",header = T)
head(muted.anno.ready)
rownames(muted.anno.ready)=muted.anno.ready$Sample.ID
#######
head(kaps.td)
#
sh.id.mu=intersect(rownames(kaps.td),rownames(muted.anno.ready))
colnames(muted.anno.ready)
caldata.mu=cbind(kaps.td[sh.id.mu,]$kaps.group,muted.anno.ready[sh.id.mu,-(100)])[,-2]
colnames(caldata.mu)[1]="kaps.group"
head(caldata.mu)
#####
#fish=chisq.test(caldata[,13],caldata$CRSS.subtype,correct = T)
#fish$p.value
#
dfm=NULL
for (k in 5:ncol(caldata.mu)){
  p.value=chisq.test(caldata.mu[,k],caldata.mu$kaps.group,correct = T)
  val=as.data.frame(p.value$p.value)
  rownames(val)=colnames(caldata.mu)[k]
  dfm=rbind(dfm,val)
}
head(dfm)
write.csv(dfm,"4.model/kaps/mutation/CRC.driver/chisq.results.for.95.CRC.drivers.kaps.four.csv")
summary(dfm$`p.value$p.value`)
dim(dfm)
sig.dfm=subset(dfm,dfm$`p.value$p.value`<0.05)
head(sig.dfm)
dim(sig.dfm)
write.csv(sig.dfm,"4.model/kaps/mutation/CRC.driver/chsiq.results.for.30.significant.CRC.drivers.kaps.four.csv")
###################draw all the gene mutation profiles
##########################################################final hetamp
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(circlize)
######################
head(caldata.mu)
#caldata.draw=caldata.mu[order(caldata.mu$kaps.group),]
caldata.draw=caldata.mu
######pvalue anova
pvalueanova = dfm$`p.value$p.value`
is_siganova = pvalueanova < 0.05
pchanova = rep("*", length(pvalueanova))
pchanova[!is_siganova] = NA
###
col = list(kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027"),
           MSI.status.bin=c("MSI"="#e7298a","MSS"="#66a61e")
           )
top.anno <- columnAnnotation(kaps.group=caldata.draw$kaps.group,
                             MSI.status.bin=kaps.td[rownames(caldata.draw),]$MSI.status.bin,
                             col=col)
####
mycol=colorRampPalette(c("white", "white", "black"))(50)
mu.matrix.total=as.data.frame(t(caldata.draw[,-c(1:4)]))
mu.matrix.total[1:4,1:4]
dim(mu.matrix.total)
#
####choose significant mutations to show
dim(sig.dfm)
showinggs=rownames(sig.dfm)
idx=match(showinggs,rownames(mu.matrix.total))
#
# color mapping
pvalue_col_funanova = colorRamp2(c(10e-5, 0.05, 0.1), c("red", "white","green")) 
row_ha_right = rowAnnotation(
  pvalueanova = anno_simple(pvalueanova, col = pvalue_col_funanova, pch = pchanova),
  sig.genes = anno_mark(at = idx, labels = showinggs,labels_gp = gpar(fontsize = 6),link_width = unit(5, "mm")),
  #width = unit(2, "mm"),
  annotation_name_side = "bottom")


#row_ha_right = rowAnnotation(sig.genes = anno_mark(at = idx, labels = showinggs))
####

#
ht=Heatmap(as.matrix(mu.matrix.total), name = "CRC.95.driver.mutations",
           col = mycol,border = F,
           cluster_columns = F,
           show_column_names = F,show_row_names = F,
           row_names_gp = gpar(fontsize = 6),
           column_names_gp = gpar(fontsize = 10),
           #left_annotation = row_ha.left,
           right_annotation = row_ha_right,
           top_annotation = top.anno
           
)
# now we generate two legends, one for the p-value
#
lgd_pvalue.anova = Legend(title = "anova.p-value", col = pvalue_col_funanova, at = c(10e-5, 0.05, 0.1), 
                          labels = c("10e-5", "0.05", "0.1"))
lgd_sig.anova = Legend(pch = "*", type = "points", labels = "< 0.05")


ht1=draw(ht, annotation_legend_list = list(lgd_pvalue.anova, lgd_sig.anova))

generate.PDF <- function(fig) {
  pdf("4.model/kaps/mutation/CRC.driver/ht.95genes.with.anova.p.sig.annoted.3.pdf",
      height = 5,width = 10)
  print(ht1)
  dev.off()
}
generate.PDF(fig)
#########
######count the number of mutated genes in each samples across subtypes
dim(caldata.draw)
cal.datasel=caldata.draw[,-c(2:4)]
cal.datasel[1:4,1:3]
d1=cal.datasel[,-1]
d1[1:3,1:3]
dim(d1)
mcounts=as.data.frame(rowSums(d1))
mcounts=cbind(mcounts,cal.datasel[rownames(mcounts),])
mcounts[1:4,1:4]
mcounts$proportion=mcounts$`rowSums(d1)`/95
colnames(mcounts)[1]="No.of.driver.genes.mutation"
head(mcounts)
write.csv(mcounts,"4.model/kaps/mutation/CRC.driver/boxplot.No.driver.mutations.csv")
####
library(ggplot2)
library(ggpubr)
cm.col= c("#00AFBB","#756bb1","#E69F00","#d73027")
my_comparisons <- list( c("set4", "set3"), c("set4", "set2"), c("set4", "set1"),
                        c("set3", "set2"), c("set3", "set1"), c("set2", "set1"))
p1=ggplot(mcounts, aes(x = kaps.group, y = No.of.driver.genes.mutation)) +
  geom_boxplot(alpha = 0.5)+
  stat_compare_means(label = "p.signif",method = "wilcox.test",comparisons = my_comparisons)+
  geom_jitter(aes(color = kaps.group),width = 0.2,size=1)+theme_bw()+
  scale_fill_manual(values=cm.col)+ # Boxplot fill color
  scale_color_manual(values = cm.col)#+ 
#geom_hline(yintercept=8, linetype="dashed",color = "black", size=0.5)
p2=ggplot(mcounts, aes(x = kaps.group, y = No.of.driver.genes.mutation)) +
  geom_boxplot(alpha = 0.5)+
  stat_compare_means(label = "p.signif",method = "kruskal.test")+
  geom_jitter(aes(color = kaps.group),width = 0.2,size=1)+theme_bw()+
  scale_fill_manual(values=cm.col)+ # Boxplot fill color
  scale_color_manual(values = cm.col)
generate.PDF <- function(fig) {
  pdf("4.model/kaps/mutation/CRC.driver/boxplot.No.driver.mutations.pdf",width = 4,height = 4)
  print(p1)
  print(p2)
  dev.off()
}
generate.PDF(fig)










