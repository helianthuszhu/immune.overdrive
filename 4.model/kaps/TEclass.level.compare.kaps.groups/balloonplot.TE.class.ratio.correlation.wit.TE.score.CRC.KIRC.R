#########
#########ratioCounts of TE
#########
geneCountdata=readRDS("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.from.CGS/output.new/GENE_1_raw_counts.RDS")
geneCountdataMatrix=as.data.frame(geneCountdata$counts)
geneCountdataMatrix[1:5,1:4]
dim(geneCountdataMatrix)
########
teCountdata=readRDS("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.from.CGS/output.new/RE_intergenic_1_raw_counts.RDS")
teCountdataMatrix=as.data.frame(teCountdata$counts)
teCountdataMatrix[1:4,1:4]
dim(teCountdataMatrix)
#TE selected
head(cdrep)
dim(cdrep)
teCountdataMatrix.sel=teCountdataMatrix[rownames(cdrep),]
teCountdataMatrix.sel[1:4,1:4]
########
head(cdrep)
dim(cdrep)
dim(cdrep.sel)
################filtering TE
stat.genecounts=as.data.frame(colSums(geneCountdataMatrix))
colnames(stat.genecounts)="geneCount"
stat.genecounts$id=rownames(stat.genecounts)
head(stat.genecounts)
#
stat.tecounts=as.data.frame(colSums(teCountdataMatrix.sel))
colnames(stat.tecounts)="teCount"
stat.tecounts$id=rownames(stat.tecounts)
head(stat.tecounts)
#
length(intersect(rownames(stat.genecounts), rownames(stat.tecounts)))
####
ratiocountdata=merge(stat.genecounts, stat.tecounts, by="id", all=FALSE)
ratiocountdata$id=gsub("_1_","",ratiocountdata$id)
ratiocountdata$platform=substr(ratiocountdata$id,nchar(as.character(ratiocountdata$id))-4,nchar(as.character(ratiocountdata$id))-3)
ratiocountdata$histology=substr(ratiocountdata$id,nchar(as.character(ratiocountdata$id))-25,nchar(as.character(ratiocountdata$id))-24)
ratiocountdata$id2=substr(ratiocountdata$id,1,nchar(as.character(ratiocountdata$id))-23)
ratiocountdata$id3=substr(ratiocountdata$id,16,16)
ratiocountdata$id4=substr(ratiocountdata$id,1,15)
head(ratiocountdata)
dim(ratiocountdata)
table(ratiocountdata$id3)
length(unique(ratiocountdata$id))
########
ratiocountdata.agg= ratiocountdata[,c(2,3,8)] %>% group_by(id4) %>% summarise_all(mean)
ratiocountdata.agg=as.data.frame(ratiocountdata.agg)
rownames(ratiocountdata.agg)=ratiocountdata.agg$id4
ratiocountdata.agg$totalCount=rowSums(ratiocountdata.agg[,2:3])
ratiocountdata.agg$ratio.teCount.vs.total=ratiocountdata.agg$teCount/ratiocountdata.agg$totalCount
ratiocountdata.agg$ratio.teCount.vs.geneCount=ratiocountdata.agg$teCount/ratiocountdata.agg$geneCount
head(ratiocountdata.agg)
#
length(intersect(rownames(ratiocountdata.agg), rownames(kaps.td)))
colnames(kaps.td)
#####
ratiodata.draw=cbind(kaps.td[,c(17,38)], ratiocountdata.agg[rownames(kaps.td),])
ratiodata.draw.sel=subset(ratiodata.draw, ratio.teCount.vs.total<=0.02)
head(ratiodata.draw)
##########
##########
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/ratio.correlation.pdf",width = 6,height = 7)
ggscatter(ratiodata.draw, x = "z.of.mean.exp", y = "ratio.teCount.vs.total", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "ratio.teCount.vs.total")+
  ggtitle("TE.score.vs.ratioCount normalized by total counts")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(ratiodata.draw, x = "z.of.mean.exp", y = "ratio.teCount.vs.geneCount", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "ratio.teCount.vs.geneCount")+
  ggtitle("TE.score.vs.ratio.teCount vs geneCount")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(ratiodata.draw.sel, x = "z.of.mean.exp", y = "ratio.teCount.vs.total", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "ratio.teCount.vs.total")+
  ggtitle("TE.score.vs.ratioCount normalized by total counts")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")

ggscatter(ratiodata.draw.sel, x = "z.of.mean.exp", y = "ratio.teCount.vs.geneCount", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "z.of.mean.exp", ylab = "ratio.teCount.vs.geneCount")+
  ggtitle("TE.score.vs.ratio.teCount vs geneCount")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
################
################
head(kaps.te.class.compare)
head(ratiodata.draw)
#
ratiodata.draw.ball=cbind(ratiodata.draw, kaps.te.class.compare[rownames(ratiodata.draw),c(3,5,7,11,13,15)])
head(ratiodata.draw.ball)
#calculate the correlation value
datalist=list()
for (j in c(7,9:14)) {
  res <- cor.test(ratiodata.draw.ball[,1], ratiodata.draw.ball[,j],  method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(ratiodata.draw.ball)[j]
  datalist[[j]] <- res
}
res.cor.ratiocount=do.call(rbind, datalist)
colnames(res.cor.ratiocount)=paste("CRC",colnames(res.cor.ratiocount),sep = "_")
res.cor.ratiocount$id=rownames(res.cor.ratiocount)
#############
res.cor.ratiocount.ball.draw=as.data.frame(res.cor.ratiocount[,2])
rownames(res.cor.ratiocount.ball.draw)=rownames(res.cor.ratiocount)
colnames(res.cor.ratiocount.ball.draw)="coefs with TE.score"
res.cor.ratiocount.ball.draw
############calculate correlation in KIRC
############
head(kaps.te.class.compare.kirc)
head(kaps.te.class.compare.kirc[,c(4,6,8,12,14,16)])
datalist=list()
for (j in c(4,6,8,12,14,16)) {
  res <- cor.test(kaps.te.class.compare.kirc[,1], kaps.te.class.compare.kirc[,j],  method = "spearman")
  pvalue=res$p.value
  cor=res$estimate
  res=data.frame(pvalue, cor)
  rownames(res)=colnames(kaps.te.class.compare.kirc)[j]
  datalist[[j]] <- res
}
res.cor.ratiocount.kirc=do.call(rbind, datalist)
colnames(res.cor.ratiocount.kirc)=paste("KIRC",colnames(res.cor.ratiocount.kirc),sep = "_")
res.cor.ratiocount.kirc$id=rownames(res.cor.ratiocount.kirc)
###################
res.cor.ratiocount.ball.draw=merge(res.cor.ratiocount,res.cor.ratiocount.kirc,by="id",all=TRUE)
rownames(res.cor.ratiocount.ball.draw)=res.cor.ratiocount.ball.draw$id
res.cor.ratiocount.ball.draw=res.cor.ratiocount.ball.draw[c(1:3,5:7,4),c(3,5)]
write.csv(res.cor.ratiocount.ball.draw,"4.model/kaps/TEclass.level.compare.kaps.groups/balloonplot.TE.class.ratio.correlation.wit.TE.score.CRC.KIRC.csv" )
############
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/balloonplot.TE.class.ratio.correlation.wit.TE.score.CRC.KIRC.pdf",width = 4,height = 5)
ggballoonplot(res.cor.ratiocount.ball.draw, fill = "value",size.range = c(7, 10))+
  scale_fill_gradientn(colors = c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
                                  "#E16462FF", "#FCA636FF", "#F0F921FF"))+
  ggtitle(paste0("correlation with TE score"))+theme_bw()+xlab("TE score")+ylab("TE class")
dev.off()

#
save(res.cor.ratiocount.ball.draw,res.cor.ratiocount, res.cor.ratiocount.kirc,ratiodata.draw,file="4.model/kaps/TEclass.level.compare.kaps.groups/balloonplot.TE.class.ratio.correlation.wit.TE.score.CRC.KIRC.RData")
####################
####################
###################
#######circle visualization
#######
library(circlize)
#
res.cor.ratiocount.ball.draw

circleCor=as.data.frame(t(res.cor.ratiocount.ball.draw))
rownames(circleCor)=c("Cor.with.TE.score.CRC","Cor.with.TE.score.KIRC")

#
circleCor.list = data.frame(from = rep(rownames(circleCor), times = ncol(circleCor)),
                to = rep(colnames(circleCor), each = nrow(circleCor)),
                value = as.vector(circleCor),
                stringsAsFactors = FALSE)

grid.col = c(Cor.with.TE.score.CRC = "yellow",Cor.with.TE.score.KIRC = "#fb6a4a",
             DNA="#1f78b4",LINE="#d95f02",LTR="#7570b3",Retroposon="#e7298a",Satellite="#66a61e",SINE="#e6ab02",ratio.teCount.vs.total="#e31a1c")

#chordDiagram(as.matrix(mat),grid.col = grid.col,transparency = 0.5)

col_circle = matrix(c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#e31a1c",
                      "#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#e31a1c"),ncol = 7,byrow = T)
dim(col_circle) = dim(circleCor)  # to make sure it is a matrix
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/circle.plot.correlation.TEscore.vs.class.level.rationCount.CRC.KIRC.pdf",width = 5,height = 5)
chordDiagram(as.matrix(circleCor), grid.col = grid.col, col = col_circle)
dev.off()

########### delete ratioCounts

circleCor=as.data.frame(t(res.cor.ratiocount.ball.draw))
rownames(circleCor)=c("Cor.with.TE.score.CRC","Cor.with.TE.score.KIRC")
circleCor=circleCor[,-7]
#
circleCor.list = data.frame(from = rep(rownames(circleCor), times = ncol(circleCor)),
                            to = rep(colnames(circleCor), each = nrow(circleCor)),
                            value = as.vector(circleCor),
                            stringsAsFactors = FALSE)

grid.col = c(Cor.with.TE.score.CRC = "yellow",Cor.with.TE.score.KIRC = "#e31a1c",
             DNA="#1f78b4",LINE="#d95f02",LTR="#7570b3",Retroposon="#e7298a",Satellite="#66a61e",SINE="#e6ab02")

#chordDiagram(as.matrix(mat),grid.col = grid.col,transparency = 0.5)
#"#e31a1c"

col_circle = matrix(c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02",
                      "#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),ncol = 6,byrow = T)
dim(col_circle) = dim(circleCor)  # to make sure it is a matrix
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/circle.plot.correlation.TEscore.vs.six.class.level.CRC.KIRC.pdf",width = 5,height = 5)
chordDiagram(as.matrix(circleCor), grid.col = grid.col, col = col_circle)
dev.off()
###################################
#################individually draw
###CRC
#
circleCor.CRC=circleCor[1,]

circleCor.list = data.frame(from = rep(rownames(circleCor.CRC), times = ncol(circleCor.CRC)),
                            to = rep(colnames(circleCor.CRC), each = nrow(circleCor.CRC)),
                            value = as.vector(circleCor.CRC),
                            stringsAsFactors = FALSE)

grid.col = c(Cor.with.TE.score.CRC = "yellow",
             #Cor.with.TE.score.KIRC = "#e31a1c",
             DNA="#1f78b4",LINE="#d95f02",LTR="#7570b3",Retroposon="#e7298a",Satellite="#66a61e",SINE="#e6ab02")

#chordDiagram(as.matrix(mat),grid.col = grid.col,transparency = 0.5)
#"#e31a1c"

col_circle = matrix(c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"
                      ),ncol = 6,byrow = T)
dim(col_circle) = dim(circleCor.CRC)  # to make sure it is a matrix
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/circle.plot.correlation.TEscore.vs.six.class.level.CRC.only.pdf",width = 5,height = 5)
chordDiagram(as.matrix(circleCor.CRC), grid.col = grid.col, col = col_circle)
dev.off()
##############
#######KIRC
circleCor.KIRC=circleCor[2,]

circleCor.list = data.frame(from = rep(rownames(circleCor.KIRC), times = ncol(circleCor.KIRC)),
                            to = rep(colnames(circleCor.KIRC), each = nrow(circleCor.KIRC)),
                            value = as.vector(circleCor.KIRC),
                            stringsAsFactors = FALSE)

grid.col = c(#Cor.with.TE.score.CRC = "yellow",
             Cor.with.TE.score.KIRC = "#e31a1c",
             DNA="#1f78b4",LINE="#d95f02",LTR="#7570b3",Retroposon="#e7298a",Satellite="#66a61e",SINE="#e6ab02")

#chordDiagram(as.matrix(mat),grid.col = grid.col,transparency = 0.5)
#"#e31a1c"

col_circle = matrix(c("#1f78b4","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),ncol = 6,byrow = T)
dim(col_circle) = dim(circleCor.KIRC)  # to make sure it is a matrix
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/circle.plot.correlation.TEscore.vs.six.class.level.KIRC.only.pdf",width = 5,height = 5)
chordDiagram(as.matrix(circleCor.KIRC), grid.col = grid.col, col = col_circle)
dev.off()
##################
##################
#################correlation TE score with ratioCounts scatter plot
###
head(ratiodata.draw)
pdf("4.model/kaps/TEclass.level.compare.kaps.groups/cor.ratioCounts.vs.TEscore.in.CRC.pdf",width = 5,height = 5)
ggscatter(ratiodata.draw, x = "z.of.mean.exp", y = "ratio.teCount.vs.total", 
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "TE score", ylab = "ratio.teCount.vs.total")+
  ggtitle("TE score vs ratio.teCount.vs.total in CRC")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")
dev.off()
