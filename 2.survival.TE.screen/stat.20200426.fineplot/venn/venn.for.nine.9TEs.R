######draw venn diagram for 9 TEs overlapped
#####
dim(te.counts)
length(sur.total.id)
length(intersect(rownames(te.counts), sur.total.id))

venndata1=data.frame(name=sur.total.id)
venndata1$type=rep("survival.assed", times=nrow(venndata1))
#
venndata2=data.frame(name=te.counts$repName)
venndata2$type=rep("immune.assed", times=nrow(venndata2))
#
venndata=rbind(venndata1, venndata2)
write.csv(venndata, "2.survival.TE.screen/stat.20200426.fineplot/venn/venn.for.nine.9TEs.csv")
####
library(VennDiagram)
#
venn.diagram(x= list(survival.associated = sur.total.id, immune.associated=rownames(te.counts)),
             compression ="lzw",main="sig TEs overlapped among four endpoints",main.cex = 0.45,
             filename = "2.survival.TE.screen/stat.20200426.fineplot/venn/venn.for.nine.9TEs.tif", height = 450,
             width = 500,resolution =300,  col ="transparent", 
             fill =c("#66c2a5","#fc8d62"),alpha = 0.5, 
             cex = 0.45,fontfamily = "serif", fontface = "bold",
             #cat.col =c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"), 
             cat.cex = 0.45,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif", 
             rotation.degree = 0)
