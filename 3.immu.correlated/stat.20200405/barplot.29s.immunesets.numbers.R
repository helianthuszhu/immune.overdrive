########
immugene=read.delim2("~/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/Geneappend3.from.immuport.txt")
head(immugene)
setnames=unique(immugene$Category)
length(setnames)
datalist=list()
for (i in 1: length(setnames)) {
  aa=subset(immugene,Category==setnames[i])
  bb=as.data.frame(length(aa$Symbol))
  rownames(bb)=setnames[i]
  colnames(bb)="No.genes"
  datalist[[i]]=bb
  
}
sets=do.call(rbind, datalist)
#####
sets2=data.frame(No.genes=rep(1,times=11))
rownames(sets2)=c("LAG3","TNFRSF14","BTLA","CD86","CD80","CTLA4","PDCD1LG2","CD274","PDCD1","CD8A","HAVCR2")
#####
sets3=data.frame(No.genes=c(141))
rownames(sets3)="Immunescore.estimate"
######
sets.29s=rbind(sets2,sets,sets3)
sets.29s$type=rep(c("type1","type2","type3"),times=c(11,17,1))
sets.29s$setName=rownames(sets.29s)
sets.29s$count=log2(sets.29s$No.genes+1)
#############
#############
library(ggplot2)
library(ggpubr)
pdf("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/3.immu.correlated/stat.20200405/barplot.29s.immunesets.numbers.pdf",width = 5,height = 5)
ggbarplot(sets.29s, x = "setName", y = "count",
          fill = "type",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          #color = "#525252", 
          palette = c("#a6d854","#ffd92f","#e78ac3"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = TRUE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal()
)+scale_y_continuous(breaks=c(0,1,3,5,7,9),position = "right")+
  #font("x.text", size = 5, vjust = 0.5)+
  theme(axis.text.x = element_text(size=8))+
  ylab("log2(Number of immune-related genes+1)")+ggtitle("29 immune sets")+coord_flip()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="right",
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
  )
dev.off()
#######################
#######################
library(ggplot2)
library(ggpubr)
cbPalette= c("#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e")


pdf("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/3.immu.correlated/stat.20200405/barplot.14.TE.counts.pdf",width = 5,height = 3)
ggbarplot(te.counts[,-1], x = "repName", y = "count",
          fill = "repClass.new2",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          #color = "#525252", 
          palette = c("#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = TRUE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal()
)+scale_y_continuous(breaks=c(0,1,3,6, 9,12),position = "right")+
  #scale_y_continuous(position = "right")+
  #font("x.text", size = 5, vjust = 0.5)+
  theme(axis.text.x = element_text(size=8))+
  ylab(paste0("Potentially immunogenic","\n", "with number of immue sets"))+ggtitle("TEs at least postively (cor>=0.4)")+coord_flip()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="right",
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
  )
dev.off()
######
save(te.counts, sets.29s, file="3.immu.correlated/stat.20200405/barplot.TE.immunecorrealted.RData")
