#chrARMdata=read.csv("4.model/kaps/mutation/CNA.arm/clin_TCGA577_combined_adding.csv",header = T,row.names = 1)
#head(chrARMdata)
#GIdata=read.table("4.model/kaps/mutation/CNA.arm/GI.clinical.txt",header = T,row.names = 1,sep = "\t")
##
GIdata=read.table("4.model/kaps/mutation/snv/data_clinical_sample.armstatus.txt",header = T,sep = "\t",row.names = 1)
head(GIdata)
########
###################3
library(ComplexHeatmap)
library(circlize)
####################
#armdata <- chrARMdata[,c(44:83)]
#rownames(armdata)=gsub("[.]","-",rownames(armdata))
#armdata[armdata=="not_available"]=NA

#armdata[] <- lapply(armdata, function(x) as.numeric(as.character(x)))

yyy=kaps.td
rownames(yyy)=substr(rownames(yyy),1,12)
head(yyy)
yyyid=intersect(rownames(yyy), rownames(GIdata))
###
statarm=cbind(yyy[yyyid, ]$kaps.group, GIdata[yyyid,])
colnames(statarm)[1]="kaps.group"
statarm=statarm[,-c(2,3)]
statarm[statarm==""]="Not Called"
head(statarm)
table(statarm$STATUS_1P)
#

library(RColorBrewer)
#mutation_col = structure(names = c("MUT", "WT", "G34R", "G34V", "K27M"), 
 #                        c("black", "white", "#4DAF4A", "#4DAF4A", "#377EB8"))
#cnv_col = c("gain" = "#E41A1C", "loss" = "#377EB8", "amp" = "#E41A1C", 
 #           "del" = "#377EB8", "normal" = "white")
#######get the pvalue of chisrq.test
head(statarm)
datalist=list()
for (i in 2: ncol(statarm)) {
  tmppvalue <- chisq.test(statarm[,1], statarm[,i],correct = T)$p.value
  tmppvalue=as.data.frame(tmppvalue)
  rownames(tmppvalue)=colnames(statarm)[i]
  datalist[[i]]=tmppvalue
}
#
armpvalue=do.call(rbind, datalist)
armpvalue$tmppvalue=round(armpvalue$tmppvalue,digits = 4)
class(armpvalue)
save(statarm, armpvalue,file="4.model/kaps/mutation/CNA.arm/CNarm.heatmap.RData")
######
kaps_col = c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
cnv_col = c("Gained" = "#E41A1C", "Lost" = "#377EB8", "Not Called" = "white")

ha = HeatmapAnnotation(
  kaps.group=statarm$kaps.group,
  STATUS_1P=statarm$STATUS_1P,
  STATUS_1Q=statarm$STATUS_1Q,
  STATUS_2P=statarm$STATUS_2P,
  STATUS_2Q= statarm$STATUS_2Q,
  STATUS_3P= statarm$STATUS_3P,
  STATUS_3Q= statarm$STATUS_3Q,
  STATUS_4P= statarm$STATUS_4P,
  STATUS_4Q= statarm$STATUS_4Q,
  STATUS_5P= statarm$STATUS_5P,
  STATUS_5Q= statarm$STATUS_5Q,
  STATUS_6P= statarm$STATUS_6P,
  STATUS_6Q= statarm$STATUS_6Q,
  STATUS_7P= statarm$STATUS_7P,
  STATUS_7Q= statarm$STATUS_7Q,
  STATUS_8P= statarm$STATUS_8P,
  STATUS_8Q= statarm$STATUS_8Q,
  STATUS_9P= statarm$STATUS_9P,
  STATUS_9Q= statarm$STATUS_9Q,
  STATUS_10P= statarm$STATUS_10P,
  STATUS_10Q= statarm$STATUS_10Q,
  STATUS_11P= statarm$STATUS_11P,
  STATUS_11Q= statarm$STATUS_11Q,
  STATUS_12P= statarm$STATUS_12P,
  STATUS_12Q= statarm$STATUS_12Q,
  STATUS_13Q= statarm$STATUS_13_13Q,
  STATUS_14Q= statarm$STATUS_14_14Q,
  STATUS_15Q= statarm$STATUS_15_15Q,
  STATUS_16P= statarm$STATUS_16P,
  STATUS_16Q= statarm$STATUS_16Q,
  STATUS_17P= statarm$STATUS_17P,
  STATUS_17Q= statarm$STATUS_17Q,
  STATUS_18P= statarm$STATUS_18P,
  STATUS_18Q= statarm$STATUS_18Q,
  STATUS_19P=  statarm$STATUS_19P,
  STATUS_19Q= statarm$STATUS_19Q,
  STATUS_20P= statarm$STATUS_20P,
  STATUS_20Q= statarm$STATUS_20Q,
  STATUS_21Q= statarm$STATUS_21_21Q,
  STATUS_22Q= statarm$STATUS_22_22Q,
  col = list(
    kaps.group=kaps_col,
    STATUS_1P=cnv_col,
    STATUS_1Q=cnv_col,
    STATUS_2P=cnv_col,
    STATUS_2Q= cnv_col,
    STATUS_3P= cnv_col,
    STATUS_3Q= cnv_col,
    STATUS_4P= cnv_col,
    STATUS_4Q= cnv_col,
    STATUS_5P= cnv_col,
    STATUS_5Q= cnv_col,
    STATUS_6P= cnv_col,
    STATUS_6Q= cnv_col,
    STATUS_7P= cnv_col,
    STATUS_7Q= cnv_col,
    STATUS_8P= cnv_col,
    STATUS_8Q= cnv_col,
    STATUS_9P= cnv_col,
    STATUS_9Q= cnv_col,
    STATUS_10P= cnv_col,
    STATUS_10Q= cnv_col,
    STATUS_11P= cnv_col,
    STATUS_11Q= cnv_col,
    STATUS_12P= cnv_col,
    STATUS_12Q= cnv_col,
    STATUS_13Q= cnv_col,
    STATUS_14Q= cnv_col,
    STATUS_15Q= cnv_col,
    STATUS_16P= cnv_col,
    STATUS_16Q= cnv_col,
    STATUS_17P= cnv_col,
    STATUS_17Q= cnv_col,
    STATUS_18P= cnv_col,
    STATUS_18Q= cnv_col,
    STATUS_19P= cnv_col,
    STATUS_19Q= cnv_col,
    STATUS_20P= cnv_col,
    STATUS_20Q= cnv_col,
    STATUS_21Q= cnv_col,
    STATUS_22Q= cnv_col),
  na_col = "grey", border = TRUE,
  show_legend = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                  FALSE, FALSE, FALSE, FALSE, FALSE
                  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    kaps.group = list(title = "TE.cluster"),
    STATUS_1Q = list(title = "CNV")
    )
)

zero_row_mat = matrix(nrow = 0, ncol = nrow(statarm))
ht = Heatmap(zero_row_mat, top_annotation = ha, column_title = "CNV.arm.profile")
pdf("4.model/kaps/mutation/snv/CNA.heatmap.pdf",height = 10,width = 15)
draw(ht, padding = unit(c(22, 60, 20,20), "mm"), #bottom, left, top, right paddings,
     annotation_legend_side = "right"
     #,heatmap_legend_side = "right"
     )
#draw(ht)
annotation_titles = c(kaps.group="TE.cluster",
                      STATUS_1P=paste0("CNA_1P"," ","p=", armpvalue$tmppvalue[1]),
                      STATUS_1Q=paste0("CNA_1Q"," ","p=", armpvalue$tmppvalue[2]),
                      STATUS_2P=paste0("CNA_2P"," ","p=", armpvalue$tmppvalue[3]),
                      STATUS_2Q= paste0("CNA_2Q"," ","p=", armpvalue$tmppvalue[4]),
                      STATUS_3P= paste0("CNA_3P"," ","p=", armpvalue$tmppvalue[5]),
                      STATUS_3Q= paste0("CNA_3Q"," ","p=", armpvalue$tmppvalue[6]),
                      STATUS_4P= paste0("CNA_4P"," ","p=", armpvalue$tmppvalue[7]),
                      STATUS_4Q= paste0("CNA_4Q"," ","p=", armpvalue$tmppvalue[8]),
                      STATUS_5P= paste0("CNA_5P"," ","p=", armpvalue$tmppvalue[9]),
                      STATUS_5Q= paste0("CNA_5Q"," ","p=", armpvalue$tmppvalue[10]),
                      STATUS_6P= paste0("CNA_6P"," ","p=", armpvalue$tmppvalue[11]),
                      STATUS_6Q= paste0("CNA_6Q"," ","p=", armpvalue$tmppvalue[12]),
                      STATUS_7P= paste0("CNA_7P"," ","p=", armpvalue$tmppvalue[13]),
                      STATUS_7Q= paste0("CNA_7Q"," ","p=", armpvalue$tmppvalue[14]),
                      STATUS_8P= paste0("CNA_8P"," ","p=", armpvalue$tmppvalue[15]),
                      STATUS_8Q= paste0("CNA_8Q"," ","p=", armpvalue$tmppvalue[16]),
                      STATUS_9P= paste0("CNA_9P"," ","p=", armpvalue$tmppvalue[17]),
                      STATUS_9Q= paste0("CNA_9Q"," ","p=", armpvalue$tmppvalue[18]),
                      STATUS_10P= paste0("CNA_10P"," ","p=", armpvalue$tmppvalue[19]),
                      STATUS_10Q= paste0("CNA_10Q"," ","p=", armpvalue$tmppvalue[20]),
                      STATUS_11P= paste0("CNA_11P"," ","p=", armpvalue$tmppvalue[21]),
                      STATUS_11Q= paste0("CNA_11Q"," ","p=", armpvalue$tmppvalue[22]),
                      STATUS_12P= paste0("CNA_12P"," ","p=", armpvalue$tmppvalue[23]),
                      STATUS_12Q= paste0("CNA_12Q"," ","p=", armpvalue$tmppvalue[24]),
                      STATUS_13Q= paste0("CNA_13Q"," ","p=", armpvalue$tmppvalue[25]),
                      STATUS_14Q= paste0("CNA_14Q"," ","p=", armpvalue$tmppvalue[26]),
                      STATUS_15Q= paste0("CNA_15Q"," ","p=", armpvalue$tmppvalue[27]),
                      STATUS_16P= paste0("CNA_16P"," ","p=", armpvalue$tmppvalue[28]),
                      STATUS_16Q= paste0("CNA_16Q"," ","p=", armpvalue$tmppvalue[29]),
                      STATUS_17P= paste0("CNA_17P"," ","p=", armpvalue$tmppvalue[30]),
                      STATUS_17Q= paste0("CNA_17Q"," ","p=", armpvalue$tmppvalue[31]),
                      STATUS_18P= paste0("CNA_18P"," ","p=", armpvalue$tmppvalue[32]),
                      STATUS_18Q=paste0("CNA_18Q"," ","p=", armpvalue$tmppvalue[33]),
                      STATUS_19P= paste0("CNA_19P"," ","p=", armpvalue$tmppvalue[34]),
                      STATUS_19Q= paste0("CNA_19Q"," ","p=", armpvalue$tmppvalue[35]),
                      STATUS_20P= paste0("CNA_20P"," ","p=", armpvalue$tmppvalue[36]),
                      STATUS_20Q= paste0("CNA_20Q"," ","p=", armpvalue$tmppvalue[37]),
                      STATUS_21Q= paste0("CNA_21Q"," ","p=", armpvalue$tmppvalue[38]),
                      STATUS_22Q= paste0("CNA_22Q"," ","p=", armpvalue$tmppvalue[39]))
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-40, "mm"), just = "left")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

dev.off()








#
for (i in 47: 87) {
  pdf(paste0("4.model/kaps/mutation/CNA.arm/",colnames(statarm)[i],".kaps.group.pdf"))
  print(ggviolin(statarm,x = "kaps.group", y = colnames(statarm)[i], fill = "kaps.group",alpha = 1,size = 0.3,
                 #palette = c("#d73027","#00AFBB", "#E69F00","#1f78b4","#d95f02","#7570b3","#e7298a","#e6ab02","#66a61e"),
                 palette = c("#00AFBB","#756bb1","#E69F00","#d73027"),
                 add = "boxplot", add.params = list(fill = "white"))+labs(fill = "kaps.group.agg")+
          stat_compare_means(label = "p.signif")+
          theme_bw()+xlab("kaps.group")+ylab(colnames(statarm)[i])+ggtitle(paste0(colnames(statarm)[i]))
  )
  dev.off()
}
#
##################
col_ha.arm=columnAnnotation(kaps.group=statarm$kaps.group,
                             col=list(
                             kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")
                             ),show_annotation_name = TRUE)


parm=pheatmap(t(statarm[,c(89:127)]),show_rownames = T,border_color = NA,main = "CN.arm",
         show_colnames = F,fontsize = 4,#na_col = "black",
         #cluster_rows = hc,
         cluster_cols =F,
         annotation_col = statarm[,1:2],
         #annotation_row = cndi.rep[rownames(ppd),c(2,3)],
         #scale="row",
         #color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8), 
         color=colorRampPalette(c("#ffffcc", "#feb24c","#e31a1c"))(256),
         #annotation_legend = T,
         annotation_colors  = list( 
                                    kaps.group=c("set1"="#00AFBB","set2"="#756bb1","set3"="#E69F00","set4"="#d73027")#,
                                    #repClass=c("DNA"="#1f78b4","LTR"="#7570b3","Retroposon"="#e7298a","SINE"="#e6ab02"),
                                    #repFamily=c("TcMar-Tigger"="#6baed6","ERV1"="#bcbddc","SVA"="#fa9fb5","Alu"="#fee391")
         )
)
save_pheatmap_pdf(parm, filename = "4.model/kaps/mutation/CNA.arm/armCN.heatmap.pdf")







