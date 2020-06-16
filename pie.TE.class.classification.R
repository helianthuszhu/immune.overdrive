#######select out TE for analysis
load("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/output.plot/repinfor.sel.RData")
#
head(repinfo.sel)
cdrep=repinfo.sel
cdrep$repClass.new=gsub('[?]',"_",cdrep$repClass)
table(cdrep$repFamily)
library(DataCombine)
Replaces <- data.frame(from = unique(cdrep$repClass.new), 
                       to = c("other.repeats","other.repeats","other.repeats","Satellite","SINE","DNA",
                              "other.repeats","LINE","other.repeats","LTR","other.repeats","other.repeats",
                              "other.repeats","other.repeats","other.repeats","Retroposon","other.repeats"))
cdrep$repClass.new2<- FindReplace(data = cdrep, Var = "repClass.new", replaceData = Replaces,
                        from = "from", to = "to",vector=TRUE)


table(cdrep$repFamily, cdrep$repClass.new2)
########
tedis=data.frame(unclass(table(cdrep$repClass.new2)))
tedis$class=rownames(tedis)
colnames(tedis)[1]='count'
#
tedis= tedis %>% group_by(class) %>% summarise_all(sum)
tedis=as.data.frame(tedis)
colnames(tedis)[2]="value"
head(tedis)
#
#colors <- c("#377eb8","#4daf4a","#e41a1c",'#27408B', '#FF0000', '#2E8B57', '#CD00CD', '#EE7942','#009ACD')
#colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3", "#ff7f00","#ffff33","#a65628")
colors= c("#1f78b4","#d95f02","#7570b3","#c6dbef","#e7298a","#66a61e","#e6ab02")
dfpie=tedis
slices <- dfpie$value
lbls <- dfpie$class
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
library(ggplotify)
pdf("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/1.TE.used/pie.TE.class.classification.pdf")
pie(slices,labels = lbls, col=colors,radius = 0.7,
    main="TE class level classification n=1204")
dev.off()
######
save(cdrep, file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/1.TE.used/TE.label.used.RData")
