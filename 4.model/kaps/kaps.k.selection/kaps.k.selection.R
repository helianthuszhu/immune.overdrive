######for kaps plot
######
load("4.model/kaps/for.kaps.classification.output.k2-4.RData")
fitkaps
####
head(kaps.td)
table(kaps.td$kaps.group.three)
aa=kaps.td
aa$kaps.group.two=cut(aa$z.of.mean.exp, 
                      breaks = c(min(aa$z.of.mean.exp),0.0128,max(aa$z.of.mean.exp)),
                      labels = c("set1","set2"), include.lowest = TRUE)
write.csv(aa, "4.model/kaps/kaps.k.selection/suvival.data.csv")
#
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group, data = kaps.td )
#"set4"="#d73027","set3"="#E69F00","set2"="#756bb1","set1"="#00AFBB" 
kfour=ggsurvplot(fit2, data =kaps.td,
                   #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                   pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "kaps.four.group",               # Change legend titles
                   legend.labs = c("set1","set2","set3","set4"),  # Change legend labels
                   #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                   #palette =c("#00AFBB","#756bb1","#E69F00","#d73027"),
                   palette =c("#1b9e77","#d95f02","#7570b3","#e7298a"),
                   #palette = c("#ca0020","#0571b0","#4daf4a"),
                   risk.table = TRUE,                  # Add No at risk table
                   cumevents = TRUE,                   # Add cumulative No of events table
                   tables.height = 0.15,               # Specify tables height
                   tables.theme = theme_cleantable(),  # Clean theme for tables
                   tables.y.text = FALSE,             # Hide tables y axis text
                   #ylab="Disease specific survival probability",
                   ylab="Overall survival probability",
                   #ylab="progression-free interval probability",
                   #ylab="disease-free interval probability",
                   xlab="Time(days)"
)
#
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group.three, data = kaps.td )
kthree=ggsurvplot(fit2, data =kaps.td,
                 #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                 pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                 surv.median.line = "hv",            # Add median survival lines
                 legend.title = "kaps.three.group",               # Change legend titles
                 legend.labs = c("set1","set2","set3"),  # Change legend labels
                 #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                 #palette =c("#00AFBB","#756bb1","#E69F00"),
                 palette =c("#1b9e77","#d95f02","#7570b3"),
                 #palette = c("#ca0020","#0571b0","#4daf4a"),
                 risk.table = TRUE,                  # Add No at risk table
                 cumevents = TRUE,                   # Add cumulative No of events table
                 tables.height = 0.15,               # Specify tables height
                 tables.theme = theme_cleantable(),  # Clean theme for tables
                 tables.y.text = FALSE,             # Hide tables y axis text
                 #ylab="Disease specific survival probability",
                 ylab="Overall survival probability",
                 #ylab="progression-free interval probability",
                 #ylab="disease-free interval probability",
                 xlab="Time(days)"
)
#
fit2<- survfit(Surv(OS.time, OS) ~ kaps.group.two, data = aa )
ktwo=ggsurvplot(fit2, data =aa,
                 #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","\n","logrank","P = ", pvalue, sep=""),
                 pval = TRUE, pval.method = TRUE, # Add p-value &  method name
                 surv.median.line = "hv",            # Add median survival lines
                 legend.title = "kaps.two.group",               # Change legend titles
                 legend.labs = c("set1","set2"),  # Change legend labels
                 #palette = c("#00C5CD","tomato2"),  # Use JCO journal color palette,
                 #palette =c("#00AFBB","#756bb1"),
                 palette =c("#1b9e77","#d95f02"),
                 #palette = c("#ca0020","#0571b0","#4daf4a"),
                 risk.table = TRUE,                  # Add No at risk table
                 cumevents = TRUE,                   # Add cumulative No of events table
                 tables.height = 0.15,               # Specify tables height
                 tables.theme = theme_cleantable(),  # Clean theme for tables
                 tables.y.text = FALSE,             # Hide tables y axis text
                 #ylab="Disease specific survival probability",
                 ylab="Overall survival probability",
                 #ylab="progression-free interval probability",
                 #ylab="disease-free interval probability",
                 xlab="Time(days)"
)
generate.PDF <- function(pbmc) {
  #pdf("bar.top.200.of.781.tes.counts.geneset.pdf",width = 10,height = 5)
  pdf("4.model/kaps/kaps.k.selection/survival.plots.pdf",width = 5,height = 6)
  print(kfour)
  print(kthree)
  print(ktwo)
  dev.off()
}
generate.PDF(pbmc)
#####
#####
outputdata.kaps=as.data.frame(t(fitkaps@test.stat))
outputdata.kaps$pvaluesofXk=c(8e-04,3e-033, 3e-03)
colnames(outputdata.kaps)[2]="pvaluesofX1"
outputdata.kaps=outputdata.kaps[,-1]
outputdata.kaps$group=c("kaps.two.group","kaps.three.group","kaps.four.group")
outputdata.kaps$group=factor(outputdata.kaps$group,levels = c("kaps.two.group","kaps.three.group","kaps.four.group"))
write.csv(outputdata.kaps, "4.model/kaps/kaps.k.selection/p.adjp.value.csv")
#
pdf("4.model/kaps/kaps.k.selection/p.adjp.value.pdf",width = 5,height = 4)
ggplot(outputdata.kaps, aes(x=group, y=pvaluesofX1,group=1)) +
  geom_line(color="red",linetype = "dashed")+
  geom_point(shape=1, fill="blue", color="blue", size=3)+theme_classic()+
  geom_hline(yintercept = 0.05,size=0.5,linetype = "dashed")

ggplot(outputdata.kaps, aes(x=group, y=pvaluesofXk,group=1)) +
  geom_line(color="red",linetype = "dashed")+
  geom_point(shape=1, fill="blue", color="blue", size=3)+theme_classic()+
  geom_hline(yintercept = 0.05,size=0.5,linetype = "dashed")
dev.off()
####
aaa=fitkaps@data
aaa$status=ifelse(aaa$OS=="1","Event","Censored")
write.csv(aaa, "4.model/kaps/kaps.k.selection/survival.os.time.event.csv")
#
pdf("4.model/kaps/kaps.k.selection/survival.os.time.event.pdf",width = 5,height = 4)
ggplot(aaa, aes(x=z.of.mean.exp, y=OS.time, group=status)) +
  geom_point(aes(shape=status, color=status))+
  scale_shape_manual(values=c(3, 1))+
  scale_color_manual(values=c("red","blue"))+ylab("OS.time (days)")+xlab("TE score")+theme_classic()
dev.off()