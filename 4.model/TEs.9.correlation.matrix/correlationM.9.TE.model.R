################correlation matrix of these 9 tes
library(GGally)
dim(dat)
ggg1=ggcorr(da1,nbreaks = 4,geom = "circle",method = c("pairwise", "spearman"),name = "spearman coefs",
            min_size=4,max_size = 8,
            low = "#6baed6", mid = "white", high = "#dd1c77",legend.size = 9)
generate.PDF <- function(fig) {
  pdf("4.model/TEs.9.correlation.matrix/correlationM.9.TE.model.pdf",height  = 5,width = 5)
  print(ggg1)
  dev.off()
}
generate.PDF(fig)
save(da1, file = "4.model/TEs.9.correlation.matrix/correlationM.9.TE.model.RData")