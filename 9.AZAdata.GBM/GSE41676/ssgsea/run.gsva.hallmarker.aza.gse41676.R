#run ssgsea C7 gmt file
##############gsva
library(GSVA)
library(Biobase)
#library(readgmt)
###
#' @export read_gmt
read_gmt <- function(file, tidy = FALSE) {
  con <- file(file, "r")
  gmt_lines <- readLines(file, warn = FALSE)
  close(con)
  rlist <- purrr::map(gmt_lines, parse_gmt_lines)
  rlist_names <- purrr::map_chr(gmt_lines, get_gmt_names)
  names(rlist) <- rlist_names
  if (tidy) rlist <- tidy_gmt(rlist)
  return(rlist)
}

parse_gmt_lines <- function(gl) {
  gl <- unlist(stringr::str_split(gl, "\\\t"), recursive = FALSE)
  gl <- gl[3:length(gl)]
  return(gl)
}

get_gmt_names <- function(gl) {
  gl <- unlist(stringr::str_split(gl, "\\\t"))
  gl_name <- gl[[1]]
  return(gl_name)
}

sel_gmt=read_gmt("~/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/hall.marker/h.all.v7.0.symbols.gmt")
sets=sel_gmt
########
#panimmune=read.csv("/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/strom.fib.emt.ecm.sigs.csv",header = T)
#head(panimmune)
#setnames=unique(panimmune$SetName)
#length(setnames)
#datalist=list()
#for (i in 1: length(setnames)) {
#  aa=subset(panimmune,SetName==setnames[i])
#  datalist[[i]]=aa$Gene
#}
#names(datalist) <- setnames
#sets=datalist
#sets[1]
#run ssgsea
load("/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/9.AZAdata.GBM/GSE41676/azaexpm.gse41676.batch.removed.RData")
exprMatrix=combat_edata.azaexp.gse41676
gsva_matrix<- gsva(as.matrix(exprMatrix), sets,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix1<- t(scale(t(gsva_matrix)))
#gsva_matrix1[gsva_matrix1< -2] <- -2
#gsva_matrix1[gsva_matrix1>2] <- 2
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1)
nor_gsva_matrix1[1:4,1:4]
#e.matirx=(nor_gsva_matrix1-rowMeans(nor_gsva_matrix1))/apply(nor_gsva_matrix1,1,sd)
score.gsva=as.data.frame(t(nor_gsva_matrix1))
#stroma.score=score.gsva
#save(stroma.score,file="/home/zhuxq/nas/Xiaoqiang/opti.data/signatureVSimmune/CRC.immu.sigs.score/strom.fib.emt.ecm.sigs.RData")
genesets.hall.marker.score.aza.gse41676=score.gsva
save(genesets.hall.marker.score.aza.gse41676,file="/home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/9.AZAdata.GBM/GSE41676/ssgsea/genesets.hallmarker.score.aza.gse41676.RData")
