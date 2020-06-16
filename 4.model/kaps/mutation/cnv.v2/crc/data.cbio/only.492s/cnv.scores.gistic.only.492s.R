
#################################
library(BSgenome.Hsapiens.UCSC.hg19)
options(stringsAsFactors = FALSE)
#
### segment information
#tumor.seg.cnv <- read.table("4.model/kaps/mutation/cnv.v2/exm/gdac.broadinstitute.org_ESCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/focal_input.seg.txt", sep="\t", header=T, stringsAsFactors=F)
#####
tumor.seg.cnv <- read.table("4.model/kaps/mutation/snv/data_cna_hg19.seg", sep="\t", header=T, stringsAsFactors=F)
tail(tumor.seg.cnv)
colnames(tumor.seg.cnv)=c("Sample", "Chromosome", "Start.bp",   "End.bp", "Num.Markers","Seg.CN")
#tumor.seg.cnv$Sample <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", tumor.seg.cnv$Sample)
tail(tumor.seg.cnv)
### subset segment file based on kaps group
####only select matched on cbioportal
cbiosample=read.table("4.model/kaps/mutation/tcga.maf/maf.with.kaps.group.clean/clin.maf.with.kaps.group.492s.txt",header = T,row.names = 1)
head(cbiosample)
dim(cbiosample)
length(intersect(rownames(cbiosample), tumor.seg.cnv$Sample))
#
# Get subtype segments information
set1.seg.cnv <- tumor.seg.cnv[tumor.seg.cnv$Sample %in% rownames(subset(cbiosample, kaps.group=="set1")),]
set2.seg.cnv <- tumor.seg.cnv[tumor.seg.cnv$Sample %in% rownames(subset(cbiosample, kaps.group=="set2")),]
set3.seg.cnv <- tumor.seg.cnv[tumor.seg.cnv$Sample %in% rownames(subset(cbiosample, kaps.group=="set3")),]
set4.seg.cnv <- tumor.seg.cnv[tumor.seg.cnv$Sample %in% rownames(subset(cbiosample, kaps.group=="set4")),]
head(set1.seg.cnv)

write.table(set1.seg.cnv, file="4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/tumor.set1.seg.txt", sep="\t", row.names=F, quote = F)
write.table(set2.seg.cnv, file="4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/tumor.set2.seg.txt", sep="\t", row.names=F, quote = F)
write.table(set3.seg.cnv, file="4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/tumor.set3.seg.txt", sep="\t", row.names=F, quote = F)
write.table(set4.seg.cnv, file="4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/tumor.set4.seg.txt", sep="\t", row.names=F, quote = F)
length(unique(set4.seg.cnv$Sample))
######use segment file into GISTIC2.0 to run and calculate gistic score
######################
######################
# Create a chromosomes reference objects function
chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}

# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)
#str(chrom)
###############
################
#画四个subtype对比的gistic score

pdf("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/cnv.scores.gistic.pdf",18,12)
par(mfrow=c(4,1), mar = par()$mar + c(3,0,0,3))


### set1 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set1/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("set1 CRC copy number gistic score"," ","n=",length(unique(set1.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.05)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### set2 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set2/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("set2 CRC copy number gistic score"," ","n=",length(unique(set2.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### set3 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set3/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("set3 CRC copy number gistic score"," ","n=",length(unique(set3.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.05)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### set4 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set4/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("set4 CRC copy number gistic score"," ","n=",length(unique(set4.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.05)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
dev.off()
############################
############################
#########frequency plot
#############
pdf("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/cnv.frequence.pdf",18,12)
par(mfrow=c(4,1), mar = par()$mar + c(3,0,0,3))

### set1 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set1/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-0.1,max(scores$frequency)+0.1)
title=paste0("set1 CRC copy number frequency"," ","n=",length(unique(set1.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "Frequency", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### set2 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set2/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-0.1,max(scores$frequency)+0.1)
title=paste0("set2 CRC copy number frequency"," ","n=",length(unique(set2.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "Frequency", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### set3 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set3/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-0.1,max(scores$frequency)+0.1)
title=paste0("set3 CRC copy number frequency"," ","n=",length(unique(set3.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "Frequency", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### set4 ###
scores <- read.table("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set4/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-0.1,max(scores$frequency)+0.1)
title=paste0("set4 CRC copy number frequency"," ","n=",length(unique(set4.seg.cnv$Sample)))

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "Frequency", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
dev.off()