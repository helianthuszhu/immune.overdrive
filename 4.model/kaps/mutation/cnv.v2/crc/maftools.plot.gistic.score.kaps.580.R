library(maftools)
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/set4.gisticout/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set4.gisticout/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set4.gisticout/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/set4.gisticout/scores.gistic", isTCGA = TRUE)
######
pdf("4.model/kaps/mutation/cnv.v2/crc/set4.gisticout/mafout.set4.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()
#
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/set3.gisticout/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set3.gisticout/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set3.gisticout/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/set3.gisticout/scores.gistic", isTCGA = TRUE)
pdf("4.model/kaps/mutation/cnv.v2/crc/set3.gisticout/mafout.set3.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()
#
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/set2.gisticout/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set2.gisticout/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set2.gisticout/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/set2.gisticout/scores.gistic", isTCGA = TRUE)
pdf("4.model/kaps/mutation/cnv.v2/crc/set2.gisticout/mafout.set2.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()
#
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/set1.gisticout/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set1.gisticout/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/set1.gisticout/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/set1.gisticout/scores.gistic", isTCGA = TRUE)
pdf("4.model/kaps/mutation/cnv.v2/crc/set1.gisticout/mafout.set1.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()