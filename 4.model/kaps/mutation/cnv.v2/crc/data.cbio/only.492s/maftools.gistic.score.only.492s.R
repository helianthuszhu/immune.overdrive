library(maftools)
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set1/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set1/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set1/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set1/scores.gistic", isTCGA = TRUE)
######
pdf("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set1/mafout.set1.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()
#
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set2/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set2/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set2/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set2/scores.gistic", isTCGA = TRUE)
######
pdf("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set2/mafout.set2.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()
#
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set3/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set3/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set3/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set3/scores.gistic", isTCGA = TRUE)
######
pdf("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set3/mafout.set3.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()
#
laml.gistic = readGistic(gisticAllLesionsFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set4/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set4/amp_genes.conf_90.txt", 
                         gisticDelGenesFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set4/del_genes.conf_90.txt", 
                         gisticScoresFile = "4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set4/scores.gistic", isTCGA = TRUE)
######
pdf("4.model/kaps/mutation/cnv.v2/crc/data.cbio/only.492s/set4/mafout.set4.pdf",height = 5,width = 8)
gisticChromPlot(gistic = laml.gistic,markBands = "all",color = c("#ce1256","#2171b5"),ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 1,
                cytobandTxtSize = 0.5)
dev.off()