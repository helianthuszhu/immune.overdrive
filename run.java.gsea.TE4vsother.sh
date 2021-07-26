/usr/bin/java -cp /home/zhuxq/nas/Xiaoqiang/bio_tools/gsea/gsea-3.0.jar \
   -Xmx20000m xtools.gsea.Gsea \
   -gmx c2.all.v7.1.symbols.gmt \
   -res /home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/gsea20201103/TEcluster4VSothers/Result/Pathway_Analysis/GSEA/TE.group1.probe.gct \
   -cls /home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/gsea20201103/TEcluster4VSothers/Result/Pathway_Analysis/GSEA/TE.group1.phenotype.cls#TE_cluster4_versus_TE_cluster.other \
   -plot_top_x 1000 \
   -collapse false \
   -out /home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/4.model/kaps/gsea20201103/TEcluster4VSothers/output.TE4VSother.v3 \
   -rpt_label TEcluster4VSother
