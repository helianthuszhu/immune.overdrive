#!/bin/bash
ls /home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/9.AZAdata.GBM/GSE41586/sra.fastq.data/*_1.fastq.gz > 1
ls /home/zhuxq/nas/Xiaoqiang/opti.data/1.CGC.TE.CRC.REdiscoverTE/stats.20200311/9.AZAdata.GBM/GSE41586/sra.fastq.data/*_2.fastq.gz > 2
cut -d"/" -f 12 1 |cut -d"_" -f 1  > 0
paste 0 1 2 > config
cat config |while read id
do
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
echo $sample $fq1 $fq2
bbduk.sh -Xmx70g in1=$fq1 \
 in2=$fq2 \
 out1=./fastq.trimed/$sample.clean_1.fastq.gz \
 out2=./fastq.trimed/$sample.clean_2.fastq.gz \
 ref=/home/bioinf/tools/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20
done
