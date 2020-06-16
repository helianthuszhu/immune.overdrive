#!/bin/bash
ls sra.fastq.data/*.fastq.gz > 1
cut -d"/" -f 2 1 |cut -d"." -f 1  > 0
paste 0 1 > config
sed -n '1,20p' config > config20s
cat config20s |while read id
do
arr=($id)
fq1=${arr[1]}
sample=${arr[0]}
echo $sample $fq1
bbduk.sh -Xmx60g in=$fq1 \
 out=./fastq.trimed/$sample.clean.fastq.gz \
 ref=/home/bioinf/tools/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20
done
