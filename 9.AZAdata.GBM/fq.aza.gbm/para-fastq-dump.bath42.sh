#!/bin/bash
sed -n '27,42p' SRR_Acc_List.GBM.cell.line.txt > bath42s
for sample in `cat bath42s`
do
parallel-fastq-dump --sra-id ${sample} --threads 14 --outdir ./sra.fastq.data --gzip
done
