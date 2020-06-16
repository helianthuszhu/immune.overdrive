#!/bin/bash
sed -n '12,26p' SRR_Acc_List.GBM.cell.line.txt > bath26s
for sample in `cat bath26s`
do
parallel-fastq-dump --sra-id ${sample} --threads 14 --outdir ./sra.fastq.data --gzip
done
