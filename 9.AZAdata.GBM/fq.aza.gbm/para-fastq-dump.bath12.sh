#!/bin/bash
sed -n '1,12p' SRR_Acc_List.GBM.cell.line.txt > bath12s
for sample in `cat bath12s`
do
parallel-fastq-dump --sra-id ${sample} --threads 14 --outdir ./sra.fastq.data --gzip
done
