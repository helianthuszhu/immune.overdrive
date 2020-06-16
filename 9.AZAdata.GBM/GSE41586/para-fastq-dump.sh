#!/bin/bash

for sample in `cat SRR_Acc_List.GSE41586.txt`
do
parallel-fastq-dump --sra-id ${sample} --threads 14 --outdir ./sra.fastq.data --split-files --gzip
done
