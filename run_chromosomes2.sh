#!/bin/bash

DIRNAME='run_chromosomes_out_'`date +%Y-%m-%d`;
mkdir $DIRNAME

for i in $(seq 1 19); do
perl /home/mueller/cxgn/gbs/call_genotypes2.pl --verbose --outfile $DIRNAME/cassava2_chr$i.snps.txt --infile /home/mueller/cassava_vcf_files/cassava_mAF0.001MnSCov40_withDepth_IGD_c$i.vcf 2> $DIRNAME/cassava_chr$i.err &
done