#!/bin/bash

# script to run the imputation R script on all chromosomes.
#
DIRNAME=$1;

for i in $(seq 1 4); do
    echo "Launching analysis for chr $i...";
    Rscript /home/mueller/cxgn/gbs/ImputeWithGlmnetScript.R $DIRNAME/cassava2_chr19.$i.snps.txt 2> $DIRNAME/cassava2_chr19.$i_imputation.err &
done