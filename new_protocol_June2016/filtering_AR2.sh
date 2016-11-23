#!/bin/bash
# Script to filter a BEAGLE-imputed VCF file on AR2 field
# Syntax: bash filterBeagleImpute.sh <vcf_filename w/out .vcf.gz> <output_filename>

set -e
#TEMPDIR='filter_AR2_BEAGLE_imputed_chromosomes_'`date +%Y-%m-%d`;
TEMPDIR="FILT_AR2_TEMP_0.3"
mkdir $TEMPDIR

for i in $(seq 1 18);do
    echo "1. Extracting AR2"
    zcat /home/gjb99/TP/TP_glmode_NEW/test_multiall_filter_igdBuildWithV6_hapmap_20151008_withRef_filt_ArielPPL_TP-CYC_chr$i.NEW.edit.imputedvcf.gz.vcf.gz  | tail -n +10 |cut -f8 |cut -d';' -f1 | cut -d'=' -f2 | sed '1d' > $TEMPDIR/temp.AR2.chr$i

    echo "2. Removing VCF header"
    zcat /home/gjb99/TP/TP_glmode_NEW/test_multiall_filter_igdBuildWithV6_hapmap_20151008_withRef_filt_ArielPPL_TP-CYC_chr$i.NEW.edit.imputedvcf.gz.vcf.gz  |grep -v ^# > $TEMPDIR/temp_nohead.chr$i

    echo "3. Adding AR2 values to original dataframe"
    paste $TEMPDIR/temp.AR2.chr$i $TEMPDIR/temp_nohead.chr$i > $TEMPDIR/temp2.chr$i

    echo "delete old tempfiles"
    rm $TEMPDIR/temp.AR2.chr$i
    rm $TEMPDIR/temp_nohead.chr$i

    echo "4. Filtering out SNPs with AR2 below 0.5"
    awk '$1 >= 0.3' $TEMPDIR/temp2.chr$i |cut -f2- > $TEMPDIR/temp3.chr$i

    echo "5. Extracting original header"
    zcat /home/gjb99/TP/TP_glmode_NEW/test_multiall_filter_igdBuildWithV6_hapmap_20151008_withRef_filt_ArielPPL_TP-CYC_chr$i.NEW.edit.imputedvcf.gz.vcf.gz  |grep ^# > $TEMPDIR/temp.header.chr$i

    echo "6. Pasting into new filtered VCF file"
    cat $TEMPDIR/temp.header.chr$i $TEMPDIR/temp3.chr$i > $TEMPDIR/beagle_onTP_from_rob-ariel_filter_glmode_chr$i.AR2filt.vcf

    echo "delete old tempfile"
    rm $TEMPDIR/temp3.chr$i
    rm $TEMPDIR/temp.header.chr$i

    echo "7. Recoding using VCFTools"
    vcftools --gzvcf $TEMPDIR/beagle_onTP_from_rob-ariel_filter_glmode_chr$i.AR2filt.vcf --recode --recode-INFO-all --stdout | gzip -c > $TEMPDIR/beagle_onTP_from_rob-ariel_filter_glmode_chr$i.AR2filt.vcf.gz
    
    #echo "8. Running tabix"
    #bgzip beagle_onTP_gt_chr$i.AR2filt.vcf
    #tabix -p vcf beagle_onTP_gt_chr$i.AR2filt.vcf.gz
    
    #echo "Some cleanup"
    #rm $TEMPDIR/temp4.chr$i.vcf.gz

#echo "8. Converting to PLINK tped format using vcftools"
#vcftools --vcf $TEMPDIR/temp4.snps_only.$1.recode.vcf --plink-tped --out $TEMPDIR/temp.$2

#echo "Some cleanup"
#rm $TEMPDIR/temp4.snps_only.$1.recode.vcf

#echo "9. Converting to binary PLINK format"
#plink --noweb --file $TEMPDIR/temp.$2 --make-bed --out $2

#echo "Some cleanup"
#rm $TEMPDIR/temp.$2*

echo "Done"
done
