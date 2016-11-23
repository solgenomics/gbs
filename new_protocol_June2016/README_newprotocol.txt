Imputation protocol June 2016

gjb99@cornell.edu

#Tools:
Beagle v4.0: https://faculty.washington.edu/browning/beagle/b3.html
PLINK/Seq v101: https://atgu.mgh.harvard.edu/plinkseq/pseq.shtml
Vcftools v1.14: https://vcftools.github.io/

#Scripts:
transpose_matrix.pl
filtering_AR2.sh
load_genotypes.pl  

# 1 Imputing without a reference panel:

-a Pre-filtering on the entiredataset to keep only bi-allelic markers, remove marker and individuals with >80% missing data:

vcftools –gzvcf REFPAN.vcf.gz --min-alleles 2 --max-alleles 2 --recode --max-missing 0.2 --out REF_PANEL.filter

-b list individuals above 80% missing data:

vcftools –vcf REF_PANEL.filter.recode.vcf --missing-indv --out REF_PANEL.filter.recode

-c list individuals above 80%:
 
cat REF_PANEL.filter.recode.imiss | tail -n +2 | awk '{if($5==$5+0 && $5>0.8)print $1}' | sort -u > list_outliers.txt

-d discard outliers from dataset:

vcftools –gzvcf chr$i.vcf.gz --keep list_outliers.txt --min-alleles 2 --max-alleles 2 --recode --out REF_PANEL.final


-b Imputing with a reference panel:

java -Xmx8000m -jar beagle.r1399.jar 
gl= REF_PANEL.final.vcf ibd=false impute=true
out= IITA_CYCLE.imputed.vcf 
nthreads=24 &




# 2a Imputing cycles (CYC) using the reference panel

-a Pre-filtering on the entire dataset to keep only bi-allelic markers, remove individuals with more than 80% missing data per chromosome

vcftools –gzvcf CYCLE.vcf.gz --min-alleles 2 --max-alleles 2 --recode --max-missing 0.2 --out CYCLE.filter

java -Xmx8000m -jar beagle.r1399.jar 
ref=REFPAN.vcf.gz  
gl= CYCLE.filter.vcf.gz
ibd=false 
impute=true 
out= IITA_CYCLE.imputed
nthreads=24

# 3 filtering for allelic correlations AR2 >0.9:

bash filterBeagleImpute.sh CYCLE.imputed.vcf.gz CYCLE.imputed09

# 4 Convert to dosage data:

pseq CYCLE.imputed09.vcf.gz v-meta-matrix --name DS &> CYCLE.dosage.txt

# 5 transpose dosage matrix:

perl transpose_matrix.pl CYCLE.dosage.txt > CYCLE.dosage.tp.txt

# 6 Load data to cassavabase, ex:

perl bin/load_genotypes.pl -H db4.sgn.cornell.edu  -D devel_cassava -a  -i ~/../gjb99/testout_tp  -p "SNP genotyping 2016 Cornell Biotech" -y 2016  -g "global_population_set_2012-2016" -m protocol "GBS ApeKI Cassava genome v6_june2016" 2> test_load_st_err.log
