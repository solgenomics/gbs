#!/usr/local/bin/perl

###### HOW TO USE THIS SCRIPT #######################################################
# perl vcfFormatChange.gl <UnImputedVCFFile> <ImputedVCFfile> <outputFile>
#
# This script take the imput from unimputed VCF and imputed vcf files
# and adds the dosage information and INFO fileds are added to unimputed VCF file
#####################################################################################

$inputUnImp = $ARGV[0];
$inputImp = $ARGV[1];
$output = $ARGV[2];

$sym = "##"; @header =();
open (IUI, $inputUnImp);
while (<IUI>) {
    if ($_ =~ /^$sym/) { push @header, $_; }
    elsif ($_ =~ /!^$sym/) { last; }
}
close(IUI);

open (IMP, $inputImp);
while(<IMP>) {
    if ($_ =~ /ID=DS/) { push @header, $_; }
    elsif ($_ =~ /REF/) {
        push @header, $_;
        $samplesInfo = $_;
        last; }
}
close(IMP);

@nonImpArray = (); @imputedArray = ();
$start = "#";

open (IUI, $inputUnImp);
while (<IUI>) {
    if ($_ =~ /^$start/) { next; }
    else {
        chomp;
        push @nonImpArray, [split /\t/];
    }
}
open (IMP, $inputImp);
while (<IMP>) {
    if ($_ =~ /^$start/) { next; }
    else {
        chomp;
        push @imputedArray, [split /\t/];
    }
}

open (RES, ">$output");
for ($a=0; $a<=$#header; $a++) {
    print RES $header[$a];
}
@samples = split(/\t/, $samplesInfo);
#print "Samples = $#samples\n";

for ($b=0; $b<=$#imputedArray; $b++) {
    for ($c=0; $c<=7; $c++) {
        print RES "$imputedArray[$b][$c]\t";
    }
    print RES "GT:AD:DP:GQ:DS:PH:PL\t";
    
    for ($c=9; $c<=$#samples; $c++) {
        my @tmp1 = split(/:/, $imputedArray[$b][$c]);
        my $phase = $tmp1[0];
	my $dosage = $tmp1[1];
	
        my @tmp2 = split(/:/, $nonImpArray[$b][$c]);
        print RES "$tmp2[0]", ":", "$tmp2[1]", ":", "$tmp2[2]", ":", "$tmp2[3]", ":", "$dosage",":","phase",":", "$tmp2[4]\t";
    }
    print RES "\n";
}
