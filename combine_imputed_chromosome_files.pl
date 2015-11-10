#!/usr/bin/perl

use warnings;
use strict;

my $dosage_file = $ARGV[0];
my $vcf_accessions = $ARGV[1];

my $fixed_dosage_file = $dosage_file . ".fixed";
my $dosage_file_dir = $dosage_file . ".dir";

my $accession;
my $accession_name;
my @values; 

# remove X that Rscript prepended to accession names

`sed -i 's/^X//g' $dosage_file`;

# read dosage file line by line and make directory for temporary output

open (LINES, "<", $dosage_file) || die "Can't open dosage file $dosage_file!\n";

`mkdir $dosage_file_dir`;

while (<LINES>) {
    if ($. == 1) {
	next;
    }
    chomp( $accession = $_ );                                   # change back - and : characters that Rscript changed to .

    if ($accession =~ m/([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)/) {
        $accession =~ s/([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)/$1-$2:$3/;
    } else {
        $accession =~ s/([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)/$1:$2/;
    }

    if ($accession =~ m/O0214:250251113/) {                     # fix the 2 accessions that actually should have an X at the start of their names
        $accession =~ s/O0214:250251113/XO0214:250251113/;
    } elsif ($accession =~ m/O0214:250399878/) {
        $accession =~ s/O0214:250399878/XO0214:250399878/;
    } else {
    }

    ($accession_name , @values) = split (/\t/ , $accession);      # extract accession name to make temp file

    #print "Accession name = $accession_name \n";
    #print "First value = $values[0] \n";

    open (TEMP, ">", "$dosage_file_dir/$accession_name") || die "Can't open temp file $accession_name!\n";     # save temp file	
    print TEMP $accession . "\n";
}

# create file for fixed output including header with snp names

`head -n 1 $dosage_file > $fixed_dosage_file`;

open (NEW, ">>", $fixed_dosage_file) || die "Can't open output file $fixed_dosage_file!\n";

open (VCF, "<", $vcf_accessions) || die "Can't open vcf accessions file $vcf_accessions!\n";

# create array of NAs to fill in lines of accessions that weren't imputed

my @filler = map {$_ = 'NA'} @values;

# go through accessions in order they are found in filtered_vcf_file, and print their info to fixed dosage file

while (<VCF>) {
    chomp(my $file = $_);                         
    if (-e "$dosage_file_dir/$file") {
	`cat $dosage_file_dir/$file >> $fixed_dosage_file`;
    } else {
	print NEW $file . "\t"; 
	print NEW join ("\t", @filler);          # print filler NAs for missing values
	print NEW "\n";
    }
}
    
`rm -r $dosage_file_dir`;
print "$dosage_file fixed!\n";
