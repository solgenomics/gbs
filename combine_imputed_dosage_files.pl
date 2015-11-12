#!/usr/bin/perl                                                                                                                                                                       
use warnings;
use strict;

my $temp_dir = "dosage_file_processing";
my $common_accessions = $temp_dir . "/" . "uniquenames.txt";
my $placeholder = "SNP";
my ($dosage_file_dir, @values);


print STDERR "Making dir for temp files . . .\n";
`mkdir $temp_dir`;
 
print STDERR "Getting header with accession names from each chromosome dosage file . . .\n";
`seq 1 20 | parallel 'head -n 1 cassava_chr{}.snps.txt.imputed.txt > dosage_file_processing/headerch{}.imp.txt'`; 

print STDERR "Transposing them to have one accession name per line . . .\n";
`seq 1 20 | parallel './transpose.pl dosage_file_processing/headerch{}.imp.txt > dosage_file_processing/headerch{}.imp.txt.t'`;

print STDERR "Putting all the names from each chromosome in one file . . .\n";
`cat \$(ls dosage_file_processing/*.imp.txt.t | sort -V) >> dosage_file_processing/allnames.txt`;

print STDERR "Getting a list of the unique accessions that are present in at least one chromosome dosage file . . .\n";
`sort -u $temp_dir/allnames.txt > $common_accessions`;

print STDERR "Putting SNPS placeholder in header line so that we can transpose . . .\n";
`seq 1 20 | parallel 'head -n 1 cassava_chr{}.snps.txt.imputed.txt | grep "^SNP" -q && echo "placeholder found in chrom {} header" || sed -i "1s/^/SNP\t/" cassava_chr{}.snps.txt.imputed.txt'`;

print STDERR "Transposing to have one accession per row . . . this will take a few minutes . . .\n";
`parallel -j 6 './transpose.pl {} > dosage_file_processing/{.}.t' ::: *.imputed.txt`;

for my $dosage_file (glob("$temp_dir/*.imputed.t")) {
    print STDERR "Storing dosage file $dosage_file individual accession data in individual files . . . \n";
    &index_by_accession($dosage_file); 
    print STDERR "Constructing new dosage file for this chrom with the number and order of accessions the same as other chroms . . .\n";
    &rebuild_with_common_accessions($dosage_file, $common_accessions); 
};



print STDERR "Pasting dosage files together columnwise to create merged dosage file with snps as columns . . . this will take a few minutes . . .\n"; 
`paste -d '\t' \$(ls dosage_file_processing/*imputed.t.rebuilt | sort -V) > cassava_allchr.snps.txt.imputed.snpcolumns`;

print STDERR "Transposing dosage files back to have one snp per row . . . this will take a few minutes . . .\n"; 
`parallel -j 6 './transpose.pl {} > {.}.cat' ::: $temp_dir/*.t.rebuilt`;

print STDERR "Concatenating transposed, rebuilt dosage files to create merged dosage file with accessions as columns . . .\n";
for $i (1..20) {`cat dosage_file_processing/cassava_chr$i.snps.txt.imputed.t.cat >> cassava_allchr.snps.txt.imputed.accessioncolumns`};

print STDERR "Finished!\n";


sub index_by_accession() {
    
    my $dosage_file = shift;
    
    open (LINES, "<", $dosage_file) || die "Can't open dosage file $dosage_file!\n";

    $dosage_file_dir = $dosage_file . "_dir";

    `mkdir $dosage_file_dir`;

    while (<LINES>) {
	chomp(my $accession = $_ );

	my($accession_name, @values) = split (/\t/ , $accession);      # extract accession name to make temp file                                                                     
	
	$accession_name =~ s/^.+:(\d+)$/$1/;  #use tassel identifier number for file name rather than file accession name in case it contains special characters

	open (TEMP, ">", "$dosage_file_dir/$accession_name") || die "Can't open temp file $accession_name!\n";
	
	if ($dosage_file =~ /chr1.snps/) {
	    print TEMP $accession . "\n"
	} else {
	    print TEMP join("\t", @values); 
	    print TEMP "\n";
    }
    }
}

sub rebuild_with_common_accessions() {

    my $dosage_file = shift;
    my $common_accessions = shift;
    my $rebuilt_dosage_file = $dosage_file . ".rebuilt";

    open (NEW, ">>", $rebuilt_dosage_file) || die "Can't open output file $rebuilt_dosage_file!\n";

    open (NAMES, "<", $common_accessions) || die "Can't open common accessions file $common_accessions!\n";

# create array of NAs to fill in lines of accessions that weren't imputed                                                                                                               
    my @filler = map {$_ = 'NA'} @values;

# go through accessions in order, and print their info to rebuilt dosage file                                                                         

    `cat $dosage_file_dir/$placeholder >> $rebuilt_dosage_file`;

    while (<NAMES>) {
	chomp(my $name = $_);
	my $file = $name;
	$file =~ s/^.+:(\d+)$/$1/;
	if (-e "$dosage_file_dir/$file") {
	    `cat $dosage_file_dir/$file >> $rebuilt_dosage_file`;
	} else {
	    if ($dosage_file =~ /chr1.snps/) {print NEW $name . "\t"};
	    print NEW join ("\t", @filler);          # print filler NAs for missing values                                                                                                  
	    print NEW "\n";
	}
    }

    `rm -r $dosage_file_dir`;
    print "$dosage_file rebuilt!\n";    
}
