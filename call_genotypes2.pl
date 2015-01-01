

=head1 NAME

call_genotypes2.pl - a script to filter genotypes from GBS experiments

=head1 SYNOPSYS

perl call_genotypes2.pl --genotype_min_good_scores 0.2 --verbose --min_maf 0.001  min_scored_marker_fraction 0.2 --max_chi 30 --min_heterozygote_count 3  --infile snp_file.vcf --outfile output.dosage

All parameters are optional, except the infile and the outfile parameters.

=head1 DESCRIPTION

The scripts implements the following checks.

First, the vcf file is parsed on a 'per row' basis and statistics about the fraction of good snp calls per genotype is calculated. If that value is below genotype_min_good_scores (default 0.2), the corresponding genotype is excluded from further analysis.

Then the vcf file is parsed on a 'per line' basis and checks on a snp population basis are applied. 

=over 3

=item * 

If the minor allele frequency is below min_maf (default 0.001), the SNP is excluded

=item *

If the SNP is monomorphic, it is excluded

=item * 

If the SNP does not correspond a Hardy Weinberg distribution with a chi square value of more than max_chi (default 30), it is excluded.

=item *

If the SNP 

=back

The SNP dosage is calculated based on the specific SNP counts and an assumed error rate.

=head2 OUTPUT FILE

The results are output in the following tab delimited format:

 header line: 
 column 1 : SNP ID
 columns 2-N: Genotype IDs

 Other lines:
 column 1 : the SNP ids
 columns 2-N: The dosage values

A stats file is also output, with the following columns:

 SNP ID
 monomorphic or polymorphic
 scored_marker_fraction
 heterozygote_count
 chi square value (Hardy Weinberg) 


=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut

use strict;

use Data::Dumper;
use Getopt::Long;

use CXGN::Genotype;
use CXGN::GenotypeIO;
use CXGN::Genotype::SNP;
use CXGN::SNPs;
use CXGN::SNPsIO;

# command line options
#
my $genotype_min_good_scores = 0.2;
my $verbose = 0;
my $min_maf = 0.001;
my $min_scored_marker_fraction = 0.2;
my $max_chi = 30;
my $min_heterozygote_count = 3;
my $infile;
my $outfile;

GetOptions("genotype_min_good_scores=f"=> \$genotype_min_good_scores,
	   "verbose" => \$verbose,
	   "min_maf=f" => \$min_maf,
	   "min_scored_marker_fraction=f"=> \$min_scored_marker_fraction,
	   "max_chi=f" => \$max_chi,
	   "min_heterozygote_count=i" => \$min_heterozygote_count,
	   "infile=s" => \$infile,
	   "outfile=s" => \$outfile,
    );

if (!$outfile) { 
    die "Output file (option --outfile) required.\n";
}

if (!$infile) { 
    die "Input file (option --infile) is required.\n";
}

my $stats = $outfile.".stats";

open(my $OUT, ">", $outfile) || die "Can't open file for writing $outfile\n";
open(my $STATS, ">", $stats) || die "Can't open file for writing stats $stats\n";

print STDERR "\nRUN SUMMARY\n";
print STDERR "***********\n";
print STDERR "Processing file: $infile\n";
print STDERR "Output in file: $outfile\n";
print STDERR "Stats in file: $stats\n";
print STDERR "genotype_min_good_scores: $genotype_min_good_scores\n";
print STDERR "min_maf: $min_maf\n\n";

my $gtio = CXGN::GenotypeIO->new({ file => $infile, format => 'vcf' });

my $acc_count = $gtio->count();
#print STDERR "TOTAL ACC: $acc_count\n";
my %genotype_info = ();
my %genotype_problems = ();
my @valid_accessions = ();

message("Gathering genotype stats...");
my $stats = $gtio->summary_stats();

message(" Done.");
#print STDERR Dumper($stats);

message("Generating a list of valid accessions...");
foreach my $acc (keys %$stats) { 
    my $valid =1;

    if ($acc =~ /failed|empty|blank|NA\:d+/) {
	push @{$genotype_problems{$acc}}, "empty or blank";
	$valid =0;
	$acc_count--;
    }
    my $good_score_fraction = $stats->{$acc} / $acc_count;
    if ( $good_score_fraction < $genotype_min_good_scores) { 
	print STDERR "$good_score_fraction is too low ($genotype_min_good_scores)\n";
	push @{$genotype_problems{$acc}}, "failed min_good_scores";
	$valid =0;
    }
    
    if ($valid && $acc) { 
	push @valid_accessions, $acc;
    }

}
message(" Done.");

$gtio->close();

message("Parsing SNPS...\n");

my $snps_io = CXGN::SNPsIO->new( { file => $infile });
$snps_io->ignore_accessions(\%genotype_problems);	
$snps_io->valid_accessions(\@valid_accessions);

print STDERR "VALID ACCESSIONS: ".Dumper(\@valid_accessions);

my $total_clone_count = scalar(@{$snps_io->valid_accessions()});
print STDERR "Total valid accessions: $total_clone_count\n";
# print header
#
print $OUT "SNP\t";
print $OUT join ("\t", @valid_accessions);
print $OUT "\n";

while (my $snps = $snps_io->next()) { 
    my $snp_id = $snps->id();
    message("ANALYZING SNP $snp_id...\n");

    my $skip = 0;

    $snps->ignore_accessions($snps_io->ignore_accessions());
    $snps->valid_accessions(\@valid_accessions);

    #message(join ",", $snps->snp_stats());
    
    my $allele_freq = $snps->calculate_allele_frequency_using_counts();
    message("ALLELE FREQ: $allele_freq\n");
    if ($allele_freq < $min_maf || $allele_freq > (1-$min_maf)) { 
	print STDERR "Ignoring snp ".$snps->id()."\n";
	$skip = 1;
    }
    
    my $valid_clone_count = scalar(@{$snps->accessions}) - scalar(keys(%{$snps_io->ignore_accessions}));
    #print STDERR "CLONE COUNT: $valid_clone_count\n";

    my @dosages;
    foreach my $k (@valid_accessions) { 
	my $s = $snps->snps()->{$k};
	#print STDERR "ACCESSION: ".$s->accession()."\n";
	
	my $dosage = $snps->calculate_snp_dosage($s);
	#print STDERR "DOSAGE: $dosage\n";
	push @dosages, $dosage;
    }
    


    my $lowN;
    if (( $allele_freq < 0.01) || ($allele_freq > 0.99)) { 
	$lowN = $total_clone_count * 0.4;
    }
    elsif ( ( $allele_freq < 0.1) || ( $allele_freq > 0.9)) { 
	$lowN = $total_clone_count * 0.3;
    }
    else {    #if (($allele_freq > 0.1) || ($allele_freq < 0.9)) { 
	$lowN = $total_clone_count * 0.2;
    }
    
    if ($snps->depth() < $lowN) { 
	print STDERR "Skipping $snp_id because of low depth (".$snps->depth()." vs. minimum of $lowN).\n";
	$skip = 1;
    }
    else { 
	message("Keeping $snp_id becasue depth is ok (".$snps->depth()." vs. minimum of $lowN).\n");
    }

    my %score = $snps->hardy_weinberg_filter(\@dosages);
    print STDERR Dumper(\%score);
    if (exists($score{monomorphic})) { 
	message("Skipping $snp_id because it is monomorphic\n");
	$skip = 1;
    }
    if (exists($score{scored_marker_fraction})){ 
	if ($score{scored_marker_fraction} < $min_scored_marker_fraction) { 
	    message("Skipping $snp_id because it is scored marker fraction ($score{scored_marker_fraction}) is below the minimum of $min_scored_marker_fraction\n");
	    $skip=1;
	}
	if ($score{heterozygote_count} < $min_heterozygote_count) { 
	    message("Skipping $snp_id because of low heterozygote count ($score{heterozygote_count})\n");
	    $skip=1;
	}
	if ($score{chi} > $max_chi) { 
	    message("Skipping $snp_id because of chi value ($score{chi})\n");
	    $skip=1;
	}
    }
    
    printf($STATS "$snp_id\t%.3f\t%.4f\t%.1f\n", ($snp_id, $score{monomorphic} ? 'monomorphic' : 'polymorphic',$score{scored_marker_fraction}, $score{heterozygote_count}, $score{chi}));

    if ($skip) { 
	message("SKIPPING!\n");
	next();
    }
    else { 
	message("SNP $snp_id IS OK! OUTPUTTING...\n");
	print $OUT $snp_id."\t";
	foreach my $d (@dosages) { printf ($OUT "\t%.2f", $d); }
	print $OUT "\n";
    }
		
}

close($OUT);
close($STATS);

sub message { 
    my $message = shift;
    if ($verbose) { 
	print STDERR $message;
    }
}
