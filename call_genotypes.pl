
use strict;

use Getopt::Long;

use lib '/home/mueller/cxgn/gbs/lib';

use StatFunc;

our $MIN_DEPTH= 0;
our $MAF = 0.001;
our $OUTFILE = "";

GetOptions("maf=i" => \$MAF,
	   "outfile=s" => \$OUTFILE,
    );

my $filename = shift;

if (! $OUTFILE) { 
    $OUTFILE = $filename.".dosage";
}

my $stats_file = $OUTFILE.".stats";
my @clone_names;
my @clone_info; 
my @bad_calls_per_clone;

print STDERR "Detecting bad clones...\n";

# Filter bad clones. Keep a list of clones that failed in @clone_info
#
my $snps_processed = 0;

open (my $F, "<", $filename) || die "Can't open file $filename\n";
my $ignore_clones_count = 0;

while(<$F>) { 
    chomp;
    if ($snps_processed % 100 == 0) { 
	print STDERR "Processing $snps_processed ...\n";
    }
    
    if (m/\#\#/) { 
	next;
    }
    
    if (/^#CHROM/ || /^CHROM/ ) { 
	my @header = split /\t/;
	@clone_names = @header[9.. $#header];	
    }

    my ($chr, $position, $snp_id, $ref_allele, $alt_allele, $qual, $filter, $info, $format,  @snps) = split /\t/;
    
    my $depth = 0;

    my @counts;
    my $ignore_clones_count = 0;

    for (my $n = 0; $n < @clone_names; $n++) { 
	
	if ($clone_names[$n] =~ /blank/i || $clone_names[$n] =~ /failed/i || $clone_names[$n] =~/empty/i || $clone_names[$n] =~ /NA\:\d+/i) { 
	    $clone_info[$n] = "B";  # remember this as a blank (do not count in stats)
	    $ignore_clones_count++;
	}
	else { 
	    $clone_info[$n] = "";
	    
	    # find allele frequencies by adding up the counts for each allele
	    #
	    my ($c1, $c2) = StatFunc::get_counts($snps[$n]);
	    if ($c1 + $c2 < 2) { 
		$bad_calls_per_clone[$n]++;
	    }
	}
    }

    $snps_processed++;
}

for( my $n = 0; $n< @bad_calls_per_clone; $n++) { 
    if ($bad_calls_per_clone[$n] >  $snps_processed * 0.95) { 
	$clone_info[$n] = 'F';
	print STDERR "Throwing away clone $clone_names[$n] (bad calls per clone = $bad_calls_per_clone[$n]\n";
    }
}

close($F);

$MIN_DEPTH = (@clone_names - $ignore_clones_count) * 1.5; 
print STDERR "Using a default of $MIN_DEPTH\n";

$snps_processed = 0;

open(my $F, "<", $filename) || die "Can't open file $filename\n";
open(my $OUT, ">", $OUTFILE) || die "Can't open outfile $OUTFILE\n";
open(my $STATS, ">", $stats_file) || die "Can't open stats file $stats_file\n";

my @clone_failed_count;

my %data;


while (<$F>) { 
    chomp;

    if ($snps_processed % 100 == 0) { 
	print STDERR "Processing $snps_processed ...\n";
    }

    $snps_processed++;

    if (m/\#\#/) { 
	next;
    }
    
    if (/^#CHROM/ || /^CHROM/ ) { 
	my @header = split /\t/;
	@clone_names = @header[9.. $#header];	
	print STDERR "Using MAF of $MAF\n";
	
	print $OUT "snps";

	print STDERR "Printing header information...\n";
	for (my $i=0; $i<@clone_names; $i++) { 
	    if (!$clone_info[$i]) { 
		print $OUT "\t".$clone_names[$i];  
	    }
	}
	print $OUT "\n";
    }
    
    

    else { 
	my ($chr, $position, $snp_id, $ref_allele, $alt_allele, $qual, $filter, $info, $format,  @snps) = split /\t/;
	
	print STDERR "SNP: $snp_id\t$ref_allele\t$alt_allele\t$qual\n";
	my $depth = 0;
	if ($info =~ /DP=(\d+)$/) { 
	    $depth = $1;
	    print STDERR "DEPTH: $depth\n";
	}
	
	my $total_c1 = 0;
	my $total_c2 = 0;
	
	if ($depth < $MIN_DEPTH) { 
	    print STDERR "Skipping $snp_id because of $depth < $MIN_DEPTH\n";
	    next();
	}

	my @counts;

	for (my $n = 0; $n < @clone_names; $n++) { 

	    # find allele frequencies by adding up the counts for each allele
	    #
	    my ($c1, $c2) = StatFunc::get_counts($snps[$n]);
	    $counts[$n] = [ $c1, $c2 ];
	    $total_c1 += $c1;
	    $total_c2 += $c2;
	    
	    # keep a sum of failed snps for each clone
	    #
	    if ($c1 == 0 && $c2 == 0) { 
		$clone_failed_count[$n]++; 
		####$ignore_clones_count++;
	    }
	}	    
	my $allele_freq = $total_c1 / ($total_c1 + $total_c2);
	
	if ($allele_freq > (1 - $MAF) || $allele_freq < $MAF) { 
	    print STDERR "Skipping ALLELE FREQ = $allele_freq\n";
	    next();
	}
	
	my $pAA = $allele_freq **2;
	my $pAB = $allele_freq * (1 - $allele_freq) * 2 ;
	my $pBB = (1 - $allele_freq) **2;
	
	my @dosages; 
	
	# calculate the genotype using the counts and some stats.
	#
	for (my $n = 0; $n < @clone_names; $n++) { 
	    if ($clone_info[$n] eq "") { # clone is ok (not monomorphic or blank)
		if ( ($counts[$n]->[0] + $counts[$n]->[1]) > 1 ) { 
		    $dosages[$n] = StatFunc::get_genotype($pAA, $pAB, $pBB, $counts[$n]->[0], $counts[$n]->[1]);
		}
		else { $dosages[$n] = undef; }
	    }
	}
	
	# hardy weinberg filter
	my %score = StatFunc::hardy_weinberg_filter(\@dosages, $ignore_clones_count);
	
	if (exists($score{invalid}) && $score{invalid} == 1) { 
	    print STDERR "Skipping monomorphic marker $snp_id\n";
	    next();
	}
	printf($STATS "$snp_id\t%.3f\t%.4f\t%.1f", $score{scored_marker_fraction}, $score{allele_freq}, $score{chi});
	
	if ($score{scored_marker_fraction} < 0.4) { 
	    print STDERR "Skipping $snp_id because of low scored markers ($score{scored_marker_fraction})\n";
	    print $STATS "\tSKIPPED\n";
		next();
	}
	#if (($score{allele_freq} > 1-$MAF) || ($score{allele_freq} < $MAF) ) { 
	#    print STDERR "Skipping $snp_id because of high allele frequency ($score{allele_freq})\n";
	#    print $STATS "\tSKIPPED\n";
	#    next();
	#}
        if ($score{heterozygote_count} < 3) { 
	    print STDERR "Skipping $snp_id because of low heterzygote count.\n";
	    print $STATS "\tSKIPPED\n";
	    next();
	}
	if ($score{chi} > 20) { 
	    print STDERR "Skipping $snp_id because of Hardy Weinberg distribution ($score{chi})\n";
	    print $STATS "\tSKIPPED\n";
	    next();
	}
	print $STATS "\tRETAINED\n";
	    
	print $OUT $snp_id;

	# output the snp
	#
	for (my $n = 0; $n < @clone_names; $n++) { 
	    if ($clone_info[$n] eq "") {
		if ( ($counts[$n]->[0] + $counts[$n]->[1]) >1 ) { 
		    print $OUT "\t";
		    printf $OUT "%.2f", $dosages[$n];
		}
		else { print $OUT "\tNA"; }
	    }
	}
	print $OUT "\n";	
    }
}

close($F);
close($OUT);
close($STATS);

open(my $F, ">", $filename.".failed_clones") || die "Can't open $filename\n";

for (my $n = 0; $n < @clone_names; $n++) { 
    if ($clone_info[$n]) { 
	print $F $clone_names[$n]."\n";
    }
}
close($F);
    
    
    
	    
