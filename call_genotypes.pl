

use strict;

use Getopt::Long;
use Math::BigInt;

our $MIN_DEPTH= 0;
our $MAF = 0.001;
our $OUTFILE = "";

GetOptions("maf=i" => \$MAF,
	   "min_depth=i" => \$MIN_DEPTH,
	   "outfile=s" => \$OUTFILE,
    );

my $filename = shift;

if (! $OUTFILE) { 
    $OUTFILE = $filename.".dosage";
}

my $stats_file = $OUTFILE.".stats";

open(my $F, "<", $filename) || die "Can't open file $filename\n";
open(my $OUT, ">", $OUTFILE) || die "Can't open outfile $OUTFILE\n";
open(my $STATS, ">", $stats_file) || die "Can't open stats file $stats_file\n";

my @clone_names;
my %data;

my $snps_processed = 0;

while (<$F>) { 
    chomp;

    if ($snps_processed % 10 == 0) { 
	print STDERR "Processing $snps_processed ...\r";
    }

    $snps_processed++;

    if (m/\#\#/) { 
	next;
    }
    
    if (/^#CHROM/ || /^CHROM/ ) { 
	my @header = split /\t/;
	@clone_names = @header[9.. $#header];	
	if (! $MIN_DEPTH) { $MIN_DEPTH = scalar(@clone_names) * 2; }
	print STDERR "Using MIN_DEPTH of $MIN_DEPTH\n";
	print STDERR "Using MAF of $MAF\n";
    }
    
    else { 
	my ($chr, $position, $snp_id, $ref_allele, $alt_allele, $qual, $filter, $info, $format, $x, @snps) = split /\t/;
	
	#print STDERR "SNP: $snp_id\t$ref_allele\t$alt_allele\t$qual\n";
	my $depth = 0;
	if ($info =~ /DP=(\d+)$/) { 
	    $depth = $1;
	    #print STDERR "DEPTH: $depth\n";
	}
	
	my $total_c1 = 0;
	my $total_c2 = 0;
	
	if ($depth > $MIN_DEPTH) { 
	    my @counts;
	    for (my $n = 0; $n < @clone_names; $n++) { 
		
		if ($clone_names[$n]!~/blank/i) { 
		    #print "PROCESSING $snp_id CLONE NAME: $clone_names[$n]\n";
		    
		    my ($c1, $c2) = get_counts($snps[$n]);
		    $counts[$n] = [ $c1, $c2 ];
		    $total_c1 += $c1;
		    $total_c2 += $c2;
		}
	    }	    
	    my $allele_freq = $total_c1 / ($total_c1 + $total_c2);

	    if ($allele_freq > (1 - $MAF) || $allele_freq < $MAF) { 
		#print STDERR "Skipping ALLELE FREQ = $allele_freq\n";
		next();
	    }

	    my $pAA = $allele_freq **2;
	    my $pAB = $allele_freq * (1 - $allele_freq) * 2 ;
	    my $pBB = (1 - $allele_freq) **2;

	    my @dosages; 

	    for (my $n = 0; $n < @clone_names; $n++) { 
		if (ref($counts[$n]) eq "ARRAY") {
		    if ( ($counts[$n]->[0] + $counts[$n]->[1]) > 1 ) { 
			$dosages[$n] = get_genotype($pAA, $pAB, $pBB, $counts[$n]->[0], $counts[$n]->[1]);
		    }
		}
	    }

	    # hardy weinberg filter
	    my %score = hardy_weinberg_filter(@dosages);

	    print $STATS "$snp_id\t$score{scored_marker_fraction}\t$score{allele_freq}\t$score{chi}\n";

	    if ($score{scored_marker_fraction} < @dosages * 0.4) { 
		print STDERR "Skipping $snp_id because of low scored markers ($score{scored_marker_fraction})\n";
		next();
	    }
	    if ($score{allele_freq} > (1-$MAF)) { 
		print STDERR "Skipping $snp_id because of high allele frequency ($score{allele_freq})\n";
		next();
	    }
	    if ($score{chi} > 20) { 
		print STDERR "Skipping $snp_id because of Hardy Weinberg distribution ($score{chi})\n";
		next();
	    }

	    
	    print $OUT $snp_id;
	    for (my $n = 0; $n < @clone_names; $n++) { 
		if (ref($counts[$n]) eq "ARRAY") {
		    if ( ($counts[$n]->[0] + $counts[$n]->[1]) >1 ) { 
			print $OUT "\t";
			printf $OUT "%.2f", $dosages[$n];
		    }
		}
	    }
	    print $OUT "\n";

	}
    }
}

close($F);
close($OUT);
close($STATS);

sub get_counts { 
    my $snp = shift;
    
    my $counts = (split /\:/, $snp)[1];

    my ($c1, $c2) = split /\,/, $counts;

    return ($c1, $c2);

}

sub get_genotype { 
    my $pAA = shift;
    my $pAB = shift;
    my $pBB = shift;
    my $c1 = shift;
    my $c2 = shift;

    my $n = $c1 + $c2;

    my $N = Math::BigInt->new($n);

    my $pDAA = ($N->bnok($c1))->numify * 0.99 ** $c1 * 0.01 ** $c2;
    my $pDAB = ($N->bnok($c1))->numify * 0.5 ** $c1 * 0.5 ** $c2;
    my $pDBB = ($N->bnok($c2))->numify * 0.01 ** $c1 * 0.99 ** $c2;

    #print STDERR "pDAA: $pDAA, pDAB $pDAB, pDBB $pDBB\n";

    my $pSAA = $pDAA * $pAA;
    my $pSAB = $pDAB * $pAB;
    my $pSBB = $pDBB * $pBB;

    my $x = 1 / ($pSAA + $pSAB + $pSBB);

    my $dosage = ($pSAB  + 2 * $pSBB) * $x;

    return $dosage;
}

sub hardy_weinberg_filter { 
    my @dosages = @_;

    my %classes = ( AA => 0, AB => 0, BB => 0, NA => 0);
    
    foreach my $d (@dosages) { 
	if (! defined($d)) { 
	    $classes{NA}++;
	}
	elsif ($d >= 0 && $d <= 0.1) { 
	    $classes{AA}++;
	}
	elsif ($d >=0.9 && $d <= 1.1) { 
	    $classes{AB}++;
	}

	elsif ($d >=1.9 && $d <= 2.0) { 
	    $classes{BB}++;
	}
	else { 
	    #print STDERR "Dosage outlier: $d\n";
	}

    }

    print STDERR "NA count: $classes{NA}\n";
 
    my $total = $classes{AA} + $classes{AB} + $classes{BB};

    my %score;
    
    if ($total < @dosages * 0.4) { 
	#print STDERR "$total too small, skipping\n";
	$score{scored_marker_fraction} = $total / @dosages;
    }
    else { 
	print STDERR "Total: $total, lines: ".scalar(@dosages)."\n";
    }

    print STDERR "AA  $classes{AA}, AB $classes{AB}, BB $classes{BB} Total: $total\n";
    my $allele_freq = (2 * $classes{AA} + $classes{AB}) / (2 * $total);

    if ($allele_freq > (1-$MAF)) { 
	#print STDERR "Skipping due to allele frequence too close to 1\n";
	$score{allele_freq} = $allele_freq;
    }

    print STDERR "Allele freq = $allele_freq       \n";

    my $expected = $allele_freq **2 * $total;
    my $x = ($classes{AA} - $expected)**2 / $expected;

    $score{chi} = $x;

#    if ($x > 20) { 
#	return 1;
#    }
#    else { 
#	return 0;
#    }

    return %score;
}
	

    
    
    
	    
