
package StatFunc;

use Math::BigInt;

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

    my $N1 = Math::BigInt->new($n);
    my $N2 = Math::BigInt->new($n);

 #   print STDERR "$N1 bnok $c1 is: ". $N1->bnok($c1)."\n";

    my $Nbnokc1 = $N1->bnok($c1)->numify();
    my $Nbnokc2 = $N2->bnok($c2)->numify();
    
#    print STDERR "NBnokc1: $Nbnokc1, NBnokc2 $Nbnokc2\n";

    my $pDAA = $Nbnokc1 * (0.999 ** $c1) * (0.001 ** $c2);
    my $pDAB = $Nbnokc1 * (0.5 ** $c1) * (0.5 ** $c2);
    my $pDBB = $Nbnokc2 * (0.999 ** $c2) * (0.001 ** $c1);

 #   print STDERR "pDAA: $pDAA pDAB $pDAB, pDBB $pDBB\n";

    my $pSAA = $pDAA * $pAA;
    my $pSAB = $pDAB * $pAB;
    my $pSBB = $pDBB * $pBB;

    if ($pSAA + $pSAB + $pSBB == 0) { #
	return "NA";
    }
    
    my $x = 1 / ($pSAA + $pSAB + $pSBB);

    my $dosage = ($pSAB  + 2 * $pSBB) * $x;

 #   print "DOSAGE: $dosage\n";
    return $dosage;
}

sub hardy_weinberg_filter { 
    my $dosages = shift;
    my $ignore_clones_count = shift;

    my %classes = ( AA => 0, AB => 0, BB => 0, NA => 0);
    
    foreach my $d (@$dosages) { 
	if (! defined($d)) { 
	    $classes{NA}++;
	}
	elsif ( ($d >= 0) && ($d <= 0.1) ) { 
	    $classes{AA}++;
	}
	elsif ( ($d >=0.9) && ($d <= 1.1) ) { 
	    $classes{AB}++;
	}

	elsif (($d >=1.9) && ($d <= 2.0)) { 
	    $classes{BB}++;
	}
	else { 
	    #print STDERR "Dosage outlier: $d\n";
	}

    }

    #print STDERR "NA count: $classes{NA}\n";
 
    if ( ( ($classes{AA} ==0) && ($classes{AB} ==0)) ||
	( ($classes{BB} == 0) && ($classes{AB} ==0)) ) { 
	    return ( invalid => 1);
    }

    my $total = $classes{AA} + $classes{AB} + $classes{BB};

    my %score = ();
    
    $score{scored_marker_fraction} = $total / (@$dosages - $ignore_clones_count);
    
    print STDERR "AA  $classes{AA}, AB $classes{AB}, BB $classes{BB} Total: $total\n";
    my $allele_freq = (2 * $classes{AA} + $classes{AB}) / (2 * $total);

    $score{heterozygote_count} = $classes{AB};

    $score{allele_freq} = $allele_freq;
    
    my $expected = $allele_freq * (1-$allele_freq) * 2 * $total;

    print STDERR "TOTAL: $total\n";

    # overrepresented heterzygocity is unlikely (underrep is not)
    # do the hardy weinberg test only if heterozygocity is higher than expected.
    #
    if ($classes{AB} > $expected) { 
	my $x = ($classes{AB} - $expected)**2 / $expected;
	$score{chi} = $x;
    }
    else { 
	$score{chi} = 1;  #essentially ignore this test.
    }

    return %score;
}
	

1;
