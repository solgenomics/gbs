#!/usr/bin/perl -w
use strict;

my $file  = shift;

open(my $F, "<", $file) || die "Can't open file $file";

my $header = <$F>;
chomp($header);

my $test_file = $file .".test";
my $mod_markers_count;

open(my $T, ">", $test_file) || die "Can't open test file $test_file";
print $T "$header\n";

while (<$F>) { 
    chomp;
    my ($marker, @scores) = split /\t/;

    # count NAs in scores...
    #
    my $na_count = 0;
    my $s;
    foreach $s (@scores) { 
	if ($s eq "NA") { 
	    $na_count++;
	}
    }
    if (($na_count < @scores * 0.10) and (rand() <0.50)) { 

	# randomly mask about half of known markers in some of the best accessions
	#
	$mod_markers_count++;
	print STDERR "Randomly masking half of known SNP values for well characterized marker $marker (Real NA count: $na_count)...\n";
	foreach $s (@scores) {
	    if ($s eq "NA") {
	    } elsif (rand() >0.50) {
		$s = "NA";
	    } else {
	    }
	}
    }
    print $T join "\t", ($marker, @scores);
    print $T "\n";
}
print "$mod_markers_count markers modified./"; 
