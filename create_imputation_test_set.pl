
use strict;

my $file  = shift;

open(my $F, "<", $file) || die "Can't open file $file";

my $header = <$F>;
chomp($header);

my $test_file = $file .".test";

open(my $T, ">", $test_file) || die "Can't open test file $test_file";
print $T "$header\n";

while (<$F>) { 
    chomp;
    my ($marker, @scores) = split /\t/;

    # count NAs in scores...
    #
    my $na_count = 0;
    foreach my $s (@scores) { 
	if ($s eq "NA") { 
	    $na_count++;
	}
    }

    if ($na_count > @scores * 0.70) { 
	# do not eliminate this line
	#
	print STDERR "Ignoring line for marker $marker (NA count: $na_count)...\n";
    }
    else { 
	if (rand() > 0.20) { 
	    print $T join "\t", ($marker, @scores);
	    print $T "\n";
	}
    }
}
