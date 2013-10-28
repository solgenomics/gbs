use strict;
use Getopt::Std;

our ($opt_f, $opt_a, $opt_b);
getopts('f:a:b:');

if (!$opt_f) {
  die "Filename required\n";
}

if (!$opt_a || !$opt_b) {
  die "Group files required\n";
}

open FILE1, "<", $opt_a or die "No group a file";
open FILE2, "<", $opt_b or die "No group a file";

my @first_group;
my @second_group;

while (<FILE1>) {
  chomp $_;
  push @first_group, $_;
}

while (<FILE2>) {
  chomp $_;
  push @second_group, $_;
}

close FILE1;
close FILE2;


my @sample_names;

my $filename=$opt_f;
open FILE, "<", $filename or die "No such file $filename";
my $firstline=1;
while (<FILE>) {
  chomp $_;
  if ($firstline==1) {
    my @firstrow =  split(/\t/, $_);
    foreach my $name (@firstrow) {
      push @sample_names, $name;
    }
    $firstline=0;
    next;
  }
  my @row =  split(/\t/, $_);

  print $row[0]."\t";
  my $first_allele;
  my $missing_count = 0;
  my %first_allele_count;
  my %second_allele_count;
  my %het_count;
  $first_allele_count{'a'} = 0;
  $first_allele_count{'b'} = 0;
  $second_allele_count{'a'} = 0;
  $second_allele_count{'b'} = 0;
  $het_count{'a'} = 0;
  $het_count{'b'} = 0;
  my $col = 0;

  my $shared;

  my %snps;

  while ($col < scalar(@row)) {
    if ($col < 11) {
      $col++;
      next;
    }

    my $group_name;

    foreach my $name (@first_group) {
      if ($sample_names[$col] eq $name) {
	$group_name = "a";
      }
    }

    foreach my $name (@second_group) {
      if ($sample_names[$col] eq $name) {
	$group_name = "b";
      }
    }

    my $snp = $row[$col];
    if ($snp eq "N") {
      $missing_count++;
    } else {

      if ($snp =~ m/[ATGC]/) {
	if ($first_allele) {
	  if ($snp eq $first_allele) {
	    $first_allele_count{$group_name}++;
	  } else {
	    $second_allele_count{$group_name}++;
	  }
	} else {
	  $first_allele = $snp;
	  $first_allele_count{$group_name}++;
	}
      } else {
	$het_count{$group_name}++;
      }
    }
  $col++;
  }


  my $group_a_total = $first_allele_count{'a'} + $het_count{'a'} + $second_allele_count{'a'};
  my $group_b_total = $first_allele_count{'b'} + $het_count{'b'} + $second_allele_count{'b'};


  my %major_allele_count;
  my %minor_allele_count;


  if ($first_allele_count{'a'} + $first_allele_count{'b'} > $second_allele_count{'a'} + $second_allele_count{'b'}) {
    $major_allele_count{'a'} = $first_allele_count{'a'};
    $major_allele_count{'b'} = $first_allele_count{'b'};
    $minor_allele_count{'a'} = $second_allele_count{'a'};
    $minor_allele_count{'b'} = $second_allele_count{'b'};
  } else {
    $major_allele_count{'a'} = $second_allele_count{'a'};
    $major_allele_count{'b'} = $second_allele_count{'b'};
    $minor_allele_count{'a'} = $first_allele_count{'a'};
    $minor_allele_count{'b'} = $first_allele_count{'b'};
  }

  my $group_a_frequency;
  my $group_b_frequency;

  if ($group_a_total > 0){
    $group_a_frequency = ($minor_allele_count{'a'} + 0.5 * $het_count{'a'})/$group_a_total;
  } else {
    $group_a_frequency = "NA";
  }


  if ($group_b_total > 0){
    $group_b_frequency = ($minor_allele_count{'b'} + 0.5 * $het_count{'b'})/$group_b_total;
  } else {
    $group_b_frequency = "NA";
  }


  my $frequence_difference;
  if ($group_a_frequency eq "NA" || $group_b_frequency eq "NA") {
    $frequence_difference = "NA";

  } else {
    $frequence_difference = abs($group_b_frequency - $group_a_frequency);
  }
  print "\t$group_a_frequency\t$group_b_frequency\t$frequence_difference\n";


}

close FILE;
