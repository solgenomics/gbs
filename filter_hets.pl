use strict;
use Getopt::Std;

our ($opt_f);
getopts('f:');

if (!$opt_f) {
  die "Filename required\n";
}


my $filename=$opt_f;
open FILE, "<", $filename or die "No such file $filename";
my $firstline=1;
my $number_of_samples = 12;
my $proportion_allowed_het = 0.5;
my $min_ratio = 0.10;
my $min_reads_per_marker = 50;
while (<FILE>) {
  chomp $_;
  if ($firstline==1) {
    $firstline=0;
    next;
  }
  my @row =  split(/\t/, $_);
  my $not_enough_reads = 0;
  my $ratio_too_low=0;
  my $het_count=0;
  my $i;
  for ($i=1; $i<=$number_of_samples; $i++) {
    my @alleles = split(/\|/,$row[$i]);
    if (($alleles[0]) > 0 && ($alleles[1] > 0)) {
      $het_count++;
      my $het_major_allele;
      my $het_minor_allele;
      if ($alleles[0] >= $alleles[1]) {
	$het_major_allele = $alleles[0];
	$het_minor_allele = $alleles[1];
      }
      else {
	$het_major_allele = $alleles[1];
	$het_minor_allele = $alleles[0];
      }
      if (($het_minor_allele/$het_major_allele) < $min_ratio) {
       $ratio_too_low = 1;
      }
    }
  }
  if ($row[3+$number_of_samples] + $row[4+$number_of_samples] < $min_reads_per_marker) {
    $not_enough_reads = 1;
  }
  # my $het_allele1 = $row[1+$number_of_samples];
  # my $het_allele2 = $row[2+$number_of_samples];
  # my $het_major_allele;
  # my $het_minor_allele;

  # if ($het_allele1 > 0 && $het_allele2 > 0) {
  #   if ($het_allele1 >= $het_allele2) {
  #     $het_major_allele = $het_allele1;
  #     $het_minor_allele = $het_allele2;
  #   }
  #   else {
  #     $het_major_allele = $het_allele2;
  #     $het_minor_allele = $het_allele1;
  #   }
  #   if (($het_minor_allele/$het_major_allele) < $min_ratio) {
  #     $ratio_too_low = 1;
  #   }
  # }
  if ($het_count/$number_of_samples > $proportion_allowed_het || $ratio_too_low || $not_enough_reads) {
    print "$row[0]\n";
  }
}


