=head1 NAME

filter_hets.pl - a script to filter potential hets from a uneak gbs pipeline result

=head1 USAGE

 perl filter_hets.pl -f HapMap.hmc.txt > filter_hets.txt
 grep -v -w -f filter_hets.txt HapMap.hmc.txt > HapMap.hmc.filtered
 grep -v -w -f filter_hets.txt HapMap.hmp.txt > HapMap.hmp.filtered

=head1 DESCRIPTION


=head1 AUTHORS

 Jeremy D. Edwards (jde22@cornell.edu)

=cut


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
my $proportion_allowed_het = 0.7;
my $min_ratio = 0.03;
#my $min_reads_per_marker = 1500;
my $min_reads_per_marker = 3500;
while (<FILE>) {
  chomp $_;
  if ($firstline==1) {
    $firstline=0;
    next;
  }
  my @row =  split(/\t/, $_);
  my $number_of_samples = 0;
  foreach my $cell (@row) {
    if ($cell =~ m/^\d+\|\d+$/){
      $number_of_samples++;
    }
  }

  my $number_of_missing = 0;

  my $not_enough_reads = 0;
  my $ratio_too_low=0;
  my $het_count=0;
  my $i;
  for ($i=1; $i<=$number_of_samples; $i++) {
    my @alleles = split(/\|/,$row[$i]);
    if (($alleles[0]) == 0 && ($alleles[1] == 0)) {
      $number_of_missing++;
    }
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
  my $number_of_nonmissing = $number_of_samples - $number_of_missing;
  # if ($number_of_missing){
  #   $number_of_nonmissing = $number_of_samples-$number_of_missing;
  # } else {
  #   $number_of_nonmissing = $number_of_samples;
  # }
  if ($number_of_nonmissing == 0){
    #print "$row[0]\n";
    #print STDERR "$row[0] $number_of_samples $number_of_missing\n"; 
    next;
  }
  if ($het_count/$number_of_nonmissing > $proportion_allowed_het || $ratio_too_low || $not_enough_reads) {
    #print "$row[0]\n";
  } else {
    print "$row[0]\n";
  }
}


