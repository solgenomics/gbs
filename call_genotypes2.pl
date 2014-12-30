use strict;
use CXGN::Genotype;
use CXGN::GenotypeIO;

my $file = shift;

my $gtio = CXGN::GenotypeIO->new({ file => $file });

print STDERR $gtio->count();
