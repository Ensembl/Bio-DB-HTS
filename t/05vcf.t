use Test::More tests => 3, 'die';

use FindBin qw( $Bin );

BEGIN { use_ok 'Bio::DB::HTS::VCF'; }

ok my $v = Bio::DB::HTS::VCF->new( filename => $Bin . "/data/test.vcf.gz" );
is $v->num_variants, 9, 'correct number of variants identified in file';

{
  my $sweep = Bio::DB::HTS::VCFSweep->new(filename => "data/test.vcf.gz");

  #At time of writing these tests with HTSlib v.1.3 there appears to be an issue
  #with the forward sweep function. Hence the tests use the reverse sweep
  my $row = $sweep->previous_row();

  #TODO add actual tests on values once functions have been written for VCF::Row object


  $sweep->close() ;
}
