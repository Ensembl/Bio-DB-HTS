use Test::More tests => 11, 'die';

use FindBin qw( $Bin );

BEGIN { use_ok 'Bio::DB::HTS::VCF'; }

ok my $v = Bio::DB::HTS::VCF->new( filename => $Bin . "/data/test.vcf.gz" );
is $v->num_variants, 9, 'correct number of variants identified in file';

{
  # Test sweep functions
  my $sweep = Bio::DB::HTS::VCFSweep->new(filename => $Bin . "/data/test.vcf.gz");
  my $h = $sweep->header ;

  #At time of writing these tests with HTSlib v.1.3 there appears to be an issue
  #with the forward sweep function. Hence the tests use the reverse sweep
  my $row = $sweep->next_row();
  $row = $sweep->previous_row();
  is $row->chromosome($h), "X", "Chromosome value read" ;
  is $row->position(), 10, "Position value read" ;
  is $row->quality(), 10, "Quality value read" ;

  $row = $sweep->previous_row();
  is $row->reference(), "T", "Quality value read" ;

  $row = $sweep->previous_row();
  is $row->id(), "microsat1", "ID value read" ;
  is $row->num_alleles(), 2, "Number of alleles value read" ;
  my $a_team = $row->get_alleles() ;
  isa_ok($a_team, 'ARRAY');
  is_deeply $a_team, ['GA', 'GAC'], 'alleles are correct';

  $sweep->close() ;
}
