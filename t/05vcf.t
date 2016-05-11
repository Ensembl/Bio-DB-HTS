use Test::More tests => 18, 'die';

use FindBin qw( $Bin );

BEGIN { use_ok 'Bio::DB::HTS::VCF'; }

{
  # Test sweep functions
  my $sweep = Bio::DB::HTS::VCF::Sweep->new(filename => $Bin . "/data/test.vcf.gz");
  my $h = $sweep->header ;

  my $row = $sweep->next_row();
  is $row->chromosome($h), "19", "Chromosome value read" ;
  #This should one day be fixed to 2 - but the HTSlib API always returns 0 at the moment
  #is $row->num_filters(), "2", "Number of filters read" ;
  #Once that is corrected, also add a test for the filters themselves

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

{
  # Test standard functions
  ok my $v = Bio::DB::HTS::VCF->new( filename => $Bin . "/data/test.vcf.gz" );
  is $v->num_variants(), 9, 'correct number of variants identified in file';

  my $h = $v->header();
  ok my $row = $v->next();
  is $row->chromosome($h), "19", "Chromosome value read" ;
  is $row->position(), "111", "Position value read" ;
  is $row->id(), "testid", "ID value read" ;
  is $row->num_filters(), 2, "Num Filters OK" ;
  is $row->has_filter($h,"DP50"), 1, "Actual Filter present" ;
  is $row->has_filter($h,"."), 0, "PASS filter absent" ;

  ok $row = $v->next();
  is $row->chromosome($h), "19", "Chromosome value read" ;
  is $row->position(), "112", "Position value read" ;
  is $row->quality(), "10", "Quality value read" ;
  is $row->reference(), "A", "Reference value read" ;
  is $row->num_alleles(), 1, "Num Alleles" ;
  is $row->is_snp(), 1, "Is SNP" ;
  my $a_team = $row->get_alleles() ;
  isa_ok($a_team, 'ARRAY');
  is_deeply $a_team, ['G'], 'alleles are correct';
  is $row->num_filters(), 1, "Num Filters OK" ;
  is $row->has_filter($h,"PASS"), 1, "PASS Filter present" ;
  is $row->has_filter($h,"DP50"), 0, "Actual Filter absent" ;
  is $row->has_filter($h,"sdkjsdf"), -1, "Made up filter not existing" ;

  $v->close();
}
