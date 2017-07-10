# Copyright [2015-2017] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
use Test::More tests => 113, 'die';

use FindBin qw( $Bin );

BEGIN {
  use_ok 'Bio::DB::HTS::VCF';
  use_ok 'Bio::DB::HTS::VCF::Row';
  use_ok 'Bio::DB::HTS::VCF::RowPtr';
  use_ok 'Bio::DB::HTS::VCF::Header';
  use_ok 'Bio::DB::HTS::VCF::HeaderPtr';    
}

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
  # Test standard functions on a VCF file
  ok my $v = Bio::DB::HTS::VCF->new( filename => $Bin . "/data/test.vcf.gz" ), "VCF file open";
  is $v->num_variants(), 9, 'correct number of variants identified in file';

  my $h = $v->header();
  is $h->version(), "VCFv4.0", "VCF Header version matches" ;
  is $h->num_samples(), 3, "Number of samples" ;
  is_deeply $h->get_sample_names(), ['NA00001','NA00002','NA00003'], "sample names correct" ;
  is $h->num_seqnames(), 3, "Number of seqnames" ;
  ok my $seqnames = $h->get_seqnames() ;
  is_deeply $seqnames, ['19','20','X'], "sequence names correct" ;

  ok my $row = $v->next(), "Next row";
  is $row->chromosome($h), "19", "Chromosome value read" ;
  is $row->position(), "111", "Position value read" ;
  is $row->id(), "testid", "ID value read" ;
  is $row->num_filters(), 2, "Num Filters OK" ;
  is $row->has_filter($h,"DP50"), 1, "Actual Filter present" ;
  is $row->has_filter($h,"."), 0, "PASS filter absent" ;
  is $row->get_variant_type(1),1, "Variant type matches" ;

  #info related tests
  my $info_result ;
  $info_result = $row->get_info($h,"DB") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, [1], 'info flag read correctly';
  is $row->get_info_type($h,"DB"), "Flag", "info flag type correct" ;

  $info_result = $row->get_info($h,"AF") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, [0.5], 'info float read correctly';
  is $row->get_info_type($h,"AF"), "Float", "info float type correct" ;

  $info_result = $row->get_info($h,"TT") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, ["TESTSTRING"], 'info string read correctly';
  is $row->get_info_type($h,"TT"), "String", "info String type correct" ;

  $info_result = $row->get_info($h,"NS") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, [3], 'info ints read correctly';
  is $row->get_info_type($h,"NS"), "Integer", "info int type correct" ;

  ok $row = $v->next(), "Next row";
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

  #Format and genotype tests
  ok $row = $v->next(), "Next row";
  is $row->chromosome($h), "20", "Chromosome value read" ;
  is $row->get_format_type($h,"DP"), "Integer", "int format type correct" ;
  my $fmt_result = $row->get_format($h,"DP") ;
  isa_ok($fmt_result, 'ARRAY');
  is_deeply $fmt_result, [1,8,5], 'format ints read correctly' ;
  $fmt_result = $row->get_genotypes($h) ;
  isa_ok($fmt_result, 'ARRAY');
  #TODO resolve how these translate to the strings in htslib
  is_deeply $fmt_result, [2,3,4,3,4,4], 'genotypes read correctly' ;
  is $row->get_format($h,"IDONTEXIST"), 'ID_NOT_FOUND', 'format id not found ok';
  is $row->get_info($h,"IDONTEXIST"), 'ID_NOT_FOUND', 'info id not found ok';

  $v->close();
}


{
  # Test standard functions on a BCF file
  ok my $v = Bio::DB::HTS::VCF->new( filename => $Bin . "/data/test.bcf" ), "BCF file open";
  is $v->num_variants(), 9, 'correct number of variants identified in file';

  my $h = $v->header();
  is $h->version(), "VCFv4.0", "VCF Header version matches" ;
  is $h->num_samples(), 3, "Number of samples" ;
  is_deeply $h->get_sample_names(), ['NA00001','NA00002','NA00003'], "sample names correct" ;
  is $h->num_seqnames(), 3, "Number of seqnames" ;
  is_deeply $h->get_seqnames(), ['19','20','X'], "sequence names correct" ;

  ok my $row = $v->next(), "Next row";
  is $row->chromosome($h), "19", "Chromosome value read" ;
  is $row->position(), "111", "Position value read" ;
  is $row->id(), "testid", "ID value read" ;
  is $row->num_filters(), 2, "Num Filters OK" ;
  is $row->has_filter($h,"DP50"), 1, "Actual Filter present" ;
  is $row->has_filter($h,"."), 0, "PASS filter absent" ;
  is $row->get_variant_type(1),1, "Variant type matches" ;

  #info related tests
  my $info_result ;
  $info_result = $row->get_info($h,"DB") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, [1], 'info flag read correctly';
  is $row->get_info_type($h,"DB"), "Flag", "info flag type correct" ;

  $info_result = $row->get_info($h,"AF") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, [0.5], 'info float read correctly';
  is $row->get_info_type($h,"AF"), "Float", "info float type correct" ;

  $info_result = $row->get_info($h,"TT") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, ["TESTSTRING"], 'info string read correctly';
  is $row->get_info_type($h,"TT"), "String", "info String type correct" ;

  $info_result = $row->get_info($h,"NS") ;
  isa_ok($info_result, 'ARRAY');
  is_deeply $info_result, [3], 'info ints read correctly';
  is $row->get_info_type($h,"NS"), "Integer", "info int type correct" ;

  ok $row = $v->next(), "Next row";
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

  #Format and genotype tests
  ok $row = $v->next(), "Next row";
  is $row->chromosome($h), "20", "Chromosome value read" ;
  is $row->get_format_type($h,"DP"), "Integer", "int format type correct" ;
  my $fmt_result = $row->get_format($h,"DP") ;
  isa_ok($fmt_result, 'ARRAY');
  is_deeply $fmt_result, [1,8,5], 'format ints read correctly' ;
  $fmt_result = $row->get_genotypes($h) ;
  isa_ok($fmt_result, 'ARRAY');
  #TODO resolve how these translate to the strings in htslib
  is_deeply $fmt_result, [2,3,4,3,4,4], 'genotypes read correctly' ;
  is $row->get_format($h,"IDONTEXIST"), 'ID_NOT_FOUND', 'format id not found ok';
  is $row->get_info($h,"IDONTEXIST"), 'ID_NOT_FOUND', 'info id not found ok';

  $v->close();
}
