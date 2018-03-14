#-*-Perl-*-
# Copyright [2015-2018] EMBL-European Bioinformatics Institute
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

use ExtUtils::MakeMaker;
use File::Temp qw(tempfile);
use FindBin '$Bin';

use Test::More;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";
use Bio::DB::HTS;

my $samfile = "$Bin/data/test.sam";

# low-level tests
my $hts_file = Bio::DB::HTSfile->open($samfile);
ok($hts_file);

#
# https://github.com/Ensembl/Bio-DB-HTS/issues/54
# reading header again should seek at start of file for sam as well
# and should not exit with samtools error:
# [E::sam_parse1] missing SAM header
# [W::sam_read1] parse error at line 4
#
# WARNING
#
# Couldn't get it to working for htslib < 1.5 (segfault)
# Skip the tests in these cases
#
SKIP: {
  skip "Cannot test header reading for htslib < 1.5", 10
    if Bio::DB::HTS->htslib_version() lt "1.5";
  
  my $header  = $hts_file->header_read();
  my $targets = $header->n_targets;
  is( $targets, 1, "Number of reference sequences" );

  $header  = $hts_file->header_read();
  $targets = $header->n_targets;
  is( $targets, 1, "num reference sequences" );

  my $target_names = $header->target_name;
  ok($target_names);
  is( scalar @$target_names, 1, "num reference sequence names" );
  is( $target_names->[0], 'ref', "reference sequence name" );

  my $target_lens = $header->target_len;
  ok($target_lens);
  is( scalar @$target_lens, 1 );
  # ok( $target_lens->[0],    100 );

  my $text = $header->text;
  ok( length $text > 0 );

  my $c = "\@CO\tThis is a comment\n";
  $header->text($c);
  is( $header->text, $c );

  my $count;
  while ( my $b = $hts_file->read1($header) ) {
    $count++;
  }
  is( $count, 2, "num reads mapped" );
}

# TODO
# parse region

# high-level tests
my $hts = Bio::DB::HTS->new( -bam  => $samfile );
ok( $hts );
is( $hts->n_targets, 1 );

#
# https://github.com/Ensembl/Bio-DB-HTS/issues/54 again
# should read header correctly and not exit with error
#
SKIP: {
  skip "Cannot test iterating over alignments for htslib < 1.5", 4
    if Bio::DB::HTS->htslib_version() lt "1.5";
  
  my $alignment_iterator = $hts->features(-iterator => 1);
  ok($alignment_iterator->next_seq->qname eq "r001");
  ok($alignment_iterator->next_seq->qname eq "r002");
  ok(!$alignment_iterator->next_seq);
}

is( $hts->length('ref'), 100 );
is( join('', $hts->seq_ids), 'ref' );

my $seg = $hts->segment('ref');
ok($seg);
is( $seg->length, 100 );
my $seq = $seg->seq;
isa_ok( $seq, 'Bio::PrimarySeq' );
is( length $seq->seq, 100 );

# TODO
# get_features_by_location
#   seq_id, qscore
#   get_tag_values
#   attributes
#   query
#   target
#   features

done_testing();
