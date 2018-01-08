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


#########################

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";
use Test::More tests => 15 ;

#########################

BEGIN { use_ok('Bio::DB::HTS') } ;
ok(1) ;

BEGIN { use_ok('Bio::DB::HTS::Faidx') } ;
ok(1) ;

my $fasta = "$Bin/data/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz" ;
my $location = "I:1-100" ;
my $index = Bio::DB::HTS::Faidx->new($fasta);
ok($index) ;

my $seq = "" ;
my $length = 0 ;
($seq, $length) = $index->get_sequence($location);
ok($seq eq
  'CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG');
ok($length eq 100);


my @seq_ids = $index->get_all_sequence_ids();
ok($seq_ids[0] eq 'I') ;

$seq = $index->get_sequence_no_length($location);
ok($seq eq
  'CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG');


$seq = "" ;
$length = 0 ;
($seq, $length) = $index->get_sequence2("I",1,99);
ok($seq eq
  'CACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG');
ok($length eq 99);


$seq = "" ;
$length = 0 ;
$seq = $index->get_sequence2_no_length("I",1,99);
ok($seq eq
  'CACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG');



$length = $index->length('I');
note $length ;
ok($length eq 230218, "length of I");

my $has_seq = $index->has_sequence('I');
ok($has_seq eq 1, "has sequence I");

$has_seq = $index->has_sequence('II');
ok($has_seq eq 0, "Does not have sequence II");
