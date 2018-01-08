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

use FindBin '$Bin';
use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch", "$Bin/../t";

use Test::More;
use Bio::DB::HTS::Kseq;

my %example_files = (
    bug2335               => {
        'variant'       => 'sanger',
        'seq'           => 'TTGGAATGTTGCAAATGGGAGGCAGTTTGAAATACTGAATAGGCCTCATC'.
                           'GAGAATGTGAAGTTTCAGTAAAGACTTGAGGAAGTTGAATGAGCTGATGA'.
                           'ATGGATATATG',
        'qual'          => '@8A8@7<=A9:8?:#6B>*2B<<<:B>*=A<(<6;6<<2;8@8A9<<=<='.
                           '=<<@88=<<A8<D?-=<<:B>+<==B:<<@8C<<A9<?79=9<=<;=<=A'.
                           '9=B:8<<=<=;',
        'phred'         => '31 23 32 23 31 22 27 28 32 24 25 23 30 25 2 21 33 '.
                           '29 9 17 33 27 27 27 25 33 29 9 28 32 27 7 27 21 '.
                           '26 21 27 27 17 26 23 31 23 32 24 27 27 28 27 28 '.
                           '28 27 27 31 23 23 28 27 27 32 23 27 35 30 12 28 '.
                           '27 27 25 33 29 10 27 28 28 33 25 27 27 31 23 34 '.
                           '27 27 32 24 27 30 22 24 28 24 27 28 27 26 28 27 '.
                           '28 32 24 28 33 25 23 27 27 28 27 28 26',
        'name'    => 'DS6BPQV01A2G0A',
        'desc'          => undef,
        'count'         => 1
        },
    test1_sanger            => {
        'variant'       => 'sanger',
        'seq'           => 'TATTGACAATTGTAAGACCACTAAGGATTTTTGGGCGGCAGCGACTTGGA'.
                           'GCTCTTGTAAAAGCGCACTGCGTTCCTTTTCTTTATTCTTTTGATCTTGA'.
                           'GAATCTTCTAAAAATGCCGAAAAGAAATGTTGGGAAGAGAGCGTAATCAG'.
                           'TTTAGAAATGCTCTTGATGGTAGCTTTATGTTGATCCATTCTTCTGCCTC'.
                           'CTTTACGAATAAAATAGAATAAAACTCAAATGACTAATTACCTGTATTTT'.
                           'ACCTAATTTTGTGATAAAATTCAAGAAAATATGTTCGCCTTCAATAATTA'.
                           'TG',
        'qual'          => 'FFFFFFFFFFFIGIIFFFHHIHHHHHFBBBBBHHC>==GHHHHHHFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGG>>>CGFFBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGGGFFFFFFFFFFFFFFFFFFFFFFFFCCCGGFFFFFFFFFFIIIIGGGGIIIGGGIIIIIIIIIIGGGGGA::::??@AA@@@@@@@@@4444477@@@@@@@@AA@A@@@@@@:::==?????@@A',
        'phred'         => '37 37 37 37 37 37 37 37 37 37 37 40 38 40 40 37 '.
                           '37 37 39 39 40 39 39 39 39 39 37 33 33 33 33 33 '.
                           '39 39 34 29 28 28 38 39 39 39 39 39 39 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 38 '.
                           '38 29 29 29 34 38 37 37 33 33 33 33 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 38 38 '.
                           '38 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 34 34 34 38 38 37 37 '.
                           '37 37 37 37 37 37 37 37 40 40 40 40 38 38 38 38 '.
                           '40 40 40 38 38 38 40 40 40 40 40 40 40 40 40 40 '.
                           '38 38 38 38 38 32 25 25 25 25 30 30 31 32 32 31 '.
                           '31 31 31 31 31 31 31 31 19 19 19 19 19 22 22 31 '.
                           '31 31 31 31 31 31 31 32 32 31 32 31 31 31 31 31 '.
                           '31 25 25 25 28 28 30 30 30 30 30 31 31 32',
        'name'    => 'SRR005406.250',
        'desc'          => 'FB9GE3J10F6I2T length=302',
        'count'         => 250
                                },
    test2_solexa            => {
        'variant'       => 'solexa',
        'seq'           => 'GTATTATTTAATGGCATACACTCAA',
        'qual'          => 'YYYYYYYYYYWYYYYWYWWUWWWQQ',
        'phred'         => '25 25 25 25 25 25 25 25 25 25 23 25 25 25 25 23 '.
                           '25 23 23 21 23 23 23 17 17',
        'name'    => 'SLXA-B3_649_FC8437_R1_1_1_183_714',
        'desc'          => undef,
        'count'         => 5
                                },
    test3_illumina          => {
        'variant'       => 'illumina',
        'seq'           => 'CCAAATCTTGAATTGTAGCTCCCCT',
        'qual'          => 'OSXOQXXXXXSXXUXXTXXXXTRMS',
        'phred'         => '15 19 24 15 17 24 24 24 24 24 19 24 24 21 24 24 '.
                           '20 24 24 24 24 20 18 13 19',
        'name'    => 'FC12044_91407_8_200_285_136',
        'desc'          => undef,
        'count'         => 25
                                },
    example                 => {
        'variant'       => 'sanger', # TODO: guessing on the format here...
        'seq'           => 'GTTGCTTCTGGCGTGGGTGGGGGGG',
        'qual'          => ';;;;;;;;;;;9;7;;.7;393333',
        'phred'         => '26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 '.
                           '13 22 26 18 24 18 18 18 18',
        'name'    => 'EAS54_6_R1_2_1_443_348',
        'desc'          => 'baz',
        'count'         => 3
                                },
    illumina_faked          => {
        'variant'       => 'illumina',
        'seq'           => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN',
        'qual'          => 'hgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@',
        'phred'         => '40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 '.
                           '21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0',
        'name'    => 'Test',
        'desc'          => 'PHRED qualities from 40 to 0 inclusive',
        'count'         => 1
                                },
    sanger_93               => {
        'variant'       => 'sanger',
        'seq'           => 'ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC'.
                           'TGACTGACTGACTGACTGACTGACTGACTGACTGAN',
        'qual'          => '~}|{zyxwvutsrqponmlkjihgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFE'.
                           'DCBA@?>=<;:9876543210/.-,+*)(\'&%$#"!',
        'phred'         => '93 92 91 90 89 88 87 86 85 84 83 82 81 80 79 78 77 76 75 '.
                           '74 73 72 71 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56 '.
                           '55 54 53 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 '.
                           '36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 '.
                           '17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0',
        'name'    => 'Test',
        'desc'          => 'PHRED qualities from 93 to 0 inclusive',
        'count'         => 1
                                },
    sanger_faked            => {
        'variant'       => 'sanger',
        'seq'           => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN',
        'qual'          => 'IHGFEDCBA@?>=<;:9876543210/.-,+*)(\'&%$#"!',
        'phred'         => '40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 '.
                           '21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0',
        'name'    => 'Test',
        'desc'          => 'PHRED qualities from 40 to 0 inclusive',
        'count'         => 1
                                },
    solexa_example          => {
        'variant'       => 'solexa',
        'seq'           => 'GTATTATTTAATGGCATACACTCAA',
        'qual'          => 'YYYYYYYYYYWYYYYWYWWUWWWQQ',
        'phred'         => '25 25 25 25 25 25 25 25 25 25 23 25 25 25 25 23 '.
                           '25 23 23 21 23 23 23 17 17',
        'name'    => 'SLXA-B3_649_FC8437_R1_1_1_183_714',
        'desc'          => undef,
        'count'         => 5
                                },
    solexa_faked            => {
        'variant'       => 'solexa',
        'seq'           => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN',
        'qual'          => 'hgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@?>=<;',
        'phred'         => '40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 '.
                           '24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 10 9 '.
                           '8 7 6 5 5 4 4 3 3 2 2 1 1',
        'name'    => 'slxa_0001_1_0001_01',
        'desc'          => undef,
        'count'         => 1
                                },
    tricky                  => {
        'variant'       => 'sanger', # TODO: guessing on the format here...
        'seq'           => 'TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA',
        'qual'          => 'IIIIIII.IIIIII1@44@-7.%<&+/$/%4(++(%',
        'phred'         => '40 40 40 40 40 40 40 13 40 40 40 40 40 40 16 31 '.
                           '19 19 31 12 22 13 4 27 5 10 14 3 14 4 19 7 10 10 '.
                           '7 4',
        'name'    => '071113_EAS56_0053:1:3:990:501',
        'desc'          => undef,
        'count'         => 4
                                },
    evil_wrapping           => {
        'variant'       => 'sanger', # TODO: guessing on the format here...
        'seq'           => 'AACCCGTCCCATCAAAGATTTTGGTTGGAACCCGAAAGGGTTTTGAATTC'.
                           'AAACCCCTTTCGGTTCCAACTATTCAATTGTTTAACTTTTTTTAAATTGA'.
                           'TGGTCTGTTGGACCATTTGTAATAATCCCCATCGGAATTTCTTT',
        'qual'          => 'A;@;%75?:#<<9EA1;=EA3%B;B;A;B;@;%9EA1EA1EA3%<B;A;8EA0D@3$EA1=B;A;B;B;:=:B;:B:A9:EA0A9<FA81+&"D?-B;4<::/<;=:A98-5?6=C>+8<<3;=4:DA3%<;=8-9.A=):B=*',
        'phred'         => '32 26 31 26 4 22 20 30 25 2 27 27 24 36 32 16 '.
                           '26 28 36 32 18 4 33 26 33 26 32 26 33 26 31 26 '.
                           '4 24 36 32 16 36 32 16 36 32 18 4 27 33 26 32 26 '.
                           '23 36 32 15 35 31 18 3 36 32 16 28 33 26 32 26 33 '.
                           '26 33 26 25 28 25 33 26 25 33 25 32 24 25 36 32 '.
                           '15 32 24 27 37 32 23 16 10 5 1 35 30 12 33 26 19 '.
                           '27 25 25 14 27 26 28 25 32 24 23 12 20 30 21 28 '.
                           '34 29 10 23 27 27 18 26 28 19 25 35 32 18 4 27 26 '.
                           '28 23 12 24 13 32 28 8 25 33 28 9',
        'name'    => 'SRR014849.203935',
        'desc'          => 'EIXKN4201B4HU6 length=144',
        'count'         => 3
                                },
    zero_qual           => {
        'variant'       => 'sanger', # TODO: guessing on the format here...
        'seq'           => 'AG',
        'qual'          => 'DD',
        'phred'         => '0 0',
        'name'          => 'someID2',
        'desc'          => undef,
        'count'         => 2
    }
);

for my $example (sort keys %example_files) {
    my $target= "$example.fastq";
    my $file = File::Spec->catfile($Bin, 'kseq_data', $target);
    my $kseq = Bio::DB::HTS::Kseq->new($file);
    SKIP: {
      skip 'Cannot go further. No kseq file found for '.$target, 8 if ! defined $kseq;
      my $it = $kseq->iterator;
      my $ct = 0;
      my $sample_seq;
      while (my $seq = $it->next_seq) {
          $ct++;
          $sample_seq = $seq; # always grab the last seq
      }
      is($ct, $example_files{$example}->{count}, "correct num. seqs in $example");
      ok(defined($sample_seq), 'sample sequence obtained');
      if ($sample_seq) {
          for my $key (qw(seq desc name qual)) {
              is($sample_seq->{$key},
                 $example_files{$example}->{$key},
                 "$key matches $example");
          }
          is(length($sample_seq->{seq}), length($sample_seq->{qual}), 'Checking qual and seq lengths are the same in '.$example);
          is(length($sample_seq->seq()), length($sample_seq->qual()), 'Checking qual and seq lengths are the same using object methods in '.$example);
      }
    }
}

done_testing();