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
use constant TEST_COUNT => 22;
use Data::Dumper;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if ($@) {
        use lib 't';
    }
    use Test;
    plan test => TEST_COUNT;
}

use Bio::DB::HTS;

{
    my $bamfile = "$Bin/data/ex1.bam";
    my $hts_file     = Bio::DB::HTSfile->open($bamfile);
    ok($hts_file);

    my $header  = $hts_file->header_read();

    my $cramfile = "$Bin/data/write_test_07.cram";
    my $hts_file2 = Bio::DB::HTSfile->open( $cramfile, 'wc' );
    ok($hts_file2);

    $hts_file2->header_write($header, "$Bin/data/ex1.fa");
    my $count;
    while ( my $b = $hts_file->read1($header) ) {
        $hts_file2->write1($header, $b);
        $count++;
    }
    $hts_file2 = undef;
    ok( $count, 3307 );

    $hts_file = Bio::DB::HTSfile->open( $cramfile );
    ok($hts_file);

    $header  = $hts_file->header_read();
    my $targets = $header->n_targets;
    ok( $targets, 2 );

    my $target_names = $header->target_name;
    ok($target_names);
    ok( scalar @$target_names, 2 );
    ok( $target_names->[0],    'seq1' );

    my $target_lens = $header->target_len;
    ok($target_lens);
    ok( scalar @$target_lens, 2 );
    ok( $target_lens->[0],    1575 );

    my $text = $header->text;
    ok( length $text > 0 );

    $count = 0;
    while ( my $b = $hts_file->read1($header) ) {
        $count++;
    }
    ok( $count, 3307 );

    my @result = $header->parse_region('seq2:51-1000');
    ok( $result[0], 1 );
    ok( $result[1], 50 );
    @result = $header->parse_region('seq_invalid:51-1000');
    ok( scalar @result, 0 );

    Bio::DB::HTSfile->index_build($cramfile);
    my $index = Bio::DB::HTSfile->index_load($hts_file);
    ok($index);

    my @a;
    my $print_region = sub {
        my ( $alignment, $data ) = @_;
        push @a, $alignment;
        return;
    };

    $index->fetch( $hts_file, $header->parse_region('seq2'),
                   $print_region, "foobar" );
    ok( scalar @a > 1 );

    my $fai = Bio::DB::HTS::Fai->open("$Bin/data/ex1.fa");

    my %matches;
    my $fetch_back = sub {
        my ( $tid, $pos, $p ) = @_;
        my $p1 = $pos + 1;
        my $r  = $fai->fetch( $header->target_name->[$tid] . ":$p1-$p1" );
        for my $pileup (@$p) {
            my $b    = $pileup->b;
            my $qpos = $pileup->qpos;
            my $base =
              $pileup->indel == 0 ? substr( $b->qseq, $qpos, 1 ) :
              $pileup->indel > 0 ? '*' :
              '-';
            $matches{matched}++ if $r eq $base;
            $matches{total}++;
        }
    };

    $index->pileup( $hts_file, $header->parse_region('seq2:1-100'), $fetch_back );
    ok( $matches{matched}/$matches{total} > 0.99 );

    # try to get coverage
    my $coverage =
      $index->coverage( $hts_file, $header->parse_region('seq2'), 100, 9000 );
    ok( scalar @$coverage, 100 );
    my @c = sort { $a <=> $b } @$coverage;
    ok( $c[0] >= 0 );
    ok( $c[-1] < 1000 );

}

exit 0;

sub reversec {
    my $dna = shift;
    $dna =~ tr/gatcGATC/ctagCTAG/;
    return scalar reverse $dna;
}

__END__
