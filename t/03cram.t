
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
use File::Temp qw(tempfile);
use FindBin '$Bin';
use constant TEST_COUNT => 139;

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
use Bio::DB::HTS::AlignWrapper;

my $cramfile = "$Bin/data/yeast.sorted.cram";
my $fastafile = "$Bin/data/yeast.fasta";

# low level tests (defined in lib/Bio/DB/HTS.pm)
{
    my $hts_file = Bio::DB::HTSfile->open($cramfile);
    ok($hts_file);

    my $header  = $hts_file->header_read();
    my $targets = $header->n_targets;
    ok( $targets, 17, "Targets Found" );

    my $target_names = $header->target_name;
    ok( $target_names );
    ok( scalar @$target_names );
    ok( $target_names->[10], 'X' );

    my $target_lens = $header->target_len;
    ok( $target_lens ) ;
    ok( scalar @$target_lens, 17 );
    ok( $target_lens->[0], 230218 );

    my $text = $header->text;
    ok( length $text > 0 );

    my $c = "\@CO\tThis is a comment\n";
    $header->text($c);
    ok( $header->text, $c );

    my $region = 'X:51-1000' ;
    my $fai = Bio::DB::HTS::Fai->open("$Bin/data/yeast.fasta");
    my $seq = $fai->fetch($region);
    ok( length $seq, 950 );

    my $count;
    while ( my $b = $hts_file->read1($header) )
    {
        $count++;
    }
    ok( $count, 50000 );

    my @result = $header->parse_region($region);
    ok( $result[0], 10 );
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

    $index->fetch( $hts_file, $header->parse_region('X'),
                   $print_region, "foobar" );
    ok( scalar @a > 1 );

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

    $index->pileup( $hts_file, $header->parse_region($region), $fetch_back );
    ok( $matches{matched}, 37 );
    ok( $matches{total}, 96 );

}

# high level tests (defined in lib/Bio/DB/HTS.pm)

#
# Note that the high level API does not reset the CRAM file pointer to the start
# of the file as the method to do so is (at time or writing) not easily accessible.
# Therefore some of these tests create a new HTS object.
#
my $dummy = eval
{
  Bio::DB::HTS->new( -fasta => "invalid_path.txt",
                     -bam   => "invalid_path.txt" );
};
ok( $dummy, undef );
ok( $@ =~ /does not exist/ );


for my $use_fasta (0,1)
{
  my $hts = Bio::DB::HTS->new(
                              -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    ok($hts);
    ok( $hts->n_targets, 17 );

    ok( $hts->length('X'), 745751 );

    my $seg = $hts->segment('I');
    ok($seg);
    ok( $seg->length, 230218 );
    my $seq = $seg->seq;
    ok( $seq->isa('Bio::PrimarySeq') );
    ok( length $seq->seq, 230218 );

    my @alignments =
      $hts->get_features_by_location( -seq_id => 'I',
                                      -start  => 500,
                                      -end    => 20000 );
    ok( scalar @alignments, 71 );
    ok( $alignments[0]->seq_id, 'I' );

    ok( scalar @{ $alignments[0]->qscore }, length $alignments[0]->dna );

    my @keys = $alignments[0]->get_all_tags;
    ok( scalar @keys,                         17 );
    ok( $alignments[0]->get_tag_values('AS'), 36 );

    my %att = $alignments[0]->attributes;
    ok( scalar( keys %att ),       17 );
    ok( $alignments[0]->cigar_str, '36M' );

    my $query = $alignments[0]->query;
    ok($query);
    ok( $query->start,          1 );
    ok( $query->end,            36 );
    ok( $query->length,         36 );
    ok( $query->dna,            $alignments[0]->dna );
    ok( $alignments[0]->strand, -1 );
    ok( $query->strand,         -1 );

    my $target = $alignments[0]->target;
    ok($target);
    ok( $target->start,  36 );
    ok( $target->end,    1 );
    ok( $target->length, 36 );
    ok( $target->dna,    reversec( $alignments[0]->dna ) );

    ok( $alignments[0]->get_tag_values('FLAGS'), $alignments[0]->flag_str );

    my @pads = $alignments[0]->padded_alignment;
    ok( scalar @pads, 3 );
    ok( $pads[0], $pads[2] );
    ok( $pads[1] =~ tr/|/|/, length( $pads[0] ) );

   # There is an issue to be resolved still with fetching features by name for CRAM files
   # so opening a new object each time to ensure the tests take place.
   $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );

    my @f = $hts->features( -name => 'SRR507778.24982' );
    ok( scalar @f, 2 );

    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    @f = $hts->features(
        -filter => sub {
            my $a = shift;
            return 1 if $a->display_name eq 'SRR507778.24982';
        } );
    ok( scalar @f, 2 );

    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    @f = $hts->features(
        -seq_id => 'XIII',
        -filter => sub {
            my $a = shift;
            return 1 if $a->display_name =~ /^SRR507778/;
        } );
    ok( scalar @f, 3296 );

    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    @f = $hts->features(
        -filter => sub {
            my $a = shift;
            return 1 if $a->display_name =~ /^SRR507778/;
        } );
    ok( scalar @f, 50000 );

    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    @f = $hts->features( -name => 'SRR507778.24982',
                         -tags => { FIRST_MATE => 1 } );
    ok( scalar @f, 1 );

    # try iteration
    my $i =
      $hts->get_seq_stream( -seq_id => 'XIII', -start => 200, -end => 10000 );
    ok($i);
    my $count = 0;
    while ( $i->next_seq ) { $count++ }
    ok( $count, 65 );

    # try tam fh retrieval
    my $fh = $hts->get_seq_fh( -seq_id => 'XIII', -start => 200, -end => 10000, );
    $count = 0;
    $count++ while <$fh>;
    ok( $count, 65 );
    $fh->close;

    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    $i = $hts->get_seq_stream();    # all features!
    ok($i);
    $count = 0;
    while ( $i->next_seq ) { $count++ }
    ok( $count, 50000 );

    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    $i = $hts->get_seq_stream( -max_features => 200, -seq_id => 'XIII' );
    ok($i);
    $count = 0;
    while ( $i->next_seq ) { $count++ }
    ok( $count,                   200 );
    ok( $hts->last_feature_count, 3296 );

    # try the read_pair aggregation
    my @pairs = $hts->features( -type => 'read_pair', -seq_id => 'XIII' );
    ok( scalar @pairs, 1682 );

    # test high level API version of pileup
    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    my @coverage = $hts->features( -type => 'coverage', -seq_id => 'XIII' );
    my ($c) = $coverage[0]->get_tag_values('coverage');
    ok( $coverage[0]->type, "coverage:924431" ) ;
    my %matches;
    my $fetch_back = sub {
        my ( $seqid, $pos, $p ) = @_;
        my $r = $hts->segment( $seqid, $pos, $pos )->dna;
        for my $pileup (@$p) {
            my $a    = $pileup->alignment;
            my $qpos = $pileup->qpos;
            my $dna  = $a->query->dna;
            my $base =
              $pileup->indel == 0 ? substr( $dna, $qpos, 1 ) :
              $pileup->indel > 0 ? '*' :
              '-';
            $matches{matched}++ if $r eq $base;
            $matches{total}++;
        }
    };
    $hts->pileup( 'XIII:1-1000', $fetch_back );
    ok( $matches{matched}, 115 );
    ok( $matches{total}, 211 );

    # test filtered coverage (no filtering)
    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    @coverage = $hts->features( -type => 'coverage', -seq_id => 'XIII',
                                -filter => sub { return 1; } );
    ($c) = $coverage[0]->get_tag_values('coverage');
    ok( $c );
    ok( $c->[30],     3 );
    ok( $c->[33148],  3 );
    ok( $c->[924333], 5 );
    ok( $c->[924351], 4 );
    ok( $coverage[0]->type, "coverage:924431" );

    # test filtered coverage (really filtering)
    $hts = Bio::DB::HTS->new(
                                 -fasta => $fastafile,
                                 -bam          => $cramfile,
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    @coverage = $hts->features( -type => 'coverage', -seq_id => 'XIII',
                                -filter => sub { my $a = shift;
                           return defined $a->start && $a->start < 100000; } );
    ($c) = $coverage[0]->get_tag_values('coverage');
    ok( $c );
    ok( $c->[30],     3 );
    ok( $c->[33148],  3 );
    ok( $c->[924333], 0 );
    ok( $c->[924351], 0 );
    ok( $coverage[0]->type, "coverage:924431" );

} ## end for my $use_fasta ( 0, ...)

exit 0;

sub reversec {
    my $dna = shift;
    $dna =~ tr/gatcGATC/ctagCTAG/;
    return scalar reverse $dna;
}

__END__
