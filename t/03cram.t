#-*-Perl-*-

use strict;
use ExtUtils::MakeMaker;
use File::Temp qw(tempfile);
use FindBin '$Bin';
use constant TEST_COUNT => 135;

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

# low level tests (defined in lib/Bio/DB/HTS.pm)
{
    my $cramfile = "$Bin/data/yeast.sorted.cram";
    my $hts_file     = Bio::DB::HTSfile->open($cramfile);
    ok($hts_file);

    my $header  = $hts_file->header_read();
    my $targets = $header->n_targets;
    ok( $targets, 17, "Targets Found" );

    my $target_names = $header->target_name;
    ok($target_names);
    ok( scalar @$target_names );
    ok( $target_names->[0],    'II' );

    my $target_lens = $header->target_len;
    ok($target_lens);
    ok( scalar @$target_lens, 2 );
    ok( $target_lens->[0],    1575 );

    my $text = $header->text;
    ok( length $text > 0 );

    my $c = "\@CO\tThis is a comment\n";
    $header->text($c);
    ok( $header->text, $c );

    my $fai = Bio::DB::HTS::Fai->open("$Bin/data/ex1.fa");
    my $seq = $fai->fetch('seq2:51-1000');
    ok( length $seq, 950 );

    my $count;
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

# high level tests (defined in lib/Bio/DB/HTS.pm)
for my $use_fasta ( 0, 1 ) {
    my $hts = Bio::DB::HTS->new( -fasta        => "$Bin/data/ex1.fa",
                                 -bam          => "$Bin/data/ex1.bam",
                                 -expand_flags => 1,
                                 -autoindex    => 1,
                                 -force_refseq => $use_fasta, );
    ok($hts);
    ok( $hts->n_targets, 2 );

    ok( $hts->length('seq1'), 1575 );
    ok( join $hts->seq_ids,   'seq1 seq2' );

    my $seg = $hts->segment('seq1');
    ok($seg);
    ok( $seg->length, 1575 );
    my $seq = $seg->seq;
    ok( $seq->isa('Bio::PrimarySeq') );
    ok( length $seq->seq, 1575 );

    my $dummy = eval {
        Bio::DB::HTS->new( -fasta => "invalid_path.txt",
                           -bam   => "invalid_path.txt" );
    };
    ok( $dummy, undef );
    ok( $@ =~ /does not exist/ );

    my @alignments =
      $hts->get_features_by_location( -seq_id => 'seq2',
                                      -start  => 500,
                                      -end    => 800 );
    ok( scalar @alignments,     442 );
    ok( $alignments[0]->seq_id, 'seq2' );

    ok( scalar @{ $alignments[0]->qscore }, length $alignments[0]->dna );

    my @keys = $alignments[0]->get_all_tags;
    ok( scalar @keys,                         18 );
    ok( $alignments[0]->get_tag_values('MF'), 18 );

    my %att = $alignments[0]->attributes;
    ok( scalar( keys %att ),       18 );
    ok( $alignments[0]->cigar_str, '35M' );

    $hts->expand_flags(0);
    @keys = $alignments[0]->get_all_tags;
    ok( scalar @keys, 7 );

    my $query = $alignments[0]->query;
    ok($query);
    ok( $query->start,          1 );
    ok( $query->end,            35 );
    ok( $query->length,         35 );
    ok( $query->dna,            $alignments[0]->dna );
    ok( $alignments[0]->strand, -1 );
    ok( $query->strand,         -1 );

    my $target = $alignments[0]->target;
    ok($target);
    ok( $target->start,  35 );
    ok( $target->end,    1 );
    ok( $target->length, 35 );
    ok( $target->dna,    reversec( $alignments[0]->dna ) );

    ok( $alignments[0]->get_tag_values('FLAGS'), $alignments[0]->flag_str );

    my @pads = $alignments[0]->padded_alignment;
    ok( @pads,               3 );
    ok( $pads[0],            $pads[2] );
    ok( $pads[1] =~ tr/|/|/, length( $pads[0] ) );

    my @f = $hts->features( -name => 'EAS114_45:2:1:1140:1206' );
    ok( scalar @f, 2 );

    @f = $hts->features(
        -filter => sub {
            my $a = shift;
            return 1 if $a->display_name eq 'EAS114_45:2:1:1140:1206';
        } );
    ok( scalar @f, 2 );

    @f = $hts->features(
        -seq_id => 'seq2',
        -filter => sub {
            my $a = shift;
            return 1 if $a->display_name =~ /^EAS114/;
        } );
    ok( scalar @f, 306 );
    @f = $hts->features(
        -filter => sub {
            my $a = shift;
            return 1 if $a->display_name =~ /^EAS114/;
        } );
    ok( scalar @f, 534 );

    @f = $hts->features( -name => 'EAS114_45:2:1:1140:1206',
                         -tags => { FIRST_MATE => 1 } );
    ok( scalar @f, 1 );

    # try iteration
    my $i =
      $hts->get_seq_stream( -seq_id => 'seq2', -start => 500, -end => 800 );
    ok($i);
    my $count = 0;
    while ( $i->next_seq ) { $count++ }
    ok( $count, 442 );

    # try tam fh retrieval
    my $fh = $hts->get_seq_fh( -seq_id => 'seq2', -start => 500, -end => 800, );
    $count = 0;
    $count++ while <$fh>;
    ok( $count, 442 );
    $fh->close;

    $i = $hts->get_seq_stream();    # all features!
    ok($i);
    $count = 0;
    while ( $i->next_seq ) { $count++ }
    ok( $count, 3307 );

    $i = $hts->get_seq_stream( -max_features => 200, -seq_id => 'seq1' );
    ok($i);
    $count = 0;
    while ( $i->next_seq ) { $count++ }
    ok( $count,                   200 );
    ok( $hts->last_feature_count, 1482 );

    # try the read_pair aggregation
    my @pairs = $hts->features( -type => 'read_pair', -seq_id => 'seq2' );
    ok( scalar @pairs, 939 );

    # try coverage
    my @coverage = $hts->features( -type => 'coverage', -seq_id => 'seq2' );
    ok( scalar @coverage, 1 );
    my ($c) = $coverage[0]->get_tag_values('coverage');
    ok($c);
    ok( $c->[0],            3 );
    ok( $c->[1],            4 );
    ok( $coverage[0]->type, "coverage:1584" );

    # test high level API version of pileup
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

    $hts->pileup( 'seq2:1-100', $fetch_back );
    ok( $matches{matched}/$matches{total} > 0.99 );
} ## end for my $use_fasta ( 0, ...)

exit 0;

sub reversec {
    my $dna = shift;
    $dna =~ tr/gatcGATC/ctagCTAG/;
    return scalar reverse $dna;
}

__END__
