
=head1 LICENSE

Copyright [2015-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 AUTHOR

Rishi Nag E<lt>rishi@ebi.ac.uk<gt>

=head1 NAME

Bio::DB::HTS::Query -- Object representing the query portion of a BAM/SAM alignment

=head1 SYNOPSIS

Given an alignment retrieved from a Bio::DB::HTS database,

 my $query = $alignment->query;

 my $name   = $query->display_name;
 my $start  = $query->start;
 my $end    = $query->end;
 my $dna    = $query->dna;    # dna string
 my $seq    = $query->seq;    # Bio::PrimarySeq object
 my @scores = $query->qscore; # quality score

=head1 DESCRIPTION

This is a simple Bio::SeqFeatureI object that represents the query
part of a SAM alignment.

=head2 Methods

=over 4

=cut

package Bio::DB::HTS::Query;
$Bio::DB::HTS::Query::VERSION = '2.11';

use strict;
use warnings;

use Bio::DB::HTS;
use Bio::DB::HTS::Constants
  qw(CIGAR_SYMBOLS BAM_CREF_SKIP BAM_CSOFT_CLIP BAM_CHARD_CLIP);

use constant CIGAR_SKIP => { CIGAR_SYMBOLS->[BAM_CREF_SKIP]  => 1,
                             CIGAR_SYMBOLS->[BAM_CSOFT_CLIP] => 1,
                             CIGAR_SYMBOLS->[BAM_CHARD_CLIP] => 1 };

sub new {
    my $self      = shift;
    my $alignment = shift;
    bless \$alignment, ref $self || $self;
}

=item $seqid = $query->seq_id

The name of the read.

=cut

sub seq_id {
    my $self = shift;
    $$self->qname;
}

=item $name = $query->name

The read name (same as seq_id in this case).

=cut

sub name {
    my $self = shift;
    $$self->qname;
}

=item $name = $query->display_name

The read display_name (same as seq_id in this case).

=cut

sub display_name { shift->name }

=item $tag = $query->primary_tag

The string "match".

=cut

sub primary_tag { ${ shift() }->primary_tag }

=item $tag = $query->source_tag

The string "sam/bam".

=cut

sub source_tag { ${ shift() }->source_tag }

=item $start = $query->start

The start of the match in read coordinates.

=cut

sub start {
    my $self = shift;
    return $self->low;
}

=item $end = $query->end

The end of the match in read coordinates;

=cut

sub end {
    my $self = shift;
    return $self->high;
}

sub low {
    my $self       = shift;
    my $cigar_arry = $$self->cigar_array;
    my $start      = 1;
    for my $c (@$cigar_arry) {
        next if CIGAR_SYMBOLS->[BAM_CHARD_CLIP] eq $c->[0];
        last unless CIGAR_SKIP->{ $c->[0] };
        $start += $c->[1];
    }
    $start;
}

sub high {
    my $self       = shift;
    my $len        = $$self->cigar2qlen;
    my $cigar_arry = $$self->cigar_array;

    # alignment stops at first non-clip CIGAR position
    my $i = $len - 1;
    for my $c ( reverse @$cigar_arry ) {
        next if CIGAR_SYMBOLS->[BAM_CHARD_CLIP] eq $c->[0];
        last unless CIGAR_SKIP->{ $c->[0] };
        $len -= $c->[1];
    }
    return $len;
}

=item $len = $query->length

The length of the read.

=cut

sub length {
    my $self = shift;
    $self->high - $self->low + 1;
    #    $$self->cigar2qlen;
}

=item $seq = $query->seq

A Bio::PrimarySeq representing the read sequence in REFERENCE
orientation.

=cut

sub seq {
    my $self = shift;
    my $dna  = $self->dna;
    return Bio::PrimarySeq->new( -seq => $dna, -id => $$self->qname );
}

=item $scores = $query->qscore

The read quality scores. In a list context, a list of integers equal
in length to the read sequence length. In a scalar context, an array
ref. The qscores are in REFERENCE sequence orientation.

=cut

sub qscore {
    my $self   = shift;
    my @qscore = $$self->qscore;
    return wantarray ? @qscore : \@qscore;
}

=item $dna = $query->dna

The DNA string in reference sequence orientation.

=cut

sub dna {
    my $self = shift;
    return $$self->qseq || ( 'N' x $self->length );
}

=item $strand = $query->strand

If the query was reversed to align it, -1. Otherwise +1.

=cut

sub strand {
    my $self = shift;
    return $$self->reversed ? -1 : 1;
}

=item $seq = $query->subseq($start,$end)

Return a Bio::PrimarySeq object representing the requested subsequence
on the read.

=cut

sub subseq {
    my $self = shift;
    my ( $start, $end ) = @_;
    $start = 1           if $start < 1;
    $end   = $self->high if $end > $self->high;
    ( $end, $start ) = ( $start, $end ) if $start > $end;
    return Bio::PrimarySeq->new(
                  -seq => substr( $self->dna, $start - 1, $end - $start + 1 ) );
}

1;

=back

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::HTS>, L<Bio::DB::HTS::Alignment>, L<Bio::DB::HTS::Constants>

=cut
