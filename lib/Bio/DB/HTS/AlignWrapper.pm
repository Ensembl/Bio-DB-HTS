
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

=head1 NAME

Bio::DB::HTS::AlignWrapper -- Add high-level methods to Bio::DB::HTS::Alignment

=head1 SYNOPSIS

See L<Bio::DB::HTS::Alignment>.

=head1 DESCRIPTION

This is a wrapper around Bio::DB::HTS::Alignment that adds the
following high-level methods. These are described in detail in
L<Bio::DB::HTS::Alignment/High-level Bio::DB::HTS::Alignment methods>.

 add_segment()         add a new subfeature to split alignments
 get_SeqFeatures()     fetch subfeatures from split alignments
 split_splices()       process cigar strings to produce split alignments
 expand_flags()        return true if flags should be expanded into tags
 seq_id()              return human-readable reference sequence name
 seq()                 return Bio::PrimarySeq object for reference sequence
 subseq()              return a subsequence across the indicated range
 mate_seq_id()         return human-readable mate reference sequence name
 dna()                 return the DNA of the reference sequence
 tam_line()            return the text representation of the alignment
 attributes()          synonym for get_tag_values()
 get_all_tags()        return all tag names
 get_tag_values()      return the values of the given tag
 has_tag()             return true if the given tag is defined

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::HTS>, L<Bio::DB::HTS::Constants>

=head1 AUTHOR

Rishi Nag E<lt>rishi@ebi.ac.uk<gt>


=cut

package Bio::DB::HTS::AlignWrapper;
$Bio::DB::HTS::AlignWrapper::VERSION = '2.11';

use strict;
use warnings;

use Bio::DB::HTS::Constants;

our $AUTOLOAD;
use Carp 'croak';

sub new {
    my $package = shift;
    my ( $align, $sam ) = @_;

    my $self = bless { sam => $sam, align => $align }, ref $package || $package;

    $self->add_segment( $self->split_splices )
      if $sam->split_splices && $align->cigar_str =~ /N/;

    return $self;
}

sub AUTOLOAD {
    my ( $pack, $func_name ) = $AUTOLOAD =~ /(.+)::([^:]+)$/;
    return if $func_name eq 'DESTROY';

    no strict 'refs';
    $_[0] or die "autoload called for non-object symbol $func_name";
    croak qq(Can't locate object method "$func_name" via package "$pack")
      unless $_[0]->{align}->can($func_name);

    *{"${pack}::${func_name}"} = sub { shift->{align}->$func_name(@_) };

    shift->$func_name(@_);
}

sub score { shift->{align}->qual }

sub can {
    my $self = shift;
    return 1 if $self->SUPER::can(@_);
    return ref $self && $self->{align}->can(@_);
}

sub add_segment {
    my $self    = shift;
    my @subfeat = @_;
    $self->{segments} ||= [];
    push @{ $self->{segments} }, @subfeat;
}

sub get_SeqFeatures {
    my $self = shift;
    return unless $self->{segments};
    return @{ $self->{segments} };
}

sub split_splices {
    my $self  = shift;
    my $cigar = $self->cigar_array;
    my @results;

    my $start         = 0;
    my $end           = 0;
    my $skip          = 0;
    my $partial_cigar = '';

    # in case sequence is missing?
    my $qseq = $self->qseq;
    $qseq ||= 'N' x $self->length;

    for my $op ( @$cigar, [ 'N', 0 ] ) {
        my ( $operation, $count ) = @$op;

        if ( $operation eq 'N' ) {
            my $s = $self->start + $start + $skip;
            my $e = $self->start + $end - 1 + $skip;
            my $f = Bio::DB::HTS::SplitAlignmentPart->new(
                -name   => $self->display_name,
                -start  => $s,
                -end    => $e,
                -seq_id => $self->seq_id,
                -strand => $self->strand,
                -seq    => [ $self, $start + $skip, $end - $start ]
                ,    # deferred rendering
                -type => $self->type );

            $f->hit( -name   => $self->display_name,
                     -seq_id => $self->display_name,
                     -start  => $start + 1,
                     -end    => $end,
                     -strand => $self->strand,
                     -seq    => substr( $qseq, $start, $end - $start ), );
            $f->cigar_str($partial_cigar);
            $partial_cigar = '';

            push @results, $f;
            $start += $end - $start;
        } ## end if ( $operation eq 'N')
        else {
            $partial_cigar .= "$operation$count";
        }
        $end  += $count if $operation =~ /^[MDSHP=X]/i;
        $skip += $count if $operation eq 'N';
        if ( $operation eq 'H' and $start == 0 ) {
            $qseq = 'N' x $count . $qseq;
        }
    } ## end for my $op ( @$cigar, [...])
    return @results;
} ## end sub split_splices

sub expand_flags {
    shift->{sam}->expand_flags(@_);
}

sub query {
    my $self = shift;
    return Bio::DB::HTS::Query->new($self);
}

sub target {
    my $self = shift;
    return Bio::DB::HTS::Target->new($self);
}

sub hit { shift->target(@_); }

sub seq_id {
    my $self = shift;
    my $tid  = $self->tid;
    return if $tid < 0 || $self->unmapped;
    $self->{sam}->target_name($tid);
}

sub mate_seq_id {
    my $self = shift;
    my $tid  = $self->mtid;
    return if $tid < 0 || $self->munmapped;
    $self->{sam}->target_name($tid);
}

sub abs_ref   { shift->seq_id }
sub abs_start { shift->start }
sub abs_end   { shift->end }
sub low       { shift->start }
sub high      { shift->end }
sub type      { shift->primary_tag }
sub method    { shift->primary_tag }
sub source    { return shift->source_tag; }
sub name      { shift->qname }
sub class     { shift->primary_tag }

sub seq {
    my $self = shift;
    my $dna  = $self->dna;
    return Bio::PrimarySeq->new( -seq => $dna, -id => $self->seq_id );
}

sub subseq {
    my $self = shift;
    my ( $start, $end ) = @_;
    $start = 1           if $start < 1;
    $end   = $self->high if $end > $self->high;
    my $dna = $self->dna;
    return Bio::PrimarySeq->new(
                        -seq => substr( $dna, $start - 1, $end - $start + 1 ) );
}

sub padded_alignment {
    my $self     = shift;
    my $cigar    = $self->cigar_array;
    my $real_ref = 0;
    $real_ref = 1 if ( $self->{sam}->force_refseq || !$self->has_tag('MD') );

    my $sdna = $self->dna;
    my $tdna = $self->query->dna;

    my ( $pad_source, $pad_target, $pad_match, $char_source, $char_target );
    for my $event (@$cigar) {
        my ( $op, $count ) = @$event;
        if ( $op eq 'I' || $op eq 'S' ) {
            $pad_source .= '-' x $count;
            $pad_target .= substr( $tdna, 0, $count, '' );
            $pad_match  .= ' ' x $count;
        }
        elsif ( $op eq 'D' ) {
            $pad_source .= substr( $sdna, 0, $count, '' );
            $pad_target .= '-' x $count;
            $pad_match  .= ' ' x $count;
        }
        elsif ( $op eq 'N' ) {
            if ($real_ref) {
                $pad_source .= substr( $sdna, 0, $count, '' );
            }
            else {
                $pad_source .= '-' x $count;
            }
            $pad_target .= '-' x $count;
            $pad_match  .= ' ' x $count;
        }
        elsif ( $op eq 'P' ) {
            $pad_source .= '*' x $count;
            $pad_target .= '*' x $count;
            $pad_match  .= ' ' x $count;
        }
        elsif ( $op eq 'H' ) {
            # nothing needs to be done in this case
        }
        else {    # everything else is assumed to be a match
            while ( $count-- > 0 ) {
                $char_source = substr( $sdna, 0, 1, '' );
                $char_target = substr( $tdna, 0, 1, '' );
                $pad_source .= $char_source;
                $pad_target .= $char_target;
                $pad_match  .= $char_source eq $char_target ? '|' : ' ';
            }
        }
    } ## end for my $event (@$cigar)
    return ( $pad_source, $pad_match, $pad_target );
} ## end sub padded_alignment

sub dna {
    my $self = shift;

    my $sam   = $self->{sam};
    my $force = $sam->force_refseq;
    if ( !$force && ( my $md = $self->get_tag_values('MD') ) )
    {    # try to use MD string
        my $qseq = $self->qseq;

        #preprocess qseq using cigar array
        my $cigar = $self->cigar_array;
        my $seq   = '';
        for my $op (@$cigar) {
            my ( $operation, $count ) = @$op;
            if ( $operation eq 'M' || $operation eq '=' || $operation eq 'X' ) {
                $seq .= substr( $qseq, 0, $count, '' ); # include these residues
            }
            elsif ( $operation eq 'S' or $operation eq 'I' ) {
                substr( $qseq, 0, $count, '' )
                  ;    # skip soft clipped and inserted residues
            }
        }

        my $start = 0;
        my $result;
        while ( $md =~ /(\d+)|\^([gatcn]+)|([gatcn]+)/ig ) {
            if ( defined $1 ) {
                $result .= substr( $seq, $start, $1 );
                $start += $1;
            }
            elsif ($2) {
                $result .= $2;
            }
            elsif ($3) {
                $result .= $3;
                $start += length $3;
            }
        }
        return $result;
    } ## end if ( !$force && ( my $md...))

    else {
        return $self->{sam}->seq( $self->seq_id, $self->start, $self->end );
    }
} ## end sub dna

sub tseq {
    shift->dna(@_);
}

sub attributes {
    my $self = shift;
    my $tag  = shift;
    if ( defined $tag ) {
        return $self->get_tag_values($tag);
    }
    else {
        return map { $_ => $self->get_tag_values($_) } $self->get_all_tags;
    }
}

sub get_all_tags {
    my $self = shift;
    return $self->{align}->get_all_tags(@_) if $self->expand_flags;
    return ( $self->aux_keys, 'FLAGS' );
}

sub get_tag_values {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;

    return $self->{align}->get_tag_values($tag) if $self->expand_flags;
    if ( $tag eq 'FLAGS' ) {
        $self->flag_str;
    }
    else {
        $self->aux_get($tag);
    }
}

sub has_tag {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    $self->{align}->get_tag_values($tag) if $self->expand_flags;
    if ( $tag eq 'FLAGS' ) {
        return 1;
    }
    else {
        my %keys = map { $_ => 1 } $self->aux_keys;
        return exists $keys{ uc $tag };
    }
}

sub gff_string { shift->gff3_string(@_) }

sub gff3_string {
    my $self      = shift;
    my $recurse   = shift;
    my $parent_id = shift;

    my $group = $self->format_attributes($parent_id);
    my $name  = $self->name;
    my $id    = $self->primary_id;

    my $class  = $self->class;
    my $strand = ( '-', '.', '+' )[ $self->strand + 1 ];
    my $p      = join( "\t",
                  $self->seq_id || '.',
                  $self->source || '.',
                  $self->method || '.',
                  $self->start  || '.',
                  $self->stop   || '.',
                  defined( $self->score ) ? $self->score : '.',
                  $strand || '.',
                  defined( $self->phase ) ? $self->phase : '.',
                  $group || '' );
    my @rsf = $self->get_SeqFeatures;
    return join( "\n", $p, map { $_->gff3_string($id) } @rsf );
}

sub phase { return }

sub escape {
    my $self     = shift;
    my $toencode = shift;
    $toencode =~
      s/([^a-zA-Z0-9_.:?^*\(\)\[\]@!+-])/uc sprintf("%%%02x",ord($1))/eg;
    $toencode;
}

sub format_attributes {
    my $self      = shift;
    my $parent_id = shift;

    my @tags = $self->get_all_tags;
    my @result;
    for my $t (@tags) {
        my @values = $self->each_tag_value($t);
        push @result, join '=', $self->escape($t),
          join( ',', map { $self->escape($_) } @values )
          if @values;
    }
    my $id = $self->escape( $self->primary_id );

    my $name = $self->display_name;
    unshift @result, "ID=" . $id                    if defined $id;
    unshift @result, "Parent=" . $parent_id         if defined $parent_id;
    unshift @result, "Name=" . $self->escape($name) if defined $name;
    return join ';', @result;
}

sub tam_line {
    my $self = shift;
    return
      join( "\t",
            $self->qname,
            $self->flag,
            $self->tid >= 0 ? $self->{sam}->target_name( $self->tid ) : '*',
            $self->pos + 1,
            $self->qual,
            $self->cigar_str || '*',
            $self->mtid >= 0 ? ( $self->mtid == $self->tid ? '=' :
                                 $self->{sam}->target_name( $self->mtid ) ) :
              '*',
            $self->mpos + 1,
            $self->isize,
            $self->qseq,
            join( '', map { chr( $_ + 33 ) } $self->qscore ),
            $self->aux || () );
}

package Bio::DB::HTS::SplitAlignmentPart;

$Bio::DB::HTS::SplitAlignmentPart::VERSION = '2.11';

use base 'Bio::SeqFeature::Lite';

sub hit {
    my $self = shift;
    my $d    = $self->{hit};
    $self->{hit} = Bio::SeqFeature::Lite->new(@_) if @_;
    return $d;
}

sub seq {
    my $self = shift;
    my $seq  = $self->{seq};
    return $self->SUPER::seq() unless ref $seq;
    return substr( $seq->[0]->dna, $seq->[1], $seq->[2] );
}

sub Bio::SeqFeature::Lite::subseq {
    my $self = shift;
    my ( $start, $end ) = @_;
    $start = 1           if $start < 1;
    $end   = $self->high if $end > $self->high;
    return Bio::PrimarySeq->new(
                  -seq => substr( $self->dna, $start - 1, $end - $start + 1 ) );
}

sub cigar_str {
    my $self = shift;
    my $d    = $self->{cigar_str};
    $self->{cigar_str} = shift if @_;
    $d;
}

1;
