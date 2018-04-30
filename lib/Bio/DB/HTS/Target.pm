
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

Bio::DB::HTS::Target -- Object representing the query portion of a BAM/SAM alignment in NATIVE alignment

=head1 SYNOPSIS

This is identical to Bio::DB::HTS::Query, except that the dna, qscores
and start and end positions are all given in the orientation in which
the read was sequenced, not in the oreintation in which it was
aligned.

=head1 AUTHOR

Rishi Nag E<lt>rishi@ebi.ac.uk<gt>

=cut

package Bio::DB::HTS::Target;
$Bio::DB::HTS::Target::VERSION = '2.11';

use strict;
use warnings;

use base 'Bio::DB::HTS::Query';

sub dna {
    my $self = shift;
    my $qseq = $self->SUPER::dna;
    return $$self->strand > 0 ? $qseq : reversec($qseq);
}

sub qscore {
    my $self   = shift;
    my @qscore = $$self->qscore;
    @qscore = reverse @qscore if $$self->strand < 0;
    return wantarray ? @qscore : \@qscore;
}

sub start {
    my $self = shift;
    return $self->strand > 0 ? $self->low : $self->high;
}

sub end {
    my $self = shift;
    return $self->strand > 0 ? $self->high : $self->low;
}

# sub strand { 1 }

sub reversec {
    my $dna = shift;
    $dna =~ tr/gatcGATC/ctagCTAG/;
    return scalar reverse $dna;
}

1;
