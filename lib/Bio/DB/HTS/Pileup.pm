
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

Bio::DB::HTS::Pileup -- Object passed to pileup() callback

=head1 SYNOPSIS

See L<Bio::DB::HTS/The generic fetch() and pileup() methods> for how
this object is passed to pileup callbacks.

=head1 DESCRIPTION

A Bio::DB::HTS::Pileup object (or a Bio::DB::HTS::PileupWrapper
object) is passed to the callback passed to the Bio::DB::HTS->pileup()
method for each column in a sequence alignment. The only difference
between the two is that the latter returns the more convenient
Bio::DB::HTS::AlignWrapper objects in response to the alignment()
method, at the cost of some performance loss.

=head2 Methods

=over 4

=item $alignment = $pileup->alignment

Return the Bio::DB::HTS::Alignment or Bio::DB::HTS::AlignWrapper
object representing the aligned read.

=item $alignment = $pileup->b

This method is an alias for alignment(). It is available for
compatibility with the C API.

=item $qpos = $pileup->qpos

Return the position of this aligned column in read coordinates, using
zero-based coordinates.

=item $pos  = $pileup->pos

Return the position of this aligned column in read coordinates, using
1-based coordinates.

=item $indel = $pileup->indel

If this column is an indel, return a positive integer for an insertion
relative to the reference, a negative integer for a deletion relative
to the reference, or 0 for no indel at this column.

=item $is_del = $pileup->is_del

True if the base on the padded read is a deletion.

=item $level  = $pileup->level

If pileup() or fast_pileup() was invoked with the "keep_level" flag,
then this method will return a positive integer indicating the level
of the read in a printed multiple alignment.

=item $pileup->is_head

=item $pileup->is_tail

These fields are defined in bam.h but their interpretation is obscure.

=back

=head1 AUTHOR

Rishi Nag E<lt>rishi@ebi.ac.uk<gt>


=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::HTS>, L<Bio::DB::HTS::Alignment>, L<Bio::DB::HTS::Constants>

=cut

# documentation only

package Bio::DB::HTS::Pileup;

use strict;
use warnings;

$Bio::DB::HTS::Pileup::VERSION = '2.11';

1;
