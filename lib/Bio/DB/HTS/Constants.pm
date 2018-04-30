
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

Bio::DB::HTS::Constants -- Constants for use with SAM/BAM

=head1 SYNOPSIS

 use Bio::DB::HTS::Constants;
 my $pad_flag = BAM_CPAD;

=head1 DESCRIPTION

This module exports several constants for use with the SAM/BAM
module. See the SAM documentation for their interpretation.

=over 4

=item Cigar operations

  BAM_CIGAR_SHIFT
  BAM_CIGAR_MASK
  BAM_CMATCH
  BAM_CINS
  BAM_CDEL
  BAM_CREF_SKIP
  BAM_CSOFT_CLIP
  BAM_CHARD_CLIP
  BAM_CPAD

=item FLAGS

A hashref that maps flag values to human-readable names. For example:

 FLAGS->{0x0008} == 'M_UNMAPPED'

=item RFLAGS

The reverse of FLAGS:

 RFLAGS->{M_UNMAPPED} == 0x0008

=back

=head1 AUTHOR

Rishi Nag E<lt>rishi@ebi.ac.uk<gt>

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::HTS>, L<Bio::DB::Bam::Alignment>

=cut

package Bio::DB::HTS::Constants;
$Bio::DB::HTS::Constants::VERSION = '2.11';

use strict;
use warnings;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(CIGAR_SYMBOLS BAM_CIGAR_SHIFT BAM_CIGAR_MASK
  BAM_CMATCH BAM_CINS BAM_CDEL BAM_CREF_SKIP
  BAM_CSOFT_CLIP BAM_CHARD_CLIP BAM_CPAD FLAGS RFLAGS);
our @EXPORT_OK = @EXPORT;

use constant CIGAR_SYMBOLS   => [qw(M I D N S H P = X)];
use constant BAM_CIGAR_SHIFT => 4;
use constant BAM_CIGAR_MASK  => ( 1 << BAM_CIGAR_SHIFT ) - 1;
use constant BAM_CMATCH      => 0;
use constant BAM_CINS        => 1;
use constant BAM_CDEL        => 2;
use constant BAM_CREF_SKIP   => 3;
use constant BAM_CSOFT_CLIP  => 4;
use constant BAM_CHARD_CLIP  => 5;
use constant BAM_CPAD        => 6;

use constant FLAGS => { 0x0001 => 'PAIRED',
                        0x0002 => 'MAP_PAIR',
                        0x0004 => 'UNMAPPED',
                        0x0008 => 'M_UNMAPPED',
                        0x0010 => 'REVERSED',
                        0x0020 => 'M_REVERSED',
                        0x0040 => 'FIRST_MATE',
                        0x0080 => 'SECOND_MATE',
                        0x0100 => 'NOT_PRIMARY',
                        0x0200 => 'QC_FAILED',
                        0x0400 => 'DUPLICATE',
                        0x0800 => 'SUPPLEMENTARY', };
use constant RFLAGS => { reverse %{ FLAGS() } };

1;
