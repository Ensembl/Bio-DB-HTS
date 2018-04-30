
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

=cut

=head1 NAME

Bio::DB::HTS::Kseq - Bindings to Kseq

=head1 DESCRIPTION

Bindings to the Kseq library for iterating through a locally held FASTA/FASTQ file very quickly. Supports compressed and uncompressed files alongside filehandles.

=head1 SYNOPSIS

  my $kseq = Bio::DB::HTS::Kseq->new('path/to/file');
  my $iter = $kseq->iterator();
  while(my $r = $iter->next_seq()) {
    say $r->name;
    say $r->desc;
    say $r->seq;
    say $r->qual;
  }

  # Allowing the object to go out of scope will close down all file handles

=head2 METHODS

=over 8

=item C<new>

  my $kseq = Bio::DB::HTS::Kseq->new('path/to/file');

Returns an instance of this object from a file path.

=item C<newfh>

  open my $fh, '<', 'path' or die "Cannot open path: $!";
  binmode $fh;
  my $kseq = Bio::DB::HTS::Kseq->newfh($fh);

Returns an instance of this object from an opened file handle. This supports any known Perl glob/filehandle type

=item C<iterator>

  my $iter = $kseq->iterator();

Returns the kseq iterator object

=back

=cut

package Bio::DB::HTS::Kseq;

use Bio::DB::HTS; #load the XS
use Bio::DB::HTS::Kseq::Record;

$Bio::DB::HTS::Kseq::VERSION = '2.11';

use strict;
use warnings;

1;
