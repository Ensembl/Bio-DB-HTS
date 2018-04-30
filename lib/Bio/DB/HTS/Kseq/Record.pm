package Bio::DB::HTS::Kseq::Record;

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

Bio::DB::HTS::Kseq::Record -- Entry from a Kseq iterator

=head1 SYNOPSIS

  while(my $r = $iter->next_seq()) {
    say $r->name;
    say $r->desc;
    say $r->seq;
    say $r->qual;
  }

=head2 METHODS

=over 8

=item C<name>

The name of a FASTA/Q record

=item C<desc>

The description from a FASTA/Q record

=item C<seq>

The sequence from a FASTA/Q record

=item C<qual>

The quality string from a FASTA/Q record

=back

=cut

use strict;
use warnings;

$Bio::DB::HTS::Kseq::Record::VERSION = '2.11';

sub name {
  return $_[0]->{name};
}

sub desc {
  return $_[0]->{desc};
}

sub seq {
  return $_[0]->{seq};
}

sub qual {
  return $_[0]->{qual};
}

1;
