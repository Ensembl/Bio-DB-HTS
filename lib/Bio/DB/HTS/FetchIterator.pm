
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

=cut

package Bio::DB::HTS::FetchIterator;
$Bio::DB::HTS::FetchIterator::VERSION = '2.11';

use strict;
use warnings;

sub new {
    my $self  = shift;
    my $list  = shift;
    my $total = shift;
    $total ||= @$list;
    return bless { list => $list, total => $total }, ref $self || $self;
}

sub next_seq {
    my $self = shift;
    return shift @{ $self->{list} };
}

sub total {
    shift->{total};
}

1;
