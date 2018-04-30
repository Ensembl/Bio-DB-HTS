
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

package Bio::DB::HTS::ReadIterator;
$Bio::DB::HTS::ReadIterator::VERSION = '2.11';

use strict;
use warnings;

sub new {
    my $self = shift;
    my ( $sam, $hts_file, $filter, $header ) = @_;
    return
      bless { sam      => $sam,
              hts_file => $hts_file,
              filter   => $filter,
              header   => $header, },
      ref $self || $self;
}

sub next_seq {
    my $self   = shift;
    my $header = $self->{header};
    while ( my $b = $self->{hts_file}->read1($header) ) {
        return Bio::DB::HTS::AlignWrapper->new( $b, $self->{sam} )
          if $self->{filter}->($b);
    }
    return;
}

1;
