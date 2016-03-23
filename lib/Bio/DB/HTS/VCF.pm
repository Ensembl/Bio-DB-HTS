=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Alex Hodgkins
Rishi Nag E<lt>rishi@ebi.ac.ukE<gt>

=cut

package Bio::DB::HTS::VCF;
$Bio::DB::HTS::VCF::VERSION = '1.11';

use Bio::DB::HTS; #load XS
use strict;
use warnings;

sub new {
  my $class         = shift;
  my (%args) = @_;
  my $filename = $args{filename};
  my $warnings = $args{warnings};

  my $reader = bcf_sr_open($filename);
  die "Error getting reader" unless $reader;

  my $self = bless {
                    bcf_reader => $reader,
                    filename => $filename,
                   }, ref $class || $class;
  return $self;
}


sub num_variants {
    my $self = shift;
    return bcf_num_variants($self->{bcf_reader});
}

sub DEMOLISH {
    my $self = shift;

    if ( $self->{bcf_reader} ) {
        bcf_sr_close($self->{bcf_reader});
    }
}

1;
