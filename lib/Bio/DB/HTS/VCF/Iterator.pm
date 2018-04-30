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

package Bio::DB::HTS::VCF::Iterator;

use Bio::DB::HTS; #load the XS
$Bio::DB::HTS::VCF::Iterator::VERSION = '2.11';

use strict;
use warnings;
use Scalar::Util qw/reftype/;

sub new {
  my $class         = shift;
  my (%args) = @_;
  my $iter = $args{iter}; 
  my $htsfile = $args{file}; # an open htsFile pointer
  my $header = $args{header};
  my $index = $args{index}; # either tabix (VCF) or hts index (BCF), depending on the file type

  my $self = bless {
                    _iter => $iter,
                    _htsfile => $htsfile,
		    _header  => $header,
                    _index => $index,
                   }, ref $class || $class;

  return $self;

}

sub next {
    my $self = shift;

    # sometimes *_query doesn't return an iterator, just NULL
    # so we have to allow a null iterator
    return unless defined $self->{_iter} and defined $self->{_header};

    return iter_next($self->{_iter}, $self->{_htsfile}, $self->{_header}, $self->{_index});
}

sub close {
  my $self = shift;

  return unless defined $self->{_iter};
  
  iter_free($self->{_iter});
  delete $self->{_iter}; # delete once you've removed it. Prevents bad re-issuing of code
}

sub DESTROY {
    my $self = shift;
    return if reftype($self) ne 'HASH';
    $self->close();
    return;
}

1;

__END__

=head1 NAME

Bio::DB::HTS::VCF::Iterator - XS module wrapping around a hts_itr_t

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 Methods

=over 12

=item C<next>

Returns a string with the line from the iterator

=back

=head1 AUTHOR

Alessandro Vullo E<lt>avullo@cpan.orgE<gt>

=cut
