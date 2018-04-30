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

package Bio::DB::HTS::Tabix::Iterator;

use Bio::DB::HTS; #load the XS
$Bio::DB::HTS::Tabix::Iterator::VERSION = '2.11';

use strict;
use warnings;
use Scalar::Util qw/reftype/;

#this class is just a wrapper around the tabix_iter_next method,
#all the attributes it needs come from the main Tabix method

sub new {
  my $class         = shift;
  my (%args) = @_;
  my $tabix_iter = $args{_tabix_iter}; #a hts_itr_t pointer which is returned from Tabix::query
  my $htsfile = $args{_htsfile}; #an open htsFile pointer
  my $tabix_index = $args{_tabix_index};

  my $self = bless {
                    _tabix_iter => $tabix_iter,
                    _htsfile => $htsfile,
                    _tabix_index => $tabix_index,
                   }, ref $class || $class;

  return $self;

}


sub next {
    my $self = shift;

    #sometimes tabix_query doesn't return an iterator, just NULL so we have to allow
    #a null iterator
    return unless defined $self->{_tabix_iter};

    #this is an xs method
    return tbx_iter_next($self->{_tabix_iter}, $self->{_htsfile}, $self->{_tabix_index});
}

sub close {
    my $self = shift;

    #xs method
    if ( defined $self->{_tabix_iter} ) {
        tbx_iter_free($self->{_tabix_iter});
        delete $self->{_tabix_iter}; # delete once you've removed it. Prevents bad re-issuing of code
    }
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

Bio::DB::HTS::Tabix::Iterator - XS module wrapping around a tabix hts_itr_t

=head1 SYNOPSIS

You shouldn't be instantiating one of these manually it needs a load of pointers.
Usage would be through L<Bio::DB::HTS::Tabix>:

    use feature qw( say );
    use Bio::DB::HTS::Tabix;

    my $tabix = Bio::DB::HTS::Tabix->new( filename => "gerp_plus_plus_31July2014.gz" );

    say $tabix->header;
    my $iter = $tabix->query("1:4000005-4000009");

    while ( my $n = $iter->next ) {
        say $n;
    }

=head1 DESCRIPTION

This is returned from L<Bio::DB::HTS::Tabix>, the only method you need to care about is 'next'.

Don't go importing this and calling new on it if you value your sanity, it won't work.

=head2 Methods

=over 12

=item C<next>

Returns a string with the line from the tabix iterator

=back

=head1 AUTHOR

Alex Hodgkins
Rishi Nag E<lt>rishi@ebi.ac.ukE<gt>

=cut
