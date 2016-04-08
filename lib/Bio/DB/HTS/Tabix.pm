
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

=cut

package Bio::DB::HTS::Tabix;

use Bio::DB::HTS; #load the XS
use Bio::DB::HTS::Tabix::Iterator;
$Bio::DB::HTS::Tabix::VERSION = '1.12';
use strict;
use warnings;


sub new {
  my $class         = shift;
  my (%args) = @_;
  my $filename = $args{filename};
  my $warnings = $args{warnings};

  # filename checks
  die "Tabix region lookup requires a gzipped, tabixed bedfile" unless $filename =~ /gz$/;

  my $htsfile = Bio::DB::HTSfile->open($filename);
  my $tabix_index = tbx_open($filename);
  die "Couldn't find index for file " . $filename unless $tabix_index;
  my $header = tbx_header($htsfile, $tabix_index);
  if( $header )
  {
    $header = join "\n", @{ $header } ;
  }

  my $self = bless {
                    _htsfile => $htsfile,
                    _filename => $filename,
                    _tabix_index => $tabix_index,
                    _header=> $header,
                   }, ref $class || $class;

  return $self;
}

sub query {
    my ( $self, $region ) = @_;

    die "Please provide a region" unless defined $region;
    my ( $chr, $start, $end ) = $region =~ /^([^:]+)(?::(\d+))?(?:-(\d+))?$/;
    unless ( defined $chr ) {
        die "You must specify a region in the format chr, chr:start or chr:start-end";
    }

    if ( defined $end ) {
        die "End in $region is less than the start" if $end < $start;
    }
    else {
        if ( defined $start ) {
            #$self->log->warn("You have not specified an end, which actually means chr:start-end_of_chromosome");
        }
    }

    my $iter = tbx_query( $self->{_tabix_index}, $region );

    unless ( $iter ) {
      #this likely means the chromosome wasn't found in the tabix index, or it couldn't parse the provided region.
      my $seqnames_hash = { map { $_ => 1 } @{ $self->seqnames } };
      if ( not exists $seqnames_hash->{ $chr } ) {
        #$self->log->warn("Specified chromosome '$chr' does not exist in file " . $self->_filename);
      }
      else {
        die "Unable to get iterator for region '$region' in file ". $self->{_filename} . " -- htslib couldn't parse your region string";
        }

    }

    return Bio::DB::HTS::Tabix::Iterator->new( _tabix_iter => $iter, _htsfile => $self->{_htsfile}, _tabix_index => $self->{_tabix_index} );
}

sub seqnames {
    my $self = shift;
    return tbx_seqnames($self->{_tabix_index});
}

sub header {
    my $self = shift;
    return $self->{_header};
}

sub header_array {
    my $self = shift;
    my @lines = split(/\n/,$self->{_header});
    return @lines ;
}


#free up memory allocated in XS code
sub close {
    my $self = shift;

    if ( $self->{_htsfile} ) {
        Bio::DB::HTSfile::close($self->{_htsfile});
    }

    if ( $self->{_tabix_index} ) {
        tbx_close($self->{_tabix_index});
    }
}

1;

__END__

=head1 NAME

Bio::DB::HTS::Tabix - Object oriented access to the underlying tbx C methods

=head1 SYNOPSIS

    use feature qw( say );
    use Bio::DB::HTS::Tabix;

    my $tabix = Bio::DB::HTS::Tabix->new( filename => "gerp_plus_plus_31July2014.gz" );

    say $tabix->header;
    my $iter = $tabix->query("1:4000005-4000009");

    while ( my $n = $iter->next ) {
        say $n;
    }

=head1 DESCRIPTION

A high level object oriented interface to the htslib tabix (tbx.h) api. Currently it only supports
retrieving regions from a tabixed file, because that's all I needed it for.

=head2 Attributes

=over 12

=item C<filename>

The gzipped file you want to query. Must have a filename.tbi (the index is not created automatically)

=back

=head2 Methods

=over 12

=item C<header>

Returns all the header lines as a single scalar from the tabixed file

=item C<query>

Takes a single region like: '1:4000005-4000009' or '12:5000000'
Note: this works exactly the same way as the tabix executable,
so '12:5000000' actually means get all results from position 5,000,000
up to the very end of the chromosome. To get results only at position
5,000,000 you should do '12:5000000-5000001'

Returns a L<Bio::DB::HTS::Tabix::Iterator> for the specified region

=item C<seqnames>

Returns an array ref of chromosomes that are in the indexed file

=back

=head1 AUTHOR

Alex Hodgkins
Rishi Nag E<lt>rishi@ebi.ac.ukE<gt>

=cut
