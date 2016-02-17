package Bio::DB::HTS::Tabix;

use Log::Log4perl qw( :easy );

use Bio::DB::HTS; #load the XS
use Bio::DB::HTS::Tabix::Iterator;


sub new {
  my (%args) = @_;
  my $filename = $args{filename} ;

    printf("TABIX.pm filename:".$filename."\n") ;
    # filename checks
    die "Tabix region lookup requires a gzipped, tabixed bedfile" unless $filename =~ /gz$/;
    die "Filename " . $filename . " does not exist" unless -e $filename;

    $self->_htsfile = Bio::DB::HTSfile->open($filename);
    $self->_tabix_index = tbx_open($filename);
    die "Couldn't find index for file " . $filename unless $index;
    my $header = tbx_header($self->_htsfile, $self->_tabix_index);
    if( $header ) {
      $self->_header = join "", @{ $header } ;
    }

    $self->seqnames_hash = { map { $_ => 1 } @{ $self->seqnames } };

    if ( not $self->warnings ) {
        my $logger = Log::Log4perl::get_logger(__PACKAGE__);
        $logger->level($TRACE);
    }
    return;
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
            $self->log->warn("You have not specified an end, which actually means chr:start-end_of_chromosome");
        }
    }

    $self->log->trace("Fetching region $region from " . $self->filename);

    my $iter = tbx_query( $self->_tabix_index, $region );

    unless ( $iter ) {
      #this likely means the chromosome wasn't found in the tabix index, or it couldn't parse the provided region.      
      if ( not exists $self->seqnames_hash->{ $chr } ) {
        $self->log->warn("Specified chromosome '$chr' does not exist in file " . $self->filename);
      }
      else {
        die "Unable to get iterator for region '$region' in file ". $self->filename . " -- htslib couldn't parse your region string";
        }

    }

    return Bio::DB::HTS::Tabix::Iterator->new( _tabix_iter => $iter, _htsfile => $self->_htsfile, _tabix_index => $self->_tabix_index );
}

sub seqnames {
    my $self = shift;
    return tbx_seqnames($self->_tabix_index);
}

#free up memory allocated in XS code
sub DEMOLISH {
    my $self = shift;

    if ( $self->_htsfile ) {
        Bio::DB::HTSfile->close($self->_htsfile);
    }

    if ( $self->_has_tabix_index ) {
        tbx_close($self->_tabix_index);
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

=head1 COPYRIGHT

Copyright 2015 Congenica Ltd.

=head1 AUTHOR

Alex Hodgkins

=cut
