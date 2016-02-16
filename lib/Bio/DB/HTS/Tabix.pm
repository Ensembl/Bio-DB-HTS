package Bio::DB::HTS::Tabix;

use Mouse;
use Log::Log4perl qw( :easy );

use Bio::DB::HTS; #load the XS
use Bio::DB::HTS::Tabix::Iterator;

with 'Bio::DB::HTS::Logger';

has 'filename' => (
    is       => 'ro',
    isa      => 'Str',
    required => 1,
);

#pointer to a htsFile
has '_htsfile' => (
    is        => 'ro',
    isa       => 'Bio::DB::HTSfile',
    builder   => '_build__htsfile',
    predicate => '_has_htsfile',
    lazy      => 1,
);

sub _build__htsfile {
    my $self = shift;

    die "Filename " . $self->filename . " does not exist" unless -e $self->filename;

    return Bio::DB::HTSfile->open($self->filename);
}

has '_tabix_index' => (
    is        => 'ro',
    isa       => 'tbx_tPtr',
    builder   => '_build__tabix_index',
    predicate => '_has_tabix_index',
    lazy      => 1,
);

sub _build__tabix_index {
    my $self = shift;

    #make sure the htsfile is instantiated
    $self->_htsfile;

    my $index = tbx_open($self->filename);

    die "Couldn't find index for file " . $self->filename unless $index;

    return $index;
}

has 'header' => (
    is      => 'ro',
    isa     => 'Str',
    builder => '_build_header',
    lazy    => 1,
);

sub _build_header {
    my $self = shift;

    #returns an arrayref
    my $header = tbx_header($self->_htsfile, $self->_tabix_index);

    return unless $header;

    return join "", @{ $header };
}

has 'seqnames_hash' => (
    is      => 'ro',
    isa     => 'HashRef',
    builder => '_build_seqnames_hash',
    lazy    => 1,
);

sub _build_seqnames_hash {
    my $self = shift;

    return { map { $_ => 1 } @{ $self->seqnames } };
}

#set to false to only output error messages
has 'warnings' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 1,
);

sub BUILD {
    my ( $self ) = @_;

    die "Tabix region lookup requires a gzipped, tabixed bedfile"
        unless $self->filename =~ /gz$/;

    #fetch the header of the file, which will in turn open the tabix index and the file
    $self->header;

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

    if ( $self->_has_htsfile ) {
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

=item C<warnings>

Set to 0 to turn off all the warnings. Default is on

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
