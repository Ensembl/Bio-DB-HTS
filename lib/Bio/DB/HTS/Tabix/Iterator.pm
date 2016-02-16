package Bio::DB::HTS::Tabix::Iterator;

use Mouse;
use Bio::DB::HTS; #load the XS

#this class is just a wrapper around the tabix_iter_next method,
#all the attributes it needs come from the main Tabix method

#a hts_itr_t pointer which is returned from Tabix::query
has "_tabix_iter" => (
    is  => 'ro',
    isa => 'Maybe[hts_itr_tPtr]',
);

#an open htsFile pointer
has '_htsfile' => (
    is       => 'ro',
    isa      => 'Bio::DB::HTSfile',
    required => 1,
);

has '_tabix_index' => (
    is       => 'ro',
    isa      => 'tbx_tPtr',
    required => 1,
);

sub next {
    my $self = shift;

    #sometimes tabix_query doesn't return an iterator, just NULL so we have to allow
    #a null iterator
    return unless defined $self->_tabix_iter;

    #this is an xs method
    return tbx_iter_next($self->_tabix_iter, $self->_htsfile, $self->_tabix_index);
}

sub DEMOLISH {
    my $self = shift;

    #xs method
    if ( defined $self->_tabix_iter ) {
        tbx_iter_free($self->_tabix_iter);
    }
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

=head1 COPYRIGHT

Copyright 2015 Congenica Ltd.

=head1 AUTHOR

Alex Hodgkins

=cut
