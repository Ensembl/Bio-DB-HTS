package Bio::DB::HTS::VCF;

use Mouse;
use Log::Log4perl qw( :easy );

use Bio::DB::HTS; #load XS
with 'Bio::DB::HTS::Logger';

has 'filename' => (
    is       => 'ro',
    isa      => 'Str',
    required => 1,
);

has '_bcf_reader' => (
    is        => 'ro',
    isa       => 'bcf_srs_tPtr',
    builder   => '_build__bcf_reader',
    predicate => '_has_bcf_reader',
    lazy      => 1,
);

sub _build__bcf_reader {
    my $self = shift;

    #we use the bcf synced reader because i feel like it
    my $reader = bcf_sr_open($self->filename);

    die "Error getting reader" unless $reader;

    return $reader;
}

sub num_variants {
    my $self = shift;

    return bcf_num_variants($self->_bcf_reader);
}

sub DEMOLISH {
    my $self = shift;

    if ( $self->_has_bcf_reader ) {
        bcf_sr_close($self->_bcf_reader);
    }
}

1;
