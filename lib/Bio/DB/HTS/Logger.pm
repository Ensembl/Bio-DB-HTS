package Bio::DB::HTS::Logger;

use Mouse::Role;
use Log::Log4perl;

has 'logger' => (
    is      => 'rw',
    isa     => 'Log::Log4perl::Logger',
    lazy    => 1,
    default => sub { return Log::Log4perl->get_logger(ref($_[0])) }
);

sub log {
    my $self = shift;
    my $cat = shift;
    if ($cat && $cat =~ m/^(\.|::)/) {
        return Log::Log4perl->get_logger(ref($self) . $cat);
    } elsif($cat)  {
        return Log::Log4perl->get_logger($cat);
    } else {
        return $self->logger;
    }
}

1;

__END__
